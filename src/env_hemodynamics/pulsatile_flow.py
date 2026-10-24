"""
Stage 1 - Step 07: Pulsatile Non-Newtonian Microvascular Blood Flow Solver.
Implements the In-Vivo Pries-Secomb Fahraeus-Lindqvist Viscosity Law,
Cerebral Perfusion Pressure (CPP = MAP - ICP) driving dynamics, and
Darcy-Brinkman hydraulic resistance across an occluded penetrating arteriole.

Part of the Q-DIAPEDESIS Framework (OSF: https://osf.io/htswk | DOI: 10.5281/zenodo.23599020).
"""

from dataclasses import dataclass
from pathlib import Path
import numpy as np
import yaml

from .vessel_geometry import PenetratingArterioleGeometry


MMHG_TO_PA = 133.322368


@dataclass
class HemodynamicSample:
    """Instantaneous hemodynamic field variables along the penetrating arteriole."""
    time_s: float
    cardiac_phase_rad: float
    is_systole: bool
    instantaneous_cpp_mmHg: float
    volumetric_flow_um3_s: float
    axial_velocity_mean_um_s: np.ndarray
    wall_shear_stress_Pa: np.ndarray
    effective_viscosity_mPa_s: np.ndarray
    collateral_skimming_fraction: float


class PriesSecombPulsatileFlow:
    """
    Solves pulsatile non-Newtonian blood flow in cerebral penetrating arterioles (15-50 um)
    coupled with porous thrombus hydraulic resistance and collateral branch skimming.
    """

    def __init__(
        self,
        geometry: PenetratingArterioleGeometry | None = None,
        config_path: Path | str | None = None,
    ):
        if config_path is None:
            root_dir = Path(__file__).resolve().parent.parent.parent
            config_path = root_dir / "configs" / "vessel_mra_config.yaml"

        with open(config_path, "r", encoding="utf-8") as f:
            self.cfg = yaml.safe_load(f)

        self.geom = geometry if geometry is not None else PenetratingArterioleGeometry(config_path)

        hemo = self.cfg["hemodynamics_pries_secomb"]
        self.map_mmHg = float(hemo["map_mmHg"])
        self.icp_mmHg = float(hemo["icp_mmHg"])
        self.cpp_baseline_mmHg = self.map_mmHg - self.icp_mmHg
        self.heart_rate_bpm = float(hemo["heart_rate_bpm"])
        self.omega_cardiac = 2.0 * np.pi * (self.heart_rate_bpm / 60.0)
        self.pulse_amp_ratio = float(hemo["pulsatile_amplitude_ratio"])
        self.hct = float(hemo["discharge_hematocrit"])
        self.mu_plasma_mPa_s = float(hemo["plasma_viscosity_mPa_s"])
        self.blood_density = float(hemo["blood_density_kg_m3"])

    def compute_in_vivo_viscosity_mPa_s(self, diameter_um: np.ndarray) -> np.ndarray:
        """
        Computes the apparent in-vivo blood viscosity using the empirical Pries-Secomb
        Fahraeus-Lindqvist microvascular formulation (Circ. Res. 1994; Am. J. Physiol. 2005).
        """
        D = np.clip(diameter_um, 4.0, 300.0)

        # In-vivo viscosity at Hct = 0.45 (Pries et al., Circ. Res. 1994). The in-vitro
        # term 220*exp(-1.3D) was used here before H1 (ELN #002); the in-vivo law uses 6*exp(-0.085D).
        mu_045 = 6.0 * np.exp(-0.085 * D) + 3.2 - 2.44 * np.exp(-0.06 * (D ** 0.645))

        # Curvature parameter C(D) describing hematocrit dependence
        term_exp = np.exp(-0.075 * D)
        denom = 1.0 + 10.0 ** (-11.0) * (D ** 12.0)
        C_D = (0.8 + term_exp) * (-1.0 + 1.0 / denom) + (1.0 / denom)

        # In-vivo endothelial surface layer factor z = (D / (D - 1.1))^2 (Pries et al. 1994);
        # the denominator is floored at 1.5 um so the factor stays finite for very narrow channels.
        d_eff_ratio = (D / np.maximum(D - 1.1, 1.5)) ** 2
        relative_visc = (
            1.0
            + (mu_045 - 1.0)
            * (((1.0 - self.hct) ** C_D - 1.0) / ((1.0 - 0.45) ** C_D - 1.0))
            * d_eff_ratio
        ) * d_eff_ratio

        return self.mu_plasma_mPa_s * np.maximum(relative_visc, 1.05)

    def instantaneous_cpp_mmHg(self, time_s: float, cpp_override_mmHg: float | None = None) -> tuple[float, float, bool]:
        """
        Evaluates arterial pulsatile waveform (fundamental + dicrotic second harmonic)
        driven by Cerebral Perfusion Pressure (CPP = MAP - ICP).
        """
        base_cpp = self.cpp_baseline_mmHg if cpp_override_mmHg is None else float(cpp_override_mmHg)
        phase = (self.omega_cardiac * time_s) % (2.0 * np.pi)

        # Realistic physiological arterial pulse wave shape
        waveform = 0.78 * np.sin(phase) + 0.22 * np.sin(2.0 * phase - 0.45)
        cpp_t = base_cpp * (1.0 + self.pulse_amp_ratio * waveform)
        is_systole = bool(phase < np.pi)
        return float(cpp_t), float(phase), is_systole

    def solve_flow_step(
        self,
        time_s: float,
        cpp_override_mmHg: float | None = None,
        robot_induced_shear_Pa: float = 0.0,
    ) -> HemodynamicSample:
        """
        Computes axial velocity, wall shear stress (tau_wall), and collateral skimming
        across the penetrating arteriole at time `time_s`.
        """
        state = self.geom.state
        z_m = state.z_um * 1.0e-6
        dz_m = float(z_m[1] - z_m[0])

        r_lumen_um = np.maximum(state.r_lumen_um, 1.25)
        d_lumen_um = 2.0 * r_lumen_um
        r_lumen_m = r_lumen_um * 1.0e-6

        mu_mPa_s = self.compute_in_vivo_viscosity_mPa_s(d_lumen_um)
        mu_Pa_s = mu_mPa_s * 1.0e-3

        # Local Poiseuille + Darcy porous thrombus hydraulic resistance per axial segment
        poiseuille_res = (8.0 * mu_Pa_s * dz_m) / (np.pi * (r_lumen_m ** 4))
        area_m2 = np.pi * (r_lumen_m ** 2)
        darcy_res = (mu_Pa_s * dz_m * (state.occlusion_fraction ** 2)) / (
            state.darcy_permeability_m2 * area_m2 + 1.0e-24
        )
        total_hydraulic_res = float(np.sum(poiseuille_res + darcy_res))

        cpp_t_mmHg, phase_rad, is_systole = self.instantaneous_cpp_mmHg(time_s, cpp_override_mmHg)
        # Microvascular pressure drop across a 500 um penetrating arteriole segment (~8% of total CPP drop)
        segment_delta_P_Pa = 0.08 * cpp_t_mmHg * MMHG_TO_PA

        q_m3_s = segment_delta_P_Pa / max(total_hydraulic_res, 1.0e-12)
        q_um3_s = q_m3_s * 1.0e18

        u_mean_m_s = q_m3_s / area_m2
        u_mean_um_s = u_mean_m_s * 1.0e6

        # Hydrodynamic wall shear stress: tau_flow = 4 * mu * u_mean / R + robot rotational shear
        tau_hydro_Pa = (4.0 * mu_Pa_s * u_mean_m_s) / r_lumen_m
        tau_total_Pa = tau_hydro_Pa + float(robot_induced_shear_Pa)

        # Collateral bifurcation skimming risk (increases with peak occlusion & systolic surge)
        peak_occ = float(np.max(state.occlusion_fraction))
        systolic_factor = 1.25 if is_systole else 0.75
        collateral_skimming = float(np.clip(0.12 * peak_occ * systolic_factor, 0.005, 0.45))

        return HemodynamicSample(
            time_s=float(time_s),
            cardiac_phase_rad=phase_rad,
            is_systole=is_systole,
            instantaneous_cpp_mmHg=cpp_t_mmHg,
            volumetric_flow_um3_s=float(q_um3_s),
            axial_velocity_mean_um_s=u_mean_um_s,
            wall_shear_stress_Pa=tau_total_Pa,
            effective_viscosity_mPa_s=mu_mPa_s,
            collateral_skimming_fraction=collateral_skimming,
        )
