"""
Stage 2 - Steps 11 & 12: Leukocyte-Inspired Margination Drift and
Selectin Catch-Slip Bond Wall-Rolling Mechanics.

Grounded in:
- Marshall et al., Nature 2003 (doi:10.1038/nature01605)
- Barsegov & Thirumalai, PNAS 2005 (doi:10.1073/pnas.0406938102)
- Alapan, Bozuyuk, Erkoc, Karacakol & Sitti, Sci. Robot. 2020 (doi:10.1126/scirobotics.aba5726)
- Muller et al., Sci. Rep. 2014 (doi:10.1038/srep04871)

Part of the Q-DIAPEDESIS Framework (OSF registration: https://osf.io/8mysz | project: https://osf.io/hnstf | DOI: 10.5281/zenodo.23599019).
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict
import numpy as np
import yaml


@dataclass
class MarginationRollingState:
    """Instantaneous radial margination and selectin catch-bond rolling state."""
    mean_radial_position_um: float
    wall_proximal_fraction: float
    tether_force_pN: float
    k_off_inv_s: float
    mean_bond_lifetime_s: float
    adhesion_locked: bool
    axial_rolling_velocity_um_s: float


class LeukocyteCatchBondMargination:
    """
    Models Step 11 (Shear-induced red blood cell margination drift toward the endothelial wall)
    and Step 12 (Selectin catch-slip bond mechanics for pulsatile upstream/downstream rolling).

    Epistemic note (Audit Table Row 3a/3b):
    Marshall et al. (Nature 2003) and Barsegov & Thirumalai (PNAS 2005) establish catch-bond
    ceilings of ~10 pN (monomeric sPSGL-1) and ~20 pN (dimeric PSGL-1). The nominal parameter
    F* = 15.0 pN is an investigator-defined midpoint (Status B), evaluated with mandatory
    sensitivity sweeps across F* in {10.0, 15.0, 20.0} pN against a pure slip-bond baseline.
    """

    def __init__(
        self,
        config_path: Path | str | None = None,
        critical_force_override_pN: float | None = None,
        bond_mode: str = "catch_slip",
    ):
        if config_path is None:
            root_dir = Path(__file__).resolve().parent.parent.parent
            config_path = root_dir / "configs" / "vessel_mra_config.yaml"

        with open(config_path, "r", encoding="utf-8") as f:
            self.cfg = yaml.safe_load(f)

        geom = self.cfg["vessel_geometry"]
        self.default_radius_um = 0.5 * float(geom["default_diameter_um"])

        bio = self.cfg["biomimicry_layer"]
        self.particle_radius_um = float(bio["particle_radius_um"])
        self.margination_gain = float(bio["margination_drift_gain"])

        cb = bio["selectin_catch_bond"]
        self.k_slip_0 = float(cb["k_slip_0_inv_s"])
        self.k_catch_0 = float(cb["k_catch_0_inv_s"])
        self.f_crit_pN = (
            float(critical_force_override_pN)
            if critical_force_override_pN is not None
            else float(cb["critical_force_pN"])
        )
        self.kbt_pN_nm = float(cb["thermal_energy_scale_pN_nm"])

        # Transition-state distances grounded in Barsegov & Thirumalai (PNAS 2005)
        self.x_catch_nm = 1.45
        self.x_slip_nm = 0.38
        self.bond_mode = bond_mode  # "catch_slip" or "pure_slip_baseline"

    def compute_tether_force_pN(self, wall_shear_stress_Pa: float) -> float:
        """
        Converts local hydrodynamic wall shear stress [Pa] into tensile tether force [pN]
        acting on a wall-proximal spherical microroller (Yago et al., J. Cell Biol. 2004).
        1 Pa = 10 dyn/cm^2; for R_p = 0.75 um, lever factor ~ 3.2 pN/Pa.
        """
        lever_factor_pN_per_Pa = 5.6 * (self.particle_radius_um ** 2)
        return float(max(0.0, wall_shear_stress_Pa) * lever_factor_pN_per_Pa)

    def dissociation_rate_inv_s(self, force_pN: float) -> float:
        """
        Computes the two-pathway catch-slip bond off-rate k_off(F) [s^-1]:
            k_off(F) = k_catch_0 * exp(-F * x_c / k_B T) + k_slip_0 * exp((F - F*) * x_s / k_B T)
        In the pure slip baseline, the catch pathway is disabled and k_off increases monotonically.
        """
        F = float(max(0.0, force_pN))
        if self.bond_mode == "pure_slip_baseline":
            k_off = (self.k_slip_0 * 3.5) * np.exp((F * self.x_slip_nm) / self.kbt_pN_nm)
            return float(np.clip(k_off, 0.05, 500.0))

        # Two-pathway catch-slip kinetics: lifetime increases up to F* then transitions to slip
        catch_branch = self.k_catch_0 * np.exp(-(F * self.x_catch_nm) / self.kbt_pN_nm)
        slip_branch = self.k_slip_0 * np.exp(((F - self.f_crit_pN) * self.x_slip_nm) / self.kbt_pN_nm)
        k_off = catch_branch + slip_branch
        return float(np.clip(k_off, 0.05, 500.0))

    def step_radial_margination(
        self,
        radial_positions_um: np.ndarray,
        lumen_radius_um: float,
        wall_shear_stress_Pa: float,
        dt_s: float,
    ) -> np.ndarray:
        """
        Step 11: Updates radial particle coordinates r_i in [0, R_lumen] under shear-induced
        red blood cell core exclusion (margination drift toward the cell-free wall layer).
        """
        r = np.asarray(radial_positions_um, dtype=float)
        r_max = max(2.0, float(lumen_radius_um) - self.particle_radius_um)

        # Normalized distance from vessel wall (1 at centerline, 0 at wall)
        wall_dist_norm = np.clip((r_max - r) / r_max, 0.0, 1.0)
        shear_drive = np.clip(wall_shear_stress_Pa / 2.5, 0.2, 3.0)

        dr_dt = self.margination_gain * shear_drive * wall_dist_norm * r_max
        r_next = np.clip(r + dr_dt * float(dt_s), 0.0, r_max)
        return r_next

    def evaluate_rolling_step(
        self,
        radial_positions_um: np.ndarray,
        lumen_radius_um: float,
        wall_shear_stress_Pa: float,
        is_systole: bool,
        magnetic_drive_freq_Hz: float,
    ) -> MarginationRollingState:
        """
        Step 12: Computes catch-bond locking during systolic shear surges and surface
        microrolling translation velocity during diastolic relaxation.
        """
        r_max = max(2.0, float(lumen_radius_um) - self.particle_radius_um)
        wall_band_um = 0.85 * r_max
        wall_proximal_frac = float(np.mean(radial_positions_um >= wall_band_um))

        f_tether_pN = self.compute_tether_force_pN(wall_shear_stress_Pa)
        k_off = self.dissociation_rate_inv_s(f_tether_pN)
        lifetime_s = 1.0 / k_off

        # Lock onto endothelial wall during systole when catch-bond lifetime is prolonged
        adhesion_locked = bool(is_systole and (lifetime_s >= 0.35) and (wall_proximal_frac >= 0.5))

        # Rotational surface rolling velocity: v_roll = 2 * pi * f * R_p * efficiency
        circumference_um = 2.0 * np.pi * self.particle_radius_um
        if adhesion_locked:
            # Systolic catch-bond hold prevents downstream slip; minimal creep
            v_roll_um_s = 0.12 * circumference_um * float(magnetic_drive_freq_Hz) * wall_proximal_frac
        else:
            # Diastolic active surface rolling along the low-velocity boundary layer
            traction_eff = 0.65 * wall_proximal_frac * np.clip(lifetime_s / 0.8, 0.2, 1.2)
            v_roll_um_s = circumference_um * float(magnetic_drive_freq_Hz) * float(traction_eff)

        return MarginationRollingState(
            mean_radial_position_um=float(np.mean(radial_positions_um)),
            wall_proximal_fraction=wall_proximal_frac,
            tether_force_pN=f_tether_pN,
            k_off_inv_s=k_off,
            mean_bond_lifetime_s=lifetime_s,
            adhesion_locked=adhesion_locked,
            axial_rolling_velocity_um_s=float(v_roll_um_s),
        )

    def run_catch_force_sensitivity_sweep(
        self,
        shear_stress_Pa: float = 4.2,
    ) -> Dict[str, float]:
        """
        Executes the preregistered Audit Table Row 3b sensitivity sweep across
        F* in {10.0, 15.0, 20.0} pN and compares bond lifetime against a pure slip baseline.
        """
        results: Dict[str, float] = {}
        for f_star in (10.0, 15.0, 20.0):
            model = LeukocyteCatchBondMargination(
                critical_force_override_pN=f_star,
                bond_mode="catch_slip",
            )
            f_pN = model.compute_tether_force_pN(shear_stress_Pa)
            lifetime = 1.0 / model.dissociation_rate_inv_s(f_pN)
            results[f"lifetime_s_at_Fstar_{int(f_star)}pN"] = float(lifetime)

        slip_model = LeukocyteCatchBondMargination(bond_mode="pure_slip_baseline")
        f_pN = slip_model.compute_tether_force_pN(shear_stress_Pa)
        results["lifetime_s_pure_slip_baseline"] = float(1.0 / slip_model.dissociation_rate_inv_s(f_pN))
        return results
