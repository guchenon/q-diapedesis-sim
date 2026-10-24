"""
Stage 1 - Steps 08 & 09: Endothelial Glycocalyx Integrity Dynamics,
Mechanosensitive Neutrophil NETosis Risk (Piezo1 Shear Activation),
Pulsed Micro-Dialysis Acidotic pH Buffering, and Pennes Bioheat Model.

Part of the Q-DIAPEDESIS Framework (OSF: https://osf.io/htswk | DOI: 10.5281/zenodo.23599020).
"""

from dataclasses import dataclass
from pathlib import Path
import numpy as np
import yaml

from .shear_dose_response import hill_fraction, load_dose_response_constants


@dataclass
class PhysiologicalState:
    """Instantaneous neurocritical-care physiological safety state."""
    time_s: float
    glycocalyx_integrity_ratio: float     # I_glyco in [0.0, 1.0] (Target >= 0.92)
    mean_glycocalyx_thickness_um: float   # Thickness [um] (Baseline 0.50 um)
    peak_wall_shear_Pa: float             # Max total wall shear stress [Pa] (Ceiling <= 5.50 Pa)
    netosis_activation_index: float       # Cumulative shear-induced NETs risk [0.0, 1.0]
    distal_efflux_pH: float               # Intravascular pH in post-occlusive zone (Safe >= 7.15)
    brain_temp_elevation_C: float         # Magnetothermal Delta T [deg C] (Ceiling < 0.30 C)
    cpp_surge_rate_mmHg_s: float          # Rate of perfusion pressure rise [mmHg/s]


class GlycocalyxNetosisDynamics:
    """
    Simulates microvascular endothelial glycocalyx shedding (Step 08) and
    Module B3/B7 physiological hazards (Step 09):
    1. Shear-induced heparan sulfate / syndecan-1 glycocalyx degradation when tau_wall > tau_crit.
    2. Piezo1-mediated neutrophil extracellular trap (NETosis) activation under excessive shear.
    3. Pulsed micro-dialysis buffering of stagnant lactic acidosis (pH 6.65 -> >= 7.15).
    4. Pennes bioheat perfusion thermal balance under bimodal alternating magnetic fields.
    """

    def __init__(
        self,
        cbf_config_path: Path | str | None = None,
        vessel_config_path: Path | str | None = None,
    ):
        root_dir = Path(__file__).resolve().parent.parent.parent
        if cbf_config_path is None:
            cbf_config_path = root_dir / "configs" / "neurocritical_cbf_config.yaml"
        if vessel_config_path is None:
            vessel_config_path = root_dir / "configs" / "vessel_mra_config.yaml"

        with open(cbf_config_path, "r", encoding="utf-8") as f:
            self.cbf_cfg = yaml.safe_load(f)
        with open(vessel_config_path, "r", encoding="utf-8") as f:
            self.vessel_cfg = yaml.safe_load(f)

        thresh = self.cbf_cfg["cbf_safety_thresholds"]
        self.tau_crit_Pa = float(thresh["max_wall_shear_stress_Pa"])
        self.min_glyco_ratio = float(thresh["min_glycocalyx_integrity_ratio"])
        self.max_delta_T_C = float(thresh["max_brain_temp_elevation_C"])

        staged = self.cbf_cfg["staged_cpp_recanalization"]
        self.max_cpp_surge = float(staged["max_pressure_surge_mmHg_per_sec"])
        self.pulsed_duty_cycle = float(staged["pulsed_pH_buffering_duty_cycle"])
        self.min_safe_pH = float(staged["min_safe_efflux_pH"])

        geom = self.vessel_cfg["vessel_geometry"]
        self.glyco_baseline_um = float(geom["glycocalyx_baseline_thickness_um"])

        thrombus = self.vessel_cfg["thrombus_and_microenvironment"]
        self.pH_initial = float(thrombus["post_occlusive_pH_initial"])
        self.pH_target = float(thrombus["physiological_pH_target"])

        # Continuous shear dose-response constants (D3, ELN #003; implementation_constants)
        dose = load_dose_response_constants(cbf_config_path)
        self.hill_tau50_Pa = float(dose["hill_tau50_Pa"])
        self.hill_n = float(dose["hill_n"])
        self.k_glyco_mech_max_inv_s = float(dose["glyco_k_mech_max_inv_s"])
        self.k_glyco_osc_inv_s = float(dose["glyco_k_oscillatory_inv_s"])
        self.k_netosis_max_inv_s = float(dose["netosis_k_max_inv_s"])
        self.k_netosis_clearance_inv_s = float(dose["netosis_clearance_inv_s"])

        self.reset()

    def reset(self) -> PhysiologicalState:
        """Resets dynamic physiological variables to baseline pre-recanalization values."""
        self.current_time_s = 0.0
        self.glyco_ratio = 1.0
        self.netosis_index = 0.0
        self.current_pH = self.pH_initial
        self.delta_T_C = 0.0
        return self.get_state(peak_wall_shear_Pa=1.2, cpp_surge_rate_mmHg_s=0.0)

    def step_physiology(
        self,
        dt_s: float,
        wall_shear_profile_Pa: np.ndarray,
        b_field_amplitude_mT: float,
        drive_frequency_Hz: float,
        blood_perfusion_velocity_um_s: float,
        pulsed_buffering_active: bool = True,
        cpp_surge_rate_mmHg_s: float = 0.0,
    ) -> PhysiologicalState:
        """
        Advances the endothelial glycocalyx, NETosis, pH buffering, and bioheat state by `dt_s`.
        """
        peak_shear_Pa = float(np.max(wall_shear_profile_Pa))
        mean_shear_Pa = float(np.mean(wall_shear_profile_Pa))

        # 1. Step 08: Endothelial glycocalyx shedding (continuous Hill dose-response, D3)
        h_peak = float(hill_fraction(peak_shear_Pa, self.hill_tau50_Pa, self.hill_n))
        sub_threshold_wear = self.k_glyco_osc_inv_s * (mean_shear_Pa / self.hill_tau50_Pa) ** 2
        d_glyco_dt = -(self.k_glyco_mech_max_inv_s * h_peak + sub_threshold_wear) * self.glyco_ratio
        self.glyco_ratio = float(np.clip(self.glyco_ratio + d_glyco_dt * dt_s, 0.0, 1.0))

        # 2. Step 09 (Module B3a): Piezo1 NETosis, continuous in shear (Baratchi et al. 2024), with clearance
        d_net_dt = (
            self.k_netosis_max_inv_s * h_peak * (1.0 - self.netosis_index)
            - self.k_netosis_clearance_inv_s * self.netosis_index
        )
        self.netosis_index = float(np.clip(self.netosis_index + d_net_dt * dt_s, 0.0, 1.0))

        # 3. Step 09 (Module B3b): Pulsed Micro-Dialysis Acidotic pH Buffering
        # Gradually neutralizes post-occlusive lactic acidosis (pH 6.65 -> 7.35)
        buffering_efficiency = self.pulsed_duty_cycle if pulsed_buffering_active else 0.05
        flow_factor = np.clip(blood_perfusion_velocity_um_s / 500.0, 0.1, 2.0)
        tau_pH_s = 2.5 / (buffering_efficiency * flow_factor)
        d_pH_dt = (self.pH_target - self.current_pH) / tau_pH_s
        self.current_pH = float(np.clip(self.current_pH + d_pH_dt * dt_s, self.pH_initial, self.pH_target))

        # 4. Step 09 (Module B7): Pennes Bioheat Equation for Magnetothermal Brain Safety
        # Low-frequency AMF power density scales with f * B^2; blood perfusion convects heat away
        q_mag_mW_cm3 = 1.8e-6 * drive_frequency_Hz * (b_field_amplitude_mT ** 2)
        perfusion_cooling_inv_s = 0.08 * (1.0 + blood_perfusion_velocity_um_s / 1000.0)
        d_temp_dt = q_mag_mW_cm3 - perfusion_cooling_inv_s * self.delta_T_C
        self.delta_T_C = float(np.clip(self.delta_T_C + d_temp_dt * dt_s, 0.0, 5.0))

        self.current_time_s += float(dt_s)
        return self.get_state(peak_wall_shear_Pa=peak_shear_Pa, cpp_surge_rate_mmHg_s=cpp_surge_rate_mmHg_s)

    def get_state(self, peak_wall_shear_Pa: float, cpp_surge_rate_mmHg_s: float) -> PhysiologicalState:
        """Packages current physiological metrics into a structured state object."""
        return PhysiologicalState(
            time_s=self.current_time_s,
            glycocalyx_integrity_ratio=self.glyco_ratio,
            mean_glycocalyx_thickness_um=self.glyco_ratio * self.glyco_baseline_um,
            peak_wall_shear_Pa=float(peak_wall_shear_Pa),
            netosis_activation_index=self.netosis_index,
            distal_efflux_pH=self.current_pH,
            brain_temp_elevation_C=self.delta_T_C,
            cpp_surge_rate_mmHg_s=float(cpp_surge_rate_mmHg_s),
        )
