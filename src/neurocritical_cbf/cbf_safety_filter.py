"""
Stage 1 - Step 10: Zeroing Control Barrier Function Quadratic Program (ZCBF-QP)
Physiological Safety Filter and 4-Stage CPP Recanalization Controller.

Grounded in:
- Ames, Xu, Grizzle & Tabuada, IEEE Trans. Autom. Control 2017 (doi:10.1109/TAC.2016.2638961)
- Ames et al., 18th European Control Conference (ECC) 2019 (doi:10.23919/ECC.2019.8796030)
- Baratchi et al., Nature Communications 2024 (doi:10.1038/s41467-024-51211-1)
- Yang et al. (ENCHANTED2/MT), The Lancet 2022 (doi:10.1016/S0140-6736(22)01882-7)

Part of the Q-DIAPEDESIS Framework (OSF: https://osf.io/htswk | DOI: 10.5281/zenodo.23599020).
"""

from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import yaml


@dataclass
class CBFQPResult:
    """Output of the real-time ZCBF-QP physiological safety filter."""
    safe_b_amplitude_mT: float
    safe_drive_freq_Hz: float
    safe_recanalization_rate_inv_s: float
    h_shear_Pa: float
    h_thermal_C: float
    h_cpp_mmHg_s: float
    slack_delta: float
    intervention_active: bool
    stage_index: int
    target_cpp_fraction: float


class NeurocriticalCBFSafetyFilter:
    """
    Enforces forward invariance of the neurocritical-care physiological safety set C = {x : h_i(x) >= 0}
    via a minimally invasive Quadratic Program (ZCBF-QP) with a relaxation slack variable delta
    to guarantee feasibility when multiple constraints are simultaneously active (Ames et al. 2017, Remark 12):

        min_{u, delta}  0.5 * ||u - u_nom||^2 + 0.5 * p_relax * delta^2
        subject to      L_f h_i(x) + L_g h_i(x) u + kappa_i * h_i(x) + delta >= 0

    Protected Barrier Functions (Status B Falsifiable Design Targets):
    1. h_shear(x)   = tau_crit (5.50 Pa) - tau_wall_total(x, u) >= 0
    2. h_thermal(x) = DeltaT_max (0.30 C) - DeltaT_brain(x) >= 0
    3. h_cpp(x)     = dCPP_max (2.50 mmHg/s) - dCPP/dt(x, u) >= 0
    """

    def __init__(self, config_path: Path | str | None = None):
        if config_path is None:
            root_dir = Path(__file__).resolve().parent.parent.parent
            config_path = root_dir / "configs" / "neurocritical_cbf_config.yaml"

        with open(config_path, "r", encoding="utf-8") as f:
            self.cfg = yaml.safe_load(f)

        thresh = self.cfg["cbf_safety_thresholds"]
        self.tau_max_Pa = float(thresh["max_wall_shear_stress_Pa"])
        self.glyco_min_ratio = float(thresh["min_glycocalyx_integrity_ratio"])
        self.temp_max_C = float(thresh["max_brain_temp_elevation_C"])

        staged = self.cfg["staged_cpp_recanalization"]
        self.enable_staged = bool(staged["enable_staged_stepping"])
        self.num_stages = int(staged["num_recanalization_stages"])
        self.max_cpp_surge_mmHg_s = float(staged["max_pressure_surge_mmHg_per_sec"])
        self.min_safe_pH = float(staged["min_safe_efflux_pH"])

        qp = self.cfg["cbf_qp_solver"]
        self.kappa_shear = float(qp["barrier_class_kappa_shear"])
        self.kappa_thermal = float(qp["barrier_class_kappa_thermal"])
        self.kappa_cpp = float(qp["barrier_class_kappa_cpp"])
        self.relaxation_penalty = float(qp["control_relaxation_penalty"])

        # Coupling gains mapping control actions to shear, thermal power, and pressure surge
        self.shear_per_mT_Hz = 0.0042       # [Pa / (mT * Hz)]
        self.Surge_mmHg_per_recanal_rate = 18.0

    def compute_staged_cpp_target(
        self,
        recanalized_fraction: float,
        distal_pH: float,
    ) -> tuple[int, float, float]:
        """
        Determines the current recanalization stage (1..4) and maximum allowable
        recanalization stepping rate to prevent post-occlusive acidotic washout and CPP surge.
        """
        stage_idx = min(
            self.num_stages,
            int(np.floor(recanalized_fraction * self.num_stages)) + 1,
        )
        target_cpp_fraction = stage_idx / float(self.num_stages)

        # Gate recanalization rate if distal pH has not yet been buffered above min_safe_pH (7.15)
        if self.enable_staged and distal_pH < self.min_safe_pH and recanalized_fraction > 0.25:
            pH_gate = np.clip((distal_pH - 6.65) / (self.min_safe_pH - 6.65), 0.15, 1.0)
        else:
            pH_gate = 1.0

        max_safe_rate = (self.max_cpp_surge_mmHg_s / self.Surge_mmHg_per_recanal_rate) * float(pH_gate)
        return stage_idx, float(target_cpp_fraction), float(max_safe_rate)

    def filter_control(
        self,
        nominal_b_mT: float,
        nominal_freq_Hz: float,
        nominal_recanal_rate_inv_s: float,
        hydrodynamic_wall_shear_Pa: float,
        current_delta_T_C: float,
        current_glycocalyx_ratio: float,
        recanalized_fraction: float,
        distal_pH: float,
    ) -> CBFQPResult:
        """
        Solves the real-time ZCBF-QP to filter nominal swarm actuation commands into safe actions.
        """
        stage_idx, target_cpp_frac, max_recanal_rate = self.compute_staged_cpp_target(
            recanalized_fraction=recanalized_fraction,
            distal_pH=distal_pH,
        )

        # Tighten effective shear ceiling if glycocalyx approaches the 92% preservation floor
        glyco_margin = max(0.0, current_glycocalyx_ratio - self.glyco_min_ratio)
        glyco_tightening = np.clip(glyco_margin / 0.08, 0.75, 1.0)
        effective_tau_max = self.tau_max_Pa * float(glyco_tightening)

        # Decision vector z = [scale_actuation, recanal_rate, slack_delta]
        # where scale_actuation in [0.05, 1.0] scales (B_mT * freq_Hz)
        u_nom = np.array([1.0, float(nominal_recanal_rate_inv_s), 0.0], dtype=float)
        nom_mag_product = float(nominal_b_mT * nominal_freq_Hz)
        nom_robot_shear = self.shear_per_mT_Hz * nom_mag_product

        # Evaluate un-filtered barrier values
        h_thermal_0 = self.temp_max_C - float(current_delta_T_C)

        def objective(z: np.ndarray) -> float:
            du_scale = z[0] - 1.0
            du_recanal = (z[1] - u_nom[1]) / max(max_recanal_rate, 1e-3)
            slack = z[2]
            return 0.5 * (du_scale ** 2 + du_recanal ** 2) + 0.5 * self.relaxation_penalty * (slack ** 2)

        def cbf_shear_constraint(z: np.ndarray) -> float:
            # h_shear = effective_tau_max - (tau_hydro + z[0] * tau_robot_nom) + slack >= 0
            total_shear = float(hydrodynamic_wall_shear_Pa) + z[0] * nom_robot_shear
            h_s = effective_tau_max - total_shear
            return self.kappa_shear * h_s + z[2]

        def cbf_thermal_constraint(z: np.ndarray) -> float:
            # L_f h_thermal + L_g h_thermal * z[0] + kappa_thermal * h_thermal + slack >= 0
            heating_rate = 1.8e-6 * float(nominal_freq_Hz) * ((z[0] * float(nominal_b_mT)) ** 2)
            cooling_rate = 0.08 * float(current_delta_T_C)
            dh_dt = cooling_rate - heating_rate
            return dh_dt + self.kappa_thermal * h_thermal_0 + z[2]

        def cbf_cpp_constraint(z: np.ndarray) -> float:
            # h_cpp = max_cpp_surge - Surge_gain * z[1] >= 0
            cpp_surge = self.Surge_mmHg_per_recanal_rate * z[1]
            h_c = self.max_cpp_surge_mmHg_s - cpp_surge
            still_acidic_penalty = max(0.0, self.min_safe_pH - distal_pH) * 2.0
            return self.kappa_cpp * (max_recanal_rate - z[1]) - still_acidic_penalty * z[1] + z[2]

        bounds = [
            (0.05, 1.0),
            (0.0, max(max_recanal_rate, float(nominal_recanal_rate_inv_s))),
            (0.0, 2.0),
        ]
        constraints = [
            {"type": "ineq", "fun": cbf_shear_constraint},
            {"type": "ineq", "fun": cbf_thermal_constraint},
            {"type": "ineq", "fun": cbf_cpp_constraint},
        ]

        res = minimize(
            objective,
            x0=np.array([min(1.0, max(0.1, (effective_tau_max - hydrodynamic_wall_shear_Pa) / max(nom_robot_shear, 1e-6))),
                         min(float(nominal_recanal_rate_inv_s), max_recanal_rate),
                         0.0]),
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"ftol": 1e-7, "maxiter": 50},
        )

        if res.success:
            scale_opt, recanal_opt, slack_opt = float(res.x[0]), float(res.x[1]), float(res.x[2])
        else:
            # Closed-form analytical safety projection fallback
            avail_shear = max(0.0, effective_tau_max - float(hydrodynamic_wall_shear_Pa))
            scale_opt = float(np.clip(avail_shear / max(nom_robot_shear, 1e-6), 0.05, 1.0))
            if h_thermal_0 <= 0.02:
                scale_opt = min(scale_opt, 0.25)
            recanal_opt = float(min(nominal_recanal_rate_inv_s, max_recanal_rate))
            slack_opt = 0.0

        # Preserve drive frequency (80 Hz / 5 Hz Floquet resonance) while scaling B-field amplitude
        safe_b_mT = float(nominal_b_mT) * scale_opt
        safe_freq_Hz = float(nominal_freq_Hz)

        final_shear_Pa = float(hydrodynamic_wall_shear_Pa) + scale_opt * nom_robot_shear
        h_shear_final = effective_tau_max - final_shear_Pa
        h_cpp_final = self.max_cpp_surge_mmHg_s - self.Surge_mmHg_per_recanal_rate * recanal_opt

        intervened = bool(
            (abs(scale_opt - 1.0) > 1e-3)
            or (abs(recanal_opt - float(nominal_recanal_rate_inv_s)) > 1e-4)
        )

        return CBFQPResult(
            safe_b_amplitude_mT=safe_b_mT,
            safe_drive_freq_Hz=safe_freq_Hz,
            safe_recanalization_rate_inv_s=recanal_opt,
            h_shear_Pa=float(h_shear_final),
            h_thermal_C=float(h_thermal_0),
            h_cpp_mmHg_s=float(h_cpp_final),
            slack_delta=slack_opt,
            intervention_active=intervened,
            stage_index=stage_idx,
            target_cpp_fraction=target_cpp_frac,
        )
