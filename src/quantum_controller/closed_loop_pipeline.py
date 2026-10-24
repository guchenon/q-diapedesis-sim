"""
Stage 3 - Step 20: Closed-Loop 4-Layer Q-DIAPEDESIS Simulation Pipeline.
Couples:
  - Layer 1: Cerebral Penetrating Arteriole Hemodynamics (Pries-Secomb Viscosity & CPP)
  - Layer 2: Leukocyte Selectin Catch-Bond Margination & Sakaguchi-Kuramoto Soliton Swarm
  - Layer 3: Bimodal Floquet-Lindblad TTN Controller & PINO Sub-Voxel Observer
  - Layer 4: Neurocritical-Care Zeroing Control Barrier Function (ZCBF-QP) Safety Filter
  - Phase 4: Exponential Degaussing & Retrograde Wall-Rolling Retrieval

Part of the Q-DIAPEDESIS Framework (OSF: https://osf.io/htswk | DOI: 10.5281/zenodo.23599020).
"""

from dataclasses import dataclass
from typing import Dict, List
import numpy as np

from src.env_hemodynamics import PenetratingArterioleGeometry, PriesSecombPulsatileFlow
from src.biomimicry_mechanics import LeukocyteCatchBondMargination, HydroMagneticSolitonSwarm
from src.neurocritical_cbf import (
    EndothelialGlycocalyxModel,
    GlycocalyxNetosisDynamics,
    NeurocriticalCBFSafetyFilter,
)
from src.observer_retrieval.pino_observer_retrieval import (
    PINOSubVoxelObserver,
    DegaussingRetrogradeRetrieval,
)
from .floquet_hamiltonian import BimodalFloquetController
from .lindblad_ttn import LindbladEquivariantTTN


@dataclass
class ClosedLoopSummary:
    """End-to-end falsifiable hypothesis evaluation metrics for a complete simulation run."""
    controller_mode: str
    final_recanalized_fraction: float
    dual_task_floquet_fidelity: float
    transient_washout_percent: float
    glycocalyx_retained_ratio: float
    peak_wall_shear_Pa: float
    peak_brain_temp_rise_C: float
    final_distal_pH: float
    netosis_activation_index: float
    pino_p95_tracking_error_um: float
    post_task_retrieval_fraction: float
    all_hypotheses_passed: bool
    peak_cpp_surge_mmHg_s: float = 0.0
    mean_nominal_b_mT: float = 0.0
    mean_applied_b_mT: float = 0.0
    elapsed_time_s: float = 0.0


class QDiapedesisClosedLoopPipeline:
    """Executes the integrated 4-layer closed-loop cerebral microvascular recanalization protocol."""

    def __init__(self, vessel_diameter_um: float | None = None):
        self.geom = PenetratingArterioleGeometry(n_axial_nodes=120)
        if vessel_diameter_um is not None:
            self.geom.diameter_um = float(vessel_diameter_um)
            self.geom.radius_um = 0.5 * float(vessel_diameter_um)
            self.geom.state = self.geom._initialize_grid()

        self.flow = PriesSecombPulsatileFlow(geometry=self.geom)
        self.margination = LeukocyteCatchBondMargination()
        self.swarm = HydroMagneticSolitonSwarm()
        self.floquet = BimodalFloquetController()
        self.lindblad_ttn = LindbladEquivariantTTN()
        self.glyco = EndothelialGlycocalyxModel(n_axial_nodes=120)
        self.physio = GlycocalyxNetosisDynamics()
        self.cbf = NeurocriticalCBFSafetyFilter()
        self.observer = PINOSubVoxelObserver()
        self.retrieval = DegaussingRetrogradeRetrieval()

    def run_episode(
        self,
        total_time_s: float | None = None,
        dt_s: float = 0.05,
        enable_cbf_filter: bool = True,
        enable_bimodal_floquet: bool = True,
        enable_biomimetic_catch_bond: bool = True,
        enable_degaussing_retrieval: bool = True,
        enable_pulsed_pH_buffering: bool = True,
        stop_at_recanalized_fraction: float | None = None,
    ) -> ClosedLoopSummary:
        """
        Runs a recanalization + retrieval episode and returns all preregistered metrics.

        Ablation fairness (D4, ELN #003): every configuration receives the same nominal command;
        each `enable_*` switch toggles exactly one module. `total_time_s` defaults to
        `simulation_meta.total_time_seconds`. `stop_at_recanalized_fraction` ends the episode at an
        iso-recanalization checkpoint for equal-progress comparisons.
        """
        if total_time_s is None:
            total_time_s = float(self.geom.cfg["simulation_meta"]["total_time_seconds"])
        self.geom.state = self.geom._initialize_grid()
        self.glyco.reset()
        self.physio.reset()
        self.lindblad_ttn.reset()
        swarm_state = self.swarm.reset()

        r_coords = swarm_state.radial_positions_um.copy()
        recanalized_frac = 0.0
        n_steps = int(round(total_time_s / dt_s))

        peak_shear_history: List[float] = []
        peak_temp_history: List[float] = []
        fidelity_history: List[float] = []
        washout_history: List[float] = []
        cpp_surge_history: List[float] = []
        nominal_b_history: List[float] = []
        applied_b_history: List[float] = []
        elapsed_s = 0.0

        for step in range(n_steps):
            t_s = step * dt_s
            hemo_0 = self.flow.solve_flow_step(time_s=t_s, robot_induced_shear_Pa=0.0)
            hydro_peak_shear = float(np.max(hemo_0.wall_shear_stress_Pa))

            # Layer 3: Evaluate Bimodal Floquet Hamiltonian & Lindblad TTN state
            f_beta = 5.0 if enable_bimodal_floquet else 80.0
            floquet_nom = self.floquet.evaluate_control_step(
                time_s=t_s, b_scale_factor=1.0, f_alpha_override_Hz=80.0, f_beta_override_Hz=f_beta
            )
            ttn_state = self.lindblad_ttn.step_lindblad_density_matrix(
                dt_s=dt_s,
                h_floquet_2x2=floquet_nom.effective_hamiltonian_2x2,
                local_viscosity_mPa_s=float(np.mean(hemo_0.effective_viscosity_mPa_s)),
            )

            nom_b_mT = floquet_nom.b_alpha_mT * ttn_state.optimal_field_correction_scale
            nom_freq_Hz = floquet_nom.omega_alpha_Hz
            nom_recanal_rate = 0.32

            # Layer 4: Apply Neurocritical ZCBF-QP Safety Filter
            if enable_cbf_filter:
                cbf_out = self.cbf.filter_control(
                    nominal_b_mT=nom_b_mT,
                    nominal_freq_Hz=nom_freq_Hz,
                    nominal_recanal_rate_inv_s=nom_recanal_rate,
                    hydrodynamic_wall_shear_Pa=hydro_peak_shear,
                    current_delta_T_C=self.physio.delta_T_C,
                    current_glycocalyx_ratio=self.glyco.get_state().mean_integrity_ratio,
                    recanalized_fraction=recanalized_frac,
                    distal_pH=self.physio.current_pH,
                )
                act_b_mT = cbf_out.safe_b_amplitude_mT
                act_recanal_rate = cbf_out.safe_recanalization_rate_inv_s
            else:
                # CBF ablation: the same nominal command is applied unfiltered (no multipliers)
                act_b_mT = nom_b_mT
                act_recanal_rate = nom_recanal_rate

            robot_shear_Pa = self.cbf.shear_per_mT_Hz * act_b_mT * nom_freq_Hz
            hemo_step = self.flow.solve_flow_step(time_s=t_s, robot_induced_shear_Pa=robot_shear_Pa)

            # Layer 2: Step Biomimetic Margination, Catch-Bond Rolling, and Soliton Synchronization
            self.margination.bond_mode = "catch_slip" if enable_biomimetic_catch_bond else "pure_slip_baseline"
            if enable_biomimetic_catch_bond:
                r_coords = self.margination.step_radial_margination(
                    radial_positions_um=r_coords,
                    lumen_radius_um=self.geom.radius_um,
                    wall_shear_stress_Pa=float(np.mean(hemo_step.wall_shear_stress_Pa)),
                    dt_s=dt_s,
                )
            roll_state = self.margination.evaluate_rolling_step(
                radial_positions_um=r_coords,
                lumen_radius_um=self.geom.radius_um,
                wall_shear_stress_Pa=float(np.mean(hemo_step.wall_shear_stress_Pa)),
                is_systole=hemo_step.is_systole,
                magnetic_drive_freq_Hz=nom_freq_Hz,
            )
            swarm_state = self.swarm.step_soliton_swarm(
                dt_s=dt_s,
                axial_rolling_velocity_um_s=roll_state.axial_rolling_velocity_um_s,
                wall_proximal_fraction=roll_state.wall_proximal_fraction if enable_biomimetic_catch_bond else 0.25,
                adhesion_locked=roll_state.adhesion_locked if enable_biomimetic_catch_bond else False,
                collateral_skimming_fraction=hemo_step.collateral_skimming_fraction,
                enable_soliton_coupling=enable_biomimetic_catch_bond,
            )

            # Update vessel recanalization and physiological state
            recanalized_frac = min(1.0, recanalized_frac + act_recanal_rate * dt_s)
            buffer_mix = min(1.0, (self.physio.current_pH - 6.65) / (7.35 - 6.65))
            self.geom.apply_staged_recanalization(recanalized_fraction=recanalized_frac, buffer_mix_ratio=buffer_mix)

            glyco_state = self.glyco.step_spatial_glycocalyx(
                dt_s=dt_s,
                wall_shear_profile_Pa=hemo_step.wall_shear_stress_Pa,
                pH_profile=self.geom.state.pH_profile,
            )
            physio_state = self.physio.step_physiology(
                dt_s=dt_s,
                wall_shear_profile_Pa=hemo_step.wall_shear_stress_Pa,
                b_field_amplitude_mT=act_b_mT,
                drive_frequency_Hz=nom_freq_Hz,
                blood_perfusion_velocity_um_s=float(np.mean(hemo_step.axial_velocity_mean_um_s)),
                pulsed_buffering_active=enable_pulsed_pH_buffering,
            )

            peak_shear_history.append(float(np.max(hemo_step.wall_shear_stress_Pa)))
            peak_temp_history.append(physio_state.brain_temp_elevation_C)
            fidelity_history.append(floquet_nom.dual_task_fidelity)
            washout_history.append(swarm_state.transient_washout_percent)
            cpp_surge_history.append(self.cbf.Surge_mmHg_per_recanal_rate * act_recanal_rate)
            nominal_b_history.append(float(nom_b_mT))
            applied_b_history.append(float(act_b_mT))
            elapsed_s = (step + 1) * dt_s

            if stop_at_recanalized_fraction is not None and recanalized_frac >= float(stop_at_recanalized_fraction):
                break

        # Evaluate Module B2 (PINO observer) and Module B5 (Degaussing + Retrograde Retrieval)
        obs_report = self.observer.evaluate_blind_window_tracking()
        mean_washout = float(np.mean(washout_history))
        ret_report = self.retrieval.execute_post_task_retrieval(
            transient_washout_percent=mean_washout,
            enable_degaussing_override=enable_degaussing_retrieval,
        )

        mean_fidelity = float(np.mean(fidelity_history))
        max_shear = float(np.max(peak_shear_history))
        max_temp = float(np.max(peak_temp_history))
        max_cpp_surge = float(np.max(cpp_surge_history))

        all_passed = bool(
            (mean_fidelity >= 0.92)
            and (mean_washout < 2.0)
            and (glyco_state.mean_integrity_ratio >= 0.92)
            and (max_shear <= 5.50 + 1e-3)
            and (max_temp < 0.30)
            and (obs_report.pino_p95_error_um < 8.0)
            and (ret_report.final_recovered_fraction >= 0.99)
            and (max_cpp_surge <= self.cbf.max_cpp_surge_mmHg_s + 1e-9)
        )

        switches = {
            "cbf": enable_cbf_filter,
            "floquet": enable_bimodal_floquet,
            "catch_bond": enable_biomimetic_catch_bond,
            "degaussing": enable_degaussing_retrieval,
            "pH_buffering": enable_pulsed_pH_buffering,
        }
        disabled = [name for name, on in switches.items() if not on]
        mode_label = "Q-DIAPEDESIS_Full" if not disabled else "Ablated_no_" + "+".join(disabled)

        return ClosedLoopSummary(
            controller_mode=mode_label,
            final_recanalized_fraction=float(recanalized_frac),
            dual_task_floquet_fidelity=mean_fidelity,
            transient_washout_percent=mean_washout,
            glycocalyx_retained_ratio=glyco_state.mean_integrity_ratio,
            peak_wall_shear_Pa=max_shear,
            peak_brain_temp_rise_C=max_temp,
            final_distal_pH=physio_state.distal_efflux_pH,
            netosis_activation_index=physio_state.netosis_activation_index,
            pino_p95_tracking_error_um=obs_report.pino_p95_error_um,
            post_task_retrieval_fraction=ret_report.final_recovered_fraction,
            all_hypotheses_passed=all_passed,
            peak_cpp_surge_mmHg_s=max_cpp_surge,
            mean_nominal_b_mT=float(np.mean(nominal_b_history)),
            mean_applied_b_mT=float(np.mean(applied_b_history)),
            elapsed_time_s=float(elapsed_s),
        )
