"""
End-to-End Integration Tests for Stage 1 (Steps 06-10):
Coupled Cerebral Penetrating Arteriole Hemodynamics (Layer 1) and
Neurocritical-Care Control Barrier Function Safety Filter (Layer 4).
"""

import numpy as np
from src.env_hemodynamics import PenetratingArterioleGeometry, PriesSecombPulsatileFlow
from src.neurocritical_cbf import (
    EndothelialGlycocalyxModel,
    GlycocalyxNetosisDynamics,
    NeurocriticalCBFSafetyFilter,
)


def test_stage1_coupled_hemodynamics_and_cbf_protection():
    geom = PenetratingArterioleGeometry(n_axial_nodes=120)
    flow_solver = PriesSecombPulsatileFlow(geometry=geom)
    glyco_model = EndothelialGlycocalyxModel(n_axial_nodes=120)
    physio_dyn = GlycocalyxNetosisDynamics()
    cbf_filter = NeurocriticalCBFSafetyFilter()

    dt_s = 0.05
    recanalized_frac = 0.0

    # Simulate 40 steps (2.0 seconds) under an overly aggressive nominal magnetic command
    # (e.g., B = 18 mT at 80 Hz -> excessive un-filtered robot shear > 6 Pa)
    nominal_b_mT = 18.0
    nominal_freq_Hz = 80.0
    nominal_recanal_rate = 0.35

    for step in range(40):
        t_s = step * dt_s
        hemo_baseline = flow_solver.solve_flow_step(time_s=t_s, robot_induced_shear_Pa=0.0)
        peak_hydro_shear_Pa = float(np.max(hemo_baseline.wall_shear_stress_Pa))

        cbf_res = cbf_filter.filter_control(
            nominal_b_mT=nominal_b_mT,
            nominal_freq_Hz=nominal_freq_Hz,
            nominal_recanal_rate_inv_s=nominal_recanal_rate,
            hydrodynamic_wall_shear_Pa=peak_hydro_shear_Pa,
            current_delta_T_C=physio_dyn.delta_T_C,
            current_glycocalyx_ratio=glyco_model.get_state().mean_integrity_ratio,
            recanalized_fraction=recanalized_frac,
            distal_pH=physio_dyn.current_pH,
        )

        # Verify CBF intervened and capped total wall shear <= 5.50 Pa
        assert cbf_res.intervention_active is True
        robot_shear_Pa = cbf_filter.shear_per_mT_Hz * cbf_res.safe_b_amplitude_mT * cbf_res.safe_drive_freq_Hz
        hemo_safe = flow_solver.solve_flow_step(time_s=t_s, robot_induced_shear_Pa=robot_shear_Pa)

        assert float(np.max(hemo_safe.wall_shear_stress_Pa)) <= 5.50 + 1e-4

        # Advance recanalization and spatial/scalar physiology under safe CBF commands
        recanalized_frac = min(1.0, recanalized_frac + cbf_res.safe_recanalization_rate_inv_s * dt_s)
        buffer_mix = min(1.0, (physio_dyn.current_pH - 6.65) / (7.35 - 6.65))
        geom.apply_staged_recanalization(recanalized_fraction=recanalized_frac, buffer_mix_ratio=buffer_mix)

        glyco_state = glyco_model.step_spatial_glycocalyx(
            dt_s=dt_s,
            wall_shear_profile_Pa=hemo_safe.wall_shear_stress_Pa,
            pH_profile=geom.state.pH_profile,
        )
        physio_state = physio_dyn.step_physiology(
            dt_s=dt_s,
            wall_shear_profile_Pa=hemo_safe.wall_shear_stress_Pa,
            b_field_amplitude_mT=cbf_res.safe_b_amplitude_mT,
            drive_frequency_Hz=cbf_res.safe_drive_freq_Hz,
            blood_perfusion_velocity_um_s=float(np.mean(hemo_safe.axial_velocity_mean_um_s)),
            pulsed_buffering_active=True,
        )

    # Verify all preregistered Stage 1 physiological safety targets hold
    assert glyco_state.passes_cbf_threshold is True
    assert glyco_state.mean_integrity_ratio >= 0.92
    assert physio_state.brain_temp_elevation_C < 0.30
    assert physio_state.netosis_activation_index < 0.05


def test_glycocalyx_baseline_sensitivity_sweep():
    glyco_model = EndothelialGlycocalyxModel(n_axial_nodes=50)
    safe_shear_profile = np.full(50, 3.8, dtype=float)
    sweep_res = glyco_model.run_baseline_sensitivity_sweep(safe_shear_profile, duration_s=1.0, dt_s=0.05)

    assert "retained_um_at_0.40um_baseline" in sweep_res
    assert "retained_um_at_0.50um_baseline" in sweep_res
    assert "retained_um_at_1.00um_baseline" in sweep_res
    assert sweep_res["integrity_ratio_at_0.50um_baseline"] >= 0.92
