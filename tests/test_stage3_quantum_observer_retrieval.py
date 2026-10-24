"""
End-to-End Verification & Baseline Comparison Tests for Stage 3 (Steps 16-21):
- Step 16 (Module B1): Bimodal Floquet Decoupling (80 Hz vs 5 Hz, Fidelity >= 92%)
- Step 17 (Module B4): Lindblad CPTP Attractor & Equivariant TTN (Absence of Barren Plateaus)
- Step 18 (Module B2): PINO Sub-Voxel Tracking (< 8.0 um error across 50-200 ms MRI blind window)
- Step 19 (Module B5): Exponential Degaussing & Retrograde Retrieval (>= 99.0% mass balance)
- Step 20 (Closed Loop): Integrated 4-Layer Pipeline + Diameter Sensitivity Sweep (15, 30, 50 um)
"""

import pytest

from src.quantum_controller import (
    BimodalFloquetController,
    LindbladEquivariantTTN,
    QDiapedesisClosedLoopPipeline,
)
from src.observer_retrieval import (
    PINOSubVoxelObserver,
    DegaussingRetrogradeRetrieval,
)


def test_step16_bimodal_floquet_and_baselines():
    floquet = BimodalFloquetController()
    state = floquet.evaluate_control_step(time_s=0.5, b_scale_factor=1.0)

    assert state.omega_alpha_Hz == 80.0
    assert state.omega_beta_Hz == 5.0
    assert state.dual_task_fidelity >= 0.92
    assert state.passes_fidelity_target is True

    # Verify superiority over single-frequency and step-out baselines (Audit Row 6b)
    benchmarks = floquet.compare_with_baselines()
    assert benchmarks["bimodal_floquet_fidelity"] > benchmarks["single_80Hz_baseline_fidelity"]
    assert benchmarks["bimodal_floquet_fidelity"] > benchmarks["single_5Hz_baseline_fidelity"]
    assert benchmarks["bimodal_floquet_fidelity"] > benchmarks["overdriven_140Hz_stepout_fidelity"]


@pytest.mark.xfail(
    strict=True,
    reason=(
        "C2 (ELN #002): gradient variance 9.5e-4 < 1e-3 threshold and 6.6 < 10 orders; "
        "kept as a reported failure until D2-b replaces the closed-form TTN variance "
        "with a real computation. Thresholds unchanged."
    ),
)
def test_step17_lindblad_dissipation_and_equivariant_ttn_barren_plateau():
    floquet = BimodalFloquetController()
    h_eff = floquet.compute_floquet_magnus_hamiltonian(b_alpha_mT=8.5, b_beta_mT=6.0)
    ttn = LindbladEquivariantTTN()

    for _ in range(25):
        state = ttn.step_lindblad_density_matrix(dt_s=0.05, h_floquet_2x2=h_eff, local_viscosity_mPa_s=3.8)

    assert state.dissipative_attractor_overlap > 0.80
    assert state.barren_plateau_free_verified is True

    bp_bench = ttn.benchmark_barren_plateau_scaling()
    assert bp_bench["equivariant_ttn_local_grad_var"] > bp_bench["unstructured_global_grad_var"]
    assert bp_bench["variance_improvement_orders_of_magnitude"] >= 10.0


@pytest.mark.xfail(
    strict=True,
    reason=(
        "C2 (ELN #002): PINO p95 loses to dead reckoning at 50 ms (1.21 vs 0.69 um) "
        "and 100 ms (1.65 vs 1.47 um); reported as a negative result until D2-b. "
        "Thresholds unchanged."
    ),
)
def test_step18_pino_subvoxel_observer_across_blind_windows():
    observer = PINOSubVoxelObserver()
    for win_ms in (50.0, 100.0, 200.0):
        rep = observer.evaluate_blind_window_tracking(blind_window_ms=win_ms)
        assert rep.pino_p95_error_um < 8.0
        assert rep.passes_subvoxel_target is True
        # PINO must outperform Kalman, Dead Reckoning, and Last-Hold baselines (Audit Row 8)
        assert rep.pino_p95_error_um < rep.kalman_p95_error_um
        assert rep.pino_p95_error_um < rep.dead_reckoning_p95_error_um
        assert rep.pino_p95_error_um < rep.last_hold_p95_error_um


def test_step19_degaussing_retrograde_retrieval_mass_balance():
    retrieval = DegaussingRetrogradeRetrieval()
    rep_with_degauss = retrieval.execute_post_task_retrieval(
        transient_washout_percent=1.10, enable_degaussing_override=True
    )
    rep_no_degauss = retrieval.execute_post_task_retrieval(
        transient_washout_percent=1.10, enable_degaussing_override=False
    )

    assert rep_with_degauss.final_recovered_fraction >= 0.99
    assert rep_with_degauss.passes_retrieval_target is True
    assert rep_with_degauss.final_recovered_fraction > rep_no_degauss.final_recovered_fraction


def test_step20_closed_loop_pipeline_and_diameter_sweep():
    # Audit Table Row 1a/1b: verify closed-loop pipeline at 15, 30, and 50 um diameters
    for d_um in (15.0, 30.0, 50.0):
        pipeline = QDiapedesisClosedLoopPipeline(vessel_diameter_um=d_um)
        # Full configured horizon (simulation_meta.total_time_seconds = 10 s), D4 / ELN #003
        summary = pipeline.run_episode(dt_s=0.05)
        assert summary.all_hypotheses_passed is True
        assert summary.peak_wall_shear_Pa <= 5.50 + 1e-3
        assert summary.glycocalyx_retained_ratio >= 0.92
