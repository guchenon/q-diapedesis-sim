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
        "Two independent failures, thresholds unchanged. "
        "(1) C2 (ELN #002): gradient variance 9.5e-4 < 1e-3 and 6.6 < 10 orders; kept as a "
        "reported failure until D2-b replaces the closed-form TTN variance with a real "
        "computation. "
        "(2) D3 Stage 3 (ELN #003): with the blend toward `rho_attractor` removed and the GKLS "
        "equation integrated as written (RK4), `dissipative_attractor_overlap` is 0.550 at the "
        "25-step mark and only 0.622 at steady state, against the preregistered > 0.80. The "
        "preregistered dissipation is 60x weaker than the Floquet Hamiltonian "
        "(gamma_eff = 0.19 /s vs ||H_eff|| = 11.4 rad/s), so coherent precession dominates and "
        "the dissipative fixed point is never approached. Reaching 0.80 at steady state needs "
        "viscous_dissipation_gamma ~ 30 instead of 0.15, a 200x increase. The earlier pass came "
        "entirely from the blend, not from the dissipator."
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


def test_step18_pino_subvoxel_observer_across_blind_windows():
    # xfail removed under D2-b (ELN #003). The earlier failure - PINO losing to dead reckoning
    # at 50 and 100 ms - was measured against an "observer" that was ground truth plus noise and
    # baselines that read the ground-truth velocity, one of them with an added error handicap.
    # With real estimators that see only the initial condition, the physics prior and noisy
    # phase-contrast measurements, this now passes at all three windows on its own merits.
    # Read the findings in ELN #003 before quoting the margin: the advantage comes from the
    # cardiac-waveform prior, and nothing in this module is a neural network.
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


@pytest.mark.xfail(
    strict=True,
    reason=(
        "D3 Stage 2 (ELN #003): `all_hypotheses_passed` now fails on the transient-washout "
        "ceiling alone. Washout is an affine function of the Kuramoto order parameter, and with "
        "the honest Stage 2 model (incoherent start, frequency dispersion, rotational noise) "
        "r drops from ~0.95 to ~0.37, giving mean washout 2.50% vs the preregistered < 2.0% at "
        "d = 30 um (2.49% at 15 um, 2.83% at 50 um). Every other gated component still passes "
        "(peak shear 3.58 Pa, glycocalyx 0.987, fidelity 0.993, retrieval 0.9916, CPP surge 2.02). "
        "Same root cause as the Stage 2 xfail; no threshold was changed."
    ),
)
def test_step20_closed_loop_pipeline_and_diameter_sweep():
    # Audit Table Row 1a/1b: verify closed-loop pipeline at 15, 30, and 50 um diameters
    for d_um in (15.0, 30.0, 50.0):
        pipeline = QDiapedesisClosedLoopPipeline(vessel_diameter_um=d_um)
        # Full configured horizon (simulation_meta.total_time_seconds = 10 s), D4 / ELN #003
        summary = pipeline.run_episode(dt_s=0.05)
        assert summary.all_hypotheses_passed is True
        assert summary.peak_wall_shear_Pa <= 5.50 + 1e-3
        assert summary.glycocalyx_retained_ratio >= 0.92
