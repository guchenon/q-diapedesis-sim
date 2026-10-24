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


# D9 (ELN #004): step17 is SPLIT. It previously asserted two independent claims in one test -
# the dissipative attractor overlap and the absence of barren plateaus - so a single xfail hid
# whichever one happened to be passing. They are now judged separately.


def test_step17a_equivariant_ttn_is_free_of_barren_plateaus():
    """
    D9: the criterion is the measured POLYNOMIAL DECAY EXPONENT and the margin over the
    unstructured global-cost baseline, which is what OSF Hypothesis 4 registers ("without
    barren-plateau gradient vanishing"). The old `Var >= 1e-3` floor was not preregistered - it
    appears only in test code, and the test file is not in the v0.1.0-alpha tag - and it tested
    the wrong property: any polynomially decaying quantity falls below any fixed constant once
    N is large enough, here from about N = 24. The xfail is removed; this now passes on the
    criterion the cited theory actually bounds.
    """
    ttn = LindbladEquivariantTTN()
    assert ttn.barren_plateau_free is True

    bp = ttn.benchmark_barren_plateau_scaling()
    assert 0.0 < bp["equivariant_polynomial_decay_exponent_p"] < ttn.max_polynomial_decay_exponent
    assert bp["equivariant_log_log_fit_max_residual_log2"] <= ttn.max_log_log_fit_residual
    assert bp["equivariant_ttn_local_grad_var"] > bp["unstructured_global_grad_var"]
    assert bp["variance_improvement_orders_of_magnitude"] >= 10.0
    # The N = 200 baseline is an extrapolation of the measured exponential and must stay labelled.
    assert bp["unstructured_global_grad_var_is_extrapolated"] is True


@pytest.mark.xfail(
    strict=True,
    reason=(
        "D8 (ELN #004): with the attractor blend removed (D3) and the dissipation rate DERIVED "
        "from the particle physics rather than chosen (D8), `dissipative_attractor_overlap` is "
        "0.663 at steady state against the > 0.80 in this test. The derived physical ratio is "
        "gamma/||H_eff|| = 0.747, up from 0.0167 for the old free value 0.15, and that is as far "
        "as the physics goes: reaching 0.80 needs a ratio near 3.3, about 4.4x more dissipation "
        "than a 0.75 um rotor in plasma provides. The 0.80 figure is not preregistered - it is "
        "in test code only, and the test file is not in the v0.1.0-alpha tag - so it is a "
        "candidate for correction, but it has not been changed here."
    ),
)
def test_step17b_lindblad_dissipative_attractor_overlap():
    floquet = BimodalFloquetController()
    h_eff = floquet.compute_floquet_magnus_hamiltonian(b_alpha_mT=8.5, b_beta_mT=6.0)
    ttn = LindbladEquivariantTTN()
    for _ in range(25):
        state = ttn.step_lindblad_density_matrix(dt_s=0.05, h_floquet_2x2=h_eff, local_viscosity_mPa_s=3.8)
    assert state.dissipative_attractor_overlap > 0.80


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
        # D11: the harmonic-aware Kalman filter must also clear the ceiling, so the registered
        # endpoint does not rest on the non-causal batch form of the estimator.
        assert rep.harmonic_kalman_p95_error_um < 8.0


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
        "D7 (ELN #004): `all_hypotheses_passed` fails on the transient-washout ceiling alone. "
        "Washout is an affine function of the Kuramoto order parameter, and with the coupling "
        "derived from physics (negative - see the Stage 2 xfail) r collapses to ~0.07, so mean "
        "washout is ~3.4% against the preregistered < 2.0%. Every other gated component still "
        "passes (peak shear 3.58 Pa, glycocalyx 0.989, fidelity 0.942 over the 10 s horizon, "
        "retrieval 0.99, CPP surge 2.02). Same root cause as the Stage 2 xfail - one physical "
        "finding, not several - and no threshold was changed. Note that since D10 the fidelity "
        "margin is horizon-dependent: 0.9064 over 4 s at d = 15 um, which is below 0.92."
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
