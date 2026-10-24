"""
Tests for D2-b / D3 (ELN #003): the sub-voxel observer must be an estimator.

Before this change the "PINO" estimate was literally the ground truth plus noise:

    pino_z_est = z_gt + rng.normal(0, sigma, n) * (1 + 1.5 * t)

so it was handed the answer and could never be wrong by more than the noise it added;
nothing was predicted. The "Kalman filter" was `cumsum(0.32 v_mean + 0.55 (v_true - v0) dt)`,
which also reads the ground-truth velocity `v_true`, and its error then had an explicit
handicap `+ 3.2 * (t / T)` added so that the PINO would win the comparison.

These tests require that no estimator sees the ground-truth trajectory or the ground-truth
velocity; that each one uses only the initial condition, the stated physics prior and the
noisy measurements; that no baseline carries an added penalty; and that the comparison
degrades in the way a real estimator must when the prior stops matching reality.
"""

import inspect
import numpy as np
import pytest
from src.observer_retrieval import PINOSubVoxelObserver
from src.observer_retrieval import pino_observer_retrieval as mod


def _executable_source(obj) -> str:
    """Source with comments and string literals removed, so prose cannot satisfy or trip a scan."""
    import io
    import tokenize

    out = []
    for tok in tokenize.generate_tokens(io.StringIO(inspect.getsource(obj)).readline):
        if tok.type in (tokenize.COMMENT, tokenize.STRING):
            continue
        out.append(tok.string)
    return " ".join(out)


def test_no_estimator_reads_the_ground_truth_or_carries_a_handicap():
    code = _executable_source(mod.PINOSubVoxelObserver)
    # The ground-truth velocity may be DEFINED (the simulator needs it) and MEASURED through
    # noise, but must not feed any estimator. Remove those two legitimate uses, then scan.
    scanned = code.replace("v_true_z = ", " GT_DEF ", 1).replace("v_true_z + rng", " GT_MEAS ", 1)
    for forbidden in ("v_true", "z_true", "r_true"):
        for estimator in ("z_pino", "z_kf", "z_dr", "r_pino", "v_fit", "coeffs"):
            assert f"{estimator} = {forbidden}" not in scanned, (
                f"{estimator} is assigned from the ground truth {forbidden}"
            )
    # No estimator may be built by perturbing the truth, which is what the old PINO did.
    assert "z_true + rng" not in scanned and "r_true + rng" not in scanned
    # And no baseline may carry an added error penalty.
    assert "3.2 *" not in scanned, "the hand-added Kalman error handicap is still present"


def test_estimators_only_receive_initial_condition_and_noisy_measurements():
    obs = PINOSubVoxelObserver()
    rep = obs.evaluate_blind_window_tracking(blind_window_ms=100.0, seed=7, return_traces=True)
    # The measurements the observers are allowed to see must differ from the truth they estimate.
    assert not np.allclose(rep["velocity_measurements_um_s"], rep["true_velocity_um_s"])
    for name in ("pino", "kalman", "dead_reckoning", "last_hold"):
        assert not np.allclose(rep[f"{name}_z_um"], rep["true_z_um"]), f"{name} tracks truth exactly"


def test_pino_error_grows_when_measurement_noise_grows():
    obs = PINOSubVoxelObserver()
    quiet = obs.evaluate_blind_window_tracking(blind_window_ms=100.0, seed=3, n_shots_override=4096)
    noisy = obs.evaluate_blind_window_tracking(blind_window_ms=100.0, seed=3, n_shots_override=16)
    assert noisy.pino_p95_error_um > quiet.pino_p95_error_um


def test_kalman_is_a_real_filter_and_is_never_handicapped_below_dead_reckoning():
    # The KF sees strictly more than dead reckoning - the same prior plus the measurements - so
    # on average it must not be WORSE. At the a-priori SNR the measurements are nearly
    # uninformative (sigma_v ~ 398 um/s against a ~90 um/s mean velocity), so the two end up
    # close; the point of this test is that the KF is not pushed below the baseline, which is
    # what the removed `+ 3.2 * (t/T)` handicap did.
    obs = PINOSubVoxelObserver()
    pairs = [
        (r.kalman_p95_error_um, r.dead_reckoning_p95_error_um)
        for r in (obs.evaluate_blind_window_tracking(blind_window_ms=200.0, seed=s) for s in range(12))
    ]
    mean_kf = np.mean([kf for kf, _ in pairs])
    mean_dr = np.mean([dr for _, dr in pairs])
    assert mean_kf <= mean_dr * 1.02, f"Kalman {mean_kf:.3f} is worse than dead reckoning {mean_dr:.3f}"


def test_kalman_clearly_beats_dead_reckoning_once_measurements_are_informative():
    # With plenty of shots the velocity measurements carry real information, and a correct KF
    # must then pull well clear of the open-loop baseline that ignores them.
    obs = PINOSubVoxelObserver()
    rs = [
        obs.evaluate_blind_window_tracking(blind_window_ms=200.0, seed=s, n_shots_override=32768)
        for s in range(8)
    ]
    mean_kf = np.mean([r.kalman_p95_error_um for r in rs])
    mean_dr = np.mean([r.dead_reckoning_p95_error_um for r in rs])
    assert mean_kf < 0.7 * mean_dr, f"Kalman {mean_kf:.3f} vs dead reckoning {mean_dr:.3f}"


def test_last_hold_is_the_worst_baseline_over_a_long_window():
    obs = PINOSubVoxelObserver()
    r = obs.evaluate_blind_window_tracking(blind_window_ms=200.0, seed=1)
    assert r.last_hold_p95_error_um > r.dead_reckoning_p95_error_um


@pytest.mark.parametrize("window_ms", [50.0, 100.0, 200.0])
def test_results_are_reported_for_every_preregistered_window(window_ms):
    r = PINOSubVoxelObserver().evaluate_blind_window_tracking(blind_window_ms=window_ms, seed=1)
    assert r.blind_window_ms == window_ms
    assert r.pino_p95_error_um > 0.0


def test_model_mismatch_is_simulated_and_degrades_the_physics_informed_estimator():
    # Inverse crime guard: the PINO's spectral basis is the cardiac harmonics, which also
    # generate the truth. If the truth contains a component OUTSIDE that basis the estimator
    # must get worse, otherwise the comparison is only testing that the basis is exact.
    obs = PINOSubVoxelObserver()
    matched = obs.evaluate_blind_window_tracking(blind_window_ms=200.0, seed=2, model_mismatch_um_s=0.0)
    mismatched = obs.evaluate_blind_window_tracking(blind_window_ms=200.0, seed=2, model_mismatch_um_s=40.0)
    assert mismatched.pino_p95_error_um > matched.pino_p95_error_um
