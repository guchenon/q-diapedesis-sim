"""
Tests for D11 (ELN #003): name the observer honestly, and separate the value of the physics
prior from the value of the estimator.

The class is a five-mode regularized Fourier (spectral) MAP estimator with a haemodynamic
prior. There is no neural network in it and no operator is learned, so the "PINO" name and the
neural-operator citations overstate what the code does. The original name is kept so imports do
not break; `PhysicsInformedSpectralObserver` is the accurate alias.

The HARMONIC-AWARE KALMAN FILTER baseline exists to answer one question that the previous
comparison could not: the spectral observer is handed the full cardiac waveform while plain
dead reckoning gets only the mean, so its margin could be entirely down to the prior. Giving a
Kalman filter the same harmonic model settles it.
"""

import numpy as np
import pytest
from src.observer_retrieval import PINOSubVoxelObserver, PhysicsInformedSpectralObserver


def test_accurate_alias_refers_to_the_same_class():
    assert PhysicsInformedSpectralObserver is PINOSubVoxelObserver


def test_module_states_plainly_that_it_is_not_a_neural_network():
    doc = (PINOSubVoxelObserver.__doc__ or "").lower()
    flat = " ".join(doc.split())
    assert "no neural network" in flat
    assert "physics-informed spectral map observer" in flat
    assert "no operator is learned" in flat


def test_harmonic_kalman_baseline_is_reported_and_is_a_real_filter():
    rep = PINOSubVoxelObserver().evaluate_blind_window_tracking(blind_window_ms=100.0, seed=1)
    assert rep.harmonic_kalman_p95_error_um > 0.0
    # Carrying the harmonic model must beat the constant-velocity filter that lacks it.
    assert rep.harmonic_kalman_p95_error_um < rep.kalman_p95_error_um


@pytest.mark.parametrize("window_ms", [50.0, 100.0, 200.0])
def test_most_of_the_advantage_comes_from_the_physics_prior(window_ms):
    """
    The prior must account for the bulk of the spectral observer's margin over a plain Kalman
    filter. This is the finding, not a defect: it means the result is evidence that knowing the
    cardiac waveform helps, and NOT evidence that a neural-operator-style estimator is needed.
    """
    obs = PINOSubVoxelObserver()
    reps = [obs.evaluate_blind_window_tracking(blind_window_ms=window_ms, seed=s) for s in range(12)]
    spectral = np.mean([r.pino_p95_error_um for r in reps])
    harmonic = np.mean([r.harmonic_kalman_p95_error_um for r in reps])
    plain = np.mean([r.kalman_p95_error_um for r in reps])

    total_gap = plain - spectral
    prior_share = (plain - harmonic) / total_gap
    assert total_gap > 0.0
    assert prior_share > 0.5, (
        f"the prior should explain most of the gap; it explains {prior_share:.0%} "
        f"(spectral {spectral:.3f}, harmonic KF {harmonic:.3f}, plain KF {plain:.3f})"
    )


@pytest.mark.parametrize("window_ms", [50.0, 100.0, 200.0])
def test_residual_edge_over_the_causal_filter_is_small_and_is_smoothing(window_ms):
    """
    A residual edge remains for the batch MAP fit over the sequential filter, and it must be
    reported for what it is: the batch fit uses the whole window, including samples later than
    the instant it estimates, so it is a SMOOTHER, not a causal filter. For real-time use inside
    the blind window the harmonic-aware Kalman filter is the fair comparator.
    """
    obs = PINOSubVoxelObserver()
    reps = [obs.evaluate_blind_window_tracking(blind_window_ms=window_ms, seed=s) for s in range(12)]
    spectral = np.mean([r.pino_p95_error_um for r in reps])
    harmonic = np.mean([r.harmonic_kalman_p95_error_um for r in reps])
    assert spectral < harmonic                       # the smoother is better, as expected
    assert harmonic < 1.6 * spectral                 # but not by a large factor
