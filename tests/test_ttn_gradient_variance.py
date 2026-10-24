"""
Tests for D2-b / D3 (ELN #003): the barren-plateau claim must rest on gradients that were
actually computed, not on a closed-form expression.

Before this change `compute_ttn_gradient_variance` returned
`Var[sin(theta_0) * cos(theta_last) / sqrt(N)]` over 64 random angle pairs for the
equivariant branch, and the literal constant `0.25 * 2**-min(N, 30)` for the
unstructured-global baseline. Neither is a gradient of any cost function, no circuit was
built, and the baseline was never evaluated at all. The claimed
"orders of magnitude improvement" was therefore a ratio of two hand-written numbers.

These tests require a real differentiable model: a parameterized circuit on the symmetric
(Dicke) subspace with a local cost observable, exact analytic gradients, and a global-cost
baseline that is genuinely evaluated at the sizes where it can be.
"""

import numpy as np
import pytest
from src.quantum_controller.ttn_gradient_variance import (
    equivariant_local_cost_gradient_variance,
    global_cost_gradient_variance,
    gradient_variance_scaling,
)


def test_equivariant_gradients_match_finite_differences():
    # The analytic gradient must agree with a central finite difference of the same cost.
    res = equivariant_local_cost_gradient_variance(n_agents=16, n_samples=4, seed=5, return_check=True)
    assert res["max_finite_difference_error"] < 1e-6, res["max_finite_difference_error"]


def test_equivariant_variance_is_a_variance_of_many_distinct_gradients():
    res = equivariant_local_cost_gradient_variance(n_agents=32, n_samples=128, seed=1, return_check=True)
    grads = np.asarray(res["gradient_samples"])
    assert grads.size >= 128
    assert np.unique(np.round(grads, 12)).size > 0.9 * grads.size, "gradients must not be degenerate"
    assert res["gradient_variance"] > 0.0


def test_global_cost_baseline_decays_exponentially_in_measured_sizes():
    # Measured, not asserted: fit log2(Var) against n over sizes we can actually simulate.
    sizes = [4, 6, 8, 10]
    variances = [global_cost_gradient_variance(n_qubits=n, n_samples=200, seed=3) for n in sizes]
    assert all(v > 0.0 for v in variances)
    slope = np.polyfit(sizes, np.log2(variances), 1)[0]
    assert slope < -0.5, f"global cost should decay at least ~2^-n/2; got slope {slope}"


def test_equivariant_variance_decays_far_more_slowly_than_global():
    sizes = [4, 6, 8, 10]
    eq = [equivariant_local_cost_gradient_variance(n_agents=n, n_samples=200, seed=3) for n in sizes]
    gl = [global_cost_gradient_variance(n_qubits=n, n_samples=200, seed=3) for n in sizes]
    eq_slope = np.polyfit(np.log2(sizes), np.log2(eq), 1)[0]      # polynomial -> finite log-log slope
    gl_slope = np.polyfit(sizes, np.log2(gl), 1)[0]               # exponential -> slope in n
    assert eq_slope > -4.0, f"equivariant branch should be polynomial in N; log-log slope {eq_slope}"
    assert gl_slope < -0.5
    assert eq[-1] > gl[-1]


def test_scaling_report_separates_measurement_from_extrapolation():
    report = gradient_variance_scaling(n_agents=200)
    assert report["equivariant_ttn_local_grad_var"] > 0.0
    # The N = 200 global-cost figure cannot be simulated directly; it must be labelled.
    assert report["unstructured_global_grad_var_is_extrapolated"] is True
    assert "global_cost_fit_sizes" in report
    assert report["equivariant_ttn_local_grad_var"] > report["unstructured_global_grad_var"]


@pytest.mark.parametrize("n_agents", [8, 64, 200])
def test_equivariant_variance_is_seed_stable(n_agents):
    a = equivariant_local_cost_gradient_variance(n_agents=n_agents, n_samples=256, seed=1)
    b = equivariant_local_cost_gradient_variance(n_agents=n_agents, n_samples=256, seed=2)
    assert a > 0.0 and b > 0.0
    assert abs(np.log10(a) - np.log10(b)) < 0.5, "variance estimate must not swing with the seed"
