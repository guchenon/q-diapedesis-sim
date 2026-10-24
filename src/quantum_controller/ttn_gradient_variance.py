"""
Stage 3 - Step 17 (Module B4), D2-b / D3: gradient-variance measurement for the
barren-plateau claim.

This module replaces a closed-form expression with an actual computation. Previously
the "equivariant TTN" variance was `Var[sin(a) cos(b) / sqrt(N)]` over random angle pairs
and the "unstructured global" baseline was the literal constant `0.25 * 2**-min(N, 30)`:
neither was the gradient of any cost function and no circuit was ever built.

What is computed here:

1. `equivariant_local_cost_gradient_variance` - an S_n-permutation-equivariant ansatz. Because
   every gate is a function of the collective spin operators, the dynamics stay inside the
   symmetric (Dicke) subspace, which for N agents has dimension N + 1 rather than 2^N. The
   circuit is `L = ceil(log2 N)` layers of `exp(-i theta J_z / j) exp(-i phi J_x / j)` with
   `j = N/2`, matching the logarithmic tree depth the cited results assume. The cost is the
   LOCAL observable `J_z / j`, which on the symmetric subspace is exactly the single-site
   Pauli-Z expectation. Gradients are exact: for `|psi> = U_L...U_1|psi_0>`,
   `dC/dtheta_k = 2 Re<d_k psi| O |psi>` with `|d_k psi> = U_{>k} (-i G_k) U_{<=k} |psi_0>`.
   The finite-difference cross-check is available via `return_check=True`.

2. `global_cost_gradient_variance` - the barren-plateau baseline: a deep unstructured
   hardware-efficient ansatz on the full 2^n Hilbert space with the GLOBAL cost
   `O = |0...0><0...0|`, differentiated by the parameter-shift rule (exact for the
   single-qubit Pauli rotations used here). This is simulated directly, not asserted.

Honest limitation, and it is the reason `gradient_variance_scaling` reports a flag:
the global-cost baseline at N = 200 cannot be simulated - 2^200 amplitudes do not exist.
It is obtained by fitting log2(Var) against n over the sizes that CAN be simulated and
extrapolating. The returned dictionary labels that number as an extrapolation and lists the
fit sizes, so no caller can mistake it for a measurement.

Grounded in:
- McClean, Boixo, Smelyanskiy, Babbush & Neven, Nat. Commun. 2018 (doi:10.1038/s41467-018-07090-4)
- Cerezo, Sone, Volkoff, Cincio & Coles, Nat. Commun. 2021 (doi:10.1038/s41467-021-21728-w)
- Schatzki, Larocca, Nguyen, Sauvage & Cerezo, npj Quantum Inf. 2024 (doi:10.1038/s41534-024-00804-1)
- Barthel & Miao, Commun. Math. Phys. 2025 (doi:10.1007/s00220-024-05217-x)
"""

from functools import lru_cache
from typing import Dict, List
import numpy as np
from scipy.linalg import expm

# These measurements are deterministic in (size, n_samples, seed) and expensive - the
# global-cost baseline simulates a full 2^n state vector - so results are memoized.
# Caching changes no number; it only stops the same measurement being repeated.

# Sizes at which the full 2^n global-cost baseline is simulated directly.
GLOBAL_COST_FIT_SIZES: List[int] = [4, 6, 8, 10, 12]


# ---------------------------------------------------------------------------------------
# 1. Equivariant ansatz on the symmetric (Dicke) subspace
# ---------------------------------------------------------------------------------------

def _collective_spin_operators(n_agents: int):
    """J_x, J_z in the spin-j irrep with j = N/2, i.e. the symmetric subspace (dim N + 1)."""
    j = 0.5 * float(n_agents)
    m = np.arange(j, -j - 1.0, -1.0)                      # m = j, j-1, ..., -j
    dim = m.size
    J_z = np.diag(m).astype(complex)
    # <j, m+1| J_+ |j, m> = sqrt(j(j+1) - m(m+1))
    off = np.sqrt(j * (j + 1.0) - m[1:] * (m[1:] + 1.0))
    J_plus = np.zeros((dim, dim), dtype=complex)
    J_plus[np.arange(dim - 1), np.arange(1, dim)] = off
    J_x = 0.5 * (J_plus + J_plus.conj().T)
    return J_x, J_z, j


def _equivariant_layer_generators(n_agents: int):
    """Normalized generators G so that ||G|| = O(1) independent of N."""
    J_x, J_z, j = _collective_spin_operators(n_agents)
    return [J_z / j, J_x / j], J_z / j


@lru_cache(maxsize=256)
def _equivariant_variance_cached(n_agents: int, n_samples: int, seed: int) -> float:
    return float(
        _equivariant_gradient_variance_impl(n_agents, n_samples, seed, return_check=False)
    )


def equivariant_local_cost_gradient_variance(
    n_agents: int,
    n_samples: int = 256,
    seed: int = 42,
    return_check: bool = False,
):
    """
    Variance of dC/dtheta over random parameter settings, for the local cost C = <J_z>/j
    on the symmetric subspace. Exact analytic gradients.
    """
    if not return_check:
        return _equivariant_variance_cached(int(n_agents), int(n_samples), int(seed))
    return _equivariant_gradient_variance_impl(n_agents, n_samples, seed, return_check=True)


def _equivariant_gradient_variance_impl(
    n_agents: int,
    n_samples: int = 256,
    seed: int = 42,
    return_check: bool = False,
):
    N = max(2, int(n_agents))
    generators, observable = _equivariant_layer_generators(N)
    n_layers = max(1, int(np.ceil(np.log2(N))))
    gens = [g for _ in range(n_layers) for g in generators]     # 2 parameters per layer
    n_params = len(gens)
    dim = N + 1

    # Symmetric initial state: the fully polarized Dicke state |j, j>, i.e. all agents aligned.
    psi_0 = np.zeros(dim, dtype=complex)
    psi_0[0] = 1.0

    # Each generator is Hermitian, so diagonalize once and build exp(-i theta G) by rephasing
    # the eigenvalues. Exact, and far cheaper than a matrix exponential per sample.
    eig = [np.linalg.eigh(g) for g in gens]

    def unitary(k: int, theta: float) -> np.ndarray:
        evals, evecs = eig[k]
        return (evecs * np.exp(-1.0j * float(theta) * evals)) @ evecs.conj().T

    rng = np.random.default_rng(seed)
    grads: List[float] = []
    fd_error = 0.0

    def cost(thetas: np.ndarray) -> float:
        psi = psi_0.copy()
        for k, t in enumerate(thetas):
            psi = unitary(k, float(t)) @ psi
        return float(np.real(psi.conj() @ (observable @ psi)))

    for _ in range(int(n_samples)):
        thetas = rng.uniform(0.0, 2.0 * np.pi, size=n_params)
        unitaries = [unitary(k, float(t)) for k, t in enumerate(thetas)]

        # Forward states: forward[k] = U_k ... U_1 |psi_0>
        forward = [psi_0.copy()]
        for u in unitaries:
            forward.append(u @ forward[-1])
        psi = forward[-1]
        o_psi = observable @ psi

        # Suffix propagation of the observable-side vector, so each gradient is O(dim^2).
        # G_k commutes with U_k, so d_k|psi> = U_{>=k} (-i G_k) forward[k]; chi must therefore
        # already carry U_k^dagger when the gradient at index k is taken.
        chi = o_psi.copy()
        sample_grads = np.empty(n_params, dtype=float)
        for k in range(n_params - 1, -1, -1):
            chi = unitaries[k].conj().T @ chi               # chi = U_{>=k}^dagger O |psi>
            d_psi_k = -1.0j * (gens[k] @ forward[k])        # (-i G_k) U_{<k} |psi_0>
            sample_grads[k] = 2.0 * float(np.real(np.vdot(chi, d_psi_k)))
        grads.extend(sample_grads.tolist())

        if return_check:
            h = 1e-5
            for k in range(n_params):
                plus, minus = thetas.copy(), thetas.copy()
                plus[k] += h
                minus[k] -= h
                fd = (cost(plus) - cost(minus)) / (2.0 * h)
                fd_error = max(fd_error, abs(fd - sample_grads[k]))

    variance = float(np.var(np.asarray(grads, dtype=float)))
    if not return_check:
        return variance
    return {
        "gradient_variance": variance,
        "gradient_samples": grads,
        "max_finite_difference_error": float(fd_error),
        "n_layers": n_layers,
        "n_parameters": n_params,
        "symmetric_subspace_dim": dim,
    }


# ---------------------------------------------------------------------------------------
# 2. Unstructured deep ansatz with a global cost - the barren-plateau baseline
# ---------------------------------------------------------------------------------------

_PAULI = {
    "x": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    "y": np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex),
    "z": np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
}


def _apply_single_qubit(psi: np.ndarray, gate: np.ndarray, qubit: int, n_qubits: int) -> np.ndarray:
    psi = psi.reshape([2] * n_qubits)
    psi = np.tensordot(gate, psi, axes=([1], [qubit]))
    psi = np.moveaxis(psi, 0, qubit)
    return psi.reshape(-1)


def _apply_cz_ring(psi: np.ndarray, n_qubits: int) -> np.ndarray:
    """Entangling layer: CZ on every neighbouring pair plus the wrap-around pair."""
    idx = np.arange(psi.size)
    bits = ((idx[:, None] >> np.arange(n_qubits - 1, -1, -1)[None, :]) & 1)
    pairs = [(q, (q + 1) % n_qubits) for q in range(n_qubits)] if n_qubits > 2 else [(0, 1)]
    phase = np.ones(psi.size, dtype=complex)
    for a, b in pairs:
        phase = phase * np.where((bits[:, a] == 1) & (bits[:, b] == 1), -1.0, 1.0)
    return psi * phase


def _global_cost_circuit(thetas: np.ndarray, axes: List[str], qubits: List[int], n_qubits: int,
                         n_layers: int) -> np.ndarray:
    psi = np.zeros(2**n_qubits, dtype=complex)
    psi[0] = 1.0
    p = 0
    for _ in range(n_layers):
        for _q in range(n_qubits):
            gate = expm(-0.5j * float(thetas[p]) * _PAULI[axes[p]])
            psi = _apply_single_qubit(psi, gate, qubits[p], n_qubits)
            p += 1
        psi = _apply_cz_ring(psi, n_qubits)
    return psi


@lru_cache(maxsize=64)
def global_cost_gradient_variance(n_qubits: int, n_samples: int = 200, seed: int = 42) -> float:
    """
    Variance of dC/dtheta_0 for the GLOBAL cost C = |<0...0|psi(theta)>|^2 on a deep
    unstructured hardware-efficient ansatz. Gradients by the parameter-shift rule, which is
    exact for these single-qubit Pauli rotations. Simulated on the full 2^n state vector.
    """
    n = int(n_qubits)
    n_layers = max(2, n)                      # deep: depth grows with width
    n_params = n_layers * n
    rng = np.random.default_rng(seed + n)
    qubits = [i % n for i in range(n_params)]

    grads = np.empty(int(n_samples), dtype=float)
    for s in range(int(n_samples)):
        # The rotation axis is drawn afresh for every gate in every sample. A fixed
        # deterministic axis pattern such as i % 3 aliases against the width whenever n is a
        # multiple of 3 - every qubit then rotates about one axis in every layer, the circuit
        # stops scrambling, and the measured variance jumps by orders of magnitude at n = 6
        # and n = 12. Drawing the axes randomly removes that artifact.
        axes = [["x", "y", "z"][i] for i in rng.integers(0, 3, size=n_params)]
        thetas = rng.uniform(0.0, 2.0 * np.pi, size=n_params)
        plus, minus = thetas.copy(), thetas.copy()
        plus[0] += 0.5 * np.pi
        minus[0] -= 0.5 * np.pi
        grads[s] = 0.5 * (
            abs(_global_cost_circuit(plus, axes, qubits, n, n_layers)[0]) ** 2
            - abs(_global_cost_circuit(minus, axes, qubits, n, n_layers)[0]) ** 2
        )
    return float(np.var(grads))


# ---------------------------------------------------------------------------------------
# 3. Reporting, with measurement and extrapolation kept apart
# ---------------------------------------------------------------------------------------

EQUIVARIANT_FIT_SIZES: List[int] = [8, 16, 32, 64, 128]


def equivariant_scaling_exponent(n_samples: int = 256, seed: int = 42) -> Dict:
    """
    Fits Var[dC/dtheta] ~ N^(-p) for the equivariant local-cost ansatz.

    This is the quantity the cited barren-plateau results actually speak about: freedom from
    barren plateaus means the variance decays only POLYNOMIALLY in N, i.e. a finite p, as
    opposed to exponentially. It is not the claim that the variance exceeds any fixed constant.
    """
    sizes = np.array(EQUIVARIANT_FIT_SIZES, dtype=float)
    variances = np.array(
        [equivariant_local_cost_gradient_variance(n_agents=int(n), n_samples=n_samples, seed=seed)
         for n in sizes]
    )
    slope, intercept = np.polyfit(np.log2(sizes), np.log2(variances), 1)
    residuals = np.log2(variances) - (slope * np.log2(sizes) + intercept)
    return {
        "polynomial_decay_exponent_p": float(-slope),
        "log_log_fit_max_residual_log2": float(np.max(np.abs(residuals))),
        "fit_sizes": [int(n) for n in sizes],
        "measured_variances": {int(n): float(v) for n, v in zip(sizes, variances)},
    }


def gradient_variance_scaling(n_agents: int, n_samples: int = 256, seed: int = 42) -> Dict:
    """
    Compares the equivariant local-cost variance at `n_agents` against the unstructured
    global-cost barren plateau.

    The global-cost figure at `n_agents` is an EXTRAPOLATION, flagged as such: the full
    2^N state vector is only simulable at the sizes in `GLOBAL_COST_FIT_SIZES`, and
    log2(Var) is linear in n there, so the value at N is read off that fit.
    """
    var_equivariant = equivariant_local_cost_gradient_variance(
        n_agents=n_agents, n_samples=n_samples, seed=seed
    )

    measured = {n: global_cost_gradient_variance(n_qubits=n, n_samples=1500, seed=seed)
                for n in GLOBAL_COST_FIT_SIZES}
    sizes = np.array(sorted(measured), dtype=float)
    log2_vars = np.array([np.log2(measured[int(n)]) for n in sizes])
    slope, intercept = np.polyfit(sizes, log2_vars, 1)
    fit_residual = float(np.max(np.abs(log2_vars - (slope * sizes + intercept))))
    extrapolated_log2 = slope * float(n_agents) + intercept
    var_global = float(2.0 ** extrapolated_log2)

    return {
        "equivariant_ttn_local_grad_var": var_equivariant,
        "unstructured_global_grad_var": var_global,
        "unstructured_global_grad_var_is_extrapolated": True,
        "global_cost_fit_sizes": [int(n) for n in sizes],
        "global_cost_measured_variances": {int(n): measured[int(n)] for n in sizes},
        "global_cost_log2_slope_per_qubit": float(slope),
        "global_cost_fit_max_residual_log2": fit_residual,
        "variance_improvement_orders_of_magnitude": float(
            np.log10(max(var_equivariant, 1e-300)) - np.log10(max(var_global, 1e-300))
        ),
    }
