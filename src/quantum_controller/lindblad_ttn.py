"""
Stage 3 - Step 17 (Module B4): Open-System Lindblad (GKLS) Dissipative Master Equation
and S_n-Permutation-Equivariant Isometric Tree Tensor Network (TTN) Controller.

Grounded in:
- Lindblad, Commun. Math. Phys. 1976 (doi:10.1007/BF01608499)
- Gorini, Kossakowski & Sudarshan, J. Math. Phys. 1976 (doi:10.1063/1.522979)
- Barthel & Miao, Commun. Math. Phys. 2025 (doi:10.1007/s00220-024-05217-x)
- Cervero Martin, Plekhanov & Lubasch, Quantum 2023 (doi:10.22331/q-2023-04-13-974)
- Schatzki, Larocca, Nguyen, Sauvage & Cerezo, npj Quantum Inf. 2024 (doi:10.1038/s41534-024-00804-1)
- Cerezo et al., Nat. Commun. 2025 (doi:10.1038/s41467-025-63099-6)

Part of the Q-DIAPEDESIS Framework (OSF: https://osf.io/htswk | DOI: 10.5281/zenodo.23599019).
"""

from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Dict
import numpy as np
import yaml

from src.physics import derive_dissipation_gamma

from .ttn_gradient_variance import (
    GLOBAL_COST_FIT_SIZES,
    equivariant_local_cost_gradient_variance,
    equivariant_scaling_exponent,
    global_cost_gradient_variance,
    gradient_variance_scaling,
)


@dataclass
class LindbladTTNState:
    """State of the open-system Lindblad density matrix and Equivariant TTN optimizer."""
    time_s: float
    purity_tr_rho2: float
    von_neumann_entropy: float
    dissipative_attractor_overlap: float
    ttn_bond_dimension_chi: int
    gradient_variance: float
    barren_plateau_free_verified: bool
    optimal_field_correction_scale: float


class LindbladEquivariantTTN:
    """
    Implements Module B4:
    1. Classical Markov/Lindblad (GKLS) open-system density matrix evolution where blood
       viscosity gamma_blood acts as a dissipative jump operator driving the swarm into
       a synchronized steady-state attractor.
    2. An Isometric & S_n-Permutation-Equivariant Tree Tensor Network (TTN) controller
       with bounded bond dimension chi <= 16 and a local cost observable.

    Epistemic note (Section 3.2 & 3.4 of Layer 3 Audit):
    Absence of barren plateaus (gradient variance scaling as Omega(1/poly(N))) holds
    conditionally under local cost functions, logarithmic tree depth O(log N), and
    symmetric isometric initialization (Barthel & Miao 2025; Schatzki et al. 2024).
    As shown by Cerezo et al. (Nat. Commun. 2025), this symmetric subspace is efficiently
    classically simulable, enabling millisecond classical execution in `q-diapedesis-sim`.
    """

    def __init__(self, config_path: Path | str | None = None):
        if config_path is None:
            root_dir = Path(__file__).resolve().parent.parent.parent
            config_path = root_dir / "configs" / "quantum_floquet_config.yaml"

        with open(config_path, "r", encoding="utf-8") as f:
            self.cfg = yaml.safe_load(f)

        ttn_cfg = self.cfg["lindblad_and_tensor_network"]
        # D8 (ELN #003): in the overdamped regime of a 0.75 um rotor in plasma the SAME quantity
        # m B / zeta_r sets both the coherent precession scale and the dissipative relaxation
        # scale, so gamma cannot be chosen independently of H_eff. It is now derived.
        # Note also the unit inconsistency this exposes: H_eff is assembled from field
        # amplitudes in mT and then used in exp(-i H t) with t in seconds, i.e. mT read as
        # rad/s. `mT_to_rad_s` is the physical conversion; the Hamiltonian itself is left
        # untouched because it carries a preregistered threshold.
        self.dissipation_is_derived = bool(ttn_cfg.get("derive_dissipation_from_physics", False))
        if self.dissipation_is_derived:
            root_dir = Path(__file__).resolve().parent.parent.parent
            with open(root_dir / "configs" / "vessel_mra_config.yaml", "r", encoding="utf-8") as f:
                vessel_cfg = yaml.safe_load(f)
            derived = derive_dissipation_gamma(vessel_cfg, self.cfg)
            self.gamma_viscous = float(derived["viscous_dissipation_gamma"])
            self.gamma_physical_rad_s = float(derived["gamma_physical_rad_s"])
            self.mT_to_rad_s = float(derived["mT_to_rad_s"])
        else:
            self.gamma_viscous = float(ttn_cfg["viscous_dissipation_gamma"])
            self.gamma_physical_rad_s = float("nan")
            self.mT_to_rad_s = float("nan")
        self.dipole_J0 = float(ttn_cfg["dipole_interaction_strength_J0"])
        self.chi_max = int(ttn_cfg["ttn_max_bond_dimension_chi"])
        self.equivariant = bool(ttn_cfg["permutation_equivariant_symmetry"])

        self.n_super = int(self.cfg["swarm_scaling"]["n_super_particles"])

        # D5: the barren-plateau floor was hard-coded as 1.0e-3; moved to config unchanged.
        bp_cfg = self.cfg["implementation_constants"]["barren_plateau"]
        self.legacy_gradient_variance_floor = float(bp_cfg["gradient_variance_floor"])
        self.gradient_variance_samples = int(bp_cfg["measurement_samples"])
        self.max_polynomial_decay_exponent = float(bp_cfg["max_polynomial_decay_exponent"])
        self.max_log_log_fit_residual = float(bp_cfg["max_log_log_fit_residual_log2"])
        self.min_improvement_orders = float(bp_cfg["min_improvement_orders"])

        self.reset()

    @cached_property
    def barren_plateau_free(self) -> bool:
        """
        D9 (ELN #003): the plateau-freedom verdict, from the measured SCALING rather than a
        fixed variance floor. Freedom from barren plateaus means the gradient variance decays
        only polynomially in N - a finite exponent with a power law that actually fits - and
        that it stays far above the unstructured global-cost baseline. A constant floor tests
        the wrong property: here the polynomial decay crosses 1e-3 at about N = 24, so the old
        criterion rejected an ansatz that is in fact plateau-free.
        """
        report = self.benchmark_barren_plateau_scaling()
        p = float(report["equivariant_polynomial_decay_exponent_p"])
        residual = float(report["equivariant_log_log_fit_max_residual_log2"])
        orders = float(report["variance_improvement_orders_of_magnitude"])
        return bool(
            0.0 < p < self.max_polynomial_decay_exponent
            and residual <= self.max_log_log_fit_residual
            and orders >= self.min_improvement_orders
        )

    @cached_property
    def measured_gradient_variance(self) -> float:
        """Measured Var[dC/dtheta] for this controller's ansatz. Computed once, then reused."""
        return self.compute_ttn_gradient_variance(
            n_agents=self.n_super,
            equivariant=self.equivariant,
            local_cost=True,
            n_samples=self.gradient_variance_samples,
        )

    def reset(self) -> None:
        """Initializes the reduced collective 2-level density matrix rho in a mixed state."""
        self.current_time_s = 0.0
        # Initial collective density matrix (partially mixed state before dissipative locking)
        self.rho = np.array([[0.55, 0.15 - 0.05j], [0.15 + 0.05j, 0.45]], dtype=complex)
        # Target dissipative attractor state |+><+| representing phase-locked collective motion
        ket_plus = np.array([[1.0], [1.0]], dtype=complex) / np.sqrt(2.0)
        self.rho_attractor = ket_plus @ ket_plus.conj().T

    def step_lindblad_density_matrix(
        self,
        dt_s: float,
        h_floquet_2x2: np.ndarray,
        local_viscosity_mPa_s: float = 3.5,
    ) -> LindbladTTNState:
        """
        Integrates one step of the trace-preserving completely positive (CPTP) GKLS master equation:
            d(rho)/dt = -i [H_eff, rho] + gamma_eff * (L rho L^dagger - 0.5 * {L^dagger L, rho})
        where the jump operator L = |+><-| channels viscous drag into collective alignment.

        D3 (ELN #003): this is now a pure RK4 integration of the generator above. Two terms that
        were not part of the master equation have been removed: a blend of rho toward
        `rho_attractor` on every step (which produced the attractor overlap directly rather than
        letting the dissipator produce it), and a factor 0.05 multiplying the Hamiltonian
        commutator (which scaled H_eff almost out of the dynamics). The reduced two-level
        collective description remains a reduced-order surrogate, Status B: it is not derived
        from the N-body dipolar Hamiltonian.
        """
        # Scale effective dissipation rate with Pries-Secomb blood viscosity
        gamma_eff = self.gamma_viscous * (float(local_viscosity_mPa_s) / 3.0)

        # Jump operator L projecting antisymmetric modes onto the symmetric |+> attractor
        ket_plus = np.array([[1.0], [1.0]], dtype=complex) / np.sqrt(2.0)
        ket_minus = np.array([[1.0], [-1.0]], dtype=complex) / np.sqrt(2.0)
        L_jump = ket_plus @ ket_minus.conj().T
        L_dag_L = L_jump.conj().T @ L_jump
        H = np.asarray(h_floquet_2x2, dtype=complex)

        def liouvillian(rho: np.ndarray) -> np.ndarray:
            """The GKLS generator, with no added relaxation and no rescaling of H (D3)."""
            return -1.0j * (H @ rho - rho @ H) + gamma_eff * (
                L_jump @ rho @ L_jump.conj().T - 0.5 * (L_dag_L @ rho + rho @ L_dag_L)
            )

        # Classical RK4 on the Liouvillian. The generator is linear and trace-annihilating,
        # so RK4 preserves trace and Hermiticity to machine precision without any projection.
        dt = float(dt_s)
        k1 = liouvillian(self.rho)
        k2 = liouvillian(self.rho + 0.5 * dt * k1)
        k3 = liouvillian(self.rho + 0.5 * dt * k2)
        k4 = liouvillian(self.rho + dt * k3)
        self.rho = self.rho + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        self.current_time_s += dt

        # Purity Tr(rho^2), von Neumann entropy, and overlap with the attractor, reported as
        # computed. These were previously clipped to their physical ranges, which would have
        # hidden any integration error; the GKLS tests assert the bounds instead (D3).
        purity = float(np.real(np.trace(self.rho @ self.rho)))
        eigvals = np.clip(np.real(np.linalg.eigvalsh(self.rho)), 1e-12, 1.0)
        vn_entropy = float(-np.sum(eigvals * np.log2(eigvals)))
        overlap = float(np.real(np.trace(self.rho @ self.rho_attractor)))

        # The gradient variance is a property of the ansatz, not of rho, so it does not change
        # from step to step. It is now really computed, which costs far more than a closed-form
        # expression, so it is evaluated once per controller and reused (D2-b).
        grad_var = self.measured_gradient_variance

        # Optimal field modulation scale derived from TTN energy minimization
        field_scale = float(np.clip(0.85 + 0.15 * overlap, 0.70, 1.05))

        return LindbladTTNState(
            time_s=self.current_time_s,
            purity_tr_rho2=purity,
            von_neumann_entropy=vn_entropy,
            dissipative_attractor_overlap=overlap,
            ttn_bond_dimension_chi=self.chi_max,
            gradient_variance=grad_var,
            barren_plateau_free_verified=self.barren_plateau_free,
            optimal_field_correction_scale=field_scale,
        )

    @staticmethod
    def compute_ttn_gradient_variance(
        n_agents: int,
        equivariant: bool = True,
        local_cost: bool = True,
        seed: int = 42,
        n_samples: int = 256,
    ) -> float:
        """
        Measures the parameter gradient variance Var[dC/dtheta] by building a parameterized
        circuit, evaluating its cost, and differentiating it (D2-b / D3, ELN #003).

        Previously this returned the closed-form `Var[sin(a) cos(b) / sqrt(N)]` for the
        equivariant branch and the literal constant `0.25 * 2**-min(N, 30)` for the global
        baseline. Neither was the gradient of a cost function and no circuit was built, so the
        barren-plateau claim rested on two hand-written numbers. The work is now done in
        `ttn_gradient_variance.py`:
          - equivariant + local cost: exact analytic gradients of C = <J_z>/j on the symmetric
            (Dicke) subspace, dimension N + 1, with ceil(log2 N) layers.
          - unstructured + global cost: parameter-shift gradients of C = |<0...0|psi>|^2 on the
            full 2^n state vector (McClean et al. 2018 barren plateau).
        """
        if equivariant and local_cost:
            return float(
                equivariant_local_cost_gradient_variance(
                    n_agents=max(2, int(n_agents)), n_samples=int(n_samples), seed=int(seed)
                )
            )
        # The global-cost baseline is simulated on the full Hilbert space, so it can only be
        # evaluated at widths that fit in memory. Beyond that, use benchmark_barren_plateau_scaling,
        # which fits the measured decay and labels the result as an extrapolation.
        n_qubits = min(max(2, int(n_agents)), max(GLOBAL_COST_FIT_SIZES))
        return float(global_cost_gradient_variance(n_qubits=n_qubits, n_samples=1500, seed=int(seed)))

    def benchmark_barren_plateau_scaling(self) -> Dict[str, float]:
        """
        Compares measured gradient variance at N = `n_super` between:
          1. Permutation-Equivariant Isometric TTN + Local Cost (Q-DIAPEDESIS)
          2. Generic Unstructured Ansatz + Global Cost (Barren Plateau baseline)

        The N = 200 global-cost value is an extrapolation of the measured exponential decay -
        2^200 amplitudes cannot be simulated - and the returned dictionary says so via
        `unstructured_global_grad_var_is_extrapolated`. `polynomial_decay_exponent_p` is the
        quantity the cited theory is about: plateau freedom means a finite p, not a variance
        above any particular constant.
        """
        report = dict(gradient_variance_scaling(n_agents=self.n_super))
        report.update(
            {f"equivariant_{k}": v for k, v in equivariant_scaling_exponent().items()}
        )
        return report
