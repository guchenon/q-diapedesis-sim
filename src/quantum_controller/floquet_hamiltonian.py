"""
Stage 3 - Step 16 (Module B1): Bimodal Floquet-Magnus Hamiltonian Engineering
and Step-Out Frequency-Selective Sub-Swarm Control under a Single Global Magnetic Field.

Grounded in:
- Mahoney, Nelson, Peyer, Nelson & Abbott, Appl. Phys. Lett. 2014 (doi:10.1063/1.4870768)
- Bukov, D'Alessio & Polkovnikov, Adv. Phys. 2015 (doi:10.1080/00018732.2015.1055918)
- Higashikawa, Fujita & Sato, arXiv:1810.01103 (Classical Floquet-Magnus expansion)
- Ceron, Gardi, Petersen & Sitti, PNAS 2023 (doi:10.1073/pnas.2221913120)

Part of the Q-DIAPEDESIS Framework (OSF: https://osf.io/htswk | DOI: 10.5281/zenodo.23599020).
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict
import numpy as np
import yaml
from scipy.linalg import expm


@dataclass
class BimodalFloquetState:
    """Instantaneous and stroboscopically averaged bimodal Floquet control metrics."""
    time_s: float
    b_alpha_mT: float
    b_beta_mT: float
    omega_alpha_Hz: float
    omega_beta_Hz: float
    alpha_rotation_rate_Hz: float
    beta_mesh_stiffness_norm: float
    crosstalk_leakage_ratio: float
    dual_task_fidelity: float
    passes_fidelity_target: bool
    effective_hamiltonian_2x2: np.ndarray


class BimodalFloquetController:
    """
    Synthesizes a bimodal periodic magnetic drive B_global(t) in a single external coil
    to simultaneously actuate two heterogeneous sub-swarms in 15-50 um arterioles:
      - Group alpha (60%): high step-out frequency (f_so_alpha = 95 Hz), driven at f_alpha = 80 Hz
        for proximal helical clot drilling.
      - Group beta (40%): low step-out frequency (f_so_beta = 12 Hz), driven at f_beta = 5 Hz
        for distal embolic protection mesh stabilization; adiabatically averages out the 80 Hz mode.

    Epistemic note (Audit Table Row 6a/6b):
    Uses the classical spin/Bloch torque representation of rotating magnetic dipoles in viscous flow.
    No hardware quantum coherence at 37 C is claimed. The >= 92% dual-task fidelity is a Status B
    falsifiable design target benchmarked against single-frequency and uniform-field baselines.
    """

    def __init__(self, config_path: Path | str | None = None):
        if config_path is None:
            root_dir = Path(__file__).resolve().parent.parent.parent
            config_path = root_dir / "configs" / "quantum_floquet_config.yaml"

        with open(config_path, "r", encoding="utf-8") as f:
            self.cfg = yaml.safe_load(f)

        drive = self.cfg["bimodal_floquet_drive"]
        self.b_static_mT = float(drive["global_b_field_static_mT"])
        self.f_alpha_Hz = float(drive["alpha_drill_frequency_Hz"])
        self.b_alpha_nom_mT = float(drive["alpha_drive_amplitude_mT"])
        self.f_beta_Hz = float(drive["beta_trap_frequency_Hz"])
        self.b_beta_nom_mT = float(drive["beta_drive_amplitude_mT"])
        self.truncation_order = int(drive["floquet_truncation_order"])
        self.target_fidelity = float(drive["target_dual_task_fidelity_min"])

        # Step-out frequencies engineered via magnetic moment / hydrodynamic drag ratio (Mahoney 2014)
        self.f_stepout_alpha_Hz = 95.0
        self.f_stepout_beta_Hz = 12.0

        # Pauli matrices used as generators of classical SO(3)/SU(2) dipole rotation
        self.sigma_x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
        self.sigma_y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
        self.sigma_z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)

    @staticmethod
    def stepout_rotation_rate_Hz(drive_freq_Hz: float, stepout_freq_Hz: float) -> float:
        """
        Mahoney-Abbott step-out frequency law (Appl. Phys. Lett. 2014, Eq. 5):
          - Below step-out (f <= f_so): synchronous rotation omega_bar = f
          - Above step-out (f > f_so):  asynchronous slip omega_bar = f * (1 - sqrt(1 - (f_so / f)^2))
        """
        f = float(max(0.0, drive_freq_Hz))
        f_so = float(max(1e-3, stepout_freq_Hz))
        if f <= f_so:
            return f
        ratio_sq = (f_so / f) ** 2
        return float(f * (1.0 - np.sqrt(max(0.0, 1.0 - ratio_sq))))

    def compute_floquet_magnus_hamiltonian(
        self,
        b_alpha_mT: float,
        b_beta_mT: float,
        f_alpha_Hz: float | None = None,
        f_beta_Hz: float | None = None,
    ) -> np.ndarray:
        """
        Computes the 2nd-order high-frequency Floquet-Magnus effective Hamiltonian H_F
        (Bukov et al., Adv. Phys. 2015):
            H_F = H^{(0)} + H^{(1)} + H^{(2)} + O(omega_alpha^{-3})
        """
        fa = self.f_alpha_Hz if f_alpha_Hz is None else float(f_alpha_Hz)
        fb = self.f_beta_Hz if f_beta_Hz is None else float(f_beta_Hz)
        omega_a = 2.0 * np.pi * max(fa, 1.0)

        # Zeroth-order time-averaged dipolar + static locking term
        h_0 = 0.5 * self.b_static_mT * self.sigma_z + 0.5 * float(b_beta_mT) * self.sigma_x

        # First-order Magnus commutator term [H_{+1}, H_{-1}] / omega_alpha
        h_plus = 0.25 * float(b_alpha_mT) * (self.sigma_x - 1.0j * self.sigma_y)
        h_minus = 0.25 * float(b_alpha_mT) * (self.sigma_x + 1.0j * self.sigma_y)
        comm_1 = h_plus @ h_minus - h_minus @ h_plus
        h_1 = comm_1 / omega_a

        if self.truncation_order < 2:
            return h_0 + h_1

        # Second-order frequency-mixing correction scaling as (f_beta / f_alpha)^2
        freq_ratio_sq = (fb / max(fa, 1.0)) ** 2
        h_2 = 0.15 * freq_ratio_sq * float(b_alpha_mT * b_beta_mT) * self.sigma_z / omega_a
        return h_0 + h_1 + h_2

    def bimodal_drive_hamiltonian(
        self,
        time_s: float,
        b_alpha_mT: float,
        b_beta_mT: float,
        f_alpha_Hz: float | None = None,
    ) -> np.ndarray:
        """
        The time-dependent bimodal drive itself, written so that its Fourier components are
        exactly the terms `compute_floquet_magnus_hamiltonian` expands (D3, ELN #003):

            H(t) = (B_static/2) sigma_z + (B_beta/2) sigma_x
                   + (B_alpha/2) [cos(omega_a t) sigma_x + sin(omega_a t) sigma_y]

        so H_0 = (B_static/2) sigma_z + (B_beta/2) sigma_x and
        H_{+-1} = (B_alpha/4)(sigma_x -+ i sigma_y). The slow beta mode is taken as static over
        one alpha period, which is the adiabatic separation the class docstring assumes
        (f_beta = 5 Hz vs f_alpha = 80 Hz).
        """
        fa = self.f_alpha_Hz if f_alpha_Hz is None else float(f_alpha_Hz)
        omega_a = 2.0 * np.pi * max(fa, 1.0)
        phase = omega_a * float(time_s)
        return (
            0.5 * self.b_static_mT * self.sigma_z
            + 0.5 * float(b_beta_mT) * self.sigma_x
            + 0.5 * float(b_alpha_mT) * (np.cos(phase) * self.sigma_x + np.sin(phase) * self.sigma_y)
        )

    def exact_one_period_propagator(
        self,
        b_alpha_mT: float,
        b_beta_mT: float,
        f_alpha_Hz: float | None = None,
        n_substeps: int = 2048,
    ) -> np.ndarray:
        """
        Numerically exact stroboscopic propagator U(T) for one alpha drive period, obtained by
        a midpoint product of short-time exponentials of the true H(t) above. Convergence in
        `n_substeps` is asserted in tests/test_floquet_propagator_fidelity.py.
        """
        fa = self.f_alpha_Hz if f_alpha_Hz is None else float(f_alpha_Hz)
        period_s = 1.0 / max(fa, 1.0)
        dt = period_s / int(n_substeps)
        u = np.eye(2, dtype=complex)
        for k in range(int(n_substeps)):
            t_mid = (k + 0.5) * dt
            h = self.bimodal_drive_hamiltonian(t_mid, b_alpha_mT, b_beta_mT, fa)
            u = expm(-1.0j * h * dt) @ u
        return u

    def magnus_propagator_fidelity(
        self,
        b_alpha_mT: float,
        b_beta_mT: float,
        f_alpha_Hz: float | None = None,
        f_beta_Hz: float | None = None,
        truncation_order: int | None = None,
    ) -> float:
        """
        Average gate fidelity between the exact one-period propagator and the propagator the
        truncated Floquet-Magnus effective Hamiltonian predicts:

            F = |Tr(U_magnus^dagger U_exact)|^2 / d^2,   d = 2

        This is a real, computed fidelity obtained by evolving a state, and it is the quantity
        that tells you whether the Floquet engineering in `compute_floquet_magnus_hamiltonian`
        is accurate. It is NOT the preregistered `dual_task_fidelity`, which is a weighted
        design score (see `evaluate_control_step`); the two measure different things and are
        reported separately on purpose (D3, ELN #003).
        """
        fa = self.f_alpha_Hz if f_alpha_Hz is None else float(f_alpha_Hz)
        order_before = self.truncation_order
        if truncation_order is not None:
            self.truncation_order = int(truncation_order)
        try:
            h_eff = self.compute_floquet_magnus_hamiltonian(b_alpha_mT, b_beta_mT, fa, f_beta_Hz)
        finally:
            self.truncation_order = order_before
        period_s = 1.0 / max(fa, 1.0)
        u_magnus = expm(-1.0j * h_eff * period_s)
        u_exact = self.exact_one_period_propagator(b_alpha_mT, b_beta_mT, fa)
        overlap = np.trace(u_magnus.conj().T @ u_exact)
        return float(abs(overlap) ** 2 / 4.0)

    def evaluate_control_step(
        self,
        time_s: float,
        b_scale_factor: float = 1.0,
        f_alpha_override_Hz: float | None = None,
        f_beta_override_Hz: float | None = None,
    ) -> BimodalFloquetState:
        """
        Evaluates sub-swarm alpha drilling rate, sub-swarm beta embolic mesh stiffness,
        cross-talk leakage, and joint dual-task fidelity under the bimodal global field.
        """
        scale = float(np.clip(b_scale_factor, 0.05, 1.5))
        b_alpha = self.b_alpha_nom_mT * scale
        b_beta = self.b_beta_nom_mT * scale

        fa = self.f_alpha_Hz if f_alpha_override_Hz is None else float(f_alpha_override_Hz)
        fb = self.f_beta_Hz if f_beta_override_Hz is None else float(f_beta_override_Hz)

        # Effective step-out scales linearly with applied B-field amplitude
        f_so_alpha_eff = self.f_stepout_alpha_Hz * np.clip(scale, 0.85, 1.25)
        f_so_beta_eff = self.f_stepout_beta_Hz * np.clip(scale, 0.85, 1.25)

        # Group alpha rotation at fa (80 Hz): synchronous if fa <= f_so_alpha_eff
        omega_rot_alpha = self.stepout_rotation_rate_Hz(fa, f_so_alpha_eff)
        alpha_sync_fidelity = float(np.clip(omega_rot_alpha / max(fa, 1.0), 0.0, 1.0))

        # Group beta response to high-frequency fa (80 Hz >> f_so_beta = 12 Hz):
        # Above step-out, parasitic rotation rate decays as ~ 0.5 * (f_so_beta / fa)^2
        beta_parasitic_rot_Hz = self.stepout_rotation_rate_Hz(fa, f_so_beta_eff)
        crosstalk_leakage = float(np.clip(beta_parasitic_rot_Hz / max(fa, 1.0), 0.0, 1.0))

        # Group beta synchronous locking to the slow 5 Hz trapping mode
        omega_rot_beta = self.stepout_rotation_rate_Hz(fb, f_so_beta_eff)
        beta_lock_efficiency = float(np.clip(omega_rot_beta / max(fb, 1.0), 0.0, 1.0))
        beta_mesh_stiffness = float(np.clip(beta_lock_efficiency * (1.0 - 1.2 * crosstalk_leakage), 0.0, 1.0))

        # Joint dual-task fidelity F_dual (Audit Table Row 6b: target >= 0.92)
        h_eff = self.compute_floquet_magnus_hamiltonian(b_alpha, b_beta, fa, fb)
        floquet_decoupling_bonus = 1.0 - 0.5 * ((fb / max(fa, 1.0)) ** 2)
        dual_fidelity = float(
            np.clip(0.52 * alpha_sync_fidelity + 0.48 * beta_mesh_stiffness * floquet_decoupling_bonus, 0.0, 0.995)
        )

        return BimodalFloquetState(
            time_s=float(time_s),
            b_alpha_mT=b_alpha,
            b_beta_mT=b_beta,
            omega_alpha_Hz=fa,
            omega_beta_Hz=fb,
            alpha_rotation_rate_Hz=omega_rot_alpha,
            beta_mesh_stiffness_norm=beta_mesh_stiffness,
            crosstalk_leakage_ratio=crosstalk_leakage,
            dual_task_fidelity=dual_fidelity,
            passes_fidelity_target=bool(dual_fidelity >= self.target_fidelity),
            effective_hamiltonian_2x2=h_eff,
        )

    def compare_with_baselines(self) -> Dict[str, float]:
        """
        Benchmarks bimodal Floquet decoupling (80 Hz + 5 Hz) against:
          1. Single-frequency high-speed drive (80 Hz only -> beta mesh steps out and collapses)
          2. Single-frequency low-speed drive (5 Hz only -> alpha drill rate too slow for lysis)
          3. Open-loop excessive frequency drive (140 Hz -> alpha population steps out)
        """
        bimodal = self.evaluate_control_step(time_s=1.0, b_scale_factor=1.0)

        # Baseline 1: Single 80 Hz drive applied to both tasks (fb = 80 Hz -> beta steps out)
        single_80Hz = self.evaluate_control_step(time_s=1.0, f_alpha_override_Hz=80.0, f_beta_override_Hz=80.0)

        # Baseline 2: Single 5 Hz drive applied to both tasks (alpha drilling loses 93.75% speed)
        single_5Hz_drill_penalty = 5.0 / 80.0
        single_5Hz_fidelity = 0.5 * single_5Hz_drill_penalty + 0.5 * bimodal.beta_mesh_stiffness_norm

        # Baseline 3: Over-driven 140 Hz (> f_stepout_alpha = 95 Hz -> alpha steps out!)
        stepout_140Hz = self.evaluate_control_step(time_s=1.0, f_alpha_override_Hz=140.0, f_beta_override_Hz=5.0)

        return {
            "bimodal_floquet_fidelity": bimodal.dual_task_fidelity,
            "single_80Hz_baseline_fidelity": single_80Hz.dual_task_fidelity,
            "single_5Hz_baseline_fidelity": float(single_5Hz_fidelity),
            "overdriven_140Hz_stepout_fidelity": stepout_140Hz.dual_task_fidelity,
            "bimodal_crosstalk_leakage": bimodal.crosstalk_leakage_ratio,
        }
