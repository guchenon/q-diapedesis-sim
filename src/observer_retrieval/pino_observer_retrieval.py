"""
Stage 3 - Steps 18 & 19 (Modules B2 & B5):
- Step 18 (Module B2): Physics-Informed Neural Operator (PINO) Sub-Voxel State Observer
  for Swarm Tracking across 50-200 ms MRI Actuation-Imaging Duty-Cycle Blind Windows.
- Step 19 (Module B5): Post-Task Exponential Degaussing & Retrograde Catch-Bond Wall-Rolling
  Retrieval with Unified Mass-Balance Accounting.

Grounded in:
- Martel et al., Comput. Aided Surg. 2008 (doi:10.3109/10929080802551274) & IJRR 2009 (doi:10.1177/0278364908104855)
- He, Feng, Yang G-Z et al., Nat. Commun. 2023 (doi:10.1038/s41467-023-43966-w)
- Li, Kovachki, Anandkumar et al., ICLR 2021 (arXiv:2010.08895) & ACM/IMS J. Data Sci. 2024 (doi:10.1145/3648506)
- Wang B et al., Sci. Adv. 2024 (doi:10.1126/sciadv.adk8970)
- Alapan, Bozuyuk, Erkoc, Karacakol & Sitti, Sci. Robot. 2020 (doi:10.1126/scirobotics.aba5726)

Part of the Q-DIAPEDESIS Framework (OSF: https://osf.io/htswk | DOI: 10.5281/zenodo.23599019).
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict
import numpy as np
import yaml


@dataclass
class SubVoxelObserverReport:
    """Performance metrics of the PINO sub-voxel observer across an MRI blind window (Step 18)."""
    blind_window_ms: float
    update_interval_ms: float
    pino_p95_error_um: float
    pino_median_error_um: float
    kalman_p95_error_um: float
    harmonic_kalman_p95_error_um: float
    dead_reckoning_p95_error_um: float
    last_hold_p95_error_um: float
    passes_subvoxel_target: bool


@dataclass
class RetrievalMassBalanceReport:
    """Post-task degaussing and retrograde wall-rolling retrieval metrics (Step 19)."""
    degaussing_tau_ms: float
    residual_dipole_cohesion_ratio: float
    retrograde_velocity_um_s: float
    clearance_time_s: float
    transient_washout_events_pct: float
    recovered_transient_pct: float
    final_recovered_fraction: float
    unrecovered_loss_fraction: float
    passes_retrieval_target: bool


class PINOSubVoxelObserver:
    """
    Step 18 (Module B2): Predicts the sub-voxel axial-radial swarm centroid and density field
    every 2.0 ms during the 50-200 ms MRI imaging-actuation blind window by coupling spectral
    Fourier modes of the pulsatile advection equation with RF phase-contrast tomography priors.

    WHAT THIS CLASS ACTUALLY IS (D11, ELN #003). Despite the name, there is no neural network
    here and no operator is learned. The estimator is a **Physics-Informed Spectral MAP
    Observer**: a five-mode regularized Fourier fit of the cerebral haemodynamic waveform
    (mean, first and second cardiac harmonics in sine and cosine), solved as a
    maximum-a-posteriori problem whose weighting lambda = sigma_v^2 / sigma_c^2 follows from the
    measurement and prior uncertainties. `PhysicsInformedSpectralObserver` is provided as an
    alias and is the accurate name; the `PINOSubVoxelObserver` name is kept so imports and the
    registered module layout do not break. The Li et al. (ICLR 2021) neural-operator citations
    describe a method this code does not use.

    WHERE THE ACCURACY COMES FROM. The spectral observer is given the full cardiac waveform as
    a prior; plain dead reckoning gets only the mean, and the plain Kalman filter gets the mean
    plus a constant-velocity model. To separate "does the physics prior help?" from "is this
    estimator better than a Kalman filter?", a HARMONIC-AWARE KALMAN FILTER is also reported:
    the same harmonic model, carried in the filter state, processed sequentially. If that
    baseline matches the spectral observer, the advantage is the prior and not the architecture.

    Epistemic note (Audit Table Row 8):
    The 50-200 ms blind window (100 ms default) and < 8.0 um 95th-percentile tracking error
    are Status B falsifiable design targets benchmarked against Kalman filter, harmonic-aware
    Kalman filter, dead-reckoning, and last-observation-hold baselines.
    """

    def __init__(self, config_path: Path | str | None = None):
        if config_path is None:
            root_dir = Path(__file__).resolve().parent.parent.parent
            config_path = root_dir / "configs" / "quantum_floquet_config.yaml"

        with open(config_path, "r", encoding="utf-8") as f:
            self.cfg = yaml.safe_load(f)

        obs_cfg = self.cfg["mri_sub_voxel_observer"]
        self.default_blind_ms = float(obs_cfg["mri_blind_window_ms"])
        self.dt_obs_ms = float(obs_cfg["pino_update_interval_ms"])
        self.n_shots = int(obs_cfg["tomography_measurement_shots"])
        self.error_ceiling_um = 8.0

    def evaluate_blind_window_tracking(
        self,
        blind_window_ms: float | None = None,
        mean_flow_velocity_um_s: float = 280.0,
        pulsatile_amplitude_ratio: float = 0.25,
        seed: int = 42,
        n_shots_override: int | None = None,
        model_mismatch_um_s: float = 40.0,
        return_traces: bool = False,
    ):
        """
        Simulates ground-truth pulsatile swarm advection during an MRI blind window and compares
        four estimators that each see ONLY the last pre-window observation, the stated physics
        prior (cardiac frequency and mean flow), and noisy RF phase-contrast VELOCITY
        measurements at `pino_update_interval_ms`.

        D2-b / D3 (ELN #003). Previously the "PINO" estimate was `z_gt + noise` - it was handed
        the ground-truth trajectory - the "Kalman filter" read the ground-truth velocity
        `v_true` with a 0.55 gain, and the Kalman error had `+ 3.2 * (t/T)` added to it so that
        the PINO would win. None of that is present now. The estimators are:

          1. PINO: a spectral, physics-informed operator. The cardiac fundamental is known, so
             the velocity is modelled as a truncated Fourier series in that frequency and its
             second harmonic, the coefficients are fitted by least squares to the noisy velocity
             measurements over the window, and the fitted waveform is integrated to position.
             Averaging over the window is what suppresses shot noise, which is also why the
             advantage is expected to appear only for long enough windows.
          2. Kalman: a textbook linear KF on [z, v] with a constant-velocity process model,
             consuming the same velocity measurements. No handicap.
          2b. Harmonic-aware Kalman (D11): the same KF but carrying the cardiac harmonic
             coefficients in its state, so it has the SAME physics prior as the spectral
             observer. Reported to separate the value of the prior from the value of the method.
          3. Dead reckoning: the pre-window velocity estimate, held constant. Ignores the
             measurements.
          4. Last-observation hold: frozen at the last pre-window position.

        INVERSE-CRIME NOTE: the PINO's basis is the cardiac harmonics, which also generate the
        truth, so a matched-basis comparison would only be testing that the basis is exact.
        `model_mismatch_um_s` therefore adds a component to the truth that lies OUTSIDE that
        basis (an exponential catch-bond retardation transient), and it is non-zero by default.
        Set it to 0.0 to see the optimistic matched-basis case.
        """
        win_ms = self.default_blind_ms if blind_window_ms is None else float(blind_window_ms)
        n_shots = self.n_shots if n_shots_override is None else int(n_shots_override)
        n_steps = max(5, int(round(win_ms / self.dt_obs_ms)))
        t_s = np.linspace(0.0, win_ms * 1.0e-3, n_steps)
        dt_s = float(t_s[1] - t_s[0])

        rng = np.random.default_rng(seed)
        omega_c = 2.0 * np.pi * 1.25  # 75 bpm cardiac fundamental, known to every estimator
        wall_coupling = 0.32          # fraction of bulk flow seen by a wall-rolling agent

        # ----- ground truth: known to the simulator only -----
        v_pulsatile = mean_flow_velocity_um_s * (
            1.0
            + pulsatile_amplitude_ratio
            * (0.78 * np.sin(omega_c * t_s + 0.4) + 0.22 * np.sin(2.0 * omega_c * t_s))
        )
        # Component deliberately outside the estimator's Fourier basis (see inverse-crime note).
        catch_bond_transient = float(model_mismatch_um_s) * np.exp(-t_s / 0.045)
        v_true_z = wall_coupling * v_pulsatile + catch_bond_transient + rng.normal(0.0, 4.5, size=n_steps)
        z_true = 250.0 + np.cumsum(v_true_z * dt_s)
        r_true = 13.2 + 0.6 * np.sin(omega_c * t_s) + rng.normal(0.0, 0.25, size=n_steps)

        # ----- what the observers are allowed to see -----
        z_0 = float(z_true[0])
        r_0 = float(r_true[0])
        v_prior = wall_coupling * mean_flow_velocity_um_s        # physics prior, no truth
        # Phase-contrast velocity noise, stated a priori: the per-sample position resolution of
        # the tomography is sigma_pos = 18 um / sqrt(shots) (the scale the previous code used as
        # a position noise), and a single 2 ms sample resolves velocity only to sigma_pos / dt.
        # At 512 shots and dt = 2 ms that is ~398 um/s, several times the mean wall velocity, so
        # a single sample is nearly useless and averaging over the window is what buys accuracy.
        sigma_pos_um = 18.0 / np.sqrt(float(n_shots))
        sigma_v = sigma_pos_um / max(dt_s, 1e-9)
        v_meas = v_true_z + rng.normal(0.0, sigma_v, size=n_steps)

        # 1. Last-observation hold
        err_hold = np.hypot(z_true - z_0, r_true - r_0)

        # 2. Constant-velocity dead reckoning from the prior (ignores the measurements)
        z_dr = z_0 + v_prior * t_s
        err_dr = np.hypot(z_true - z_dr, r_true - r_0)

        # 3. Linear Kalman filter on [z, v], constant-velocity model, velocity measurements
        x = np.array([z_0, v_prior])
        P = np.diag([1.0, (0.25 * v_prior) ** 2])
        Q = np.diag([(2.0 * dt_s) ** 2, (900.0 * dt_s) ** 2])   # process noise on z and v
        R = np.array([[sigma_v**2]])
        H = np.array([[0.0, 1.0]])
        F = np.array([[1.0, dt_s], [0.0, 1.0]])
        z_kf = np.empty(n_steps)
        for k in range(n_steps):
            if k > 0:
                x = F @ x
                P = F @ P @ F.T + Q
            resid = np.array([v_meas[k]]) - H @ x
            S = H @ P @ H.T + R
            K = P @ H.T @ np.linalg.inv(S)
            x = x + (K @ resid).ravel()
            P = P - K @ H @ P
            z_kf[k] = x[0]
        err_kf = np.hypot(z_true - z_kf, r_true - r_0)

        # 3b. D11: HARMONIC-AWARE Kalman filter. Same measurements, same physics prior as the
        #     spectral observer - the cardiac harmonics are carried in the filter state - but
        #     processed sequentially rather than as one batch MAP fit. This isolates the value of
        #     the prior from the value of the estimator: if this matches the spectral observer,
        #     the advantage is the prior.
        basis_full = np.column_stack(
            [
                np.ones_like(t_s),
                np.sin(omega_c * t_s),
                np.cos(omega_c * t_s),
                np.sin(2.0 * omega_c * t_s),
                np.cos(2.0 * omega_c * t_s),
            ]
        )
        c_prior_hk = np.array(
            [
                wall_coupling * mean_flow_velocity_um_s,
                wall_coupling * mean_flow_velocity_um_s * pulsatile_amplitude_ratio * 0.78 * np.cos(0.4),
                wall_coupling * mean_flow_velocity_um_s * pulsatile_amplitude_ratio * 0.78 * np.sin(0.4),
                wall_coupling * mean_flow_velocity_um_s * pulsatile_amplitude_ratio * 0.22,
                0.0,
            ]
        )
        sigma_c_hk = np.maximum(0.25 * np.abs(c_prior_hk), 5.0)
        x_h = np.concatenate([[z_0], c_prior_hk])                       # [z, c0..c4]
        P_h = np.diag(np.concatenate([[1.0], sigma_c_hk**2]))
        R_h = np.array([[sigma_v**2]])
        z_hk = np.empty(n_steps)
        for k in range(n_steps):
            b_k = basis_full[k]
            if k > 0:
                F_h = np.eye(6)
                F_h[0, 1:] = b_k * dt_s                                  # z advances by basis . c
                x_h = F_h @ x_h
                P_h = F_h @ P_h @ F_h.T
            H_h = np.concatenate([[0.0], b_k]).reshape(1, 6)             # velocity is observed
            S_h = H_h @ P_h @ H_h.T + R_h
            K_h = P_h @ H_h.T @ np.linalg.inv(S_h)
            x_h = x_h + (K_h @ (np.array([v_meas[k]]) - H_h @ x_h)).ravel()
            P_h = P_h - K_h @ H_h @ P_h
            z_hk[k] = x_h[0]

        # 4. PINO: spectral least-squares fit of the velocity on the known cardiac harmonics,
        #    then integrated. The design matrix contains no information about the truth beyond
        #    the frequency, which is part of the stated prior.
        basis = np.column_stack(
            [
                np.ones_like(t_s),
                np.sin(omega_c * t_s),
                np.cos(omega_c * t_s),
                np.sin(2.0 * omega_c * t_s),
                np.cos(2.0 * omega_c * t_s),
            ]
        )
        # The fit is a MAP estimate, not an unregularized least squares. The physics prior gives
        # the expected waveform coefficients, and the weighting is DERIVED from the two
        # uncertainties - lambda = sigma_v^2 / sigma_c^2 - with no tunable knob: the estimator
        # falls back on the prior exactly as far as the measurements are uninformative. An
        # unregularized 5-parameter fit is a strawman at this SNR, losing to dead reckoning,
        # which gets the same prior for free and ignores the noise.
        c_prior = np.array(
            [
                wall_coupling * mean_flow_velocity_um_s,
                wall_coupling * mean_flow_velocity_um_s * pulsatile_amplitude_ratio * 0.78 * np.cos(0.4),
                wall_coupling * mean_flow_velocity_um_s * pulsatile_amplitude_ratio * 0.78 * np.sin(0.4),
                wall_coupling * mean_flow_velocity_um_s * pulsatile_amplitude_ratio * 0.22,
                0.0,
            ]
        )
        # A-priori prior spread: 25% of each coefficient, floored so a zero coefficient is not
        # pinned. Same 25% figure as the configured pulsatile amplitude uncertainty.
        sigma_c = np.maximum(0.25 * np.abs(c_prior), 5.0)
        lhs = basis.T @ basis / sigma_v**2 + np.diag(1.0 / sigma_c**2)
        rhs = basis.T @ v_meas / sigma_v**2 + c_prior / sigma_c**2
        coeffs = np.linalg.solve(lhs, rhs)
        v_fit = basis @ coeffs
        z_pino = z_0 + np.cumsum(v_fit * dt_s)
        # Radial position: the prior is the same cardiac oscillation about the margination band.
        r_pino = r_0 + 0.6 * (np.sin(omega_c * t_s) - np.sin(0.0))
        err_pino = np.hypot(z_true - z_pino, r_true - r_pino)
        # The harmonic-aware KF shares the spectral observer's radial prior, so the only thing
        # that differs between them is how the axial estimate is formed (D11).
        err_hk = np.hypot(z_true - z_hk, r_true - r_pino)

        p95_pino = float(np.percentile(err_pino, 95))
        report = SubVoxelObserverReport(
            blind_window_ms=win_ms,
            update_interval_ms=self.dt_obs_ms,
            pino_p95_error_um=p95_pino,
            pino_median_error_um=float(np.median(err_pino)),
            kalman_p95_error_um=float(np.percentile(err_kf, 95)),
            harmonic_kalman_p95_error_um=float(np.percentile(err_hk, 95)),
            dead_reckoning_p95_error_um=float(np.percentile(err_dr, 95)),
            last_hold_p95_error_um=float(np.percentile(err_hold, 95)),
            passes_subvoxel_target=bool(p95_pino < self.error_ceiling_um),
        )
        if not return_traces:
            return report
        return {
            "report": report,
            "true_z_um": z_true,
            "true_r_um": r_true,
            "true_velocity_um_s": v_true_z,
            "velocity_measurements_um_s": v_meas,
            "velocity_measurement_sigma_um_s": sigma_v,
            "pino_z_um": z_pino,
            "kalman_z_um": z_kf,
            "harmonic_kalman_z_um": z_hk,
            "dead_reckoning_z_um": z_dr,
            "last_hold_z_um": np.full_like(z_true, z_0),
        }


class DegaussingRetrogradeRetrieval:
    """
    Step 19 (Module B5): Executes post-recanalization exponential alternating-field
    demagnetization (degaussing) to break magnetic chains into single sub-RBC colloidal
    rollers, followed by retrograde catch-bond wall rolling to the parent microcatheter.

    Epistemic note (Audit Table Row 10 & Section 3 Cross-Parameter Consistency Check #1):
    Explicitly reconciles transient bifurcation washout (< 2.0% transient branch escape events)
    with final mass-balance retrieval (>= 99.0% of total deployed dose recovered) by tracking
    retrograde recapture of transiently escaped wall-rolling agents.
    """

    def __init__(self, cbf_config_path: Path | str | None = None):
        if cbf_config_path is None:
            root_dir = Path(__file__).resolve().parent.parent.parent
            cbf_config_path = root_dir / "configs" / "neurocritical_cbf_config.yaml"

        with open(cbf_config_path, "r", encoding="utf-8") as f:
            self.cfg = yaml.safe_load(f)

        ret_cfg = self.cfg["post_task_retrieval"]
        self.enable_degaussing = bool(ret_cfg["enable_degaussing_pulse"])
        self.tau_degauss_ms = float(ret_cfg["degaussing_decay_tau_ms"])
        self.v_retro_um_s = float(ret_cfg["retrograde_wall_velocity_um_s"])
        self.target_retrieval_min = float(ret_cfg["target_retrieval_efficiency_min"])

        # D5 (ELN #003): retrieval accounting constants moved out of the source, values unchanged.
        model = self.cfg["implementation_constants"]["post_task_retrieval_model"]
        self.branch_recapture_efficiency_max = float(model["branch_recapture_efficiency_max"])
        self.mainstem_loss_floor_pct = float(model["mainstem_plugging_loss_floor_pct"])
        self.mainstem_loss_per_cohesion_pct = float(model["mainstem_plugging_loss_per_cohesion_pct"])
        self.default_degaussing_duration_ms = float(model["default_degaussing_duration_ms"])

    def execute_post_task_retrieval(
        self,
        degaussing_duration_ms: float | None = None,
        retrieval_distance_um: float = 250.0,
        transient_washout_percent: float = 1.15,
        enable_degaussing_override: bool | None = None,
    ) -> RetrievalMassBalanceReport:
        """
        Computes residual dipole clumping after the degaussing pulse and evaluates the final
        mass-balance fraction of deployed microrobots retrieved at the parent catheter mouth.
        """
        use_degauss = self.enable_degaussing if enable_degaussing_override is None else bool(enable_degaussing_override)
        duration_ms = (
            self.default_degaussing_duration_ms
            if degaussing_duration_ms is None
            else float(degaussing_duration_ms)
        )

        if use_degauss:
            # Exponentially damped alternating pulse decays remanent dipole clumping: exp(-t / tau_d)
            residual_cohesion = float(np.exp(-duration_ms / max(self.tau_degauss_ms, 1.0)))
        else:
            # Without degaussing, dipole chains remain aggregated and risk capillary plugging
            residual_cohesion = 0.85

        # Time required for disaggregated microrollers to roll upstream along the no-slip wall
        effective_retro_v = self.v_retro_um_s * (1.0 - 0.55 * residual_cohesion)
        clearance_time_s = float(retrieval_distance_um) / max(effective_retro_v, 1.0)

        # Unified Mass-Balance Accounting (Audit Section 3 Check #1):
        # Transient washout agents (< 2.0%) that entered the proximal bifurcation collar
        # are actively recalled along the branch wall during retrograde rolling if disaggregated.
        transient_pct = float(np.clip(transient_washout_percent, 0.0, 100.0))
        branch_recapture_efficiency = self.branch_recapture_efficiency_max * (1.0 - residual_cohesion)
        recovered_from_branch_pct = transient_pct * branch_recapture_efficiency
        unrecovered_branch_pct = transient_pct - recovered_from_branch_pct

        # Main-stem capillary trapping loss scales with residual magnetic chain clumping
        mainstem_plugging_loss_pct = (
            self.mainstem_loss_floor_pct + self.mainstem_loss_per_cohesion_pct * residual_cohesion
        )

        total_unrecovered_pct = unrecovered_branch_pct + mainstem_plugging_loss_pct
        final_recovered_fraction = float(np.clip(1.0 - (total_unrecovered_pct / 100.0), 0.0, 1.0))
        unrecovered_fraction = float(1.0 - final_recovered_fraction)

        return RetrievalMassBalanceReport(
            degaussing_tau_ms=self.tau_degauss_ms,
            residual_dipole_cohesion_ratio=residual_cohesion,
            retrograde_velocity_um_s=float(effective_retro_v),
            clearance_time_s=clearance_time_s,
            transient_washout_events_pct=transient_pct,
            recovered_transient_pct=float(recovered_from_branch_pct),
            final_recovered_fraction=final_recovered_fraction,
            unrecovered_loss_fraction=unrecovered_fraction,
            passes_retrieval_target=bool(final_recovered_fraction >= self.target_retrieval_min),
        )


# D11 (ELN #003): the accurate name for what this class implements. The original name is kept so
# imports and the registered module layout do not break, but it is the alias that describes the
# method: a five-mode regularized Fourier (spectral) MAP estimator with a haemodynamic prior,
# not a trained neural operator.
PhysicsInformedSpectralObserver = PINOSubVoxelObserver
