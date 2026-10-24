"""
Stage 2 - Steps 13 & 14: Mycelial Hydro-Magnetic Soliton Phase Synchronization
(Sakaguchi-Kuramoto Model) and Autonomous Dual Sub-Swarm Partitioning (Alpha Drill vs. Beta Mesh).

Grounded in:
- Sakaguchi & Kuramoto, Prog. Theor. Phys. 1986 (doi:10.1143/PTP.76.576)
- Uchida & Golestanian, Phys. Rev. Lett. 2010 (doi:10.1103/PhysRevLett.104.178103)
- O'Keeffe, Hong & Strogatz, Nat. Commun. 2017 (doi:10.1038/s41467-017-01190-3)
- Ceron, Gardi, Petersen & Sitti, PNAS 2023 (doi:10.1073/pnas.2221913120)

Part of the Q-DIAPEDESIS Framework (OSF: https://osf.io/htswk | DOI: 10.5281/zenodo.23599020).
"""

from dataclasses import dataclass
from pathlib import Path
import numpy as np
import yaml

from src.physics import derive_swarm_coupling_K


@dataclass
class BiomimeticSwarmState:
    """State of the dual-sub-swarm biomimetic collective (Steps 13 & 14)."""
    time_s: float
    n_super_particles: int
    n_alpha_drill: int
    n_beta_mesh: int
    kuramoto_order_r: float
    alpha_cohesion_r: float
    beta_mesh_coverage_ratio: float
    transient_washout_percent: float
    passes_washout_ceiling: bool
    axial_positions_um: np.ndarray
    radial_positions_um: np.ndarray
    phases_rad: np.ndarray
    sub_swarm_labels: np.ndarray  # 0 = alpha (proximal drill), 1 = beta (distal mesh)


class HydroMagneticSolitonSwarm:
    """
    Implements Step 13 (Decentralized hydro-magnetic phase synchronization via the
    Sakaguchi-Kuramoto active rotator model) and Step 14 (Bimodal sub-swarm partitioning
    into proximal group alpha for staged clot drilling and distal group beta for embolic caging).

    Epistemic note (Audit Table Row 7):
    `transient_washout_percent` measures the percentage of agents transiently swept into
    collateral bifurcation branches per pass (< 2.0% design target), kept distinct from
    final post-task unrecovered mass balance.
    """

    def __init__(
        self,
        vessel_config_path: Path | str | None = None,
        floquet_config_path: Path | str | None = None,
        seed: int | None = None,
        lumen_diameter_um: float | None = None,
    ):
        root_dir = Path(__file__).resolve().parent.parent.parent
        if vessel_config_path is None:
            vessel_config_path = root_dir / "configs" / "vessel_mra_config.yaml"
        if floquet_config_path is None:
            floquet_config_path = root_dir / "configs" / "quantum_floquet_config.yaml"

        with open(vessel_config_path, "r", encoding="utf-8") as f:
            self.vessel_cfg = yaml.safe_load(f)
        with open(floquet_config_path, "r", encoding="utf-8") as f:
            self.floquet_cfg = yaml.safe_load(f)

        meta_seed = int(self.vessel_cfg["simulation_meta"]["random_seed"])
        self.rng = np.random.default_rng(meta_seed if seed is None else int(seed))

        geom = self.vessel_cfg["vessel_geometry"]
        diameter_um = (
            float(geom["default_diameter_um"]) if lumen_diameter_um is None else float(lumen_diameter_um)
        )
        self.vessel_length_um = float(geom["vessel_length_um"])
        # M1 (ELN #003): the closed loop can run at 15, 30 or 50 um, but this was always built
        # from the config default, so an overridden episode placed the swarm in the wrong lumen.
        self.lumen_diameter_um = diameter_um
        self.lumen_radius_um = 0.5 * diameter_um

        thrombus = self.vessel_cfg["thrombus_and_microenvironment"]
        self.thrombus_z_um = float(thrombus["thrombus_axial_position_um"])
        self.thrombus_len_um = float(thrombus["thrombus_length_um"])

        soliton = self.vessel_cfg["biomimicry_layer"]["hydro_magnetic_soliton"]
        # D7 (ELN #003): the coupling is derived from the magnetic dipole-dipole torque and the
        # rotational drag of the configured particle, in the geometry this model actually seeds.
        # The SIGN is physical and is kept: a negative K means the dipolar torques anti-align
        # the phases, which is what an axially extended near-wall band produces.
        self.coupling_is_derived = bool(soliton.get("derive_coupling_from_physics", False))
        if self.coupling_is_derived:
            detail = derive_swarm_coupling_K(
                self.vessel_cfg, self.floquet_cfg, diameter_um, return_detail=True
            )
            self.derived_coupling_K_rad_s = float(detail["coupling_K_rad_s"])
            self.derived_moment_A_m2 = float(detail["moment_A_m2"])
            self.placeable_dose_fraction = float(detail["placeable_fraction"])
            self.coupling_K = self.derived_coupling_K_rad_s
        else:
            self.derived_coupling_K_rad_s = float("nan")
            self.derived_moment_A_m2 = float("nan")
            self.placeable_dose_fraction = float("nan")
            self.coupling_K = float(soliton["kuramoto_coupling_K"])
        self.phase_lag_alpha = float(soliton["phase_lag_alpha_rad"])
        self.washout_ceiling_pct = float(soliton["washout_loss_ceiling_percent"])

        scaling = self.floquet_cfg["swarm_scaling"]
        self.n_super = int(scaling["n_super_particles"])
        self.alpha_ratio = float(scaling["sub_swarm_alpha_ratio"])
        self.beta_ratio = float(scaling["sub_swarm_beta_ratio"])

        drive = self.floquet_cfg["bimodal_floquet_drive"]
        self.omega_alpha_rad_s = 2.0 * np.pi * float(drive["alpha_drill_frequency_Hz"])
        self.omega_beta_rad_s = 2.0 * np.pi * float(drive["beta_trap_frequency_Hz"])

        # D3 Stage 2: frequency dispersion and rotational noise (implementation constants)
        disp = self.vessel_cfg["implementation_constants"]["kuramoto_dispersion"]
        self.natural_frequency_cv = float(disp["natural_frequency_cv"])
        self.temperature_K = float(disp["temperature_K"])
        self.particle_radius_um = float(soliton_layer_radius := self.vessel_cfg["biomimicry_layer"]["particle_radius_um"])
        plasma_mu_Pa_s = 1.0e-3 * float(
            self.vessel_cfg["hemodynamics_pries_secomb"]["plasma_viscosity_mPa_s"]
        )
        # Stokes-Einstein-Debye rotational diffusion of a sphere: D_r = kT / (8 pi eta a^3)
        BOLTZMANN_J_PER_K = 1.380649e-23
        radius_m = 1.0e-6 * float(soliton_layer_radius)
        self.rotational_diffusion_rad2_s = (BOLTZMANN_J_PER_K * self.temperature_K) / (
            8.0 * np.pi * plasma_mu_Pa_s * radius_m**3
        )

        self.reset()

    def reset(self) -> BiomimeticSwarmState:
        """Initializes super-particle positions, natural frequencies, and alpha/beta roles."""
        self.current_time_s = 0.0
        self.n_alpha = int(round(self.n_super * self.alpha_ratio))
        self.n_beta = self.n_super - self.n_alpha

        # Labels: 0 for proximal drilling sub-swarm alpha, 1 for distal trapping mesh beta
        self.labels = np.zeros(self.n_super, dtype=int)
        self.labels[self.n_alpha :] = 1

        # Initialize axial positions: alpha proximal to thrombus, beta deployed distally
        z_alpha = self.rng.normal(loc=self.thrombus_z_um - 45.0, scale=12.0, size=self.n_alpha)
        z_beta = self.rng.normal(loc=self.thrombus_z_um + 55.0, scale=10.0, size=self.n_beta)
        self.z_um = np.clip(np.concatenate([z_alpha, z_beta]), 10.0, self.vessel_length_um - 10.0)

        # Initialize radial positions near vessel wall (post-margination band)
        self.r_um = self.rng.uniform(0.75 * self.lumen_radius_um, 0.94 * self.lumen_radius_um, size=self.n_super)

        # D3 Stage 2: start incoherent over the whole circle. A narrow initial window made
        # the order parameter start near 0.95, so "synchronization" was the initial condition.
        self.phases_rad = self.rng.uniform(-np.pi, np.pi, size=self.n_super)

        # D3 Stage 2: per-rotor natural frequencies, dispersed about each sub-swarm's drive
        # frequency. Without this the oscillators are identical and cannot desynchronize.
        self.natural_omegas_rad_s = np.empty(self.n_super, dtype=float)
        for group_id, omega_0 in ((0, self.omega_alpha_rad_s), (1, self.omega_beta_rad_s)):
            mask = self.labels == group_id
            self.natural_omegas_rad_s[mask] = self.rng.normal(
                loc=omega_0,
                scale=self.natural_frequency_cv * abs(omega_0),
                size=int(np.count_nonzero(mask)),
            )
        return self.get_state(collateral_skimming_fraction=0.05, wall_locked=True)

    @staticmethod
    def _kuramoto_order_parameter(phases: np.ndarray) -> float:
        """Computes the Kuramoto phase coherence order parameter r in [0, 1]."""
        if phases.size == 0:
            return 0.0
        complex_mean = np.mean(np.exp(1j * phases))
        return float(np.abs(complex_mean))

    def step_soliton_swarm(
        self,
        dt_s: float,
        axial_rolling_velocity_um_s: float,
        wall_proximal_fraction: float,
        adhesion_locked: bool,
        collateral_skimming_fraction: float,
        enable_soliton_coupling: bool = True,
    ) -> BiomimeticSwarmState:
        """
        Advances the Sakaguchi-Kuramoto hydro-magnetic phase synchronization (Step 13)
        and updates the spatial positions of sub-swarms alpha and beta (Step 14).
        """
        K_eff = self.coupling_K if enable_soliton_coupling else 0.0

        # Sakaguchi-Kuramoto phase update within each frequency-decoupled sub-swarm.
        # Written in the rotating frame of each sub-swarm's drive frequency, so the detuning
        # of rotor i is its natural frequency minus that drive frequency (D3 Stage 2).
        dt = float(dt_s)
        noise_scale = np.sqrt(2.0 * self.rotational_diffusion_rad2_s * dt)
        for group_id, omega_0 in ((0, self.omega_alpha_rad_s), (1, self.omega_beta_rad_s)):
            mask = self.labels == group_id
            theta_g = self.phases_rad[mask]
            n_g = theta_g.size
            if n_g > 1 and K_eff != 0.0:
                # Pairwise phase differences: theta_j - theta_i - alpha
                diff_matrix = theta_g[np.newaxis, :] - theta_g[:, np.newaxis] - self.phase_lag_alpha
                coupling_term = (K_eff / float(n_g)) * np.sum(np.sin(diff_matrix), axis=1)
            else:
                coupling_term = np.zeros_like(theta_g)

            detuning = self.natural_omegas_rad_s[mask] - omega_0
            rotational_noise = noise_scale * self.rng.standard_normal(n_g)
            self.phases_rad[mask] = (
                theta_g + (detuning + coupling_term) * dt + rotational_noise
            ) % (2.0 * np.pi)

        # Step 14: Axial dynamics for alpha (drilling toward thrombus center) and beta (stabilizing mesh)
        alpha_mask = self.labels == 0
        beta_mask = self.labels == 1

        # Alpha sub-swarm advances toward thrombus plug face during rolling
        target_alpha_z = self.thrombus_z_um - 5.0
        dir_alpha = np.sign(target_alpha_z - self.z_um[alpha_mask])
        self.z_um[alpha_mask] += dir_alpha * float(axial_rolling_velocity_um_s) * 0.15 * float(dt_s)

        # Beta sub-swarm maintains a synchronized distal filtration mesh across the lumen
        target_beta_z = self.thrombus_z_um + 0.65 * self.thrombus_len_um
        restore_beta = 0.8 * (target_beta_z - self.z_um[beta_mask])
        self.z_um[beta_mask] += restore_beta * float(dt_s)

        self.z_um = np.clip(self.z_um, 5.0, self.vessel_length_um - 5.0)
        self.current_time_s += float(dt_s)

        return self.get_state(
            collateral_skimming_fraction=collateral_skimming_fraction,
            wall_locked=adhesion_locked,
            wall_proximal_fraction=wall_proximal_fraction,
        )

    def get_state(
        self,
        collateral_skimming_fraction: float,
        wall_locked: bool,
        wall_proximal_fraction: float = 0.88,
    ) -> BiomimeticSwarmState:
        """Computes swarm order parameters, distal mesh coverage, and transient washout."""
        alpha_phases = self.phases_rad[self.labels == 0]
        beta_phases = self.phases_rad[self.labels == 1]

        r_alpha = self._kuramoto_order_parameter(alpha_phases)
        r_beta = self._kuramoto_order_parameter(beta_phases)
        r_global = 0.6 * r_alpha + 0.4 * r_beta

        # Distal embolic mesh coverage depends on beta phase synchronization
        beta_mesh_coverage = float(np.clip(0.55 + 0.43 * r_beta, 0.0, 0.995))

        # Transient bifurcation washout percentage: suppressed by wall margination,
        # catch-bond systolic locking, and soliton phase cohesion (Target < 2.0%)
        retention_shield = (
            0.45 * float(wall_proximal_fraction)
            + 0.35 * r_global
            + (0.15 if wall_locked else 0.05)
        )
        raw_washout_pct = 100.0 * float(collateral_skimming_fraction) * max(0.02, 1.0 - retention_shield)
        transient_washout_pct = float(np.clip(raw_washout_pct, 0.05, 45.0))

        return BiomimeticSwarmState(
            time_s=self.current_time_s,
            n_super_particles=self.n_super,
            n_alpha_drill=self.n_alpha,
            n_beta_mesh=self.n_beta,
            kuramoto_order_r=r_global,
            alpha_cohesion_r=r_alpha,
            beta_mesh_coverage_ratio=beta_mesh_coverage,
            transient_washout_percent=transient_washout_pct,
            passes_washout_ceiling=bool(transient_washout_pct < self.washout_ceiling_pct),
            axial_positions_um=self.z_um.copy(),
            radial_positions_um=self.r_um.copy(),
            phases_rad=self.phases_rad.copy(),
            sub_swarm_labels=self.labels.copy(),
        )
