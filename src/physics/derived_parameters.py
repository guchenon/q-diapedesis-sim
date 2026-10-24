"""
First-principles derivations of parameters that were previously free numbers in the YAML
configs (D7 / D8, ELN #003).

Everything here is derived from quantities already present in the configs - particle radius,
plasma viscosity, field amplitudes, the Mahoney step-out frequency, the swarm dose and the
vessel geometry - so no new free parameter is introduced.

Grounded in:
- Mahoney, Nelson, Peyer, Nelson & Abbott, Appl. Phys. Lett. 2014 (doi:10.1063/1.4870768)
  for the step-out balance between magnetic torque and rotational drag.
- Jackson, Classical Electrodynamics, 3rd ed., Sec. 5.7, for the dipole-dipole interaction.
- Einstein-Debye rotational diffusion of a sphere, D_r = kT / (8 pi eta a^3).
"""

from dataclasses import dataclass
from functools import lru_cache
import numpy as np

MU_0 = 4.0e-7 * np.pi          # vacuum permeability [T m / A]
BOLTZMANN_J_PER_K = 1.380649e-23


def rotational_drag_coefficient(radius_m: float, viscosity_Pa_s: float) -> float:
    """Rotational drag of a sphere in a Newtonian fluid: zeta_r = 8 pi eta a^3 [N m s]."""
    return 8.0 * np.pi * float(viscosity_Pa_s) * float(radius_m) ** 3


def rotational_diffusion_rad2_s(radius_m: float, viscosity_Pa_s: float, temperature_K: float) -> float:
    """Einstein-Debye rotational diffusion coefficient D_r = kT / zeta_r [rad^2/s]."""
    return BOLTZMANN_J_PER_K * float(temperature_K) / rotational_drag_coefficient(radius_m, viscosity_Pa_s)


def magnetic_moment_from_stepout(
    stepout_frequency_Hz: float, drive_field_T: float, radius_m: float, viscosity_Pa_s: float
) -> float:
    """
    Magnetic moment implied by the configured step-out frequency.

    At step-out the maximum magnetic torque equals the rotational drag torque
    (Mahoney et al. 2014):  m B = zeta_r * omega_so,  so  m = 2 pi f_so zeta_r / B.

    This is how the particle's moment is pinned without introducing a new free parameter:
    f_so is already in the config, and the moment follows from it.
    """
    zeta_r = rotational_drag_coefficient(radius_m, viscosity_Pa_s)
    return 2.0 * np.pi * float(stepout_frequency_Hz) * zeta_r / float(drive_field_T)


def magnetic_relaxation_rate_inv_s(
    moment_A_m2: float, field_T: float, radius_m: float, viscosity_Pa_s: float
) -> float:
    """
    Overdamped magnetic-viscous relaxation rate of a dipole toward the field: m B / zeta_r [1/s].

    In the low-Reynolds regime of a 0.75 um particle in plasma (Re ~ 1e-6) inertia is
    irrelevant, so this rate sets BOTH the coherent precession scale and the dissipative
    relaxation scale. That is the point of D8: the two cannot be independently chosen.
    """
    zeta_r = rotational_drag_coefficient(radius_m, viscosity_Pa_s)
    return float(moment_A_m2) * float(field_T) / zeta_r


def sample_near_wall_shell(
    n_particles: int,
    lumen_radius_m: float,
    radial_band: tuple[float, float],
    axial_span_m: float,
    particle_radius_m: float,
    seed: int = 0,
    max_tries_per_particle: int = 200,
    allow_partial: bool = False,
) -> np.ndarray:
    """
    Hard-sphere sample of `n_particles` centres inside the near-wall annular shell the swarm
    model seeds into: radii in `radial_band` (fractions of the lumen radius) over `axial_span_m`.

    Returns an (n, 3) array in metres with the vessel axis along z. Raises if the requested
    dose cannot be placed without overlap, which is itself a physical result worth reporting.
    """
    rng = np.random.default_rng(seed)
    r_lo = float(radial_band[0]) * float(lumen_radius_m)
    r_hi = float(radial_band[1]) * float(lumen_radius_m)
    min_sep = 2.0 * float(particle_radius_m)
    n_target = int(n_particles)

    # Batched hard-sphere (Poisson-disk) sampling with a cKDTree, so a dose that does not fit
    # fails quickly instead of grinding through rejections one candidate at a time.
    from scipy.spatial import cKDTree

    accepted = np.empty((0, 3), dtype=float)
    batch = max(256, n_target // 8)
    stalled = 0
    for _ in range(int(max_tries_per_particle)):
        rad = np.sqrt(rng.uniform(r_lo**2, r_hi**2, size=batch))
        phi = rng.uniform(0.0, 2.0 * np.pi, size=batch)
        cand = np.column_stack([rad * np.cos(phi), rad * np.sin(phi),
                                rng.uniform(0.0, float(axial_span_m), size=batch)])
        # Thin the candidate batch against itself, then against what is already accepted.
        keep = []
        tree_c = cKDTree(cand)
        blocked = set()
        for i, j in tree_c.query_pairs(min_sep):
            if i not in blocked:
                blocked.add(j)
        order = [i for i in range(batch) if i not in blocked]
        if accepted.shape[0]:
            tree_a = cKDTree(accepted)
            order = [i for i in order if not tree_a.query_ball_point(cand[i], min_sep)]
        keep = cand[order] if order else np.empty((0, 3))
        if keep.shape[0] == 0:
            stalled += 1
            if stalled >= 8:
                break
        else:
            stalled = 0
            accepted = np.vstack([accepted, keep])
        if accepted.shape[0] >= n_target:
            return accepted[:n_target]

    if allow_partial and accepted.shape[0] >= 20:
        # The configured dose exceeds what the geometry can hold. Returning what fits, with the
        # shortfall visible to the caller, is more useful than failing - but the shortfall is a
        # physical result and must be reported, not absorbed.
        return accepted
    raise ValueError(
        f"could not place {n_target} particles of radius {particle_radius_m:.3e} m without "
        f"overlap in the shell (placed {accepted.shape[0]}); the configured dose does not fit "
        f"this geometry"
    )


@dataclass
class DipolarCouplingResult:
    """Outcome of the first-principles dipolar Kuramoto coupling calculation (D7)."""
    coupling_K_rad_s: float            # signed: positive synchronizes, negative anti-aligns
    nearest_neighbour_distance_m: float
    volume_fraction: float
    mean_geometric_factor: float       # <(1.5 sin^2 psi - 1) / r^3> summed over neighbours [1/m^3]
    n_particles_placed: int
    axial_span_m: float


def dipolar_kuramoto_coupling(
    moment_A_m2: float,
    positions_m: np.ndarray,
    radius_m: float,
    viscosity_Pa_s: float,
    rotation_axis: np.ndarray | None = None,
    n_reference: int = 200,
    seed: int = 0,
) -> DipolarCouplingResult:
    """
    Kuramoto coupling produced by magnetic dipole-dipole torques between co-rotating rotors.

    Derivation. For two dipoles of magnitude m rotating in the plane normal to the drive axis
    with phases theta_1, theta_2, separated by r with polar angle psi to that axis, the
    interaction energy (Jackson Sec. 5.7) is

        U = (mu_0 m^2 / 4 pi r^3) [ cos(theta_1 - theta_2)
                                    - 3 sin^2(psi) cos(theta_1 - phi) cos(theta_2 - phi) ]

    and the torque on rotor 1 about the drive axis is -dU/dtheta_1. Expanding the product
    splits it into a secular part in sin(theta_1 - theta_2) and a counter-rotating part in
    sin(theta_1 + theta_2 - 2 phi), which averages away over a drive period. The secular
    coefficient is

        tau = (mu_0 m^2 / 4 pi r^3) (1.5 sin^2(psi) - 1) sin(theta_2 - theta_1),

    so dividing by the rotational drag and summing over neighbours gives

        K = (mu_0 m^2) / (4 pi zeta_r) * sum_j (1.5 sin^2(psi_ij) - 1) / r_ij^3 .

    The geometric factor is what makes this interesting and is why a nearest-neighbour
    estimate is not enough: it is +0.5 for a neighbour in the rotation plane, -1 for a
    neighbour along the drive axis, zero at sin^2(psi) = 2/3, and its ISOTROPIC AVERAGE IS
    EXACTLY ZERO (<sin^2 psi> = 2/3). A dipolar cloud therefore has no net phase coupling at
    all unless its geometry is anisotropic, and the sign depends on which way.
    """
    pts = np.asarray(positions_m, dtype=float)
    n = pts.shape[0]
    axis = np.array([0.0, 0.0, 1.0]) if rotation_axis is None else np.asarray(rotation_axis, float)
    axis = axis / np.linalg.norm(axis)

    rng = np.random.default_rng(seed)
    idx = rng.choice(n, size=min(int(n_reference), n), replace=False)

    factors = np.empty(idx.size, dtype=float)
    nn_dists = np.empty(idx.size, dtype=float)
    for k, i in enumerate(idx):
        delta = np.delete(pts, i, axis=0) - pts[i]
        dist = np.linalg.norm(delta, axis=1)
        cos_psi = (delta @ axis) / dist
        sin2_psi = 1.0 - cos_psi**2
        factors[k] = float(np.sum((1.5 * sin2_psi - 1.0) / dist**3))
        nn_dists[k] = float(np.min(dist))

    zeta_r = rotational_drag_coefficient(radius_m, viscosity_Pa_s)
    prefactor = MU_0 * float(moment_A_m2) ** 2 / (4.0 * np.pi * zeta_r)
    mean_factor = float(np.mean(factors))

    span = float(pts[:, 2].max() - pts[:, 2].min())
    radii = np.linalg.norm(pts[:, :2], axis=1)
    # Diagnostic only. A flat or single-layer arrangement has no shell volume, so the volume
    # fraction is reported as nan rather than dividing by zero.
    shell_vol = np.pi * (radii.max() ** 2 - radii.min() ** 2) * span
    phi_vol = (
        n * (4.0 / 3.0) * np.pi * float(radius_m) ** 3 / shell_vol if shell_vol > 0.0 else float("nan")
    )

    return DipolarCouplingResult(
        coupling_K_rad_s=prefactor * mean_factor,
        nearest_neighbour_distance_m=float(np.mean(nn_dists)),
        volume_fraction=float(phi_vol),
        mean_geometric_factor=mean_factor,
        n_particles_placed=n,
        axial_span_m=span,
    )


@lru_cache(maxsize=32)
def _cached_shell_coupling(
    moment_A_m2: float,
    lumen_radius_m: float,
    band_lo: float,
    band_hi: float,
    axial_span_m: float,
    radius_m: float,
    viscosity_Pa_s: float,
    n_particles: int,
    n_reference: int,
    seed: int,
) -> tuple[float, int]:
    """Signed K and the particle count actually placed. Memoized: the sample is costly."""
    pts = sample_near_wall_shell(
        n_particles, lumen_radius_m, (band_lo, band_hi), axial_span_m, radius_m, seed=seed,
        allow_partial=True,
    )
    k = dipolar_kuramoto_coupling(
        moment_A_m2, pts, radius_m, viscosity_Pa_s, n_reference=n_reference, seed=seed
    ).coupling_K_rad_s
    return float(k), int(pts.shape[0])


def derive_swarm_coupling_K(
    vessel_cfg: dict, floquet_cfg: dict, lumen_diameter_um: float, return_detail: bool = False
):
    """
    Signed Kuramoto coupling for the swarm as the model actually seeds it (D7).

    Every input comes from the existing configs. Returns rad/s; a NEGATIVE value means the
    dipolar torques anti-align the phases, which is what this geometry gives.

    With `return_detail`, also reports how many of the configured dose could actually be placed
    without overlap. At 15 um only ~38% of the preregistered 10 000-particle dose fits in the
    near-wall band, which is a physical limit of the registered geometry, not a solver failure.
    """
    bio = vessel_cfg["biomimicry_layer"]
    hemo = vessel_cfg["hemodynamics_pries_secomb"]
    drive = floquet_cfg["bimodal_floquet_drive"]
    deriv = vessel_cfg["implementation_constants"]["dipolar_coupling_derivation"]

    a = 1.0e-6 * float(bio["particle_radius_um"])
    eta = 1.0e-3 * float(hemo["plasma_viscosity_mPa_s"])
    b_alpha_T = 1.0e-3 * float(drive["alpha_drive_amplitude_mT"])
    f_so = float(floquet_cfg["implementation_constants"]["step_out"]["alpha_stepout_frequency_Hz"])
    moment = magnetic_moment_from_stepout(f_so, b_alpha_T, a, eta)

    band = deriv["radial_band_fraction_of_lumen"]
    n_dose = int(floquet_cfg["swarm_scaling"]["n_full_swarm_equivalent"])
    k, n_placed = _cached_shell_coupling(
        float(moment),
        0.5e-6 * float(lumen_diameter_um),
        float(band[0]),
        float(band[1]),
        1.0e-6 * float(deriv["axial_span_um"]),
        a,
        eta,
        n_dose,
        int(deriv["reference_particles"]),
        int(deriv["sampling_seed"]),
    )
    if not return_detail:
        return k
    return {
        "coupling_K_rad_s": k,
        "moment_A_m2": float(moment),
        "configured_dose": n_dose,
        "particles_placed": n_placed,
        "placeable_fraction": n_placed / float(n_dose),
    }


def derive_dissipation_gamma(vessel_cfg: dict, floquet_cfg: dict, reference_viscosity_mPa_s: float = 3.0) -> dict:
    """
    Dissipation rate for the GKLS generator, derived rather than chosen (D8).

    In the overdamped regime the SAME quantity m B / zeta_r sets both the coherent precession
    scale and the dissipative relaxation scale, so they cannot be picked independently.

    There is also a unit inconsistency in the existing Hamiltonian: `H_eff` is assembled from
    field amplitudes in mT and then used in exp(-i H t) with t in seconds, i.e. mT is being
    read as rad/s. The physical conversion is m / zeta_r per tesla, so
    `mT_to_rad_s = (m / zeta_r) * 1e-3`. Returning it lets callers state the physical scale
    without silently rescaling the preregistered Hamiltonian.
    """
    bio = vessel_cfg["biomimicry_layer"]
    hemo = vessel_cfg["hemodynamics_pries_secomb"]
    drive = floquet_cfg["bimodal_floquet_drive"]

    a = 1.0e-6 * float(bio["particle_radius_um"])
    eta = 1.0e-3 * float(hemo["plasma_viscosity_mPa_s"])
    b_alpha_T = 1.0e-3 * float(drive["alpha_drive_amplitude_mT"])
    f_so = float(floquet_cfg["implementation_constants"]["step_out"]["alpha_stepout_frequency_Hz"])

    zeta_r = rotational_drag_coefficient(a, eta)
    moment = magnetic_moment_from_stepout(f_so, b_alpha_T, a, eta)
    gamma_phys = magnetic_relaxation_rate_inv_s(moment, b_alpha_T, a, eta)
    mT_to_rad_s = (moment / zeta_r) * 1.0e-3

    # The code's gamma_eff is gamma_viscous * (mu_local / reference). Express the physical rate
    # in the same mT-based units the Hamiltonian uses, so the RATIO the dynamics see is right.
    gamma_in_code_units = (gamma_phys / mT_to_rad_s) * (float(reference_viscosity_mPa_s) / 3.8)
    return {
        "moment_A_m2": float(moment),
        "rotational_drag_N_m_s": float(zeta_r),
        "gamma_physical_rad_s": float(gamma_phys),
        "mT_to_rad_s": float(mT_to_rad_s),
        "viscous_dissipation_gamma": float(gamma_in_code_units),
    }
