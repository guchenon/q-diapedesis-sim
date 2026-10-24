"""First-principles derivations of model parameters from the configured physical quantities."""

from .derived_parameters import (
    DipolarCouplingResult,
    derive_dissipation_gamma,
    derive_swarm_coupling_K,
    dipolar_kuramoto_coupling,
    magnetic_moment_from_stepout,
    magnetic_relaxation_rate_inv_s,
    rotational_diffusion_rad2_s,
    rotational_drag_coefficient,
    sample_near_wall_shell,
)

__all__ = [
    "DipolarCouplingResult",
    "derive_dissipation_gamma",
    "derive_swarm_coupling_K",
    "dipolar_kuramoto_coupling",
    "magnetic_moment_from_stepout",
    "magnetic_relaxation_rate_inv_s",
    "rotational_diffusion_rad2_s",
    "rotational_drag_coefficient",
    "sample_near_wall_shell",
]
