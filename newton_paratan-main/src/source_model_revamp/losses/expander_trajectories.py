"""Classify adiabatic lost ion trajectories from a mirror throat to a selected material surface

The input path begins at the throat and ends at the selected surface
Invariant energy and global Lambda are combined with B_tilde(z) and P(z) = -q * ϕ(z)
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike

@dataclass(frozen=True)
class ExpanderTrajectoryResult:
    """Classification and material deposition state for one directed expander branch
    
    trajectory_status_v_lambda has shape (N_v, N_Lambda)
    parallel_kinetic_energy_J_v_lambda_z has shape (N_v, N_Lambda, N_z)
    wall_kinetic_energy_J_v has one total kinetic energy per invariant energy bin
    Particle fractions are weighted by directed_particle_rate_v_lambda_s
    """
    direction: str
    selected_surface_id: str
    selected_surface_particle_role: str
    trajectory_status_v_lambda: np.ndarray
    parallel_kinetic_energy_J_v_lambda_z: np.ndarray
    wall_kinetic_energy_J_v: np.ndarray
    wall_particle_rate_s: float
    wall_power_W: float
    reflected_particle_fraction: float
    locally_trapped_particle_fraction: float
    unclassified_particle_fraction: float
    deposited_particle_fraction: float
    maximum_wall_energy_residual_J: float
    finite_confinement_time_assigned: bool

def trace_adiabatic_lost_ions_to_surface(*, direction: str, invariant_total_energy_J: ArrayLike, global_lambda: ArrayLike, directed_particle_rate_v_lambda_s: ArrayLike, B_tilde_path: ArrayLike, potential_drop_magnitude_path_J: ArrayLike, selected_surface_id: str, selected_surface_particle_role: str) -> ExpanderTrajectoryResult:
    """Trace Egedal invariants through one geometry selected exhaust path
    
    Parallel energy is U + P(z) - Λ * U * B_tilde(z) for P = -q * ϕ
    
    invariant_total_energy_J has shape (N_v,)
    global_lambda has shape (N_Lambda,)
    directed_particle_rate_v_lambda_s has shape (N_v, N_Lambda)
    B_tilde_path and potential_drop_magnitude_path_J have shape (N_z,)
    
    The first path sample is the mirror throat and the last sample is the selected surface
    A connected path that reaches a nonreflecting selected surface is deposited
    A path blocked before the surface is reflected unless it disconnects and reconnects which is classified as local_trapped
    """
    branch = str(direction).strip().lower()
    if branch not in {"left", "right"}:
        raise ValueError("direction must be left or right")
    energy = np.asarray(invariant_total_energy_J, dtype=float)
    Lambda = np.asarray(global_lambda, dtype=float)
    rate = np.asarray(directed_particle_rate_v_lambda_s, dtype=float)
    B_path = np.asarray(B_tilde_path, dtype=float)
    P_path = np.asarray(potential_drop_magnitude_path_J, dtype=float)
    if energy.ndim != 1 or energy.size == 0 or np.any(~np.isfinite(energy)):
        raise ValueError("invariant_total_energy_J must be a nonempty finite vector")
    if Lambda.ndim != 1 or Lambda.size == 0 or np.any(~np.isfinite(Lambda)) or np.any(Lambda < 0.0):
        raise ValueError("global_lambda must be a nonempty finite nonnegative vector")
    if rate.shape != (energy.size, Lambda.size) or np.any(~np.isfinite(rate)) or np.any(rate < 0.0):
        raise ValueError("directed_particle_rate_v_lambda_s must be finite and match energy and Lambda")
    if B_path.ndim != 1 or B_path.size < 2 or np.any(~np.isfinite(B_path)) or np.any(B_path <= 0.0):
        raise ValueError("B_tilde_path must contain at least two positive finite values")
    if P_path.shape != B_path.shape or np.any(~np.isfinite(P_path)) or np.any(P_path < 0.0):
        raise ValueError("potential_drop_magnitude_path_J must be finite, nonnegative, and match the field path")
    role = str(selected_surface_particle_role).strip().lower()
    if role not in {"absorbing", "collector", "end_ring", "reflecting", "excluded"}:
        raise ValueError("selected_surface_particle_role is unsupported")
    if role == "excluded":
        raise ValueError("an excluded surface cannot be selected as a loss surface")

    # Evaluate parallel kinetic energy for every invariant energy Lambda and path sample
    U = energy[:, None, None]
    lam = Lambda[None, :, None]
    parallel = U + P_path[None, None, :] - lam * U * B_path[None, None, :]
    energy_scale = max(float(np.max(np.abs(parallel))), 1.0)
    tolerance = 256.0 * np.finfo(float).eps * energy_scale
    parallel = np.where(np.abs(parallel) <= tolerance, 0.0, parallel)
    status = np.full(rate.shape, "unclassified", dtype="U16")
    positive_energy = energy[:, None] > 0.0
    connected = parallel >= -tolerance
    reaches_throat = connected[:, :, 0]
    interior_blocked = np.any(~connected[:, :, 1:-1], axis=2)
    reaches_wall = connected[:, :, -1]
    # A disconnected interval followed by reconnection marks a local well topology
    disconnected_seen = np.maximum.accumulate(~connected, axis=2)
    reconnects_after_disconnected_interval = np.any(disconnected_seen & connected, axis=2)
    local_well = np.broadcast_to(~positive_energy, rate.shape) | (positive_energy & reaches_throat & reconnects_after_disconnected_interval)
    reflected = (positive_energy & reaches_throat & (interior_blocked | ~reaches_wall) & ~local_well)
    open_path = (positive_energy & reaches_throat & ~interior_blocked & reaches_wall & ~local_well)
    status[local_well] = "local_trapped"
    status[reflected] = "reflected"
    # Connected trajectories deposit unless the selected surface is reflecting
    if role == "reflecting":
        status[open_path] = "reflected"
    else:
        status[open_path] = "deposited"
    total_rate = float(np.sum(rate))

    def weighted_fraction(mask: np.ndarray) -> float:
        return float(np.sum(rate[mask]) / total_rate) if total_rate > 0.0 else 0.0

    deposited = status == "deposited"
    reflected_mask = status == "reflected"
    local_mask = status == "local_trapped"
    unresolved_mask = status == "unclassified"
    wall_energy = energy + float(P_path[-1])
    deposited_rate = float(np.sum(rate[deposited]))
    wall_power = float(np.sum(rate * deposited * wall_energy[:, None]))
    direct_wall_energy = energy[:, None] + float(P_path[-1])
    residual = direct_wall_energy - wall_energy[:, None]

    return ExpanderTrajectoryResult(
        direction=branch,
        selected_surface_id=str(selected_surface_id),
        selected_surface_particle_role=role,
        trajectory_status_v_lambda=status,
        parallel_kinetic_energy_J_v_lambda_z=parallel,
        wall_kinetic_energy_J_v=wall_energy,
        wall_particle_rate_s=deposited_rate,
        wall_power_W=wall_power,
        reflected_particle_fraction=weighted_fraction(reflected_mask),
        locally_trapped_particle_fraction=weighted_fraction(local_mask),
        unclassified_particle_fraction=weighted_fraction(unresolved_mask),
        deposited_particle_fraction=weighted_fraction(deposited),
        maximum_wall_energy_residual_J=float(np.max(np.abs(residual))),
        finite_confinement_time_assigned=False,
    )