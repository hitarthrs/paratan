"""Invariant and local coordinate mapping for confined ion reconstruction"""
from __future__ import annotations
import numpy as np
from source_model_revamp.numerical_quadrature import gauss_legendre_rule
from scipy.interpolate import RegularGridInterpolator
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid
from source_model_revamp.fbis.modal.types import ModalFBISBasis
from source_model_revamp.fbis.modal.utils import _safe_interp_rows
from source_model_revamp.scattering.physical_eigenbasis import PhysicalEigenbasis

def _modal_distribution_to_eta(modal_f_j_v: np.ndarray, physical_basis: PhysicalEigenbasis) -> np.ndarray:
    """Return f(v, η) = Σ_j f_j(v) I_j(η) with shape (n_speed, n_eta)"""
    return np.einsum("jv,je->ve", modal_f_j_v, physical_basis.eigenfunctions)

def _modal_distribution_to_lambda_grid(*, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, basis: ModalFBISBasis, distribution_v_eta: np.ndarray) -> np.ndarray:
    """Sample confined f(v, η) onto magnetic Λ cells and set passing cells to zero"""
    Lambda = np.asarray(lambda_grid.centers, dtype=float)
    result = np.zeros((speed_grid.centers_m_s.size, Lambda.size), dtype=float)
    trapped = Lambda >= basis.lambda_boundary
    if np.any(trapped):
        eta = basis.eta_lambda_map.eta_of_lambda(Lambda[trapped])
        sampled = _safe_interp_rows(basis.physical_basis.eta_grid, distribution_v_eta, eta, left=0.0, right=0.0)
        result[:, trapped] = sampled

    return np.maximum(result, 0.0)

def _density_from_v_eta(speed_grid: SpeedGrid, eta_grid: np.ndarray, distribution_v_eta: np.ndarray) -> float:
    """Return density from nonnegative f(v, η) using ∫ dη and exact 4π v² dv shell measures"""
    eta_integral = np.trapezoid(np.maximum(distribution_v_eta, 0.0), eta_grid, axis=1)

    return float(np.sum(speed_grid.shell_volumes_m3_s3 * eta_integral))

def _effective_temperature_from_v_eta(speed_grid: SpeedGrid, eta_grid: np.ndarray, distribution_v_eta: np.ndarray, mass_kg: float) -> float:
    """Return T_eff = (2/3)<E_kin> in J from nonnegative Egedal f(v, η)"""
    f = np.maximum(np.asarray(distribution_v_eta, dtype=float), 0.0)
    eta_integral = np.trapezoid(f, eta_grid, axis=1)
    density = float(np.sum(speed_grid.shell_volumes_m3_s3 * eta_integral))
    if density <= 0.0:
        return 0.0
    kinetic = 0.5 * float(mass_kg) * np.asarray(speed_grid.centers_m_s, dtype=float) ** 2
    energy_density = float(np.sum(speed_grid.shell_volumes_m3_s3 * eta_integral * kinetic))

    return float((2.0 / 3.0) * energy_density / density)

def _validated_base_distribution(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, base_distribution_v_lambda: np.ndarray) -> np.ndarray:
    """Validate f(v, Λ) shape and finite values and clip only roundoff scale negatives"""
    base = np.asarray(base_distribution_v_lambda, dtype=float)
    expected = (speed_grid.centers_m_s.size, lambda_grid.centers.size)
    if base.shape != expected:
        raise ValueError(f"base_distribution_v_lambda must have shape {expected}, got {base.shape}")
    if not np.all(np.isfinite(base)):
        raise ValueError("base_distribution_v_lambda must contain only finite values")
    scale = max(float(np.max(np.abs(base))) if base.size else 0.0, np.finfo(float).tiny)
    tolerance = 64.0 * np.finfo(float).eps * scale
    if float(np.min(base)) < -tolerance:
        raise ValueError("base_distribution_v_lambda contains a significant negative value")

    return np.where(base < 0.0, 0.0, base)

def _base_distribution_interpolator(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, base_distribution_v_lambda: np.ndarray) -> RegularGridInterpolator:
    """Build linear interpolation of the magnetic reference f(v, Λ) over its stored domain"""
    base = _validated_base_distribution(speed_grid, lambda_grid, base_distribution_v_lambda)
    speed_coordinates = np.concatenate(([float(speed_grid.faces_m_s[0])], np.asarray(speed_grid.centers_m_s, dtype=float), [float(speed_grid.faces_m_s[-1])]))
    lambda_coordinates = np.concatenate(([float(lambda_grid.faces[0])], np.asarray(lambda_grid.centers, dtype=float), [float(lambda_grid.faces[-1])]))
    # Extend cell centered values to the outer faces so interpolation spans the stored domain
    face_extended = np.pad(base, ((1, 1), (1, 1)), mode="edge")

    return RegularGridInterpolator((speed_coordinates, lambda_coordinates), face_extended, method="linear", bounds_error=False, fill_value=0.0)

def _eq71_global_trapped_passing_boundary(total_energy_J: np.ndarray, mirror_ratio: float, throat_potential_drop_magnitude_J: float) -> np.ndarray:
    """
    Return Egedal Eq 71 using q * ϕ_throat = -potential_drop_magnitude

    Λ_TP(U) = (1 / R_M) * (1 + ΔΦ_throat / U)
    """
    U = np.asarray(total_energy_J, dtype=float)
    R = float(mirror_ratio)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    throat_drop = float(throat_potential_drop_magnitude_J)
    if not np.isfinite(throat_drop) or throat_drop < 0.0:
        raise ValueError("throat_potential_drop_magnitude_J must be finite and nonnegative")
    boundary = np.full_like(U, np.inf, dtype=float)
    positive = U > 0.0
    boundary[positive] = (1.0 / R) * (1.0 + throat_drop / U[positive])

    return boundary

def _eq71_closed_interval_threshold_J(mirror_ratio: float, throat_potential_drop_magnitude_J: float) -> float:
    """Return U = ΔΦ_throat / (R_M − 1) where the Eq 71 confined interval closes"""
    R = float(mirror_ratio)
    throat_drop = float(throat_potential_drop_magnitude_J)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    if not np.isfinite(throat_drop) or throat_drop < 0.0:
        raise ValueError("throat_potential_drop_magnitude_J must be finite and nonnegative")

    return throat_drop / (R - 1.0)

def _eq72_compressed_invariant(global_lambda: np.ndarray, trapped_passing_boundary: np.ndarray, mirror_ratio: float) -> np.ndarray:
    """Return the Egedal Eq 72 compressed invariant Λ* without a physical domain clip"""
    Lambda = np.asarray(global_lambda, dtype=float)
    boundary = np.asarray(trapped_passing_boundary, dtype=float)
    denominator = 1.0 - boundary
    result = np.full(np.broadcast_shapes(Lambda.shape, boundary.shape), np.nan, dtype=float)
    valid = denominator > 0.0
    np.divide((1.0 - 1.0 / float(mirror_ratio)) * (1.0 - Lambda), denominator, out=result, where=valid)

    return 1.0 - result

def _ion_invariants_from_local_coordinates(*, local_speed_m_s: np.ndarray, local_pitch: np.ndarray, B_tilde: float, local_potential_drop_magnitude_J: float, particle_mass_kg: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Return invariant U and Λ from local speed and pitch ξ = v_parallel / v

    U = K_local − ΔΦ_local
    Λ = K_local * (1 − ξ²) / (B_tilde * U)
    """
    speed, pitch = np.broadcast_arrays(np.asarray(local_speed_m_s, dtype=float), np.asarray(local_pitch, dtype=float))
    B = float(B_tilde)
    local_drop = float(local_potential_drop_magnitude_J)
    if not np.isfinite(B) or B <= 0.0:
        raise ValueError("B_tilde must be positive and finite")
    if not np.isfinite(local_drop) or local_drop < 0.0:
        raise ValueError("local_potential_drop_magnitude_J must be finite and nonnegative")
    if np.any(~np.isfinite(speed)) or np.any(speed < 0.0):
        raise ValueError("local_speed_m_s must be finite and nonnegative")
    pitch_tolerance = 128.0 * np.finfo(float).eps
    if np.any(~np.isfinite(pitch)) or np.any(np.abs(pitch) > 1.0 + pitch_tolerance):
        raise ValueError("local_pitch must be finite and lie inside [-1, 1]")
    pitch = np.clip(pitch, -1.0, 1.0)
    local_kinetic = 0.5 * float(particle_mass_kg) * speed**2
    total_energy = local_kinetic - local_drop
    global_lambda = np.divide(local_kinetic * (1.0 - pitch**2), B * total_energy, out=np.full_like(total_energy, np.nan), where=total_energy > 0.0)

    return total_energy, global_lambda

def _cell_quadrature(faces: np.ndarray, order: int) -> tuple[np.ndarray, np.ndarray]:
    """Return Gauss Legendre nodes and coordinate weights for every grid cell"""
    nodes, weights = gauss_legendre_rule(int(order))
    lower = np.asarray(faces[:-1], dtype=float)
    upper = np.asarray(faces[1:], dtype=float)
    midpoint = 0.5 * (lower + upper)
    half_width = 0.5 * (upper - lower)

    return midpoint[:, None] + half_width[:, None] * nodes[None, :], half_width[:, None] * weights[None, :]
