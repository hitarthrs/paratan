"""Orbit averaged Egedal modal basis construction

File contains the magnetic/invariant space part of FBIS
    Λ -> η(Λ), square mirror basis and physical orbit averaged Lorentz eigenbasis

This module builds the magnetic invariant mapping from Λ to η, projects the orbit averaged Lorentz operator onto the square mirror basis, and evaluates grid and retained mode convergence
Eq 42 density weighting enters the orbit averaged scattering coefficients while the η mapping remains magnetic
"""
from __future__ import annotations
from collections.abc import Callable
from dataclasses import dataclass, replace
from time import perf_counter
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.numerical_quadrature import gauss_legendre_rule
from scipy.integrate import trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import linear_sum_assignment
from source_model_revamp.fbis.modal.loss_geometry import loss_geometry_factor_from_midpoint_profile
from source_model_revamp.orbits.bounce_averages import SEPARATRIX_BOUNCE_TIME_MODEL, TURNING_POINT_QUADRATURE_MODEL, normalized_bounce_time
from source_model_revamp.orbits.eta_mapping import EtaLambdaMap, build_eta_lambda_map, eta_trapped_passing_boundary
from source_model_revamp.scattering.lorentz_operator import orbit_averaged_lorentz_coefficients_from_bounce
from source_model_revamp.scattering.physical_eigenbasis import LORENTZ_GALERKIN_FORM_CONSERVATIVE_WEAK, build_physical_eigenbasis
from source_model_revamp.scattering.square_mirror_basis import build_square_mirror_basis
from source_model_revamp.fbis.modal.types import ModalFBISBasis, ModalFBISNumerics
from source_model_revamp.fbis.modal.density_weighting import EQ42_UNIFORM_DENSITY_REFERENCE, Eq42DensityProfile
from source_model_revamp.fbis.modal.utils import _positive_int
from source_model_revamp.fbis.velocity_space_grid import _distribution_nonnegativity_metrics

@dataclass(frozen=True)
class ModalBasisConvergenceAssessment:
    """Successive grid convergence state for the physical modal basis
    
    Histories follow the supplied Λ and η refinement sequence and retain the built basis at every resolution
    """
    assessed: bool
    converged: bool
    grid_sequence: tuple[tuple[int, int], ...]
    eigenvalue_history: tuple[np.ndarray, ...]
    eigenvalue_relative_change_history: tuple[np.ndarray, ...]
    first_eigenvalue_relative_changes: np.ndarray
    eigenfunction_overlap_history: tuple[np.ndarray, ...]
    eigenfunction_mode_permutation_history: tuple[np.ndarray, ...]
    eigenfunction_sign_history: tuple[np.ndarray, ...]
    eigenpair_residual_history: tuple[np.ndarray, ...]
    boundary_residual_history: tuple[np.ndarray, ...]
    normalization_residual_history: tuple[np.ndarray, ...]
    orthogonality_residual_history: tuple[np.ndarray, ...]
    gram_condition_number_history: np.ndarray
    refinement_runtime_s: np.ndarray
    max_eigenpair_residual: float
    gram_condition_number: float
    eta_monotonic: bool
    galerkin_symmetry_relative_error: float
    gram_symmetry_relative_error: float
    max_eigenfunction_offdiagonal_overlap: float
    normalization_max_abs_error: float
    boundary_max_abs_error: float
    turning_point_quadrature_model: str
    separatrix_bounce_time_model: str
    failure_reason: str
    basis_history: tuple[ModalFBISBasis, ...]

    @property
    def finest_basis(self) -> ModalFBISBasis:
        """Return the highest resolution basis without rebuilding it"""
        return self.basis_history[-1]

@dataclass(frozen=True)
class EigenmodeMatch:
    """Matched fine grid modes on a common physical η interval
    
    Mode indices follow the coarse mode order and signs align the matched eigenfunctions
    """
    fine_mode_indices: np.ndarray
    signs: np.ndarray
    weighted_overlaps: np.ndarray
    common_eta_grid: np.ndarray

@dataclass(frozen=True)
class RetainedModeReconstructionAssessment:
    """Convergence trend for physical reconstructions from retained mode prefixes"""
    assessed: bool
    converged: bool | None
    mode_count_sequence: np.ndarray
    weighted_relative_change_history: np.ndarray
    maximum_pointwise_relative_change_history: np.ndarray
    inventory_history: np.ndarray
    inventory_relative_change_history: np.ndarray
    negative_reconstruction_fraction_history: np.ndarray
    relative_tolerance: float | None
    status: str

def assess_retained_mode_reconstruction_convergence(*, basis: ModalFBISBasis, modal_distribution_f_j_v: ArrayLike, speed_cell_weights: ArrayLike | None = None, mode_count_sequence: tuple[int, ...] | None = None, relative_tolerance: float | None = None) -> RetainedModeReconstructionAssessment:
    """Compare nested retained mode reconstructions on the physical η grid
    
    modal_distribution_f_j_v has shape (n_mode, n_speed)
    Optional speed weights define the speed integration measure used for the relative change and inventory diagnostics
    """
    modal = np.asarray(modal_distribution_f_j_v, dtype=float)
    physical = basis.physical_basis
    mode_total = physical.eigenvalues.size
    if modal.ndim != 2 or modal.shape[0] != mode_total:
        raise ValueError("modal_distribution_f_j_v must have shape (n_mode, n_speed)")
    if np.any(~np.isfinite(modal)):
        raise ValueError("modal_distribution_f_j_v must contain only finite values")
    if speed_cell_weights is None:
        speed_weights = np.ones(modal.shape[1], dtype=float)
    else:
        speed_weights = np.asarray(speed_cell_weights, dtype=float)
        if speed_weights.shape != (modal.shape[1],):
            raise ValueError("speed_cell_weights must contain one value per speed cell")
        if np.any(~np.isfinite(speed_weights)) or np.any(speed_weights <= 0.0):
            raise ValueError("speed_cell_weights must be finite and positive")
    if mode_count_sequence is None:
        counts = tuple(range(1, mode_total + 1))
    else:
        counts = tuple(_positive_int("retained mode count", value) for value in mode_count_sequence)
    if len(counts) < 2:
        raise ValueError("mode_count_sequence must contain at least two levels")
    if counts[-1] != mode_total or any(current <= previous for previous, current in zip(counts[:-1], counts[1:], strict=True)):
        raise ValueError("mode_count_sequence must increase strictly and end at the available mode count")
    if counts[0] < 1 or counts[-1] > mode_total:
        raise ValueError("retained mode counts must lie inside the available mode range")
    tolerance: float | None
    if relative_tolerance is None:
        tolerance = None
    else:
        tolerance = float(relative_tolerance)
        if not np.isfinite(tolerance) or tolerance <= 0.0:
            raise ValueError("relative_tolerance must be finite and positive")
    eta = np.asarray(physical.eta_grid, dtype=float)
    reconstructions: list[np.ndarray] = []
    inventories: list[float] = []
    negative_fractions: list[float] = []
    for count in counts:
        reconstructed = np.einsum("jv,je->ve", modal[:count], physical.eigenfunctions[:count], optimize=True)
        reconstructions.append(reconstructed)
        eta_integral = trapezoid(reconstructed, eta, axis=1)
        inventories.append(float(np.sum(speed_weights * eta_integral)))
        _, _, _, negative_fraction, _ = _distribution_nonnegativity_metrics(reconstructed, name="retained mode reconstruction")
        negative_fractions.append(negative_fraction)
    weighted_changes: list[float] = []
    pointwise_changes: list[float] = []
    inventory_changes: list[float] = []
    for previous, current, previous_inventory, current_inventory in zip(reconstructions[:-1], reconstructions[1:], inventories[:-1], inventories[1:], strict=True):
        difference = current - previous
        difference_norm = np.sqrt(float(np.sum(speed_weights * trapezoid(difference**2, eta, axis=1))))
        current_norm = np.sqrt(float(np.sum(speed_weights * trapezoid(current**2, eta, axis=1))))
        weighted_changes.append(difference_norm / max(current_norm, np.finfo(float).tiny))
        pointwise_changes.append(float(np.max(np.abs(difference))) / max(float(np.max(np.abs(current))), np.finfo(float).tiny))
        inventory_changes.append(abs(current_inventory - previous_inventory) / max(abs(current_inventory), abs(previous_inventory), np.finfo(float).tiny))
    if tolerance is None:
        converged = None
        status = "diagnostic_only_no_authoritative_retained_mode_tolerance"
    else:
        converged = bool(weighted_changes[-1] <= tolerance and inventory_changes[-1] <= tolerance)
        status = "converged" if converged else "not_converged"

    return RetainedModeReconstructionAssessment(
        assessed=True,
        converged=converged,
        mode_count_sequence=np.asarray(counts, dtype=int),
        weighted_relative_change_history=np.asarray(weighted_changes, dtype=float),
        maximum_pointwise_relative_change_history=np.asarray(pointwise_changes, dtype=float),
        inventory_history=np.asarray(inventories, dtype=float),
        inventory_relative_change_history=np.asarray(inventory_changes, dtype=float),
        negative_reconstruction_fraction_history=np.asarray(negative_fractions, dtype=float),
        relative_tolerance=tolerance,
        status=status,
    )

def clustered_lambda_grid(n_points: int, mirror_ratio: float, loss_boundary_interval_fraction: float = 0.35) -> np.ndarray:
    """Build a full Lambda grid clustered at 0, 1/R_M, and 1

    Cosine clustering increases resolution near zero, the magnetic loss boundary, and one
    """
    count = _positive_int("n_points", n_points, minimum=9)
    R = float(mirror_ratio)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and > 1")
    fraction = float(loss_boundary_interval_fraction)
    if not np.isfinite(fraction) or not 0.1 <= fraction <= 0.9:
        raise ValueError("loss_boundary_interval_fraction must satisfy 0.1 <= value <= 0.9")
    total_intervals = count - 1
    left_intervals = int(round(fraction * total_intervals))
    left_intervals = min(max(left_intervals, 4), total_intervals - 4)
    right_intervals = total_intervals - left_intervals
    lambda_boundary = 1.0 / R
    left_angle = np.linspace(0.0, np.pi, left_intervals + 1)
    right_angle = np.linspace(0.0, np.pi, right_intervals + 1)
    left = 0.5 * lambda_boundary * (1.0 - np.cos(left_angle))
    right = lambda_boundary + 0.5 * (1.0 - lambda_boundary) * (1.0 - np.cos(right_angle))
    grid = np.concatenate((left, right[1:]))
    grid[0] = 0.0
    grid[left_intervals] = lambda_boundary
    grid[-1] = 1.0

    return grid

def _next_refinement_count(count: int) -> int:
    """Return an odd intermediate grid count above the requested resolution"""
    value = int(count) + max(4, (int(count) - 1) // 3 + 2)
    if value % 2 == 0:
        value += 1

    return value

def _default_grid_sequence(n_lambda_grid: int, n_eta_grid: int) -> tuple[tuple[int, int], ...]:
    """Return requested, intermediate, and fine Λ and η resolutions"""
    lambda_count = _positive_int("n_lambda_grid", n_lambda_grid, minimum=9)
    eta_count = _positive_int("n_eta_grid", n_eta_grid, minimum=9)
    requested = (lambda_count, eta_count)
    intermediate = (_next_refinement_count(lambda_count), _next_refinement_count(eta_count),)
    fine = (2 * lambda_count - 1, 2 * eta_count - 1)
    if fine[0] <= intermediate[0] or fine[1] <= intermediate[1]:
        fine = (2 * intermediate[0] - 1, 2 * intermediate[1] - 1)

    return requested, intermediate, fine

def _validate_grid_sequence(grid_sequence: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    """Validate at least three strictly increasing Λ and η grid resolutions"""
    sequence = tuple((_positive_int("n_lambda_grid", resolution[0], minimum=9), _positive_int("n_eta_grid", resolution[1], minimum=9),) for resolution in grid_sequence)
    if len(sequence) < 3:
        raise ValueError("grid_sequence must contain at least three resolutions")
    for previous, current in zip(sequence[:-1], sequence[1:], strict=True):
        if current[0] <= previous[0] or current[1] <= previous[1]:
            raise ValueError("grid_sequence resolutions must increase in both Lambda and eta")
        
    return sequence

def match_eigenmodes_on_common_eta(coarse_basis: ModalFBISBasis, fine_basis: ModalFBISBasis, comparison_points: int = 1001) -> EigenmodeMatch:
    """Match coarse and fine eigenmodes by normalized η overlap and align their signs"""
    point_count = _positive_int("comparison_points", comparison_points, minimum=33)
    coarse = coarse_basis.physical_basis
    fine = fine_basis.physical_basis
    mode_count = min(coarse.eigenvalues.size, fine.eigenvalues.size)
    common_upper = min(float(coarse.eta_trapped_passing), float(fine.eta_trapped_passing))
    if not np.isfinite(common_upper) or common_upper <= 0.0:
        raise ValueError("physical eta comparison interval must be finite and positive")
    common_eta = np.linspace(0.0, common_upper, point_count)
    coarse_values = np.asarray([PchipInterpolator(coarse.eta_grid, row, extrapolate=False)(common_eta) for row in coarse.eigenfunctions[:mode_count]], dtype=float)
    fine_values = np.asarray([PchipInterpolator(fine.eta_grid, row, extrapolate=False)(common_eta) for row in fine.eigenfunctions[:mode_count]], dtype=float)
    coarse_norms = np.sqrt(trapezoid(coarse_values**2, common_eta, axis=1))
    fine_norms = np.sqrt(trapezoid(fine_values**2, common_eta, axis=1))
    cross = trapezoid(coarse_values[:, None, :] * fine_values[None, :, :], common_eta, axis=2)
    normalized_cross = cross / np.maximum(coarse_norms[:, None] * fine_norms[None, :], np.finfo(float).tiny)
    # Match by overlap because discrete mode order can change between grid refinements
    coarse_indices, fine_indices = linear_sum_assignment(-np.abs(normalized_cross))
    order = np.argsort(coarse_indices)
    permutation = np.asarray(fine_indices[order], dtype=int)
    signed_overlaps = normalized_cross[np.arange(mode_count), permutation]
    signs = np.where(signed_overlaps < 0.0, -1.0, 1.0)
    overlaps = np.clip(np.abs(signed_overlaps), 0.0, 1.0)

    return EigenmodeMatch(fine_mode_indices=permutation, signs=np.asarray(signs, dtype=float), weighted_overlaps=np.asarray(overlaps, dtype=float), common_eta_grid=common_eta)

def sign_aligned_eigenfunction_overlaps(coarse_basis: ModalFBISBasis, fine_basis: ModalFBISBasis, comparison_points: int = 1001) -> np.ndarray:
    """Return absolute overlaps for matched modes on a common physical η interval"""
    return match_eigenmodes_on_common_eta(coarse_basis, fine_basis, comparison_points=comparison_points).weighted_overlaps

def _eta_map_with_cell_averaged_separatrix(*, lambda_grid: np.ndarray, B_tilde_function: Callable[[float], float], mirror_ratio: float) -> EtaLambdaMap:
    """Build the magnetic η mapping with a finite cell averaged separatrix bounce value
    
    The exact separatrix point is not sampled directly because its normalized bounce time is singular
    A finite boundary value is chosen so trapezoid integration across the adjacent cell reproduces the quadrature integral
    """
    Lambda = np.asarray(lambda_grid, dtype=float)
    boundary = 1.0 / float(mirror_ratio)
    boundary_matches = np.flatnonzero(Lambda == boundary)
    if boundary_matches.size != 1:
        raise ValueError("lambda_grid must contain the exact magnetic loss boundary once")
    boundary_index = int(boundary_matches[0])
    if boundary_index == 0 or boundary_index == Lambda.size - 1:
        raise ValueError("magnetic loss boundary must be an interior Lambda node")
    tau = np.empty_like(Lambda)
    for index, Lambda_value in enumerate(Lambda):
        if index == boundary_index:
            continue
        tau[index] = normalized_bounce_time(Lambda=float(Lambda_value), B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, zeta_throat=1.0)
    # Replace the singular separatrix point by a cell average that preserves its local integral
    outer_nodes, outer_weights = gauss_legendre_rule(16)
    left = float(Lambda[boundary_index - 1])
    right = float(Lambda[boundary_index + 1])
    cell_integral = 0.0
    for interval_left, interval_right in ((left, boundary), (boundary, right)):
        mapped_lambda = 0.5 * (interval_right - interval_left) * outer_nodes + 0.5 * (interval_left + interval_right)
        mapped_tau = np.asarray([normalized_bounce_time(Lambda=float(Lambda_value), B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, zeta_throat=1.0) for Lambda_value in mapped_lambda], dtype=float)
        cell_integral += 0.5 * (interval_right - interval_left) * float(np.dot(outer_weights, mapped_tau))
    left_width = boundary - left
    right_width = right - boundary
    effective_boundary_tau = (2.0 * cell_integral - left_width * tau[boundary_index - 1] - right_width * tau[boundary_index + 1]) / (left_width + right_width)
    if not np.isfinite(effective_boundary_tau) or effective_boundary_tau <= 0.0:
        raise ValueError("failed to construct a finite positive cell averaged separatrix bounce time")
    tau[boundary_index] = effective_boundary_tau

    return build_eta_lambda_map(lambda_grid=Lambda, tau_tilde_grid=tau)

def build_modal_fbis_basis(*, mirror_ratio: float, B_tilde_function: Callable[[float], float], zeta_faces: ArrayLike, B_tilde_midpoints: ArrayLike, numerics: ModalFBISNumerics, density_ratio_function: Callable[[float], float] | None = None, eq42_density_profile: Eq42DensityProfile | None = None) -> ModalFBISBasis:
    """Build the Eq 47 η mapping and Eq 50 physical eigenbasis
    
    The magnetic mapping spans Λ from zero to one and the physical eigenbasis spans η from zero to the trapped passing boundary
    Eq 42 density weighting modifies the orbit averaged Lorentz coefficients and the coupled loss geometry when supplied
    """
    R = float(mirror_ratio)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and > 1")
    n_lambda = _positive_int("n_lambda_grid", numerics.n_lambda_grid, minimum=9)
    n_eta = _positive_int("n_eta_grid", numerics.n_eta_grid, minimum=9)
    n_square = _positive_int("n_square_basis_modes", numerics.n_square_basis_modes, minimum=1)
    n_modes = _positive_int("n_physical_modes", numerics.n_physical_modes, minimum=1)
    if n_square < n_modes:
        raise ValueError("n_square_basis_modes must be >= n_physical_modes")
    if density_ratio_function is not None and eq42_density_profile is not None:
        raise ValueError("Specify either density_ratio_function or eq42_density_profile, not both")
    active_density_ratio_function = density_ratio_function
    if eq42_density_profile is not None:
        active_density_ratio_function = eq42_density_profile.density_ratio_function
    full_lambda = clustered_lambda_grid(n_points=n_lambda, mirror_ratio=R)
    # The η mapping is magnetic while Eq 42 density enters the orbit averaged scattering coefficients
    eta_map = _eta_map_with_cell_averaged_separatrix(lambda_grid=full_lambda, B_tilde_function=B_tilde_function, mirror_ratio=R)
    eta_tp = float(eta_trapped_passing_boundary(eta_map, R))
    eta_grid = np.linspace(0.0, eta_tp, n_eta)
    square = build_square_mirror_basis(mirror_ratio=R, n_modes=n_square, scan_points=int(numerics.square_scan_points))
    lambda_on_eta = eta_map.lambda_of_eta(eta_grid)
    lorentz = orbit_averaged_lorentz_coefficients_from_bounce(lambda_values=lambda_on_eta, B_tilde_function=B_tilde_function, mirror_ratio=R, density_ratio_function=active_density_ratio_function, zeta_throat=1.0)
    physical = build_physical_eigenbasis(eta_grid=eta_grid, eta_trapped_passing=eta_tp, eta_lambda_map=eta_map, lorentz_coefficients=lorentz, square_basis=square, n_modes=n_modes, galerkin_form=LORENTZ_GALERKIN_FORM_CONSERVATIVE_WEAK)
    lambda_boundary = 1.0 / R
    d_eta_dlambda_boundary = float(eta_map.d_eta_d_lambda(lambda_boundary))
    slopes = np.asarray(physical.eigenfunction_derivatives_eta[:, -1], dtype=float) * d_eta_dlambda_boundary
    zeta_face_values = np.asarray(zeta_faces, dtype=float)
    zeta_midpoints = 0.5 * (zeta_face_values[:-1] + zeta_face_values[1:])
    if active_density_ratio_function is None:
        loss_density_ratio_midpoints = np.ones_like(zeta_midpoints)
    else:
        loss_density_ratio_midpoints = np.asarray([active_density_ratio_function(float(zeta_value)) for zeta_value in zeta_midpoints],dtype=float)
    geom_G, volume_geom = loss_geometry_factor_from_midpoint_profile(zeta_faces=zeta_face_values, B_tilde_midpoints=B_tilde_midpoints, mirror_ratio=R, density_ratio_midpoints=loss_density_ratio_midpoints)

    return ModalFBISBasis(
        eta_lambda_map=eta_map,
        square_basis=square,
        physical_basis=physical,
        eta_boundary=eta_tp,
        lambda_boundary=lambda_boundary,
        mode_slopes_dI_dlambda_at_boundary=slopes,
        geometry_factor_G=geom_G,
        volume_geometry_integral=volume_geom,
        eq61_geometry_density_weighted=bool(active_density_ratio_function is not None),
        eq61_geometry_weighting_model=("Egedal_2022_equation_62_with_equation_42_density_ratio" if active_density_ratio_function is not None else "uniform_density_reference"),
        eq61_geometry_density_ratio_minimum=float(np.min(loss_density_ratio_midpoints)),
        eq61_geometry_density_ratio_maximum=float(np.max(loss_density_ratio_midpoints)),
        eq42_density_weighting_model=(eq42_density_profile.model if eq42_density_profile is not None else (EQ42_UNIFORM_DENSITY_REFERENCE if active_density_ratio_function is None else "supplied_nonuniform_density_ratio_function")),
        eq42_density_weighting_coupled=(eq42_density_profile is not None or active_density_ratio_function is not None),
        eq42_nonuniform_density_weighting_active=(eq42_density_profile.nonuniform_weighting_active if eq42_density_profile is not None else active_density_ratio_function is not None),
        eq42_density_profile_source=(eq42_density_profile.source if eq42_density_profile is not None else ("uniform_density_reference" if active_density_ratio_function is None else "supplied_density_ratio_callable")),
        eq42_volume_average_density_m3=(None if eq42_density_profile is None else eq42_density_profile.volume_average_density_m3),
        eq42_density_ratio_zeta=(None if eq42_density_profile is None else np.asarray(eq42_density_profile.zeta_positive_nodes, dtype=float)),
        eq42_density_ratio_values=(None if eq42_density_profile is None else np.asarray(eq42_density_profile.density_ratio_positive_nodes, dtype=float)),
        eq42_density_profile_identity=(None if eq42_density_profile is None else eq42_density_profile.profile_identity),
        eq42_density_profile_symmetry_relative_error=(0.0 if eq42_density_profile is None else eq42_density_profile.symmetry_relative_error),
        eq42_density_profile_symmetry_tolerance=(numerics.eq42_density_symmetry_tolerance if eq42_density_profile is None else eq42_density_profile.symmetry_tolerance),
        eq42_density_profile_symmetric=(True if eq42_density_profile is None else eq42_density_profile.symmetric_within_tolerance),
    )

def assess_modal_fbis_basis_convergence(*, mirror_ratio: float, B_tilde_function: Callable[[float], float], zeta_faces: ArrayLike, B_tilde_midpoints: ArrayLike, numerics: ModalFBISNumerics, density_ratio_function: Callable[[float], float] | None = None, eq42_density_profile: Eq42DensityProfile | None = None, grid_sequence: tuple[tuple[int, int], ...] | None = None, eigenvalue_relative_tolerance: float = 1.0e-2, eigenfunction_minimum_overlap: float = 0.999, comparison_points: int = 1001) -> ModalBasisConvergenceAssessment:
    """Build at least three basis resolutions and assess successive eigenpair stability
    
    Modes are matched by common interval overlap before eigenvalue changes are evaluated
    The final qualification also checks residual, symmetry, normalization, boundary, and conditioning diagnostics
    """
    eigenvalue_tolerance = float(eigenvalue_relative_tolerance)
    overlap_threshold = float(eigenfunction_minimum_overlap)
    if not np.isfinite(eigenvalue_tolerance) or eigenvalue_tolerance <= 0.0:
        raise ValueError("eigenvalue_relative_tolerance must be finite and positive")
    if not np.isfinite(overlap_threshold) or not 0.0 < overlap_threshold <= 1.0:
        raise ValueError("eigenfunction_minimum_overlap must satisfy 0 < value <= 1")
    if grid_sequence is None:
        grid_sequence = _default_grid_sequence(numerics.n_lambda_grid, numerics.n_eta_grid)
    sequence = _validate_grid_sequence(grid_sequence)
    bases: list[ModalFBISBasis] = []
    runtimes: list[float] = []
    for n_lambda, n_eta in sequence:
        resolution_numerics = replace(numerics, n_lambda_grid=n_lambda, n_eta_grid=n_eta)
        start = perf_counter()
        bases.append(build_modal_fbis_basis(mirror_ratio=mirror_ratio, B_tilde_function=B_tilde_function, zeta_faces=zeta_faces, B_tilde_midpoints=B_tilde_midpoints, numerics=resolution_numerics, density_ratio_function=density_ratio_function, eq42_density_profile=eq42_density_profile))
        runtimes.append(perf_counter() - start)
    eigenvalue_history = tuple(np.asarray(basis.physical_basis.eigenvalues, dtype=float) for basis in bases)
    relative_change_history: list[np.ndarray] = []
    overlap_history: list[np.ndarray] = []
    permutation_history: list[np.ndarray] = []
    sign_history: list[np.ndarray] = []
    for coarse, fine in zip(bases[:-1], bases[1:], strict=True):
        coarse_values = coarse.physical_basis.eigenvalues
        fine_values = fine.physical_basis.eigenvalues
        mode_count = min(coarse_values.size, fine_values.size)
        match = match_eigenmodes_on_common_eta(coarse, fine, comparison_points=comparison_points)
        aligned_fine_values = fine_values[match.fine_mode_indices]
        changes = np.abs(aligned_fine_values - coarse_values[:mode_count]) / np.maximum(np.abs(aligned_fine_values), np.finfo(float).tiny)
        relative_change_history.append(np.asarray(changes, dtype=float))
        overlap_history.append(match.weighted_overlaps)
        permutation_history.append(match.fine_mode_indices)
        sign_history.append(match.signs)
    latest_changes = relative_change_history[-1]
    previous_changes = relative_change_history[-2]
    eigenvalues_converged = bool(np.all(previous_changes <= eigenvalue_tolerance) and np.all(latest_changes <= eigenvalue_tolerance))
    eigenfunctions_converged = bool(all(np.all(overlaps >= overlap_threshold) for overlaps in overlap_history))
    physical_history = tuple(basis.physical_basis for basis in bases)
    eigenpair_residual_history = tuple(np.asarray(physical.eigenpair_relative_residuals, dtype=float) for physical in physical_history)
    boundary_residual_history = tuple(np.asarray(physical.eigenfunction_boundary_residuals, dtype=float) for physical in physical_history)
    normalization_residual_history = tuple( np.asarray(physical.eigenfunction_normalization_residuals, dtype=float) for physical in physical_history)
    orthogonality_residual_history = tuple(np.asarray(physical.eigenfunction_orthogonality_residuals, dtype=float) for physical in physical_history)
    gram_condition_number_history = np.asarray([physical.gram_condition_number for physical in physical_history], dtype=float)
    max_eigenpair_residual = max(float(physical.max_eigenpair_relative_residual) for physical in physical_history)
    gram_condition_number = max(float(physical.gram_condition_number) for physical in physical_history)
    eta_monotonic = all(bool(physical.eta_monotonic) for physical in physical_history)
    finest_physical = physical_history[-1]
    galerkin_symmetry_tolerance = max(2.0e-2, 2.0 * eigenvalue_tolerance)
    orthogonality_tolerance = max(2.0e-2, 5.0 * (1.0 - overlap_threshold))
    diagnostics_valid = bool(eta_monotonic and np.isfinite(max_eigenpair_residual) and max_eigenpair_residual <= 1.0e-8 and np.isfinite(gram_condition_number) and gram_condition_number <= 1.0e12 and finest_physical.galerkin_symmetry_relative_error <= galerkin_symmetry_tolerance and finest_physical.gram_symmetry_relative_error <= 1.0e-10 and finest_physical.max_eigenfunction_offdiagonal_overlap <= orthogonality_tolerance and finest_physical.normalization_max_abs_error <= 1.0e-10 and finest_physical.boundary_max_abs_error <= 1.0e-8)
    converged = bool(eigenvalues_converged and eigenfunctions_converged and diagnostics_valid)
    failures: list[str] = []
    if not eigenvalues_converged:
        failures.append("retained eigenvalues did not reach successive grid tolerance")
    if not eigenfunctions_converged:
        failures.append("retained eigenfunction overlaps did not reach successive grid tolerance")
    if not eta_monotonic:
        failures.append("eta mapping was not monotonic")
    if not np.isfinite(max_eigenpair_residual) or max_eigenpair_residual > 1.0e-8:
        failures.append("projected eigenpair residual exceeded tolerance")
    if finest_physical.galerkin_symmetry_relative_error > galerkin_symmetry_tolerance:
        failures.append("Galerkin matrix asymmetry exceeded tolerance")
    if finest_physical.gram_symmetry_relative_error > 1.0e-10:
        failures.append("Gram matrix symmetry residual exceeded tolerance")
    if finest_physical.max_eigenfunction_offdiagonal_overlap > orthogonality_tolerance:
        failures.append("eigenfunction orthogonality residual exceeded tolerance")
    if finest_physical.normalization_max_abs_error > 1.0e-10:
        failures.append("eigenfunction eta=0 normalization residual exceeded tolerance")
    if finest_physical.boundary_max_abs_error > 1.0e-8:
        failures.append("eigenfunction loss boundary residual exceeded tolerance")

    return ModalBasisConvergenceAssessment(
        assessed=True,
        converged=converged,
        grid_sequence=sequence,
        eigenvalue_history=eigenvalue_history,
        eigenvalue_relative_change_history=tuple(relative_change_history),
        first_eigenvalue_relative_changes=np.asarray([changes[0] for changes in relative_change_history], dtype=float),
        eigenfunction_overlap_history=tuple(overlap_history),
        eigenfunction_mode_permutation_history=tuple(permutation_history),
        eigenfunction_sign_history=tuple(sign_history),
        eigenpair_residual_history=eigenpair_residual_history,
        boundary_residual_history=boundary_residual_history,
        normalization_residual_history=normalization_residual_history,
        orthogonality_residual_history=orthogonality_residual_history,
        gram_condition_number_history=gram_condition_number_history,
        refinement_runtime_s=np.asarray(runtimes, dtype=float),
        max_eigenpair_residual=float(max_eigenpair_residual),
        gram_condition_number=float(gram_condition_number),
        eta_monotonic=bool(eta_monotonic),
        galerkin_symmetry_relative_error=float(finest_physical.galerkin_symmetry_relative_error),
        gram_symmetry_relative_error=float(finest_physical.gram_symmetry_relative_error),
        max_eigenfunction_offdiagonal_overlap=float(finest_physical.max_eigenfunction_offdiagonal_overlap),
        normalization_max_abs_error=float(finest_physical.normalization_max_abs_error),
        boundary_max_abs_error=float(finest_physical.boundary_max_abs_error),
        turning_point_quadrature_model=TURNING_POINT_QUADRATURE_MODEL,
        separatrix_bounce_time_model=SEPARATRIX_BOUNCE_TIME_MODEL,
        failure_reason="; ".join(failures),
        basis_history=tuple(bases),
    )