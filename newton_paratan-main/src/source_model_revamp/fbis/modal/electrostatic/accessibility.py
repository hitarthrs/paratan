"""Eq 71 accessibility checks and low energy confined inventory diagnostics"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import PchipInterpolator
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid
from source_model_revamp.fbis.modal.local.mapping import _eq71_global_trapped_passing_boundary
from source_model_revamp.fbis.modal.utils import _EPS

from source_model_revamp.fbis.modal.electrostatic.eq70 import _energy_roundoff_tolerance

def _scale_conditioned_potential_interpolation(distance: np.ndarray, potential_energy_J: np.ndarray, sample_distance: np.ndarray) -> np.ndarray:
    """Interpolate potential energy after scaling by its largest magnitude to protect tiny J values"""
    scale_J = float(np.max(np.abs(potential_energy_J)))
    if scale_J == 0.0:
        return np.zeros_like(sample_distance, dtype=float)

    return np.asarray(PchipInterpolator(distance, potential_energy_J / scale_J)(sample_distance), dtype=float) * scale_J

def _effective_potential_throat_boundary_diagnostics(*, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, base_distribution_v_lambda: np.ndarray, zeta: np.ndarray, B_tilde: np.ndarray, potential_drop_magnitude_J: np.ndarray, mirror_ratio: float, throat_potential_left_energy_J: float, throat_potential_right_energy_J: float, particle_mass_kg: float, relative_tolerance: float = 1.0e-9) -> dict[str, float | int | bool | str | None]:
    """Check whether the Eq 71 trapped passing boundary is set by each throat

    For each populated invariant energy U this compares the throat value
    Λ_TP = (1 + ΔΦ_throat / U) / R_M with the minimum of
    (1 + ΔΦ(z) / U) / B_tilde(z) along each mirror half
    """
    z = np.asarray(zeta, dtype=float)
    B = np.asarray(B_tilde, dtype=float)
    phi = np.asarray(potential_drop_magnitude_J, dtype=float)
    base = np.asarray(base_distribution_v_lambda, dtype=float)
    expected = (speed_grid.centers_m_s.size, lambda_grid.centers.size)
    if z.ndim != 1 or B.shape != z.shape or phi.shape != z.shape:
        raise ValueError("zeta, B_tilde, and potential profile must be matching vectors")
    if base.shape != expected or np.any(~np.isfinite(base)):
        raise ValueError(f"base_distribution_v_lambda must be finite with shape {expected}")
    base_scale = max(float(np.max(np.abs(base))) if base.size else 0.0, np.finfo(float).tiny)
    negative_tolerance = 128.0 * np.finfo(float).eps * base_scale
    if float(np.min(base)) < -negative_tolerance:
        raise ValueError("base_distribution_v_lambda contains a significant negative value")
    base = np.maximum(base, 0.0)
    R = float(mirror_ratio)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    tolerance = float(relative_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("relative_tolerance must be finite and nonnegative")
    if np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise ValueError("B_tilde must contain positive finite values")
    if np.any(~np.isfinite(phi)) or np.any(phi < 0.0):
        raise ValueError("potential_drop_magnitude_J must contain finite nonnegative values")
    left_drop = float(throat_potential_left_energy_J)
    right_drop = float(throat_potential_right_energy_J)
    if not np.isfinite(left_drop) or left_drop < 0.0 or not np.isfinite(right_drop) or right_drop < 0.0:
        raise ValueError("fitted throat potential drops must be finite and nonnegative")
    both_sides_present = bool(np.any(z < 0.0) and np.any(z > 0.0))

    def resolved_side(sign: float, throat_drop: float) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
        """Build one midplane to throat profile with exact endpoint anchors for the diagnostic"""
        mask = sign * z >= 0.0
        if not np.any(mask):
            return None
        distance = np.abs(z[mask])
        side_B = B[mask]
        side_phi = phi[mask]
        order = np.argsort(distance)
        distance = distance[order]
        side_B = side_B[order]
        side_phi = side_phi[order]
        inside = distance < 1.0 - 128.0 * np.finfo(float).eps
        distance = distance[inside]
        side_B = side_B[inside]
        side_phi = side_phi[inside]
        if distance.size == 0 or distance[0] > 1.0e-14:
            distance = np.concatenate(([0.0], distance))
            side_B = np.concatenate(([float(np.min(B))], side_B))
            side_phi = np.concatenate(([0.0], side_phi))
        distance = np.concatenate((distance, [1.0]))
        side_B = np.concatenate((side_B, [R]))
        side_phi = np.concatenate((side_phi, [throat_drop]))
        unique_distance, unique_indices = np.unique(distance, return_index=True)
        unique_B = side_B[unique_indices]
        unique_phi = side_phi[unique_indices]
        if unique_distance.size < 2:
            return None
        sample_count = max(257, 32 * (unique_distance.size - 1) + 1)
        sample_distance = np.linspace(0.0, 1.0, sample_count)
        sampled_B = PchipInterpolator(unique_distance, unique_B)(sample_distance)
        sampled_phi = _scale_conditioned_potential_interpolation(unique_distance, unique_phi, sample_distance)
        if np.any(~np.isfinite(sampled_B)) or np.any(sampled_B <= 0.0):
            raise ValueError("interpolated B_tilde profile is nonfinite or nonpositive")
        interpolation_roundoff = _energy_roundoff_tolerance(throat_drop, float(np.max(np.abs(unique_phi))))
        if np.any(~np.isfinite(sampled_phi)) or np.any(sampled_phi < -interpolation_roundoff):
            raise ValueError("interpolated potential profile is nonfinite or negative")
        sampled_phi = np.maximum(sampled_phi, 0.0)

        return sign * sample_distance, sampled_B, sampled_phi

    left_profile = resolved_side(-1.0, left_drop)
    right_profile = resolved_side(1.0, right_drop)
    total_energy = 0.5 * float(particle_mass_kg) * np.asarray(speed_grid.centers_m_s, dtype=float) ** 2
    row_measure = np.asarray(speed_grid.shell_volumes_m3_s3, dtype=float) * np.sum(base * np.asarray(lambda_grid.widths, dtype=float)[None, :], axis=1)
    active_threshold = 128.0 * np.finfo(float).eps * max(float(np.max(row_measure)), np.finfo(float).tiny)
    active = row_measure > active_threshold
    average_drop = 0.5 * (left_drop + right_drop)
    average_boundary = _eq71_global_trapped_passing_boundary(total_energy, R, average_drop)
    open_interval = active & (average_boundary < 1.0)
    violating = np.zeros_like(open_interval)
    side_asymmetric = np.zeros_like(open_interval)
    relative_shortfall = np.zeros_like(total_energy)
    relative_asymmetry = np.zeros_like(total_energy)
    minimum_location = np.full_like(total_energy, np.nan)

    for index in np.flatnonzero(open_interval):
        U = float(total_energy[index])
        left_throat_boundary = (1.0 + left_drop / U) / R
        right_throat_boundary = (1.0 + right_drop / U) / R
        average = 0.5 * (left_throat_boundary + right_throat_boundary)
        relative_asymmetry[index] = abs(left_throat_boundary - right_throat_boundary) / max(abs(average), _EPS)
        side_asymmetric[index] = relative_asymmetry[index] > tolerance
        side_shortfalls: list[tuple[float, float]] = []
        if left_profile is not None:
            left_zeta, left_B, left_phi = left_profile
            left_values = (1.0 + left_phi / U) / left_B
            left_index = int(np.argmin(left_values))
            left_minimum = float(left_values[left_index])
            left_z = float(left_zeta[left_index])
            side_shortfalls.append((max(left_throat_boundary - left_minimum, 0.0), left_z))
        if right_profile is not None:
            right_zeta, right_B, right_phi = right_profile
            right_values = (1.0 + right_phi / U) / right_B
            right_index = int(np.argmin(right_values))
            right_minimum = float(right_values[right_index])
            right_z = float(right_zeta[right_index])
            side_shortfalls.append((max(right_throat_boundary - right_minimum, 0.0), right_z))
        if side_shortfalls:
            shortfall, location = max(side_shortfalls, key=lambda item: item[0])
            scale = max(left_throat_boundary, right_throat_boundary, _EPS)
            relative_shortfall[index] = shortfall / scale
            if shortfall > 0.0:
                minimum_location[index] = location
            violating[index] = relative_shortfall[index] > tolerance

    open_weight = float(np.sum(row_measure[open_interval]))
    violation_weight = float(np.sum(row_measure[violating]))
    maximum_shortfall_index = int(np.argmax(relative_shortfall)) if relative_shortfall.size else 0
    maximum_asymmetry = float(np.max(relative_asymmetry[open_interval])) if np.any(open_interval) else 0.0
    assessed = bool(both_sides_present and np.any(open_interval))
    valid = bool(assessed and not np.any(violating) and not np.any(side_asymmetric))
    if not both_sides_present:
        reason = "effective_potential_check_requires_both_mirror_halves"
    elif not np.any(active):
        reason = "effective_potential_check_has_no_populated_energy_rows"
    elif not np.any(open_interval):
        reason = "eq71_has_no_populated_open_confined_energy_interval"
    elif np.any(side_asymmetric):
        reason = "left_and_right_eq71_throat_boundaries_are_inconsistent"
    elif np.any(violating):
        reason = "interior_effective_potential_sets_passing_boundary"
    else:
        reason = None

    return {
        "assessed": assessed,
        "valid": valid,
        "failure_reason": reason,
        "active_energy_sample_count": int(np.count_nonzero(active)),
        "open_confined_energy_sample_count": int(np.count_nonzero(open_interval)),
        "interior_minimum_energy_sample_count": int(np.count_nonzero(violating)),
        "interior_minimum_weight_fraction": violation_weight / open_weight if open_weight > 0.0 else 0.0,
        "max_relative_boundary_shortfall": float(np.max(relative_shortfall)) if relative_shortfall.size else 0.0,
        "max_relative_throat_asymmetry": maximum_asymmetry,
        "worst_total_energy_J": float(total_energy[maximum_shortfall_index]) if relative_shortfall.size else 0.0,
        "worst_interior_location_zeta": (float(minimum_location[maximum_shortfall_index]) if relative_shortfall.size and np.isfinite(minimum_location[maximum_shortfall_index]) else None),
        "relative_tolerance": tolerance,
    }

def _closed_eq71_modal_inventory_fraction(*, speed_grid: SpeedGrid, mirror_ratio: float, throat_potential_drop_magnitude_J: float, invariant_energy_cell_population_weights: np.ndarray | None, particle_mass_kg: float) -> dict[str, float | int | bool | str]:
    """Measure modal inventory below U = ΔΦ_throat / (R_M − 1)

    Partial speed cells use the v² dv shell measure so the closed fraction is not center sampled
    """
    R = float(mirror_ratio)
    throat_drop = float(throat_potential_drop_magnitude_J)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    if not np.isfinite(throat_drop) or throat_drop < 0.0:
        raise ValueError("throat_potential_drop_magnitude_J must be finite and nonnegative")
    threshold = throat_drop / (R - 1.0)
    if throat_drop == 0.0:
        return {"assessed": True, "fraction": 0.0, "threshold_J": 0.0, "intersected_population_cell_count": 0, "intersects_distribution_support": False, "model": "eta_weighted_modal_inventory_with_partial_speed_cell_measure"}
    if invariant_energy_cell_population_weights is None:
        return {"assessed": False, "fraction": 0.0, "threshold_J": float(threshold), "intersected_population_cell_count": 0, "intersects_distribution_support": False, "model": "not_assessed_missing_eta_weighted_modal_speed_cell_inventory"}
    weights = np.asarray(invariant_energy_cell_population_weights, dtype=float)
    expected = np.asarray(speed_grid.centers_m_s, dtype=float).shape
    if weights.shape != expected:
        raise ValueError(f"invariant_energy_cell_population_weights must have shape {expected}")
    if np.any(~np.isfinite(weights)) or np.any(weights < 0.0):
        raise ValueError("invariant_energy_cell_population_weights must be finite and nonnegative")
    total = float(np.sum(weights))
    if total <= 0.0:
        return {"assessed": False, "fraction": 0.0, "threshold_J": float(threshold), "intersected_population_cell_count": 0, "intersects_distribution_support": False, "model": "not_assessed_zero_eta_weighted_modal_inventory",}
    threshold_speed = np.sqrt(2.0 * threshold / float(particle_mass_kg))
    lower = np.asarray(speed_grid.faces_m_s[:-1], dtype=float)
    upper = np.asarray(speed_grid.faces_m_s[1:], dtype=float)
    closed_upper = np.clip(threshold_speed, lower, upper)
    full_measure = upper**3 - lower**3
    closed_measure = np.maximum(closed_upper**3 - lower**3, 0.0)
    cell_fraction = np.divide(closed_measure, full_measure, out=np.zeros_like(closed_measure), where=full_measure > 0.0)
    closed_inventory = float(np.sum(weights * cell_fraction))
    active_tolerance = 128.0 * np.finfo(float).eps * max(float(np.max(weights)), np.finfo(float).tiny)
    intersected = (weights > active_tolerance) & (cell_fraction > 0.0)

    return {"assessed": True, "fraction": closed_inventory / total, "threshold_J": float(threshold), "intersected_population_cell_count": int(np.count_nonzero(intersected)), "intersects_distribution_support": bool(np.any(intersected)), "model": "eta_weighted_modal_inventory_with_partial_speed_cell_measure"}
