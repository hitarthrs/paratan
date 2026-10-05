"""Shared multi species Eq 70 to Eq 72 electrostatic reconstruction"""
from __future__ import annotations
from collections.abc import Callable, Mapping, MutableMapping
from dataclasses import dataclass
from hashlib import blake2b
from time import perf_counter
import numpy as np
from scipy.optimize import brentq
from source_model_revamp.constants import ELECTRON_CHARGE_C
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid
from source_model_revamp.fbis.modal.local.mapping import _base_distribution_interpolator, _cell_quadrature
from source_model_revamp.fbis.modal.local.pitch import _density_from_local_pitch, _phi_corrected_local_pitch_distribution, _prepare_local_pitch_quadrature
from source_model_revamp.fbis.modal.local.speed import derive_local_physical_speed_grid
from source_model_revamp.fbis.modal.full_device_lost_reconstruction import _G_profiles, _fixed_boundary_eq63_charge_density_exact_pitch_from_state, _prepare_fixed_boundary_eq63_charge_density_exact_pitch
from source_model_revamp.fbis.modal.types import ModalElectrostaticProfile
from source_model_revamp.fbis.modal.electrostatic.eq70 import _clip_potential_roundoff_only, _electron_density_fraction_eq70, _energy_roundoff_tolerance, _quasineutrality_residual_metrics
from source_model_revamp.fbis.modal.electrostatic.nodes import _exact_electrostatic_node_profiles, _symmetric_average, _symmetric_input, _midplane_density_value
from source_model_revamp.fbis.modal.electrostatic.accessibility import _closed_eq71_modal_inventory_fraction, _effective_potential_throat_boundary_diagnostics

@dataclass(frozen=True)
class _DirectEq70SpeciesTrialState:
    """Precomputed local quadrature state for repeated direct Eq 70 density roots

    Flattened arrays share one point ordering over local speed and pitch quadrature nodes
    """
    local_kinetic_energy_J: np.ndarray
    perpendicular_energy_numerator_J: np.ndarray
    quadrature_weight: np.ndarray
    interpolator: object
    normalization: float
    particle_mass_kg: float

@dataclass(frozen=True)
class _PreparedConfinedInterpolationAuditState:
    """Prepared regular grid arrays used only by the confined interpolation audit"""
    speed_coordinates: np.ndarray
    lambda_coordinates: np.ndarray
    values: np.ndarray
    supported: bool
    interpolator_type: str
    method: str
    bounds_error: bool
    fill_value: float | None

def _spacing_metadata(values: np.ndarray) -> dict[str, object]:
    """Summarize coordinate spacing for interpolation audit metadata"""
    coordinates = np.asarray(values, dtype=float)
    differences = np.diff(coordinates)
    if differences.size == 0:
        return {"coordinate_count": int(coordinates.size), "minimum_spacing": None, "maximum_spacing": None, "maximum_to_minimum_spacing_ratio": None, "exact_uniform_spacing": True}
    minimum = float(np.min(differences))
    maximum = float(np.max(differences))
    
    return {
        "coordinate_count": int(coordinates.size),
        "minimum_spacing": minimum,
        "maximum_spacing": maximum,
        "maximum_to_minimum_spacing_ratio": None if minimum <= 0.0 else float(maximum / minimum),
        "exact_uniform_spacing": bool(np.all(differences == differences[0])),
    }

def _prepare_confined_interpolation_audit_state(interpolator: object) -> tuple[_PreparedConfinedInterpolationAuditState | None, dict[str, object]]:
    """Extract linear regular grid data when the interpolation audit can reproduce the active interpolator"""
    grid = getattr(interpolator, "grid", None)
    values = getattr(interpolator, "values", None)
    method = str(getattr(interpolator, "method", ""))
    bounds_error = bool(getattr(interpolator, "bounds_error", True))
    raw_fill = getattr(interpolator, "fill_value", None)
    fill_value = None
    if raw_fill is not None:
        try:
            fill_value = float(raw_fill)
        except Exception:
            fill_value = None
    interpolator_type = f"{type(interpolator).__module__}.{type(interpolator).__qualname__}"
    supported = bool(isinstance(grid, tuple) and len(grid) == 2 and values is not None and method == "linear" and not bounds_error and fill_value == 0.0)
    prepared = None
    speed = np.asarray([], dtype=float)
    lambdas = np.asarray([], dtype=float)
    value_array = np.asarray([], dtype=float)
    if supported:
        speed = np.asarray(grid[0], dtype=float)
        lambdas = np.asarray(grid[1], dtype=float)
        value_array = np.asarray(values, dtype=float)
        supported = bool(
            speed.ndim == 1
            and lambdas.ndim == 1
            and speed.size >= 2
            and lambdas.size >= 2
            and np.all(np.isfinite(speed))
            and np.all(np.isfinite(lambdas))
            and np.all(np.diff(speed) > 0.0)
            and np.all(np.diff(lambdas) > 0.0)
            and value_array.shape == (speed.size, lambdas.size)
            and np.all(np.isfinite(value_array))
        )
        if supported:
            prepared = _PreparedConfinedInterpolationAuditState(
                speed_coordinates=speed,
                lambda_coordinates=lambdas,
                values=value_array,
                supported=True,
                interpolator_type=interpolator_type,
                method=method,
                bounds_error=bounds_error,
                fill_value=fill_value,
            )
    metadata = {
        "interpolator_type": interpolator_type,
        "method": method,
        "bounds_error": bounds_error,
        "fill_value": fill_value,
        "prepared_linear_candidate_supported": bool(prepared is not None),
        "speed_grid": _spacing_metadata(speed),
        "lambda_grid": _spacing_metadata(lambdas),
    }
   
    return prepared, metadata

def _diagnostic_prepared_linear_interpolation(*, state: _PreparedConfinedInterpolationAuditState, points: np.ndarray) -> tuple[np.ndarray, float, float, bytes, int]:
    """Evaluate bilinear interpolation for timing and equality checks without changing solver values"""
    query = np.asarray(points, dtype=float)
    speed = query[:, 0]
    lambdas = query[:, 1]
    locate_started = perf_counter()
    outside = ((speed < state.speed_coordinates[0]) | (speed > state.speed_coordinates[-1]) | (lambdas < state.lambda_coordinates[0]) | (lambdas > state.lambda_coordinates[-1]))
    speed_index = np.searchsorted(state.speed_coordinates, speed, side="right") - 1
    lambda_index = np.searchsorted(state.lambda_coordinates, lambdas, side="right") - 1
    speed_index = np.clip(speed_index, 0, state.speed_coordinates.size - 2)
    lambda_index = np.clip(lambda_index, 0, state.lambda_coordinates.size - 2)
    locate_s = perf_counter() - locate_started
    evaluate_started = perf_counter()
    speed_lower = state.speed_coordinates[speed_index]
    speed_upper = state.speed_coordinates[speed_index + 1]
    lambda_lower = state.lambda_coordinates[lambda_index]
    lambda_upper = state.lambda_coordinates[lambda_index + 1]
    speed_weight = (speed - speed_lower) / (speed_upper - speed_lower)
    lambda_weight = (lambdas - lambda_lower) / (lambda_upper - lambda_lower)
    values = state.values
    v00 = values[speed_index, lambda_index]
    v10 = values[speed_index + 1, lambda_index]
    v01 = values[speed_index, lambda_index + 1]
    v11 = values[speed_index + 1, lambda_index + 1]
    result = ((1.0 - speed_weight) * (1.0 - lambda_weight) * v00 + speed_weight * (1.0 - lambda_weight) * v10 + (1.0 - speed_weight) * lambda_weight * v01 + speed_weight * lambda_weight * v11)
    result = np.asarray(result, dtype=float)
    result[outside] = 0.0
    evaluate_s = perf_counter() - evaluate_started
    packed = np.column_stack((speed_index, lambda_index, outside.astype(np.int8, copy=False)))
    digest = blake2b(np.ascontiguousarray(packed).tobytes(), digest_size=10).digest()
 
    return result, float(locate_s), float(evaluate_s), digest, int(np.count_nonzero(outside))

def _direct_eq70_confined_density(*, state: _DirectEq70SpeciesTrialState, mirror_ratio: float, B_tilde: float, local_potential_drop_magnitude_J: float, throat_potential_drop_magnitude_J: float, _interpolation_audit: MutableMapping[str, object] | None = None, _interpolation_audit_state: _PreparedConfinedInterpolationAuditState | None = None, _interpolation_audit_key: str = "confined_density", _interpolation_tracker: MutableMapping[str, object] | None = None) -> float:
    """Evaluate confined charge density from the Eq 71 and Eq 72 mapped magnetic reference

    Local quadrature points map to U and Λ, then Eq 72 maps accessible points to Λ* before f(v, Λ*) interpolation
    The return value uses the local gyrotropic measure 2π v² dv dξ
    """
    R = float(mirror_ratio)
    B = float(B_tilde)
    local_drop = float(local_potential_drop_magnitude_J)
    throat_drop = float(throat_potential_drop_magnitude_J)
    total_energy = state.local_kinetic_energy_J - local_drop
    positive = total_energy > 0.0
    global_lambda = np.divide(state.perpendicular_energy_numerator_J, B * total_energy, out=np.full_like(total_energy, np.nan), where=positive)
    boundary = np.full_like(total_energy, np.inf)
    boundary[positive] = (1.0 / R) * (1.0 + throat_drop / total_energy[positive])
    denominator = 1.0 - boundary
    lambda_star = np.full_like(total_energy, np.nan)
    valid_denominator = denominator > 0.0
    np.divide((1.0 - 1.0 / R) * (1.0 - global_lambda), denominator, out=lambda_star, where=valid_denominator)
    lambda_star = 1.0 - lambda_star
    roundoff_tolerance = 128.0 * np.finfo(float).eps
    valid = (positive & (boundary < 1.0) & (global_lambda >= boundary - roundoff_tolerance) & (global_lambda <= 1.0 + roundoff_tolerance) & (lambda_star >= 1.0 / R - roundoff_tolerance) & (lambda_star <= 1.0 + roundoff_tolerance))
    if not np.any(valid):
        if _interpolation_tracker is not None:
            tracker = _interpolation_tracker.setdefault(str(_interpolation_audit_key), {"empty_calls": 0})
            if isinstance(tracker, dict):
                tracker["empty_calls"] = int(tracker.get("empty_calls", 0)) + 1
        return 0.0
    points = np.column_stack((np.sqrt(2.0 * total_energy[valid] / float(state.particle_mass_kg)), np.clip(lambda_star[valid], 1.0 / R, 1.0)))
    production_started = perf_counter() if _interpolation_audit is not None else None
    production_values = np.asarray(state.interpolator(points), dtype=float)
    if production_started is not None:
        production_s = perf_counter() - production_started
        runtime_group = _interpolation_audit.setdefault("interpolation_runtime_s", {})
        count_group = _interpolation_audit.setdefault("interpolation_counts", {})
        if not isinstance(runtime_group, dict) or not isinstance(count_group, dict):
            raise TypeError("Eq 70 interpolation audit storage must be a mapping")
        runtime_group[f"{_interpolation_audit_key}:production"] = float(runtime_group.get(f"{_interpolation_audit_key}:production", 0.0)) + float(production_s)
        count_group[f"{_interpolation_audit_key}:calls"] = int(count_group.get(f"{_interpolation_audit_key}:calls", 0)) + 1
        count_group[f"{_interpolation_audit_key}:points"] = int(count_group.get(f"{_interpolation_audit_key}:points", 0)) + int(points.shape[0])
        if _interpolation_audit_state is not None:
            candidate, locate_s, evaluate_s, index_digest, outside_count = _diagnostic_prepared_linear_interpolation(state=_interpolation_audit_state, points=points)
            runtime_group[f"{_interpolation_audit_key}:candidate_locate"] = float(runtime_group.get(f"{_interpolation_audit_key}:candidate_locate", 0.0)) + locate_s
            runtime_group[f"{_interpolation_audit_key}:candidate_evaluate"] = float(runtime_group.get(f"{_interpolation_audit_key}:candidate_evaluate", 0.0)) + evaluate_s
            count_group[f"{_interpolation_audit_key}:candidate_calls"] = int(count_group.get(f"{_interpolation_audit_key}:candidate_calls", 0)) + 1
            count_group[f"{_interpolation_audit_key}:outside_points"] = int(count_group.get(f"{_interpolation_audit_key}:outside_points", 0)) + outside_count
            difference = np.abs(candidate - production_values)
            scale = np.maximum(np.maximum(np.abs(candidate), np.abs(production_values)), np.finfo(float).tiny)
            maximum_group = _interpolation_audit.setdefault("interpolation_maximum_difference", {})
            if not isinstance(maximum_group, dict):
                raise TypeError("Eq 70 interpolation maximum difference storage must be a mapping")
            maximum_group[f"{_interpolation_audit_key}:absolute"] = max(float(maximum_group.get(f"{_interpolation_audit_key}:absolute", 0.0)), float(np.max(difference)) if difference.size else 0.0)
            maximum_group[f"{_interpolation_audit_key}:relative"] = max(float(maximum_group.get(f"{_interpolation_audit_key}:relative", 0.0)), float(np.max(difference / scale)) if difference.size else 0.0)
            if np.array_equal(candidate, production_values):
                count_group[f"{_interpolation_audit_key}:exact_equal_calls"] = int(count_group.get(f"{_interpolation_audit_key}:exact_equal_calls", 0)) + 1
            if _interpolation_tracker is not None:
                tracker = _interpolation_tracker.setdefault(str(_interpolation_audit_key), {"index_seen": set(), "index_calls": 0, "index_repeat_calls": 0, "empty_calls": 0})
                if not isinstance(tracker, dict):
                    raise TypeError("Eq 70 interpolation tracker must be a mapping")
                seen = tracker.get("index_seen")
                if not isinstance(seen, set):
                    raise TypeError("Eq 70 interpolation index tracker must be a set")
                tracker["index_calls"] = int(tracker.get("index_calls", 0)) + 1
                if index_digest in seen:
                    tracker["index_repeat_calls"] = int(tracker.get("index_repeat_calls", 0)) + 1
                else:
                    seen.add(index_digest)
    sampled = np.zeros_like(total_energy)
    sampled[valid] = production_values
    sampled = np.maximum(sampled * float(state.normalization), 0.0)
  
    return 2.0 * np.pi * float(np.sum(sampled * state.quadrature_weight))

@dataclass(frozen=True)
class _BoundedRootSearchResult:
    """Result of one bounded quasineutral root search with branch and scan diagnostics"""
    root: float
    available: bool
    scan_used: bool
    root_count: int
    evaluation_count: int

def _bounded_continuation_root(*, residual: Callable[[float], float], lower: float, upper: float, seed: float, scan_points: int, xtol: float, rtol: float) -> _BoundedRootSearchResult:
    """Find a bounded scalar root while preferring continuation from the supplied seed

    Endpoint and seed brackets are tried first, then an odd point interval scan locates interior roots
    If multiple roots exist the root nearest the seed is selected
    """
    lo = float(lower)
    hi = float(upper)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi < lo:
        raise ValueError("bounded root interval must be finite and ordered")
    points = int(scan_points)
    if points < 5 or points % 2 == 0:
        raise ValueError("scan_points must be an odd integer of at least five")
    seed_value = float(np.clip(float(seed), lo, hi))
    cache: dict[float, float] = {}

    def evaluate(value: float) -> float:
        """Evaluate and cache one finite residual value"""
        key = float(value)
        if key not in cache:
            result = float(residual(key))
            if not np.isfinite(result):
                raise ValueError("bounded root residual must be finite")
            cache[key] = result
        return cache[key]

    low = evaluate(lo)
    high = evaluate(hi)
    zero_scale = max(abs(low), abs(high), 1.0)
    zero_tolerance = 512.0 * np.finfo(float).eps * zero_scale

    def is_zero(value: float) -> bool:
        """Return whether a residual is zero to the floating point endpoint scale"""
        return abs(float(value)) <= zero_tolerance

    def opposite(left: float, right: float) -> bool:
        """Return whether two residuals strictly bracket a sign change"""
        return bool((left < 0.0 < right) or (right < 0.0 < left))

    def solve_bracket(left: float, right: float) -> float | None:
        """Solve one validated sign change bracket with Brent root finding"""
        f_left = evaluate(left)
        f_right = evaluate(right)
        if is_zero(f_left):
            return float(left)
        if is_zero(f_right):
            return float(right)
        if not opposite(f_left, f_right):
            return None
        try:
            return float(brentq(evaluate, float(left), float(right), xtol=float(xtol), rtol=float(rtol)))
        except ValueError:
            return None

    def unique_roots(values: list[float]) -> list[float]:
        """Merge numerically coincident roots using the configured absolute tolerance"""
        if not values:
            return []
        merge_tolerance = max(8.0 * float(xtol), 64.0 * np.finfo(float).eps * max(abs(lo), abs(hi), 1.0))
        result: list[float] = []
        for value in sorted(float(item) for item in values):
            if not result or abs(value - result[-1]) > merge_tolerance:
                result.append(value)
        return result

    if is_zero(low):
        return _BoundedRootSearchResult(lo, True, False, 1, len(cache))
    if is_zero(high):
        return _BoundedRootSearchResult(hi, True, False, 1, len(cache))
    if opposite(low, high):
        root = solve_bracket(lo, hi)
        if root is not None:
            return _BoundedRootSearchResult(root, True, False, 1, len(cache))

    seed_residual = evaluate(seed_value)
    if is_zero(seed_residual):
        return _BoundedRootSearchResult(seed_value, True, False, 1, len(cache))
    continuation_roots: list[float] = []
    if lo < seed_value < hi:
        left_root = solve_bracket(lo, seed_value)
        right_root = solve_bracket(seed_value, hi)
        if left_root is not None:
            continuation_roots.append(left_root)
        if right_root is not None:
            continuation_roots.append(right_root)
    continuation_roots = unique_roots(continuation_roots)
    if continuation_roots:
        selected = min(continuation_roots, key=lambda value: (abs(value - seed_value), value))
        return _BoundedRootSearchResult(float(selected), True, False, len(continuation_roots), len(cache))

    probe_values = np.unique(np.concatenate((np.linspace(lo, hi, points, dtype=float), np.asarray([seed_value], dtype=float))))
    sampled = [(float(value), evaluate(float(value))) for value in probe_values]
    roots: list[float] = [value for value, function_value in sampled if is_zero(function_value)]
    for (left, f_left), (right, f_right) in zip(sampled[:-1], sampled[1:], strict=True):
        if opposite(f_left, f_right):
            root = solve_bracket(left, right)
            if root is not None:
                roots.append(root)
    roots = unique_roots(roots)
    if roots:
        selected = min(roots, key=lambda value: (abs(value - seed_value), value))
        return _BoundedRootSearchResult(float(selected), True, True, len(roots), len(cache))
    fallback = lo if abs(low) <= abs(high) else hi
    return _BoundedRootSearchResult(float(fallback), False, True, 0, len(cache))

@dataclass(frozen=True)
class ModalElectrostaticSpeciesInput:
    """One fast ion species input to the shared electrostatic solve

    base_distribution_v_lambda has shape (n_speed, n_lambda)
    Optional Eq 63 fields activate the fixed magnetic boundary lost ion charge contribution
    """
    species_id: str
    speed_grid: SpeedGrid
    lambda_grid: LambdaGrid
    pitch_grid: PitchGrid
    base_distribution_v_lambda: np.ndarray
    particle_mass_kg: float
    eta_to_local_phase_space_normalization: float
    invariant_energy_cell_population_weights: np.ndarray | None = None
    charge_number: float = 1.0
    eq63_left_H_U: np.ndarray | None = None
    eq63_right_H_U: np.ndarray | None = None
    eq63_left_throat_rate_v_lambda_s: np.ndarray | None = None
    eq63_right_throat_rate_v_lambda_s: np.ndarray | None = None
    eq63_left_parallel_temperature_J: float | None = None
    eq63_right_parallel_temperature_J: float | None = None
    eq63_geometry_factor_G: float | None = None


@dataclass(frozen=True)
class ModalElectrostaticSystemResult:
    """Shared electrostatic profile and species resolved local reconstructions

    Local Λ and pitch arrays use shape (n_z, n_local_speed, n_coordinate)
    Density arrays use shape (n_z,)
    """
    profile: ModalElectrostaticProfile
    local_distribution_z_v_lambda_by_species: Mapping[str, np.ndarray]
    local_distribution_z_v_pitch_by_species: Mapping[str, np.ndarray]
    local_density_m3_by_species: Mapping[str, np.ndarray]
    local_speed_grid_by_species: Mapping[str, SpeedGrid]
    mapping_diagnostics_by_species: Mapping[str, Mapping[str, float | int]]

def _solve_system_phi_profile_quasineutrality(*, species_inputs: Mapping[str, ModalElectrostaticSpeciesInput], zeta: np.ndarray, B_tilde: np.ndarray, cell_volumes_m3: np.ndarray, mirror_ratio: float, wall_barrier_energy_J: float, electron_temperature_J: float, iterations: int, root_scan_points: int=33, relaxation: float, relative_tolerance: float, electron_midplane_density_m3: float | None=None, electron_collision_density_m3: float | None=None, background_positive_charge_density_m3: np.ndarray | None=None, background_midplane_positive_charge_density_m3: float | None=None, background_left_throat_positive_charge_density_m3: float=0.0, background_right_throat_positive_charge_density_m3: float=0.0, target_volume_averaged_ion_density_m3: float | None=None, local_velocity_quadrature_order: int=3, low_energy_weight_fraction_tolerance: float=0.001, initial_potential_energy_J: np.ndarray | None=None, zeta_faces: np.ndarray | None=None, _confined_interpolation_audit: MutableMapping[str, object] | None=None) -> ModalElectrostaticSystemResult:
    """Solve the shared multi species Eq 70 quasineutral potential profile

    When no electron midplane density is prescribed the Eq 70 parent Maxwellian amplitude floats with the maximum zero drop positive charge state
    Left and right throat barriers are solved together with this scalar closure, then each axial cell receives a direct bounded quasineutral root
    Confined Eq 71 and Eq 72 charge and fixed boundary Eq 63 lost ion charge both enter n_i
    """
    audit = _confined_interpolation_audit
    audit_started = perf_counter() if audit is not None else None
    if audit is not None:
        audit.clear()
        audit.update({"schema": "eq70_confined_interpolation_audit_v1"})
    inputs = dict(species_inputs)
    if not inputs:
        raise ValueError("the shared electrostatic solve requires at least one fast species")
    for key, value in inputs.items():
        if str(key) != value.species_id:
            raise ValueError("electrostatic species keys must match species identifiers")
    z = np.asarray(zeta, dtype=float)
    B = np.asarray(B_tilde, dtype=float)
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    n_z = z.size
    if n_z == 0 or B.shape != (n_z,) or volumes.shape != (n_z,):
        raise ValueError("zeta, B_tilde, and cell volumes must have matching nonempty lengths")
    if np.any(~np.isfinite(z)) or np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise ValueError("electrostatic coordinates and magnetic field must be finite and positive where required")
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("cell_volumes_m3 must contain positive finite values")
    faces = None if zeta_faces is None else np.asarray(zeta_faces, dtype=float)
    if faces is not None and (faces.ndim != 1 or faces.size != n_z + 1 or np.any(~np.isfinite(faces)) or np.any(np.diff(faces) <= 0.0)):
        raise ValueError("zeta_faces must be a finite increasing edge grid matching the electrostatic cells")
    R = float(mirror_ratio)
    wall = float(wall_barrier_energy_J)
    Te = float(electron_temperature_J)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    magnetic_tolerance = 128.0 * np.finfo(float).eps * max(R, float(np.max(B)))
    if np.any(B > R + magnetic_tolerance):
        raise ValueError("B_tilde exceeds mirror_ratio inside the fitted throat domain")
    if not np.isfinite(wall) or wall < 0.0:
        raise ValueError("wall_barrier_energy_J must be finite and nonnegative")
    if not np.isfinite(Te) or Te <= 0.0:
        raise ValueError("electron_temperature_J must be positive and finite")
    n_iter = max(int(iterations), 1)
    root_scan_count = int(root_scan_points)
    if root_scan_count < 5 or root_scan_count % 2 == 0:
        raise ValueError("root_scan_points must be an odd integer of at least five")
    relax = float(relaxation)
    tolerance = float(relative_tolerance)
    low_energy_tolerance = float(low_energy_weight_fraction_tolerance)
    if not np.isfinite(relax) or relax <= 0.0 or relax > 1.0:
        raise ValueError("relaxation must lie in the interval (0, 1]")
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("relative_tolerance must be finite and nonnegative")
    if not np.isfinite(low_energy_tolerance) or low_energy_tolerance < 0.0:
        raise ValueError("low_energy_weight_fraction_tolerance must be finite and nonnegative")
    target_volume_average = None if target_volume_averaged_ion_density_m3 is None else float(target_volume_averaged_ion_density_m3)
    if target_volume_average is not None and (not np.isfinite(target_volume_average) or target_volume_average < 0.0):
        raise ValueError("target_volume_averaged_ion_density_m3 must be finite and nonnegative")
    background = np.zeros(n_z, dtype=float) if background_positive_charge_density_m3 is None else np.asarray(background_positive_charge_density_m3, dtype=float)
    if background.shape != (n_z,) or np.any(~np.isfinite(background)) or np.any(background < 0.0):
        raise ValueError("background_positive_charge_density_m3 must be finite and nonnegative with one value per cell")
    background_left = float(background_left_throat_positive_charge_density_m3)
    background_right = float(background_right_throat_positive_charge_density_m3)
    if not np.isfinite(background_left) or background_left < 0.0 or not np.isfinite(background_right) or background_right < 0.0:
        raise ValueError("background throat densities must be finite and nonnegative")
    background_midplane_input = None if background_midplane_positive_charge_density_m3 is None else float(background_midplane_positive_charge_density_m3)
    if background_midplane_input is not None and (not np.isfinite(background_midplane_input) or background_midplane_input < 0.0):
        raise ValueError("background_midplane_positive_charge_density_m3 must be finite and nonnegative")
    electron_midplane_prescribed = None if electron_midplane_density_m3 is None else float(electron_midplane_density_m3)
    if electron_midplane_prescribed is not None and (not np.isfinite(electron_midplane_prescribed) or electron_midplane_prescribed <= 0.0):
        raise ValueError("electron_midplane_density_m3 must be positive and finite when prescribed")
    electron_collision = None if electron_collision_density_m3 is None else float(electron_collision_density_m3)
    if electron_collision is not None and (not np.isfinite(electron_collision) or electron_collision <= 0.0):
        raise ValueError("electron_collision_density_m3 must be positive and finite")
    symmetric = _symmetric_input(z, B, volumes)
    background_order = np.argsort(z)
    background_midplane = float(np.interp(0.0, z[background_order], background[background_order])) if background_midplane_input is None else background_midplane_input

    local_speed_grids: dict[str, SpeedGrid] = {}
    interpolators: dict[str, object] = {}
    local_quadrature_by_species: dict[str, object] = {}
    direct_trial_state_by_species: dict[str, _DirectEq70SpeciesTrialState] = {}
    # Prepare species local quadrature once because direct root searches reuse the same velocity nodes
    for species_id, value in inputs.items():
        mass = float(value.particle_mass_kg)
        normalization = float(value.eta_to_local_phase_space_normalization)
        if not np.isfinite(mass) or mass <= 0.0 or not np.isfinite(normalization) or normalization <= 0.0:
            raise ValueError("electrostatic species mass and normalization must be positive and finite")
        local_speed_grids[species_id], _ = derive_local_physical_speed_grid(invariant_speed_grid=value.speed_grid, maximum_potential_drop_magnitude_J=wall, particle_mass_kg=mass)
        interpolators[species_id] = _base_distribution_interpolator(value.speed_grid, value.lambda_grid, value.base_distribution_v_lambda)
        quadrature = _prepare_local_pitch_quadrature(local_speed_grids[species_id], value.pitch_grid, int(local_velocity_quadrature_order))
        local_quadrature_by_species[species_id] = quadrature
        speed, pitch = np.broadcast_arrays(quadrature.speed_nodes_m_s, quadrature.pitch_nodes)
        local_kinetic = 0.5 * mass * speed**2
        direct_trial_state_by_species[species_id] = _DirectEq70SpeciesTrialState(
            local_kinetic_energy_J=np.ascontiguousarray(local_kinetic.reshape(-1)),
            perpendicular_energy_numerator_J=np.ascontiguousarray((local_kinetic * (1.0 - pitch**2)).reshape(-1)),
            quadrature_weight=np.ascontiguousarray(np.broadcast_to(quadrature.quadrature_weight, local_kinetic.shape).reshape(-1)),
            interpolator=interpolators[species_id],
            normalization=normalization,
            particle_mass_kg=mass,
        )

    interpolation_audit_state_by_species: dict[str, _PreparedConfinedInterpolationAuditState | None] = {}
    interpolation_tracker: dict[str, object] | None = {} if audit is not None else None
    if audit is not None:
        grid_metadata: dict[str, object] = {}
        for species_id, interpolator in interpolators.items():
            prepared, metadata = _prepare_confined_interpolation_audit_state(interpolator)
            interpolation_audit_state_by_species[species_id] = prepared
            grid_metadata[species_id] = metadata
        audit["interpolator_metadata_by_species"] = grid_metadata

    eq63_species_ids: list[str] = []
    for species_id, value in inputs.items():
        fields = (value.eq63_left_H_U, value.eq63_right_H_U, value.eq63_left_parallel_temperature_J, value.eq63_right_parallel_temperature_J, value.eq63_geometry_factor_G)
        if any(item is not None for item in fields):
            if not all(item is not None for item in fields):
                raise ValueError(f"{species_id} Eq 63 electrostatic input is incomplete")
            if faces is None:
                raise ValueError("zeta_faces are required when Eq 63 central lost density is active")
            eq63_species_ids.append(species_id)
    eq63_G_by_species: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    if eq63_species_ids:
        for species_id in eq63_species_ids:
            value = inputs[species_id]
            eq63_G_by_species[species_id] = _G_profiles(
                z_edges_m=faces,
                B_tilde_centers=B,
                left_throat_z_m=-1.0,
                right_throat_z_m=1.0,
                half_length_m=1.0,
                mirror_ratio=R,
                target_G_at_throat=float(value.eq63_geometry_factor_G),
            )

    def exact_midplane_G(species_id: str) -> tuple[float, float]:
        """Split the Eq 63 geometry integral at the exact ζ = 0 location for one species"""
        value = inputs[species_id]
        if faces is None:
            return 0.0, 0.0
        widths = np.diff(faces)
        integrand = widths / B * np.sqrt(np.maximum(1.0 - B / R, 0.0))
        total_raw = float(np.sum(integrand))
        if total_raw <= 0.0:
            raise ValueError("central magnetic profile gives a nonpositive Eq 63 geometry integral")
        scale = float(value.eq63_geometry_factor_G) / total_raw
        raw_density = np.divide(np.sqrt(np.maximum(1.0 - B / R, 0.0)), B, out=np.zeros_like(B), where=B > 0.0) * scale
        right = 0.0
        for index in range(n_z):
            lo = float(faces[index])
            hi = float(faces[index + 1])
            if hi <= 0.0:
                right += (hi - lo) * float(raw_density[index])
            elif lo < 0.0 < hi:
                right += (0.0 - lo) * float(raw_density[index])
        total = float(value.eq63_geometry_factor_G)
        return max(total - right, 0.0), max(right, 0.0)

    confined_density_cache: dict[tuple[str, float, float, float], float] = {}

    def species_confined_density(species_id: str, *, field_value: float, potential_drop_J: float, throat_drop_J: float) -> float:
        """Return one species confined density from cached direct Eq 71 and Eq 72 quadrature"""
        key = (str(species_id), float(field_value), float(potential_drop_J), float(throat_drop_J))
        if key not in confined_density_cache:
            confined_density_cache[key] = _direct_eq70_confined_density(
                state=direct_trial_state_by_species[species_id],
                mirror_ratio=R,
                B_tilde=key[1],
                local_potential_drop_magnitude_J=key[2],
                throat_potential_drop_magnitude_J=key[3],
                _interpolation_audit=audit,
                _interpolation_audit_state=interpolation_audit_state_by_species.get(species_id),
                _interpolation_audit_key=f"confined_density:{species_id}",
                _interpolation_tracker=interpolation_tracker,
            )
        return confined_density_cache[key]

    def confined_charge_density(*, field_value: float, potential_drop_J: float, throat_drop_J: float) -> float:
        """Return charge number weighted confined fast ion density summed over species"""
        return float(sum(float(inputs[key].charge_number) * species_confined_density(key, field_value=field_value, potential_drop_J=potential_drop_J, throat_drop_J=throat_drop_J) for key in inputs))

    # Eq 63 lost branches contribute positive charge separately from the confined Eq 71 mapping
    eq63_density_state_by_species = {
        species_id: _prepare_fixed_boundary_eq63_charge_density_exact_pitch(
            invariant_speed_grid=inputs[species_id].speed_grid,
            local_speed_grid=local_speed_grids[species_id],
            pitch_grid=inputs[species_id].pitch_grid,
            left_H_U=np.asarray(inputs[species_id].eq63_left_H_U, dtype=float),
            right_H_U=np.asarray(inputs[species_id].eq63_right_H_U, dtype=float),
            left_parallel_temperature_J=float(inputs[species_id].eq63_left_parallel_temperature_J),
            right_parallel_temperature_J=float(inputs[species_id].eq63_right_parallel_temperature_J),
            mirror_ratio=R,
            particle_mass_kg=float(inputs[species_id].particle_mass_kg),
            charge_number=float(inputs[species_id].charge_number),
            speed_quadrature_order=int(local_velocity_quadrature_order),
        )
        for species_id in eq63_species_ids
    }
    eq63_density_cache: dict[tuple[str, float, float, float, float], float] = {}

    def eq63_species_charge_density(species_id: str, *, field_value: float, potential_drop_J: float, left_G: float, right_G: float) -> float:
        """Return cached charge number weighted Eq 63 density for one species and local state"""
        if species_id not in eq63_species_ids:
            return 0.0
        key = (str(species_id), float(field_value), float(potential_drop_J), float(left_G), float(right_G))
        if key not in eq63_density_cache:
            eq63_density_cache[key] = _fixed_boundary_eq63_charge_density_exact_pitch_from_state(
                state=eq63_density_state_by_species[species_id],
                B_tilde=key[1],
                potential_drop_magnitude_J=key[2],
                left_geometry_G=key[3],
                right_geometry_G=key[4],
            )
       
        return eq63_density_cache[key]

    def eq63_cell_charge_density(index: int, candidate: float) -> float:
        """Return total Eq 63 lost ion charge density in one axial cell"""
        total = 0.0
        for species_id in eq63_species_ids:
            left_G, right_G = eq63_G_by_species[species_id]
            total += eq63_species_charge_density(species_id, field_value=float(B[index]), potential_drop_J=float(candidate), left_G=float(left_G[index]), right_G=float(right_G[index]))
       
        return float(total)

    exact_midplane_G_by_species = {species_id: exact_midplane_G(species_id) for species_id in eq63_species_ids}

    def eq63_midplane_charge_density(candidate: float) -> float:
        """Return total Eq 63 lost ion charge density at the exact midplane"""
        total = 0.0
        for species_id in eq63_species_ids:
            left_G, right_G = exact_midplane_G_by_species[species_id]
            total += eq63_species_charge_density(species_id, field_value=1.0, potential_drop_J=float(candidate), left_G=left_G, right_G=right_G)
      
        return float(total)

    def eq63_throat_charge_density(side: str, candidate: float) -> float:
        """Return total directed Eq 63 lost ion charge density at one exact throat"""
        total = 0.0
        for species_id in eq63_species_ids:
            value = inputs[species_id]
            G = float(value.eq63_geometry_factor_G)
            total += eq63_species_charge_density(species_id, field_value=R, potential_drop_J=float(candidate), left_G=G if side == "left" else 0.0, right_G=G if side == "right" else 0.0)
       
        return float(total)


    electron_fraction_cache: dict[float, float] = {}

    def electron_fraction(candidate: float) -> float:
        """Return and cache the Eq 70 electron fraction at one potential energy coordinate"""
        key = float(candidate)
        if key not in electron_fraction_cache:
            electron_fraction_cache[key] = float(_electron_density_fraction_eq70(key, wall, Te))
       
        return electron_fraction_cache[key]

    electron_fraction_zero = max(electron_fraction(0.0), np.finfo(float).tiny)
    root_xtol = _energy_roundoff_tolerance(wall, Te)
    root_rtol = 8.0 * np.finfo(float).eps

    def normalized_root_residual(electron: float, ion: float, parent_n0: float) -> float:
        """Return |n_e − n_i| normalized by local support or the Eq 70 reference density"""
        local_scale = max(abs(float(electron)), abs(float(ion)))
        reference = max(float(parent_n0) * electron_fraction_zero, np.finfo(float).tiny)
        support_floor = max(reference * 1.0e-12, np.finfo(float).tiny)
        denominator = local_scale if local_scale >= support_floor and local_scale > 0.0 else reference
      
        return abs(float(electron) - float(ion)) / denominator

    def throat_state(side: str, candidate: float, parent_n0: float, background_value: float) -> tuple[float, float, float]:
        """Evaluate electron, confined ion, Eq 63, and background charge at one throat candidate"""
        confined = confined_charge_density(field_value=R, potential_drop_J=float(candidate), throat_drop_J=float(candidate))
        lost = eq63_throat_charge_density(side, float(candidate))
        ion = confined + lost + float(background_value)
        electron = float(parent_n0) * electron_fraction(float(candidate))
      
        return electron - ion, electron, ion

    def solve_throat(side: str, parent_n0: float, seed: float, background_value: float) -> tuple[float, bool, _BoundedRootSearchResult, float]:
        """Solve and qualify one left or right throat quasineutrality root"""
        state_cache: dict[float, tuple[float, float, float]] = {}

        def state(candidate: float) -> tuple[float, float, float]:
            """Return cached throat residual, electron density, and ion charge density"""
            key = float(candidate)
            if key not in state_cache:
                state_cache[key] = throat_state(side, key, parent_n0, background_value)
            return state_cache[key]

        def residual(candidate: float) -> float:
            """Return the cached throat quasineutrality residual"""
            return state(candidate)[0]

        search = _bounded_continuation_root(residual=residual, lower=0.0, upper=wall, seed=float(seed), scan_points=root_scan_count, xtol=root_xtol, rtol=root_rtol)
        _, electron, ion = state(float(search.root))
        quality = normalized_root_residual(electron, ion, parent_n0)
     
        return float(search.root), bool(search.available and quality <= tolerance), search, float(quality)

    initial_roundoff_clips = 0
    continuation_accepted = False
    if initial_potential_energy_J is not None:
        initial_phi = np.asarray(initial_potential_energy_J, dtype=float)
        if initial_phi.shape != (n_z,) or np.any(~np.isfinite(initial_phi)):
            raise ValueError("initial_potential_energy_J must contain one finite value per axial cell")
        continuation_tolerance = _energy_roundoff_tolerance(wall, float(np.max(np.abs(initial_phi))) if initial_phi.size else 0.0)
        continuation_accepted = bool(np.all(initial_phi >= -continuation_tolerance) and np.all(initial_phi <= wall + continuation_tolerance))
        if continuation_accepted:
            phi_seed, initial_roundoff_clips = _clip_potential_roundoff_only(initial_phi, wall, context="initial floating potential continuation profile")
    if not continuation_accepted:
        if wall <= 0.0:
            phi_seed = np.zeros(n_z, dtype=float)
        else:
            B_midplane = float(np.min(B))
            denominator = R - B_midplane
            if denominator <= magnetic_tolerance:
                raise ValueError("mirror_ratio must exceed the midplane normalized magnetic field")
            phi_seed, initial_roundoff_clips = _clip_potential_roundoff_only(wall * (B - B_midplane) / denominator, wall, context="initial floating potential profile")
    if symmetric:
        phi_seed = _symmetric_average(phi_seed)

    floating_active = electron_midplane_prescribed is None
    average_throat = wall
    if floating_active:
        initial_zero_mid = confined_charge_density(field_value=1.0, potential_drop_J=0.0, throat_drop_J=average_throat) + eq63_midplane_charge_density(0.0) + background_midplane
        parent_n0 = max(float(initial_zero_mid) / electron_fraction_zero, np.finfo(float).tiny)
    else:
        parent_n0 = float(electron_midplane_prescribed) / electron_fraction_zero
    throat_left = wall
    throat_right = wall
    scalar_history: list[Mapping[str, object]] = []
    scalar_converged = False
    zero_reference_kind = "benchmark_prescribed_midplane_density" if not floating_active else "not_evaluated"
    zero_reference_zeta: float | None = 0.0 if not floating_active else None
    zero_reference_B: float | None = 1.0 if not floating_active else None
    # Couple the Eq 70 parent Maxwellian amplitude to left and right throat barrier roots
    for iteration_index in range(1, n_iter + 1):
        old_parent = parent_n0
        old_left = throat_left
        old_right = throat_right
        if floating_active:
            zero_cell = np.asarray([
                confined_charge_density(field_value=float(B[index]), potential_drop_J=0.0, throat_drop_J=float(throat_left if z[index] < 0.0 else throat_right))
                + eq63_cell_charge_density(index, 0.0)
                + float(background[index])
                for index in range(n_z)
            ], dtype=float)
            zero_mid = confined_charge_density(field_value=1.0, potential_drop_J=0.0, throat_drop_J=0.5 * (throat_left + throat_right)) + eq63_midplane_charge_density(0.0) + background_midplane
            peak_index = int(np.argmax(zero_cell))
            if float(zero_mid) > float(zero_cell[peak_index]):
                raw_parent = float(zero_mid) / electron_fraction_zero
                zero_reference_kind = "exact_midplane"
                zero_reference_zeta = 0.0
                zero_reference_B = 1.0
                zero_reference_density = float(zero_mid)
            else:
                raw_parent = float(zero_cell[peak_index]) / electron_fraction_zero
                zero_reference_kind = "cell_center"
                zero_reference_zeta = float(z[peak_index])
                zero_reference_B = float(B[peak_index])
                zero_reference_density = float(zero_cell[peak_index])
            parent_n0 = (1.0 - relax) * old_parent + relax * raw_parent
        else:
            zero_reference_density = float(electron_midplane_prescribed)
        raw_left, left_valid, left_search, left_quality = solve_throat("left", parent_n0, min(old_left, wall), background_left)
        if symmetric and np.isclose(background_left, background_right, rtol=0.0, atol=0.0):
            raw_right = raw_left
            right_valid = left_valid
            right_search = left_search
            right_quality = left_quality
        else:
            raw_right, right_valid, right_search, right_quality = solve_throat("right", parent_n0, min(old_right, wall), background_right)
        throat_left = (1.0 - relax) * old_left + relax * raw_left
        throat_right = (1.0 - relax) * old_right + relax * raw_right
        state_change = max(
            abs(parent_n0 - old_parent) / max(abs(parent_n0), abs(old_parent), 1.0),
            abs(throat_left - old_left) / max(abs(throat_left), abs(old_left), Te),
            abs(throat_right - old_right) / max(abs(throat_right), abs(old_right), Te),
        )
        scalar_history.append({
            "iteration": int(iteration_index),
            "parent_n0_m3": float(parent_n0),
            "left_throat_potential_energy_J": float(throat_left),
            "right_throat_potential_energy_J": float(throat_right),
            "zero_reference_kind": zero_reference_kind,
            "zero_reference_zeta": zero_reference_zeta,
            "zero_reference_B_tilde": zero_reference_B,
            "zero_reference_positive_charge_density_m3": float(zero_reference_density),
            "left_throat_root_valid": bool(left_valid),
            "right_throat_root_valid": bool(right_valid),
            "left_throat_root_quality": float(left_quality),
            "right_throat_root_quality": float(right_quality),
            "left_throat_scan_used": bool(left_search.scan_used),
            "right_throat_scan_used": bool(right_search.scan_used),
            "maximum_relative_state_change": float(state_change),
        })
        scalar_converged = bool(left_valid and right_valid and state_change <= tolerance)

    # Reconcile the floating zero drop reference exactly after the relaxed scalar iteration
    final_reference_relative_correction = 0.0
    if floating_active:
        pre_reconciliation_parent_n0 = float(parent_n0)
        final_zero_cell = np.asarray([
            confined_charge_density(field_value=float(B[index]), potential_drop_J=0.0, throat_drop_J=float(throat_left if z[index] < 0.0 else throat_right))
            + eq63_cell_charge_density(index, 0.0)
            + float(background[index])
            for index in range(n_z)
        ], dtype=float)
        final_zero_mid = confined_charge_density(field_value=1.0, potential_drop_J=0.0, throat_drop_J=0.5 * (throat_left + throat_right)) + eq63_midplane_charge_density(0.0) + background_midplane
        final_peak_index = int(np.argmax(final_zero_cell))
        if float(final_zero_mid) > float(final_zero_cell[final_peak_index]):
            final_raw_parent = float(final_zero_mid) / electron_fraction_zero
            zero_reference_kind = "exact_midplane"
            zero_reference_zeta = 0.0
            zero_reference_B = 1.0
            zero_reference_density = float(final_zero_mid)
        else:
            final_raw_parent = float(final_zero_cell[final_peak_index]) / electron_fraction_zero
            zero_reference_kind = "cell_center"
            zero_reference_zeta = float(z[final_peak_index])
            zero_reference_B = float(B[final_peak_index])
            zero_reference_density = float(final_zero_cell[final_peak_index])
        parent_n0 = max(float(final_raw_parent), np.finfo(float).tiny)
        final_reference_relative_correction = abs(parent_n0 - pre_reconciliation_parent_n0) / max(abs(parent_n0), abs(pre_reconciliation_parent_n0), 1.0)
        scalar_history[-1] = {
            **scalar_history[-1],
            "pre_exact_reference_parent_n0_m3": float(pre_reconciliation_parent_n0),
            "parent_n0_m3": float(parent_n0),
            "zero_reference_kind": zero_reference_kind,
            "zero_reference_zeta": zero_reference_zeta,
            "zero_reference_B_tilde": zero_reference_B,
            "zero_reference_positive_charge_density_m3": float(zero_reference_density),
            "final_exact_reference_reconciliation_applied": True,
            "final_exact_reference_relative_correction": float(final_reference_relative_correction),
            "maximum_relative_state_change": max(float(scalar_history[-1]["maximum_relative_state_change"]), float(final_reference_relative_correction)),
        }

    throat_drop_left = float(throat_left)
    throat_drop_right = float(throat_right)
    throat_drop = 0.5 * (throat_drop_left + throat_drop_right)
    _, final_left_electron, final_left_ion = throat_state("left", throat_drop_left, parent_n0, background_left)
    final_left_quality = normalized_root_residual(final_left_electron, final_left_ion, parent_n0)
    left_throat_valid = bool(left_search.available and final_left_quality <= tolerance)
    if symmetric and np.isclose(background_left, background_right, rtol=0.0, atol=0.0):
        final_right_quality = final_left_quality
        right_throat_valid = left_throat_valid
    else:
        _, final_right_electron, final_right_ion = throat_state("right", throat_drop_right, parent_n0, background_right)
        final_right_quality = normalized_root_residual(final_right_electron, final_right_ion, parent_n0)
        right_throat_valid = bool(right_search.available and final_right_quality <= tolerance)
    scalar_history[-1] = {
        **scalar_history[-1],
        "left_throat_root_valid": bool(left_throat_valid),
        "right_throat_root_valid": bool(right_throat_valid),
        "left_throat_root_quality": float(final_left_quality),
        "right_throat_root_quality": float(final_right_quality),
    }
    scalar_converged = bool(left_throat_valid and right_throat_valid and float(scalar_history[-1]["maximum_relative_state_change"]) <= tolerance)

    def local_root_state(index: int, candidate: float) -> tuple[float, float, float]:
        """Evaluate the complete local quasineutrality state for one cell and candidate potential"""
        confined = confined_charge_density(field_value=float(B[index]), potential_drop_J=float(candidate), throat_drop_J=float(throat_drop_left if z[index] < 0.0 else throat_drop_right))
        ion = confined + eq63_cell_charge_density(index, float(candidate)) + float(background[index])
        electron = parent_n0 * electron_fraction(float(candidate))
       
        return electron - ion, electron, ion

    phi = np.asarray(phi_seed, dtype=float).copy()
    central_root_available: list[bool] = []
    scan_cells = 0
    multiple_root_cells = 0
    maximum_root_count = 0
    # Follow the continuation profile first and scan the bounded interval only when needed
    for index in range(n_z):
        state_cache: dict[float, tuple[float, float, float]] = {}

        def state(candidate: float, index: int = index) -> tuple[float, float, float]:
            """Return cached cell residual, electron density, and ion charge density"""
            key = float(candidate)
            if key not in state_cache:
                state_cache[key] = local_root_state(index, key)
        
            return state_cache[key]

        def residual(candidate: float) -> float:
            """Return the cached cell quasineutrality residual"""
            return state(candidate)[0]

        search = _bounded_continuation_root(residual=residual, lower=0.0, upper=wall, seed=float(phi[index]), scan_points=root_scan_count, xtol=root_xtol, rtol=root_rtol)
        phi[index] = float(search.root)
        _, electron, ion = state(float(search.root))
        quality_valid = bool(search.available and normalized_root_residual(electron, ion, parent_n0) <= tolerance)
        central_root_available.append(bool(quality_valid))
        scan_cells += int(search.scan_used)
        multiple_root_cells += int(search.root_count > 1)
        maximum_root_count = max(maximum_root_count, int(search.root_count))
    phi, cell_clip_count = _clip_potential_roundoff_only(phi, wall, context="floating direct coupled potential profile")
    roundoff_clip_count = int(initial_roundoff_clips) + int(cell_clip_count)
    central_root_diagnostics = {
        "no_root_cell_count": int(sum(not value for value in central_root_available)),
        "interior_scan_cell_count": int(scan_cells),
        "multiple_root_cell_count": int(multiple_root_cells),
        "maximum_root_count": int(maximum_root_count),
    }

    def exact_midplane_state(candidate: float) -> tuple[float, float, float]:
        """Evaluate quasineutrality at exact ζ = 0 independently of cell center roots"""
        confined = confined_charge_density(field_value=1.0, potential_drop_J=float(candidate), throat_drop_J=throat_drop)
        ion = confined + eq63_midplane_charge_density(float(candidate)) + background_midplane
        electron = parent_n0 * electron_fraction(float(candidate))
      
        return electron - ion, electron, ion

    midplane_state_cache: dict[float, tuple[float, float, float]] = {}

    def cached_exact_midplane_state(candidate: float) -> tuple[float, float, float]:
        """Cache exact midplane states during the bounded root search"""
        key = float(candidate)
        if key not in midplane_state_cache:
            midplane_state_cache[key] = exact_midplane_state(key)
      
        return midplane_state_cache[key]

    # The exact midplane root is solved separately because ζ = 0 need not be a cell center
    midplane_seed = float(np.interp(0.0, z[np.argsort(z)], phi[np.argsort(z)]))
    midplane_search = _bounded_continuation_root(residual=lambda value: cached_exact_midplane_state(value)[0], lower=0.0, upper=wall, seed=midplane_seed, scan_points=root_scan_count, xtol=root_xtol, rtol=root_rtol)
    exact_midplane_potential = float(midplane_search.root)
    _, mid_electron, mid_ion = cached_exact_midplane_state(exact_midplane_potential)
    direct_midplane_root_valid = bool(midplane_search.available and normalized_root_residual(mid_electron, mid_ion, parent_n0) <= tolerance)

    densities: dict[str, np.ndarray] = {}
    local_lambda: dict[str, np.ndarray] = {}
    local_pitch: dict[str, np.ndarray] = {}
    mapping: dict[str, dict[str, float | int]] = {}
    low_energy_sample_count = 0
    nonpositive_energy_sample_count = 0
    closed_boundary_sample_count = 0
    for species_id, value in inputs.items():
        local_grid = local_speed_grids[species_id]
        species_lambda = np.zeros((n_z, local_grid.centers_m_s.size, value.lambda_grid.centers.size), dtype=float)
        species_pitch = np.zeros((n_z, local_grid.centers_m_s.size, value.pitch_grid.centers.size), dtype=float)
        species_density = np.zeros(n_z, dtype=float)
        totals: dict[str, float | int] = {"roundoff_clip_count": 0, "throat_extrapolation_clip_count": 0, "low_energy_approximation_sample_count": 0, "nonpositive_total_energy_sample_count": 0, "closed_eq71_boundary_sample_count": 0}
        for index in range(n_z):
            species_lambda[index], species_pitch[index], item_mapping = _phi_corrected_local_pitch_distribution(
                speed_grid=value.speed_grid,
                local_speed_grid=local_grid,
                lambda_grid=value.lambda_grid,
                pitch_grid=value.pitch_grid,
                base_distribution_v_lambda=value.base_distribution_v_lambda,
                mirror_ratio=R,
                B_tilde=float(B[index]),
                local_potential_drop_magnitude_J=float(phi[index]),
                throat_potential_drop_magnitude_J=float(throat_drop_left if z[index] < 0.0 else throat_drop_right),
                eta_to_local_phase_space_normalization=value.eta_to_local_phase_space_normalization,
                quadrature_order=int(local_velocity_quadrature_order),
                base_distribution_interpolator=interpolators[species_id],
                quadrature_state=local_quadrature_by_species[species_id],
                particle_mass_kg=value.particle_mass_kg,
            )
            species_density[index] = _density_from_local_pitch(local_grid, value.pitch_grid, species_pitch[index])
            for key in totals:
                totals[key] = int(totals[key]) + int(item_mapping.get(key, 0))
        if symmetric:
            species_density = _symmetric_average(species_density)
            species_lambda = 0.5 * (species_lambda + species_lambda[::-1])
            species_pitch = 0.5 * (species_pitch + species_pitch[::-1])
        densities[species_id] = species_density
        local_lambda[species_id] = species_lambda
        local_pitch[species_id] = species_pitch
        mapping[species_id] = totals
        roundoff_clip_count += int(totals["roundoff_clip_count"])
        low_energy_sample_count += int(totals["low_energy_approximation_sample_count"])
        nonpositive_energy_sample_count += int(totals["nonpositive_total_energy_sample_count"])
        closed_boundary_sample_count += int(totals["closed_eq71_boundary_sample_count"])

    confined_fast_number_density = np.zeros(n_z, dtype=float)
    confined_positive_charge_density = np.zeros(n_z, dtype=float)
    for species_id, density in densities.items():
        confined_fast_number_density += density
        confined_positive_charge_density += float(inputs[species_id].charge_number) * density
    lost_positive_charge_density = np.asarray([eq63_cell_charge_density(index, float(phi[index])) for index in range(n_z)], dtype=float)
    ion_density = confined_positive_charge_density + lost_positive_charge_density + background
    electron_density = np.asarray([parent_n0 * electron_fraction(float(value)) for value in phi], dtype=float)
    residual_metrics = _quasineutrality_residual_metrics(electron_density_m3=electron_density, ion_density_m3=ion_density, zeta=z, cell_volumes_m3=volumes, reference_density_m3=max(parent_n0, np.finfo(float).tiny))
    absolute_residual = np.asarray(residual_metrics["signed_absolute_density_residual_m3"], dtype=float)
    relative_residual = np.asarray(residual_metrics["supported_relative_residual"], dtype=float)
    max_relative = float(residual_metrics["maximum_supported_relative_residual"])
    maximum_absolute_normalized = float(residual_metrics["maximum_absolute_residual_normalized_to_reference"])
    volume_integrated_absolute_fraction = float(residual_metrics["volume_integrated_absolute_particle_mismatch_fraction"])
    final_residual = max(max_relative, maximum_absolute_normalized, volume_integrated_absolute_fraction)
    residual_history = [float(final_residual)]
    n0_history = [float(item["parent_n0_m3"]) for item in scalar_history]

    exact_midplane_by_species = {key: species_confined_density(key, field_value=1.0, potential_drop_J=exact_midplane_potential, throat_drop_J=throat_drop) for key in inputs}
    exact_left_by_species = {key: species_confined_density(key, field_value=R, potential_drop_J=throat_drop_left, throat_drop_J=throat_drop_left) for key in inputs}
    exact_right_by_species = {key: species_confined_density(key, field_value=R, potential_drop_J=throat_drop_right, throat_drop_J=throat_drop_right) for key in inputs}
    exact_midplane_lost = eq63_midplane_charge_density(exact_midplane_potential)
    exact_left_lost = eq63_throat_charge_density("left", throat_drop_left)
    exact_right_lost = eq63_throat_charge_density("right", throat_drop_right)
    exact_midplane_ion = float(sum(float(inputs[key].charge_number) * exact_midplane_by_species[key] for key in inputs)) + exact_midplane_lost + background_midplane
    exact_left_ion = float(sum(float(inputs[key].charge_number) * exact_left_by_species[key] for key in inputs)) + exact_left_lost + background_left
    exact_right_ion = float(sum(float(inputs[key].charge_number) * exact_right_by_species[key] for key in inputs)) + exact_right_lost + background_right
    node_profiles = _exact_electrostatic_node_profiles(
        cell_zeta=z,
        cell_B_tilde=B,
        cell_potential_energy_J=phi,
        cell_electron_density_m3=electron_density,
        cell_ion_density_m3=ion_density,
        cell_fast_ion_density_m3=confined_fast_number_density,
        cell_background_density_m3=background,
        mirror_ratio=R,
        left_throat_potential_energy_J=throat_drop_left,
        right_throat_potential_energy_J=throat_drop_right,
        electron_parent_maxwellian_n0_m3=parent_n0,
        wall_barrier_energy_J=wall,
        electron_temperature_J=Te,
        exact_midplane_fast_density_m3=float(sum(exact_midplane_by_species.values())),
        exact_midplane_potential_energy_J=exact_midplane_potential,
        exact_left_throat_fast_density_m3=float(sum(exact_left_by_species.values())),
        exact_right_throat_fast_density_m3=float(sum(exact_right_by_species.values())),
        exact_midplane_background_density_m3=background_midplane,
        exact_left_throat_background_density_m3=background_left,
        exact_right_throat_background_density_m3=background_right,
        exact_midplane_ion_density_m3=exact_midplane_ion,
        exact_left_throat_ion_density_m3=exact_left_ion,
        exact_right_throat_ion_density_m3=exact_right_ion,
    )
    node_zeta = np.asarray(node_profiles["zeta"], dtype=float)
    left_index = int(node_profiles["left_throat_index"])
    midplane_index = int(node_profiles["midplane_index"])
    right_index = int(node_profiles["right_throat_index"])
    node_density_by_species: dict[str, np.ndarray] = {}
    order = np.argsort(z)
    for species_id, density in densities.items():
        node_density = np.interp(node_zeta, z[order], np.asarray(density, dtype=float)[order])
        node_density[[left_index, midplane_index, right_index]] = [exact_left_by_species[species_id], exact_midplane_by_species[species_id], exact_right_by_species[species_id]]
        node_density_by_species[species_id] = node_density
    node_relative = np.asarray(node_profiles["relative_residual_profile"], dtype=float)
    node_absolute = np.asarray(node_profiles["absolute_residual_profile_m3"], dtype=float)
    node_electron = np.asarray(node_profiles["electron_density_m3"], dtype=float)
    node_ion = np.asarray(node_profiles["ion_density_m3"], dtype=float)
    anchor_indices = np.asarray([left_index, midplane_index, right_index], dtype=int)
    anchor_scale = np.maximum(node_electron[anchor_indices], node_ion[anchor_indices])
    anchor_supported = anchor_scale >= float(residual_metrics["density_support_floor_m3"])
    anchor_normalized_absolute = np.abs(node_absolute[anchor_indices]) / float(residual_metrics["reference_density_m3"])
    anchor_residual = np.where(anchor_supported, np.abs(node_relative[anchor_indices]), anchor_normalized_absolute)
    exact_left_residual = float(anchor_residual[0])
    exact_midplane_residual = float(anchor_residual[1])
    exact_right_residual = float(anchor_residual[2])
    exact_throat_residual = max(exact_left_residual, exact_right_residual)
    exact_midplane_valid = bool(direct_midplane_root_valid and exact_midplane_residual <= tolerance)
    exact_throat_valid = bool(left_throat_valid and right_throat_valid and exact_throat_residual <= tolerance)

    effective_by_species = {
        species_id: _effective_potential_throat_boundary_diagnostics(
            speed_grid=value.speed_grid,
            lambda_grid=value.lambda_grid,
            base_distribution_v_lambda=value.base_distribution_v_lambda,
            zeta=z,
            B_tilde=B,
            potential_drop_magnitude_J=phi,
            mirror_ratio=R,
            throat_potential_left_energy_J=throat_drop_left,
            throat_potential_right_energy_J=throat_drop_right,
            particle_mass_kg=value.particle_mass_kg,
        )
        for species_id, value in inputs.items()
    }
    closed_by_species = {
        species_id: _closed_eq71_modal_inventory_fraction(
            speed_grid=value.speed_grid,
            mirror_ratio=R,
            throat_potential_drop_magnitude_J=throat_drop,
            invariant_energy_cell_population_weights=value.invariant_energy_cell_population_weights,
            particle_mass_kg=value.particle_mass_kg,
        )
        for species_id, value in inputs.items()
    }
    low_energy_valid = bool(all(bool(item["assessed"]) and float(item["fraction"]) <= low_energy_tolerance for item in closed_by_species.values()))
    effective_valid = bool(all(bool(item["valid"]) for item in effective_by_species.values()))
    effective_failure = None if effective_valid else ";".join(f"{key}:{value['failure_reason']}" for key, value in effective_by_species.items() if not bool(value["valid"]))

    convergence_failures: list[str] = []
    if max_relative > tolerance:
        convergence_failures.append("maximum_supported_relative_residual_above_tolerance")
    if maximum_absolute_normalized > tolerance:
        convergence_failures.append("maximum_absolute_reference_normalized_residual_above_tolerance")
    if volume_integrated_absolute_fraction > tolerance:
        convergence_failures.append("volume_integrated_absolute_residual_above_tolerance")
    if not exact_midplane_valid:
        convergence_failures.append("exact_midplane_quasineutrality_residual_above_tolerance")
    if not exact_throat_valid:
        convergence_failures.append("exact_throat_quasineutrality_residual_above_tolerance")
    if not all(central_root_available):
        convergence_failures.append("direct_central_quasineutrality_root_unavailable")
    if floating_active and not scalar_converged:
        convergence_failures.append("floating_nonnegative_reference_scalar_closure_not_converged")
    converged = not convergence_failures

    total_volume = float(np.sum(volumes))
    inventory_by_species = {species_id: float(np.sum(density * volumes)) for species_id, density in densities.items()}
    average_by_species = {species_id: value / total_volume for species_id, value in inventory_by_species.items()}
    total_inventory = float(sum(inventory_by_species.values()))
    total_average = total_inventory / total_volume
    target_absolute_error = None if target_volume_average is None else total_average - target_volume_average
    target_relative_error = None if target_volume_average is None or target_volume_average == 0.0 else target_absolute_error / target_volume_average
    actual_midplane_electron = float(node_electron[midplane_index])
    profile = ModalElectrostaticProfile(
        zeta=z,
        B_tilde=B,
        potential_energy_J=phi,
        potential_relative_to_midplane_V=-(phi - exact_midplane_potential) / ELECTRON_CHARGE_C,
        electron_density_m3=electron_density,
        ion_density_m3=ion_density,
        iterations=len(scalar_history),
        max_relative_quasineutrality_error=max_relative,
        converged=converged,
        electron_density_scale_n0_m3=float(parent_n0),
        electron_midplane_density_m3=actual_midplane_electron,
        electron_parent_maxwellian_n0_m3=float(parent_n0),
        electron_collision_density_m3=float(electron_collision if electron_collision is not None else actual_midplane_electron),
        electron_volume_average_density_m3=float(np.sum(electron_density * volumes) / total_volume),
        fast_ion_density_m3=confined_fast_number_density,
        background_positive_charge_density_m3=background,
        electron_density_input_semantics=("floating_parent_Maxwellian_from_maximum_zero_drop_positive_charge" if floating_active else "benchmark_prescribed_actual_total_electron_density_at_magnetic_midplane"),
        electron_density_closure_model=("eq70_floating_nonnegative_drop_reference" if floating_active else "benchmark_prescribed_midplane_density"),
        relative_residual_profile=relative_residual,
        absolute_residual_profile_m3=absolute_residual,
        residual_history=tuple(residual_history),
        electron_density_scale_history_m3=tuple(n0_history),
        stagnated=False,
        oscillatory=False,
        diverged=False,
        failure_reason=None if converged else ";".join(convergence_failures),
        roundoff_clip_count=int(roundoff_clip_count),
        throat_potential_energy_J=float(throat_drop),
        throat_potential_left_energy_J=float(throat_drop_left),
        throat_potential_right_energy_J=float(throat_drop_right),
        low_energy_approximation_sample_count=int(low_energy_sample_count),
        nonpositive_total_energy_sample_count=int(nonpositive_energy_sample_count),
        closed_eq71_boundary_sample_count=int(closed_boundary_sample_count),
        low_energy_approximation_max_weight_fraction=max(float(item["fraction"]) for item in closed_by_species.values()),
        low_energy_approximation_valid=low_energy_valid,
        low_energy_approximation_weight_fraction_assessed=bool(all(bool(item["assessed"]) for item in closed_by_species.values())),
        low_energy_approximation_weight_model="maximum_species_closed_Eq71_inventory_fraction",
        eq71_closed_interval_intersected_population_cell_count=int(sum(int(item["intersected_population_cell_count"]) for item in closed_by_species.values())),
        eq71_closed_interval_intersects_distribution_support=bool(any(bool(item["intersects_distribution_support"]) for item in closed_by_species.values())),
        eta_to_local_phase_space_normalization=1.0,
        phase_space_measure_conversion_factor=1.0,
        phase_space_measure_conversion_history=(1.0,),
        phase_space_measure_conversion_model="species_specific_eq45_to_eq48_eta_measure_to_local_phase_space",
        fast_ion_inventory_particles=total_inventory,
        fast_ion_volume_average_density_m3=total_average,
        target_volume_averaged_fast_ion_density_m3=target_volume_average,
        fast_ion_volume_average_density_absolute_error_m3=target_absolute_error,
        fast_ion_volume_average_density_relative_error=target_relative_error,
        symmetric_input=symmetric,
        midplane_reference_V=0.0,
        throat_extrapolation_valid=bool(left_throat_valid and right_throat_valid),
        effective_potential_check_assessed=bool(all(bool(item["assessed"]) for item in effective_by_species.values())),
        effective_potential_throat_boundary_valid=effective_valid,
        effective_potential_failure_reason=effective_failure,
        effective_potential_active_energy_sample_count=int(sum(int(item["active_energy_sample_count"]) for item in effective_by_species.values())),
        effective_potential_open_energy_sample_count=int(sum(int(item["open_confined_energy_sample_count"]) for item in effective_by_species.values())),
        effective_potential_interior_minimum_energy_sample_count=int(sum(int(item["interior_minimum_energy_sample_count"]) for item in effective_by_species.values())),
        effective_potential_interior_minimum_weight_fraction=max(float(item["interior_minimum_weight_fraction"]) for item in effective_by_species.values()),
        effective_potential_max_relative_boundary_shortfall=max(float(item["max_relative_boundary_shortfall"]) for item in effective_by_species.values()),
        effective_potential_max_relative_throat_asymmetry=max(float(item["max_relative_throat_asymmetry"]) for item in effective_by_species.values()),
        effective_potential_worst_total_energy_J=max(float(item["worst_total_energy_J"]) for item in effective_by_species.values()),
        effective_potential_relative_tolerance=max(float(item["relative_tolerance"]) for item in effective_by_species.values()),
        absolute_residual_normalized_to_reference=np.asarray(residual_metrics["absolute_residual_normalized_to_reference"], dtype=float),
        density_support_mask=np.asarray(residual_metrics["density_support_mask"], dtype=bool),
        density_support_reference_m3=float(residual_metrics["reference_density_m3"]),
        density_support_relative_floor=float(residual_metrics["density_support_relative_floor"]),
        density_support_floor_m3=float(residual_metrics["density_support_floor_m3"]),
        density_support_floor_derivation=str(residual_metrics["density_support_floor_derivation"]),
        density_support_masked_node_count=int(residual_metrics["masked_node_count"]),
        density_support_masked_volume_m3=float(residual_metrics["masked_volume_m3"]),
        density_support_masked_volume_fraction=float(residual_metrics["masked_volume_fraction"]),
        volume_integrated_electron_minus_ion_particles=float(residual_metrics["volume_integrated_electron_minus_ion_particles"]),
        volume_integrated_absolute_particle_mismatch=float(residual_metrics["volume_integrated_absolute_particle_mismatch"]),
        volume_integrated_signed_charge_mismatch_C=float(residual_metrics["volume_integrated_signed_charge_mismatch_C"]),
        volume_integrated_absolute_charge_mismatch_C=float(residual_metrics["volume_integrated_absolute_charge_mismatch_C"]),
        central_region_maximum_supported_relative_residual=float(residual_metrics["central_region_maximum_supported_relative_residual"]),
        throat_region_maximum_supported_relative_residual=float(residual_metrics["throat_region_maximum_supported_relative_residual"]),
        maximum_absolute_density_residual_m3=float(residual_metrics["maximum_absolute_density_residual_m3"]),
        volume_integrated_signed_particle_mismatch_fraction=float(residual_metrics["volume_integrated_signed_particle_mismatch_fraction"]),
        maximum_absolute_density_residual_normalized_to_reference=maximum_absolute_normalized,
        volume_integrated_absolute_particle_mismatch_fraction=volume_integrated_absolute_fraction,
        electrostatic_node_zeta=node_zeta,
        electrostatic_node_B_tilde=np.asarray(node_profiles["B_tilde"], dtype=float),
        electrostatic_node_potential_energy_J=np.asarray(node_profiles["potential_energy_J"], dtype=float),
        electrostatic_node_potential_relative_to_midplane_V=np.asarray(node_profiles["potential_relative_to_midplane_V"], dtype=float),
        electrostatic_node_electron_density_m3=node_electron,
        electrostatic_node_ion_density_m3=node_ion,
        electrostatic_node_fast_ion_density_m3=np.asarray(node_profiles["fast_ion_density_m3"], dtype=float),
        electrostatic_node_background_positive_charge_density_m3=np.asarray(node_profiles["background_positive_charge_density_m3"], dtype=float),
        electrostatic_node_absolute_residual_profile_m3=node_absolute,
        electrostatic_node_relative_residual_profile=node_relative,
        electrostatic_left_throat_index=left_index,
        electrostatic_midplane_index=midplane_index,
        electrostatic_right_throat_index=right_index,
        electrostatic_left_throat_node_exact=True,
        electrostatic_midplane_node_exact=True,
        electrostatic_right_throat_node_exact=True,
        electrostatic_midplane_gauge_exact=bool(exact_midplane_potential == 0.0),
        electrostatic_node_model=str(node_profiles["node_model"]),
        floating_nonnegative_reference_active=bool(floating_active),
        floating_zero_reference_kind=zero_reference_kind,
        floating_zero_reference_zeta=zero_reference_zeta,
        floating_zero_reference_B_tilde=zero_reference_B,
        exact_midplane_potential_energy_J=float(exact_midplane_potential),
        floating_scalar_closure_converged=bool(scalar_converged if floating_active else True),
        floating_scalar_closure_history=tuple(scalar_history),
        direct_exact_midplane_root_valid=bool(direct_midplane_root_valid),
        local_speed_grid=None if len(local_speed_grids) != 1 else next(iter(local_speed_grids.values())),
        exact_anchor_maximum_quasineutrality_residual=max(exact_midplane_residual, exact_throat_residual),
        exact_anchor_quasineutrality_valid=bool(exact_midplane_valid and exact_throat_valid),
        exact_anchor_residual_support_model="supported_relative_else_absolute_normalized_to_parent_Maxwellian_density",
        exact_midplane_quasineutrality_residual=exact_midplane_residual,
        exact_left_throat_quasineutrality_residual=exact_left_residual,
        exact_right_throat_quasineutrality_residual=exact_right_residual,
        exact_throat_maximum_quasineutrality_residual=exact_throat_residual,
        exact_midplane_quasineutrality_valid=exact_midplane_valid,
        exact_throat_quasineutrality_valid=exact_throat_valid,
        direct_left_throat_root_valid=bool(left_throat_valid),
        direct_right_throat_root_valid=bool(right_throat_valid),
        direct_throat_roots_valid=bool(left_throat_valid and right_throat_valid),
        direct_central_roots_valid=bool(all(central_root_available)),
        direct_central_no_root_cell_count=int(central_root_diagnostics["no_root_cell_count"]),
        direct_central_interior_scan_cell_count=int(central_root_diagnostics["interior_scan_cell_count"]),
        direct_central_multiple_root_cell_count=int(central_root_diagnostics["multiple_root_cell_count"]),
        direct_central_maximum_root_count=int(central_root_diagnostics["maximum_root_count"]),
        direct_central_root_scan_points=int(root_scan_count),
        fast_ion_density_by_species_m3={key: np.asarray(value, dtype=float) for key, value in densities.items()},
        electrostatic_node_fast_ion_density_by_species_m3=node_density_by_species,
        local_speed_grid_by_species=local_speed_grids,
        fast_ion_inventory_particles_by_species=inventory_by_species,
        fast_ion_volume_average_density_m3_by_species=average_by_species,
    )
    
    if audit_started is not None and audit is not None:
        audit["solver_total_s"] = float(perf_counter() - audit_started)
        tracker_summary: dict[str, object] = {}
        if interpolation_tracker is not None:
            for key, raw in interpolation_tracker.items():
                if not isinstance(raw, Mapping):
                    continue
                seen = raw.get("index_seen")
                calls = int(raw.get("index_calls", 0))
                repeats = int(raw.get("index_repeat_calls", 0))
                tracker_summary[str(key)] = {
                    "index_calls": calls,
                    "unique_index_patterns": len(seen) if isinstance(seen, set) else None,
                    "index_repeat_calls": repeats,
                    "index_repeat_fraction": None if calls <= 0 else float(repeats / calls),
                    "empty_calls": int(raw.get("empty_calls", 0)),
                }
        audit["interpolation_index_recurrence"] = tracker_summary
        audit["output"] = {
            "converged": bool(profile.converged),
            "direct_central_no_root_cell_count": int(profile.direct_central_no_root_cell_count),
            "direct_central_roots_valid": bool(profile.direct_central_roots_valid),
            "direct_exact_midplane_root_valid": bool(profile.direct_exact_midplane_root_valid),
            "direct_throat_roots_valid": bool(profile.direct_throat_roots_valid),
        }

    return ModalElectrostaticSystemResult(
        profile=profile,
        local_distribution_z_v_lambda_by_species=local_lambda,
        local_distribution_z_v_pitch_by_species=local_pitch,
        local_density_m3_by_species=densities,
        local_speed_grid_by_species=local_speed_grids,
        mapping_diagnostics_by_species=mapping,
    )
