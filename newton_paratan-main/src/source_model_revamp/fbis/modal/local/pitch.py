"""Local pitch reconstruction from the invariant confined ion distribution"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid, lambda_grid_from_faces
from source_model_revamp.fbis.mirror_pitch_remap import conservative_lambda_distribution_from_local_pitch_distribution, conservative_local_pitch_distribution_from_lambda_distribution
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid, gyrotropic_velocity_cell_volumes
from source_model_revamp.fbis.modal.local.mapping import _base_distribution_interpolator, _cell_quadrature, _eq71_closed_interval_threshold_J, _eq71_global_trapped_passing_boundary, _eq72_compressed_invariant, _ion_invariants_from_local_coordinates, _validated_base_distribution

def conservative_speed_remap_gyrotropic_distribution(*, source_speed_grid: SpeedGrid, target_speed_grid: SpeedGrid, pitch_grid: PitchGrid, source_distribution_v_pitch: np.ndarray, particle_mass_kg: float) -> tuple[np.ndarray, dict[str, float | bool | str]]:
    """
    Conservatively remap piecewise constant speed cells in v^2 dv at fixed pitch

    Particle measure is conserved when the target domain covers the source domain
    Kinetic energy is reported as a remap diagnostic and is not forced to match
    """
    source = np.asarray(source_distribution_v_pitch, dtype=float)
    expected = (source_speed_grid.centers_m_s.size, pitch_grid.centers.size)
    if source.shape != expected:
        raise ValueError(f"source_distribution_v_pitch must have shape {expected}")
    if np.any(~np.isfinite(source)):
        raise ValueError("source_distribution_v_pitch must contain only finite values")
    scale = max(float(np.max(np.abs(source))) if source.size else 0.0, np.finfo(float).tiny)
    if float(np.min(source)) < -128.0 * np.finfo(float).eps * scale:
        raise ValueError("source_distribution_v_pitch contains a significant negative value")
    source = np.maximum(source, 0.0)
    mass = float(particle_mass_kg)
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    source_faces = np.asarray(source_speed_grid.faces_m_s, dtype=float)
    target_faces = np.asarray(target_speed_grid.faces_m_s, dtype=float)
    source_lo = source_faces[:-1]
    source_hi = source_faces[1:]
    target_lo = target_faces[:-1]
    target_hi = target_faces[1:]
    overlap_lo = np.maximum(target_lo[:, None], source_lo[None, :])
    overlap_hi = np.minimum(target_hi[:, None], source_hi[None, :])
    overlap_v3 = np.maximum(overlap_hi**3 - overlap_lo**3, 0.0) / 3.0
    target_v3 = (target_hi**3 - target_lo**3) / 3.0
    numerator = overlap_v3 @ source
    target = np.divide(numerator, target_v3[:, None], out=np.zeros((target_lo.size, pitch_grid.centers.size), dtype=float), where=target_v3[:, None] > 0.0)
    pitch_measure = 2.0 * np.pi * np.asarray(pitch_grid.widths, dtype=float)
    source_particle = float(np.sum(source * ((source_hi**3 - source_lo**3) / 3.0)[:, None] * pitch_measure[None, :]))
    target_particle = float(np.sum(target * target_v3[:, None] * pitch_measure[None, :]))
    source_v5 = (source_hi**5 - source_lo**5) / 5.0
    target_v5 = (target_hi**5 - target_lo**5) / 5.0
    source_energy = float(np.sum(source * source_v5[:, None] * pitch_measure[None, :]) * 0.5 * mass)
    target_energy = float(np.sum(target * target_v5[:, None] * pitch_measure[None, :]) * 0.5 * mass)
    covered = bool(target_faces[0] <= source_faces[0] + 128.0 * np.finfo(float).eps * max(source_faces[-1], 1.0) and target_faces[-1] >= source_faces[-1] - 128.0 * np.finfo(float).eps * max(source_faces[-1], 1.0))
    diagnostics: dict[str, float | bool | str] = {"source_particle_measure": source_particle, "target_particle_measure": target_particle, "particle_inventory_relative_error": ((target_particle - source_particle) / source_particle if source_particle > 0.0 else 0.0), "source_kinetic_energy_measure_J": source_energy, "target_kinetic_energy_measure_J": target_energy, "kinetic_energy_relative_error": ((target_energy - source_energy) / source_energy if source_energy > 0.0 else 0.0), "target_domain_covers_source_domain": covered, "remap_model": "piecewise_constant_exact_overlap_in_v_cubed_measure"}

    return target, diagnostics

@dataclass(frozen=True)
class _LocalPitchQuadratureState:
    """Cached cell quadrature arrays for local speed and pitch integration"""
    speed_faces_m_s: np.ndarray
    pitch_faces: np.ndarray
    order: int
    speed_nodes_m_s: np.ndarray
    pitch_nodes: np.ndarray
    quadrature_weight: np.ndarray
    denominator: np.ndarray

def _prepare_local_pitch_quadrature(local_speed_grid: SpeedGrid, pitch_grid: PitchGrid, order: int) -> _LocalPitchQuadratureState:
    """Prepare local speed and pitch quadrature including the v² velocity measure"""
    count = int(order)
    if count < 2:
        raise ValueError("quadrature_order must be at least two")
    speed_faces = np.asarray(local_speed_grid.faces_m_s, dtype=float)
    pitch_faces = np.asarray(pitch_grid.faces, dtype=float)
    speed_nodes, speed_weights = _cell_quadrature(speed_faces, count)
    pitch_nodes, pitch_weights = _cell_quadrature(pitch_faces, count)
    local_speed = speed_nodes[:, None, :, None]
    local_pitch = pitch_nodes[None, :, None, :]
    quadrature_weight = speed_weights[:, None, :, None] * pitch_weights[None, :, None, :] * local_speed**2
    speed_measure = (speed_faces[1:] ** 3 - speed_faces[:-1] ** 3) / 3.0
    # Denominator is ∫v² dv dξ for each local velocity cell without the common 2π factor
    denominator = speed_measure[:, None] * np.asarray(pitch_grid.widths, dtype=float)[None, :]
   
    return _LocalPitchQuadratureState(
        speed_faces_m_s=speed_faces,
        pitch_faces=pitch_faces,
        order=count,
        speed_nodes_m_s=local_speed,
        pitch_nodes=local_pitch,
        quadrature_weight=quadrature_weight,
        denominator=denominator,
    )

def _phi_corrected_local_coordinates(*, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, pitch_grid: PitchGrid, base_distribution_v_lambda: np.ndarray, mirror_ratio: float, B_tilde: float, local_potential_drop_magnitude_J: float, throat_potential_drop_magnitude_J: float, eta_to_local_phase_space_normalization: float, particle_mass_kg: float, quadrature_order: int, local_speed_grid: SpeedGrid | None, base_distribution_interpolator: RegularGridInterpolator | None, quadrature_state: _LocalPitchQuadratureState | None) -> tuple[SpeedGrid, _LocalPitchQuadratureState, RegularGridInterpolator, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Map local quadrature nodes (v, ξ) to U, Λ, the Eq 71 boundary and Eq 72 Λ*

    The returned validity mask selects positive U states accessible inside the confined interval
    """
    B = float(B_tilde)
    R = float(mirror_ratio)
    local_drop = float(local_potential_drop_magnitude_J)
    throat_drop = float(throat_potential_drop_magnitude_J)
    normalization = float(eta_to_local_phase_space_normalization)
    if not np.isfinite(B) or B <= 0.0:
        raise ValueError("B_tilde must be positive and finite")
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    if not np.isfinite(local_drop) or local_drop < 0.0:
        raise ValueError("local_potential_drop_magnitude_J must be finite and nonnegative")
    if not np.isfinite(throat_drop) or throat_drop < 0.0:
        raise ValueError("throat_potential_drop_magnitude_J must be finite and nonnegative")
    if not np.isfinite(normalization) or normalization <= 0.0:
        raise ValueError("eta_to_local_phase_space_normalization must be positive and finite")
    if not np.isfinite(float(particle_mass_kg)) or float(particle_mass_kg) <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    physical_speed_grid = speed_grid if local_speed_grid is None else local_speed_grid
    quadrature = _prepare_local_pitch_quadrature(physical_speed_grid, pitch_grid, int(quadrature_order)) if quadrature_state is None else quadrature_state
    if quadrature.order != int(quadrature_order) or not np.array_equal(quadrature.speed_faces_m_s, np.asarray(physical_speed_grid.faces_m_s, dtype=float)) or not np.array_equal(quadrature.pitch_faces, np.asarray(pitch_grid.faces, dtype=float)):
        raise ValueError("local pitch quadrature state does not match the requested grids and order")
    interpolator = _base_distribution_interpolator(speed_grid, lambda_grid, base_distribution_v_lambda) if base_distribution_interpolator is None else base_distribution_interpolator
    total_energy, global_lambda = _ion_invariants_from_local_coordinates(local_speed_m_s=quadrature.speed_nodes_m_s, local_pitch=quadrature.pitch_nodes, B_tilde=B, local_potential_drop_magnitude_J=local_drop, particle_mass_kg=particle_mass_kg)
    total_energy_full = np.broadcast_to(total_energy, global_lambda.shape)
    boundary = np.broadcast_to(_eq71_global_trapped_passing_boundary(total_energy, R, throat_drop), global_lambda.shape)
    lambda_star = _eq72_compressed_invariant(global_lambda, boundary, R)
    roundoff_tolerance = 128.0 * np.finfo(float).eps
    positive_energy = total_energy_full > 0.0
    # Sample only states whose Eq 71 interval is open and whose Eq 72 image is magnetically confined
    valid = positive_energy & (boundary < 1.0) & (global_lambda >= boundary - roundoff_tolerance) & (global_lambda <= 1.0 + roundoff_tolerance) & (lambda_star >= 1.0 / R - roundoff_tolerance) & (lambda_star <= 1.0 + roundoff_tolerance)
   
    return physical_speed_grid, quadrature, interpolator, normalization, total_energy_full, global_lambda, boundary, lambda_star, valid

def _sample_phi_corrected_local_distribution(*, interpolator: RegularGridInterpolator, normalization: float, total_energy_full: np.ndarray, lambda_star: np.ndarray, valid: np.ndarray, mirror_ratio: float, particle_mass_kg: float) -> np.ndarray:
    """Sample magnetic reference f(v, Λ*) at v = sqrt(2U / m) and apply the phase space normalization"""
    sampled = np.zeros_like(lambda_star)
    if np.any(valid):
        points = np.column_stack((_invariant_speed_from_total_energy(total_energy_full[valid], float(particle_mass_kg)), np.clip(lambda_star[valid], 1.0 / float(mirror_ratio), 1.0)))
        sampled[valid] = interpolator(points)
  
    return np.maximum(sampled * float(normalization), 0.0)

def _invariant_speed_from_total_energy(total_energy_J: np.ndarray, particle_mass_kg: float) -> np.ndarray:
    """Return invariant speed v = sqrt(2U / m) for positive total energy samples"""
  
    return np.sqrt(2.0 * np.asarray(total_energy_J, dtype=float) / float(particle_mass_kg))

def _phi_corrected_local_density(*, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, pitch_grid: PitchGrid, base_distribution_v_lambda: np.ndarray, mirror_ratio: float, B_tilde: float, local_potential_drop_magnitude_J: float, throat_potential_drop_magnitude_J: float, eta_to_local_phase_space_normalization: float, particle_mass_kg: float, quadrature_order: int = 3, local_speed_grid: SpeedGrid | None = None, base_distribution_interpolator: RegularGridInterpolator | None = None, quadrature_state: _LocalPitchQuadratureState | None = None) -> float:
    """Evaluate n = 2π ∫ v² dv dξ f directly from local quadrature samples"""
    physical_speed_grid, quadrature, interpolator, normalization, total_energy_full, _, _, lambda_star, valid = _phi_corrected_local_coordinates(speed_grid=speed_grid, lambda_grid=lambda_grid, pitch_grid=pitch_grid, base_distribution_v_lambda=base_distribution_v_lambda, mirror_ratio=mirror_ratio, B_tilde=B_tilde, local_potential_drop_magnitude_J=local_potential_drop_magnitude_J, throat_potential_drop_magnitude_J=throat_potential_drop_magnitude_J, eta_to_local_phase_space_normalization=eta_to_local_phase_space_normalization, particle_mass_kg=particle_mass_kg, quadrature_order=quadrature_order, local_speed_grid=local_speed_grid, base_distribution_interpolator=base_distribution_interpolator, quadrature_state=quadrature_state)
    sampled = _sample_phi_corrected_local_distribution(interpolator=interpolator, normalization=normalization, total_energy_full=total_energy_full, lambda_star=lambda_star, valid=valid, mirror_ratio=mirror_ratio, particle_mass_kg=particle_mass_kg)
    numerator = np.sum(sampled * quadrature.quadrature_weight, axis=(2, 3))
  
    return 2.0 * np.pi * float(np.sum(numerator))

def _phi_corrected_local_pitch_distribution(*, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, pitch_grid: PitchGrid, base_distribution_v_lambda: np.ndarray, mirror_ratio: float, B_tilde: float, local_potential_drop_magnitude_J: float, throat_potential_drop_magnitude_J: float, eta_to_local_phase_space_normalization: float, particle_mass_kg: float, quadrature_order: int = 3, local_speed_grid: SpeedGrid | None = None, base_distribution_interpolator: RegularGridInterpolator | None = None, quadrature_state: _LocalPitchQuadratureState | None = None) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Map invariant f(v, Λ) to local cell averages on (v_local, ξ) and Λ grids"""
    physical_speed_grid, quadrature, interpolator, normalization, total_energy_full, global_lambda, boundary, lambda_star, valid = _phi_corrected_local_coordinates(speed_grid=speed_grid, lambda_grid=lambda_grid, pitch_grid=pitch_grid, base_distribution_v_lambda=base_distribution_v_lambda, mirror_ratio=mirror_ratio, B_tilde=B_tilde, local_potential_drop_magnitude_J=local_potential_drop_magnitude_J, throat_potential_drop_magnitude_J=throat_potential_drop_magnitude_J, eta_to_local_phase_space_normalization=eta_to_local_phase_space_normalization, particle_mass_kg=particle_mass_kg, quadrature_order=quadrature_order, local_speed_grid=local_speed_grid, base_distribution_interpolator=base_distribution_interpolator, quadrature_state=quadrature_state)
    sampled = _sample_phi_corrected_local_distribution(interpolator=interpolator, normalization=normalization, total_energy_full=total_energy_full, lambda_star=lambda_star, valid=valid, mirror_ratio=mirror_ratio, particle_mass_kg=particle_mass_kg)
    numerator = np.sum(sampled * quadrature.quadrature_weight, axis=(2, 3))
    local_pitch = np.divide(numerator, quadrature.denominator, out=np.zeros_like(numerator), where=quadrature.denominator > 0.0)
    R = float(mirror_ratio)
    throat_drop = float(throat_potential_drop_magnitude_J)
    local_drop = float(local_potential_drop_magnitude_J)
    roundoff_tolerance = 128.0 * np.finfo(float).eps
    positive_energy = total_energy_full > 0.0
    clip_mask = valid & ((global_lambda < boundary) | (global_lambda > 1.0) | (lambda_star < 1.0 / R) | (lambda_star > 1.0))
    closed_interval_threshold = _eq71_closed_interval_threshold_J(R, throat_drop)
    low_energy = positive_energy & (total_energy_full <= closed_interval_threshold)
    positive_velocity_measure = float(np.sum(quadrature.quadrature_weight * positive_energy))
    closed_velocity_measure = float(np.sum(quadrature.quadrature_weight * low_energy))
    diagnostics: dict[str, object] = {
        "roundoff_clip_count": int(np.count_nonzero(clip_mask)),
        "nonpositive_total_energy_sample_count": int(np.count_nonzero(~positive_energy)),
        "closed_eq71_boundary_sample_count": int(np.count_nonzero(positive_energy & (boundary >= 1.0))),
        "low_energy_approximation_sample_count": int(np.count_nonzero(low_energy)),
        "closed_interval_velocity_space_measure_fraction": closed_velocity_measure / positive_velocity_measure if positive_velocity_measure > 0.0 else 0.0,
        "eq71_closed_interval_threshold_J": float(closed_interval_threshold),
        "invariant_speed_max_m_s": float(speed_grid.faces_m_s[-1]),
        "required_local_physical_speed_max_m_s": float(np.sqrt(speed_grid.faces_m_s[-1] ** 2 + 2.0 * local_drop / float(particle_mass_kg))),
        "local_physical_speed_max_m_s": float(physical_speed_grid.faces_m_s[-1]),
        "local_speed_domain_sufficient": bool(physical_speed_grid.faces_m_s[-1] >= np.sqrt(speed_grid.faces_m_s[-1] ** 2 + 2.0 * local_drop / float(particle_mass_kg)) - 128.0 * np.finfo(float).eps * max(float(physical_speed_grid.faces_m_s[-1]), 1.0)),
    }
    local_lambda = conservative_lambda_distribution_from_local_pitch_distribution(pitch_grid=pitch_grid, local_distribution_v_xi=local_pitch, lambda_grid=lambda_grid, B_tilde=float(B_tilde), fill_value=0.0)
   
    return local_lambda, local_pitch, diagnostics

def reconstruct_zero_potential_magnetic_reference(*, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, pitch_grid: PitchGrid, base_distribution_v_lambda: np.ndarray, B_tilde: np.ndarray, cell_volumes_m3: np.ndarray, mirror_ratio: float, eta_to_local_phase_space_normalization: float, expected_global_inventory_particles: float | None = None, quadrature_order: int = 3) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
    """
    ϕ = 0 magnetic local reconstruction without Eqs 70 through 72

    The magnetic boundary Λ = 1 / R_M is inserted as a temporary cell face when needed
    """
    B = np.asarray(B_tilde, dtype=float)
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    if B.ndim != 1 or volumes.shape != B.shape or B.size == 0:
        raise ValueError("B_tilde and cell_volumes_m3 must be matching nonempty vectors")
    if np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise ValueError("B_tilde must contain positive finite values")
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("cell_volumes_m3 must contain positive finite values")
    R = float(mirror_ratio)
    if not np.isfinite(R) or R <= 1.0:
        raise ValueError("mirror_ratio must be finite and greater than one")
    field_tolerance = 128.0 * np.finfo(float).eps * max(R, float(np.max(B)))
    if np.any(B > R + field_tolerance):
        raise ValueError("B_tilde exceeds mirror_ratio inside the fixed magnetic domain")
    normalization = float(eta_to_local_phase_space_normalization)
    if not np.isfinite(normalization) or normalization <= 0.0:
        raise ValueError("eta_to_local_phase_space_normalization must be positive and finite")
    if int(quadrature_order) < 1:
        raise ValueError("quadrature_order must be positive when supplied")
    base = _validated_base_distribution(speed_grid, lambda_grid, base_distribution_v_lambda)
    boundary = 1.0 / R
    original_faces = np.asarray(lambda_grid.faces, dtype=float)
    face_tolerance = 128.0 * np.finfo(float).eps * max(float(original_faces[-1]), 1.0)
    # Insert Λ = 1 / R_M as a face so no remap cell straddles the magnetic boundary
    if np.any(np.isclose(original_faces, boundary, rtol=0.0, atol=face_tolerance)):
        magnetic_lambda_grid = lambda_grid
    else:
        magnetic_lambda_grid = lambda_grid_from_faces(np.sort(np.concatenate((original_faces, np.array([boundary], dtype=float)))))
    original_cell_indices = np.searchsorted(original_faces, np.asarray(magnetic_lambda_grid.centers, dtype=float), side="right") - 1
    original_cell_indices = np.clip(original_cell_indices, 0, base.shape[1] - 1)
    trapped_distribution = base[:, original_cell_indices] * normalization
    trapped_distribution[:, magnetic_lambda_grid.centers < boundary] = 0.0
    local_lambda = np.zeros((B.size, speed_grid.centers_m_s.size, lambda_grid.centers.size), dtype=float)
    local_pitch = np.zeros((B.size, speed_grid.centers_m_s.size, pitch_grid.centers.size), dtype=float)
    density = np.zeros(B.size, dtype=float)
    for index, field_value in enumerate(B):
        local_pitch[index] = conservative_local_pitch_distribution_from_lambda_distribution(lambda_grid=magnetic_lambda_grid, distribution_v_lambda=trapped_distribution, pitch_grid=pitch_grid, B_tilde=float(field_value), fill_value=0.0)
        local_lambda[index] = conservative_lambda_distribution_from_local_pitch_distribution(pitch_grid=pitch_grid, local_distribution_v_xi=local_pitch[index], lambda_grid=lambda_grid, B_tilde=float(field_value), fill_value=0.0)
        density[index] = _density_from_local_pitch(speed_grid, pitch_grid, local_pitch[index])
    local_inventory = float(np.sum(density * volumes))
    expected = None if expected_global_inventory_particles is None else float(expected_global_inventory_particles)
    if expected is not None and (not np.isfinite(expected) or expected < 0.0):
        raise ValueError("expected_global_inventory_particles must be finite and nonnegative")
    inventory_error = ((local_inventory - expected) / expected if expected is not None and expected > 0.0 else 0.0)
    diagnostics: dict[str, object] = {
        "model": "zero_potential_magnetic_reference",
        "potential_identically_zero": True,
        "wall_barrier_energy_J": 0.0,
        "eq70_evaluated": False,
        "eq71_evaluated": False,
        "eq72_evaluated": False,
        "magnetic_trapped_passing_boundary": float(boundary),
        "magnetic_boundary_inserted_as_temporary_lambda_face": bool(magnetic_lambda_grid.faces.size != lambda_grid.faces.size),
        "local_mapping_model": "direct_fixed_magnetic_moment_pitch_remap",
        "quadrature_evaluated": False,
        "local_inventory_particles": local_inventory,
        "expected_global_inventory_particles": expected,
        "local_to_global_inventory_relative_error": (None if expected is None else float(inventory_error)),
        "active_electrostatic_lost_population_available": False,
        "electron_density_profile_available": False,
        "electron_current_balance_available": False,
        "electron_wall_power_available": False,
    }

    return local_lambda, local_pitch, density, diagnostics

def _density_from_local_pitch(speed_grid: SpeedGrid, pitch_grid: PitchGrid, local_v_pitch: np.ndarray) -> float:
    """Integrate a local f(v, ξ) with exact gyrotropic velocity cell measures"""
    measures = gyrotropic_velocity_cell_volumes(speed_grid, pitch_grid)

    return float(np.sum(np.maximum(local_v_pitch, 0.0) * measures))
