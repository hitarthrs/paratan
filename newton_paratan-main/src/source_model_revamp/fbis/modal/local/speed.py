"""Physical local speed domains, moments, and refinement checks"""
from __future__ import annotations
import numpy as np
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid, speed_grid_from_faces
from source_model_revamp.fbis.modal.local.mapping import _base_distribution_interpolator
from source_model_revamp.fbis.modal.local.pitch import _phi_corrected_local_pitch_distribution

def derive_local_physical_speed_grid(*, invariant_speed_grid: SpeedGrid, maximum_potential_drop_magnitude_J: float, particle_mass_kg: float) -> tuple[SpeedGrid, dict[str, float | int | bool | str]]:
    """
    Extend the invariant speed domain to cover electrostatically accelerated ions

    v_local,max = sqrt(v_invariant,max² + 2 ΔΦ_max / m)
    """
    drop = float(maximum_potential_drop_magnitude_J)
    mass = float(particle_mass_kg)
    if not np.isfinite(drop) or drop < 0.0:
        raise ValueError("maximum_potential_drop_magnitude_J must be finite and nonnegative")
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    invariant_faces = np.asarray(invariant_speed_grid.faces_m_s, dtype=float)
    invariant_max = float(invariant_faces[-1])
    required_max = float(np.sqrt(invariant_max**2 + 2.0 * drop / mass))
    tolerance = 128.0 * np.finfo(float).eps * max(required_max, invariant_max, 1.0)
    if required_max <= invariant_max + tolerance:
        local_faces = invariant_faces.copy()
        appended = 0
    else:
        last_width = float(invariant_faces[-1] - invariant_faces[-2])
        if not np.isfinite(last_width) or last_width <= 0.0:
            raise ValueError("invariant speed grid has an invalid final cell width")
        appended = max(1, int(np.ceil((required_max - invariant_max) / last_width)))
        extension = np.linspace(invariant_max, required_max, appended + 1)[1:]
        local_faces = np.concatenate((invariant_faces, extension))
    local_grid = speed_grid_from_faces(local_faces)
    diagnostics: dict[str, float | int | bool | str] = {
        "invariant_speed_max_m_s": invariant_max,
        "maximum_potential_drop_magnitude_J": drop,
        "required_local_physical_speed_max_m_s": required_max,
        "local_physical_speed_max_m_s": float(local_grid.faces_m_s[-1]),
        "local_speed_appended_cell_count": int(appended),
        "local_speed_domain_sufficient": bool(local_grid.faces_m_s[-1] >= required_max - tolerance),
        "local_speed_grid_model": "invariant_faces_preserved_with_energy_derived_upper_extension",
        "local_speed_overflow_particle_fraction": 0.0,
        "local_speed_overflow_energy_fraction": 0.0,
        "local_speed_tail_clipped": False,
    }

    return local_grid, diagnostics

def subdivide_speed_grid_cells(speed_grid: SpeedGrid, subdivision_factor: int) -> SpeedGrid:
    """Subdivide every speed cell while preserving all original faces and domain limits"""
    factor = int(subdivision_factor)
    if factor < 1:
        raise ValueError("subdivision_factor must be at least one")
    if factor == 1:
        return speed_grid
    faces = np.asarray(speed_grid.faces_m_s, dtype=float)
    refined: list[float] = [float(faces[0])]
    for lower, upper in zip(faces[:-1], faces[1:], strict=True):
        segment = np.linspace(float(lower), float(upper), factor + 1)
        refined.extend(float(value) for value in segment[1:])

    return speed_grid_from_faces(np.asarray(refined, dtype=float))

def _local_speed_profile_moments(*, local_speed_grid: SpeedGrid, pitch_grid: PitchGrid, local_distribution_z_v_pitch: np.ndarray, cell_volumes_m3: np.ndarray, particle_mass_kg: float) -> dict[str, object]:
    """
    Return exact piecewise constant density and kinetic energy moments

    Input shape is (n_z, n_local_speed, n_pitch)
    The gyrotropic measures are evaluated as
        2π * integral v^2 dv dξ for density and
        2π * integral (m v^2/2) v^2 dv dξ for kinetic energy
    """
    distribution = np.asarray(local_distribution_z_v_pitch, dtype=float)
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    expected = (volumes.size, local_speed_grid.centers_m_s.size, pitch_grid.centers.size)
    if distribution.shape != expected:
        raise ValueError("local_distribution_z_v_pitch has an inconsistent shape")
    if np.any(~np.isfinite(distribution)):
        raise ValueError("local distribution must contain only finite values")
    scale = max(float(np.max(np.abs(distribution))) if distribution.size else 0.0, np.finfo(float).tiny)
    if float(np.min(distribution)) < -128.0 * np.finfo(float).eps * scale:
        raise ValueError("local distribution contains a significant negative value")
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("cell_volumes_m3 must be positive and finite")
    mass = float(particle_mass_kg)
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    distribution = np.maximum(distribution, 0.0)
    lower = np.asarray(local_speed_grid.faces_m_s[:-1], dtype=float)
    upper = np.asarray(local_speed_grid.faces_m_s[1:], dtype=float)
    speed_v3 = (upper**3 - lower**3) / 3.0
    speed_v5 = (upper**5 - lower**5) / 5.0
    pitch_measure = 2.0 * np.pi * np.asarray(pitch_grid.widths, dtype=float)
    density_profile = np.sum(distribution * speed_v3[None, :, None] * pitch_measure[None, None, :], axis=(1, 2))
    kinetic_energy_density_profile = 0.5 * mass * np.sum(distribution * speed_v5[None, :, None] * pitch_measure[None, None, :], axis=(1, 2))

    return {"density_profile_m3": density_profile, "kinetic_energy_density_profile_J_m3": kinetic_energy_density_profile, "particle_inventory": float(np.sum(density_profile * volumes)), "kinetic_energy_inventory_J": float(np.sum(kinetic_energy_density_profile * volumes)),}

def _volume_weighted_profile_relative_change(current: np.ndarray, previous: np.ndarray, cell_volumes_m3: np.ndarray) -> float:
    """Return the cell volume weighted relative L2 change between axial profiles"""
    difference = np.asarray(current, dtype=float) - np.asarray(previous, dtype=float)
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    numerator = np.sqrt(float(np.sum(volumes * difference**2)))
    denominator = np.sqrt(float(np.sum(volumes * np.asarray(current, dtype=float) ** 2)))

    return numerator / max(denominator, np.finfo(float).tiny)

def assess_local_physical_speed_grid_convergence(*, invariant_speed_grid: SpeedGrid, local_speed_grid: SpeedGrid, lambda_grid: LambdaGrid, pitch_grid: PitchGrid, base_distribution_v_lambda: np.ndarray, mirror_ratio: float, zeta: np.ndarray, B_tilde: np.ndarray, potential_drop_magnitude_J: np.ndarray, left_throat_potential_drop_magnitude_J: float, right_throat_potential_drop_magnitude_J: float, eta_to_local_phase_space_normalization: float, cell_volumes_m3: np.ndarray, quadrature_order: int, relative_tolerance: float, subdivision_factors: tuple[int, ...] = (1, 2, 4), reference_local_distribution_z_v_pitch: np.ndarray | None = None, particle_mass_kg: float) -> dict[str, object]:
    """
    Assess local speed discretization convergence on one fixed physical domain

    Refinement compares density and kinetic energy profiles plus integrated inventories
    """
    factors = tuple(int(value) for value in subdivision_factors)
    if len(factors) < 2 or factors[0] != 1:
        raise ValueError("subdivision_factors must contain at least two levels and start at one")
    if any(value < 1 for value in factors) or any(current <= previous for previous, current in zip(factors[:-1], factors[1:])):
        raise ValueError("subdivision_factors must be strictly increasing")
    tolerance = float(relative_tolerance)
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("relative_tolerance must be positive and finite")
    z = np.asarray(zeta, dtype=float)
    B = np.asarray(B_tilde, dtype=float)
    phi = np.asarray(potential_drop_magnitude_J, dtype=float)
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    if (B.ndim != 1 or z.shape != B.shape or phi.shape != B.shape or volumes.shape != B.shape):
        raise ValueError("zeta, B_tilde, potential_drop_magnitude_J, and cell volumes must match")
    if np.any(~np.isfinite(z)):
        raise ValueError("zeta must contain only finite values")
    if np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise ValueError("B_tilde must be positive and finite")
    if np.any(~np.isfinite(phi)) or np.any(phi < 0.0):
        raise ValueError("potential_drop_magnitude_J must be finite and nonnegative")
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("cell_volumes_m3 must be positive and finite")
    base_interpolator = _base_distribution_interpolator(invariant_speed_grid, lambda_grid, base_distribution_v_lambda,)
    histories: list[dict[str, object]] = []
    density_profiles: list[np.ndarray] = []
    energy_profiles: list[np.ndarray] = []
    particle_inventories: list[float] = []
    energy_inventories: list[float] = []
    all_domains_sufficient = True
    # Refinement changes cell resolution without changing the physical speed limits
    for level_index, factor in enumerate(factors):
        refined_grid = subdivide_speed_grid_cells(local_speed_grid, factor)
        # Reuse the supplied nominal reconstruction and recompute each finer level
        if level_index == 0 and reference_local_distribution_z_v_pitch is not None:
            local_pitch = np.asarray(reference_local_distribution_z_v_pitch, dtype=float)
            expected = (B.size, refined_grid.centers_m_s.size, pitch_grid.centers.size)
            if local_pitch.shape != expected:
                raise ValueError("reference local distribution does not match the nominal local grid")
            mapping_diagnostics: list[dict[str, object]] = []
        else:
            local_pitch = np.zeros((B.size, refined_grid.centers_m_s.size, pitch_grid.centers.size), dtype=float,)
            mapping_diagnostics = []
            for z_index in range(B.size):
                _, local_pitch[z_index], mapping = (
                    _phi_corrected_local_pitch_distribution(
                        speed_grid=invariant_speed_grid,
                        local_speed_grid=refined_grid,
                        lambda_grid=lambda_grid,
                        pitch_grid=pitch_grid,
                        base_distribution_v_lambda=base_distribution_v_lambda,
                        mirror_ratio=mirror_ratio,
                        B_tilde=float(B[z_index]),
                        local_potential_drop_magnitude_J=float(phi[z_index]),
                        throat_potential_drop_magnitude_J=float(left_throat_potential_drop_magnitude_J if z[z_index] < 0.0 else right_throat_potential_drop_magnitude_J),
                        eta_to_local_phase_space_normalization=(eta_to_local_phase_space_normalization),
                        quadrature_order=int(quadrature_order),
                        base_distribution_interpolator=base_interpolator,
                        particle_mass_kg=particle_mass_kg,
                    ))
                mapping_diagnostics.append(mapping)
        moments = _local_speed_profile_moments(local_speed_grid=refined_grid, pitch_grid=pitch_grid, local_distribution_z_v_pitch=local_pitch, cell_volumes_m3=volumes, particle_mass_kg=particle_mass_kg)
        density_profile = np.asarray(moments["density_profile_m3"], dtype=float)
        energy_profile = np.asarray(moments["kinetic_energy_density_profile_J_m3"], dtype=float)
        density_profiles.append(density_profile)
        energy_profiles.append(energy_profile)
        particle_inventories.append(float(moments["particle_inventory"]))
        energy_inventories.append(float(moments["kinetic_energy_inventory_J"]))
        required_upper = float(np.sqrt(invariant_speed_grid.faces_m_s[-1] ** 2 + 2.0 * float(np.max(phi)) / float(particle_mass_kg)))
        domain_sufficient = bool(refined_grid.faces_m_s[-1] >= required_upper - 128.0 * np.finfo(float).eps * max(required_upper, float(refined_grid.faces_m_s[-1]), 1.0))
        all_domains_sufficient = all_domains_sufficient and domain_sufficient
        histories.append({"subdivision_factor": int(factor), "speed_cell_count": int(refined_grid.centers_m_s.size), "local_physical_speed_max_m_s": float(refined_grid.faces_m_s[-1]), "required_local_physical_speed_max_m_s": required_upper, "local_speed_domain_sufficient": domain_sufficient, "overflow_particle_fraction": 0.0 if domain_sufficient else None, "overflow_energy_fraction": 0.0 if domain_sufficient else None, "particle_inventory": float(moments["particle_inventory"]), "kinetic_energy_inventory_J": float(moments["kinetic_energy_inventory_J"]), "mapping_roundoff_clip_count": int(sum(int(item.get("roundoff_clip_count", 0)) for item in mapping_diagnostics)),})
    density_changes: list[float] = []
    energy_changes: list[float] = []
    particle_inventory_changes: list[float] = []
    energy_inventory_changes: list[float] = []
    for index in range(1, len(factors)):
        density_changes.append(_volume_weighted_profile_relative_change(density_profiles[index], density_profiles[index - 1], volumes))
        energy_changes.append(_volume_weighted_profile_relative_change(energy_profiles[index], energy_profiles[index - 1], volumes))
        particle_inventory_changes.append(abs(particle_inventories[index] - particle_inventories[index - 1]) / max(abs(particle_inventories[index]), abs(particle_inventories[index - 1]), np.finfo(float).tiny))
        energy_inventory_changes.append(abs(energy_inventories[index] - energy_inventories[index - 1]) / max(abs(energy_inventories[index]), abs(energy_inventories[index - 1]), np.finfo(float).tiny))
    converged = bool(all_domains_sufficient and all(value <= tolerance for value in density_changes) and all(value <= tolerance for value in energy_changes) and all(value <= tolerance for value in particle_inventory_changes) and all(value <= tolerance for value in energy_inventory_changes))
    failure_reasons: list[str] = []
    if not all_domains_sufficient:
        failure_reasons.append("local_speed_domain_insufficient")
    if any(value > tolerance for value in density_changes):
        failure_reasons.append("local_density_profile_not_converged")
    if any(value > tolerance for value in energy_changes):
        failure_reasons.append("local_energy_density_profile_not_converged")
    if any(value > tolerance for value in particle_inventory_changes):
        failure_reasons.append("local_particle_inventory_not_converged")
    if any(value > tolerance for value in energy_inventory_changes):
        failure_reasons.append("local_energy_inventory_not_converged")

    return {
        "local_grid_refinement_assessed": True,
        "local_grid_refinement_converged": converged,
        "local_grid_refinement_status": ("converged" if converged else "not_converged"),
        "local_grid_refinement_failure_reason": ";".join(failure_reasons),
        "local_grid_refinement_relative_tolerance": tolerance,
        "local_grid_refinement_history": histories,
        "local_grid_density_relative_change_history": density_changes,
        "local_grid_energy_density_relative_change_history": energy_changes,
        "local_grid_particle_inventory_relative_change_history": (particle_inventory_changes),
        "local_grid_energy_inventory_relative_change_history": (energy_inventory_changes),
        "local_grid_final_density_relative_change": (density_changes[-1] if density_changes else 0.0),
        "local_grid_final_energy_density_relative_change": (energy_changes[-1] if energy_changes else 0.0),
        "local_grid_final_particle_inventory_relative_change": (particle_inventory_changes[-1] if particle_inventory_changes else 0.0),
        "local_grid_final_energy_inventory_relative_change": (energy_inventory_changes[-1] if energy_inventory_changes else 0.0),
        "local_speed_domain_sufficient": all_domains_sufficient,
        "local_speed_tail_clipped": False,
        "local_speed_overflow_particle_fraction": (0.0 if all_domains_sufficient else None),
        "local_speed_overflow_energy_fraction": (0.0 if all_domains_sufficient else None),
        "local_speed_refinement_model": ("fixed_energy_derived_domain_with_exact_cell_subdivision_1x_2x_4x"),
        "eq71_72_population_change_included_in_grid_error": False,
    }
