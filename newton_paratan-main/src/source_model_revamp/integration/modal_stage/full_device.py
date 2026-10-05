"""
Full device fast ion population mapping and expander trajectory support for the modal stage

The confined local distribution is mapped to the full geometry and the active fixed magnetic boundary Eq 63 branches are reconstructed before the later source connected expander closure
"""
from __future__ import annotations
from typing import Any
import numpy as np
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState
from source_model_revamp.fbis.modal.full_device_lost_reconstruction import reconstruct_full_device_fixed_boundary_eq63
from source_model_revamp.fbis.modal.lost_ion_distribution import EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY
from source_model_revamp.fbis.modal import ModalFBISResult
from source_model_revamp.integration.full_device_populations import conservative_confined_distribution_on_grid
from source_model_revamp.integration.pipeline_types import BeamEnsembleResult, GeometryStageResult
from source_model_revamp.losses.expander_trajectories import trace_adiabatic_lost_ions_to_surface

def _expander_trajectory_support_metadata(geometry: GeometryStageResult, *, modal: ModalFBISResult | None = None, material_loss_rate_s: float) -> dict[str, Any]:
    """
    Classify whether directed Eq 63 throat rates can be traced from each throat to its selected material surface
    
    The trajectory audit keeps throat crossing distinct from material deposition and reports wall, reflected, local well, and unresolved fractions when directed bin rates are available
    """
    loss_is_material = bool(np.isfinite(material_loss_rate_s) and material_loss_rate_s > 0.0)
    selection = geometry.plasma_boundaries
    full_centers = geometry.full_device_z_centers_m
    direct_field = geometry.B_T_function
    if selection is None or full_centers is None or direct_field is None:
        status = "not_available_without_full_device_boundaries_and_direct_field"
        return {
            "expander_trajectory_classifier": "source_model_revamp.losses.expander_trajectories.trace_adiabatic_lost_ions_to_surface",
            "expander_trajectory_geometry_support_available": False,
            "expander_trajectory_classification_rate_assessed": False,
            "modal_active_material_deposition_rate_available": False,
            "left_expander_trajectory_status": status,
            "right_expander_trajectory_status": status,
            "left_wall_particle_rate_s": None,
            "right_wall_particle_rate_s": None,
            "left_wall_power_W": None,
            "right_wall_power_W": None,
            "reflected_particle_fraction": None,
            "locally_trapped_expander_fraction": None,
            "unclassified_expander_fraction": None,
            "wall_surface_deposition_breakdown": None,
            "expander_unclassified_population_check_passed": not loss_is_material,
            "expander_trajectory_failure_reason": ( "material throat crossing loss exists but full device trajectory geometry is unavailable" if loss_is_material else None),
        }

    z_full = np.asarray(full_centers, dtype=float)
    side_metadata: dict[str, Any] = {}
    path_data: dict[str, tuple[np.ndarray, np.ndarray, Any]] = {}
    support_available = bool(selection.valid)
    for side, throat_z, surface in (("left", float(geometry.z_edges_m[0]), selection.left_selected_loss_surface), ("right", float(geometry.z_edges_m[-1]), selection.right_selected_loss_surface),):
        wall_z = surface.z_m
        if wall_z is None:
            support_available = False
            side_metadata[f"{side}_expander_trajectory_status"] = ("not_rate_assessed_selected_nonaxial_surface_requires_explicit_intersection_path")
            side_metadata[f"{side}_expander_path_z_m"] = []
            side_metadata[f"{side}_expander_path_B_tilde"] = []
            continue
        lower, upper = sorted((throat_z, float(wall_z)))
        interior = z_full[(z_full > lower) & (z_full < upper)]
        ordered_interior = np.sort(interior)
        if side == "left":
            ordered_interior = ordered_interior[::-1]
        path_z = np.concatenate(([throat_z], ordered_interior, [float(wall_z)]))
        path_B_tilde = np.asarray([float(direct_field(float(value))) / float(geometry.B0_T) for value in path_z], dtype=float)
        if path_z.size < 2 or np.any(~np.isfinite(path_B_tilde)) or np.any(path_B_tilde <= 0.0):
            support_available = False
            side_metadata[f"{side}_expander_trajectory_status"] = "invalid_direct_field_path"
        else:
            side_metadata[f"{side}_expander_trajectory_status"] = ("geometry_ready_rate_not_assessed_missing_directed_Eq63_bin_rates")
        side_metadata[f"{side}_expander_path_z_m"] = path_z
        side_metadata[f"{side}_expander_path_B_tilde"] = path_B_tilde
        side_metadata[f"{side}_selected_loss_surface_id"] = surface.surface_id
        side_metadata[f"{side}_selected_loss_surface_particle_role"] = surface.particle_role
        electrical_model = str(getattr(surface, "electrical_model", "floating"))
        side_metadata[f"{side}_wall_electrical_model"] = electrical_model
        if electrical_model == "floating":
            path_data[side] = (path_z, path_B_tilde, surface)
        else:
            support_available = False
            side_metadata[f"{side}_expander_trajectory_status"] = ("electrical_boundary_not_coupled_to_midplane_reference_potential")
    left_directed_rate = None if modal is None else modal.metadata.get("modal_eq63_left_throat_rate_v_lambda_s")
    right_directed_rate = None if modal is None else modal.metadata.get("modal_eq63_right_throat_rate_v_lambda_s")
    shared_directed_rate = None if modal is None else modal.metadata.get("modal_directed_eq63_throat_rate_v_lambda_s")
    if left_directed_rate is None:
        left_directed_rate = shared_directed_rate
    if right_directed_rate is None:
        right_directed_rate = shared_directed_rate
    if (left_directed_rate is None or right_directed_rate is None or len(path_data) != 2 or modal is None or modal.electrostatic_profile is None):
        failure_reason = None
        if loss_is_material:
            failure_reason = ("material throat crossing loss exists but directed Eq63 bin rates or paths are unavailable, wall deposition and unclassified population cannot be assessed")
    
        return {
            "expander_trajectory_classifier": "source_model_revamp.losses.expander_trajectories.trace_adiabatic_lost_ions_to_surface",
            "expander_trajectory_geometry_support_available": support_available,
            "expander_trajectory_classification_rate_assessed": False,
            "modal_active_material_deposition_rate_available": False,
            **side_metadata,
            "left_wall_particle_rate_s": None,
            "right_wall_particle_rate_s": None,
            "left_wall_power_W": None,
            "right_wall_power_W": None,
            "reflected_particle_fraction": None,
            "locally_trapped_expander_fraction": None,
            "unclassified_expander_fraction": None,
            "wall_surface_deposition_breakdown": None,
            "expander_unclassified_population_check_passed": not loss_is_material,
            "expander_trajectory_failure_reason": failure_reason,
            "throat_crossing_and_material_wall_deposition_are_distinct": True,
        }

    rate_matrices = { "left": np.asarray(left_directed_rate, dtype=float), "right": np.asarray(right_directed_rate, dtype=float)}
    invariant_energy = ( 0.5 * modal.species.mass_kg * np.asarray(modal.speed_grid.centers_m_s, dtype=float) ** 2)
    lambda_faces = np.asarray(modal.lambda_grid.faces, dtype=float)
    lambda_open_lower = np.maximum(lambda_faces[:-1], 0.0)
    lambda_open_upper = np.minimum( lambda_faces[1:], 1.0 / float(geometry.mirror_ratio),)
    global_lambda = 0.5 * (lambda_open_lower + lambda_open_upper)
    profile = modal.electrostatic_profile
    throat_drops = { "left": float(profile.throat_potential_left_energy_J), "right": float(profile.throat_potential_right_energy_J),}
    wall_drop = float(modal.wall_barrier_energy_J)
    results = {}
    for side in ("left", "right"):
        path_z, path_B, surface = path_data[side]
        distance = np.abs(path_z - path_z[0])
        fraction = distance / max(float(distance[-1]), np.finfo(float).tiny)
        potential_path = throat_drops[side] + fraction * (wall_drop - throat_drops[side])
        results[side] = trace_adiabatic_lost_ions_to_surface(direction=side, invariant_total_energy_J=invariant_energy, global_lambda=global_lambda, directed_particle_rate_v_lambda_s=rate_matrices[side], B_tilde_path=path_B, potential_drop_magnitude_path_J=potential_path, selected_surface_id=surface.surface_id, selected_surface_particle_role=surface.particle_role,)
        side_metadata[f"{side}_expander_trajectory_status"] = "rate_assessed_from_published_eq63_directed_throat_flux"
        side_metadata[f"{side}_expander_potential_model"] = "linear_numerical_connection_from_solved_throat_to_floating_wall_barrier"
    branch_rates = { side: float(np.sum(rate_matrices[side])) for side in ("left", "right")}
    total_classified_rate = branch_rates["left"] + branch_rates["right"]
    prompt_rate = float(modal.metadata.get("modal_prompt_loss_birth_rate_s", 0.0))
    unresolved_rate = prompt_rate + sum( results[side].unclassified_particle_fraction * branch_rates[side] for side in ("left", "right"))
    total_rate_scale = max(total_classified_rate + prompt_rate, np.finfo(float).tiny)
    unclassified_fraction = unresolved_rate / total_rate_scale
    reflected_rates = { side: results[side].reflected_particle_fraction * branch_rates[side] for side in ("left", "right")}
    local_rates = { side: results[side].locally_trapped_particle_fraction * branch_rates[side] for side in ("left", "right")}
    reflected_fraction = sum(reflected_rates.values()) / total_rate_scale
    local_fraction = sum(local_rates.values()) / total_rate_scale
    unresolved_closure_fraction = unclassified_fraction + reflected_fraction + local_fraction
    failure_reason = ( None if unresolved_closure_fraction <= 1.0e-6 else "prompt, reflected, locally trapped, or Eq63 expander population lacks a steady return/deposition closure")
 
    return {
        "expander_trajectory_classifier": "source_model_revamp.losses.expander_trajectories.trace_adiabatic_lost_ions_to_surface",
        "expander_trajectory_geometry_support_available": support_available,
        "expander_trajectory_classification_rate_assessed": True,
        "modal_active_material_deposition_rate_available": True,
        **side_metadata,
        "left_wall_particle_rate_s": float(results["left"].wall_particle_rate_s),
        "right_wall_particle_rate_s": float(results["right"].wall_particle_rate_s),
        "fixed_boundary_left_throat_crossing_rate_s": float(branch_rates["left"]),
        "fixed_boundary_right_throat_crossing_rate_s": float(branch_rates["right"]),
        "fixed_boundary_left_wall_deposition_rate_s": float(results["left"].wall_particle_rate_s),
        "fixed_boundary_right_wall_deposition_rate_s": float(results["right"].wall_particle_rate_s),
        "fixed_boundary_left_reflected_rate_s": float(reflected_rates["left"]),
        "fixed_boundary_right_reflected_rate_s": float(reflected_rates["right"]),
        "fixed_boundary_left_unclassified_rate_s": float(results["left"].unclassified_particle_fraction * branch_rates["left"]),
        "fixed_boundary_right_unclassified_rate_s": float(results["right"].unclassified_particle_fraction * branch_rates["right"]),
        "fixed_boundary_local_well_topology_incident_rate_s": float(sum(local_rates.values())),
        "fixed_boundary_local_well_population_rate_s": None,
        "left_wall_power_W": float(results["left"].wall_power_W),
        "right_wall_power_W": float(results["right"].wall_power_W),
        "reflected_particle_fraction": float(reflected_fraction),
        "locally_trapped_expander_fraction": float(local_fraction),
        "unclassified_expander_fraction": float(unclassified_fraction),
        "wall_surface_deposition_breakdown": { results["left"].selected_surface_id: float(results["left"].wall_particle_rate_s), results["right"].selected_surface_id: float(results["right"].wall_particle_rate_s),},
        "expander_unresolved_closure_fraction": float(unresolved_closure_fraction),
        "fixed_boundary_rate_accounting_relative_error": float(abs(sum( results[side].wall_particle_rate_s + reflected_rates[side] + results[side].unclassified_particle_fraction * branch_rates[side] + local_rates[side] for side in ("left", "right")) - total_classified_rate) / max(total_classified_rate, np.finfo(float).tiny)),
        "expander_unclassified_population_check_passed": bool(unresolved_closure_fraction <= 1.0e-6),
        "expander_trajectory_failure_reason": failure_reason,
        "throat_crossing_and_material_wall_deposition_are_distinct": True,
    }

def _full_device_fast_ion_population(*, geometry: GeometryStageResult, beam: BeamEnsembleResult, modal: ModalFBISResult, collision_state: FBISCollisionParameterState, quadrature_order: int) -> tuple[np.ndarray | None, dict[str, Any]]:
    """
    Map confined and fixed magnetic boundary Eq 63 fast ions onto the full device local `(z, v, ξ)` grid
    
    The confined distribution is conservatively transferred first
    For the active fixed magnetic boundary, left and right Eq 63 branches are reconstructed with local total energy mapping through the solved Eq 70 potential
    """
    del collision_state
    local = modal.local_reconstruction
    loss_model = str(modal.metadata.get("ion_loss_closure_model", modal.metadata.get("modal_ion_loss_closure_model", ""),)).strip().lower()
    fixed_boundary_selected = loss_model == EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY
    required_geometry = (geometry.full_device_z_edges_m, geometry.full_device_z_centers_m, geometry.full_device_B_tilde_centers, geometry.full_device_cell_volumes_m3,)
    metadata: dict[str, Any] = {
        "full_device_lost_ion_reconstruction_status": "not_evaluated",
        "full_device_lost_ion_reconstruction_available": False,
        "full_device_lost_ion_reconstruction_failure_reason": None,
        "left_lost_ion_distribution": None,
        "right_lost_ion_distribution": None,
        "full_device_lost_ion_distribution": None,
        "full_device_lost_ion_potential_drop_magnitude_J": None,
        "full_device_fast_ion_distribution_includes_confined": False,
        "full_device_fast_ion_distribution_includes_directed_lost": False,
        "full_device_lost_ion_reconstruction_matches_active_hot_or_electrostatic_sink": False,
        "full_device_eq63_reference_is_separate_diagnostic": not fixed_boundary_selected,
        "full_device_eq63_reference_may_substitute_for_active_population": False,
        "full_device_lost_ion_forced_normalization_applied": False,
    }
    if local is None:
        metadata.update({"full_device_confined_population_status": "unavailable_no_local_reconstruction", "full_device_lost_ion_reconstruction_status": "unavailable_no_local_reconstruction", "full_device_lost_ion_reconstruction_failure_reason": "local_confined_distribution_is_unavailable",})
        return None, metadata
    if any(value is None for value in required_geometry):
        metadata.update({"full_device_confined_population_status": "unavailable_no_full_device_grid", "full_device_lost_ion_reconstruction_status": "unavailable_no_full_device_grid", "full_device_lost_ion_reconstruction_failure_reason": "geometry_does_not_supply_a_full_device_grid",})
        return None, metadata

    full_edges = np.asarray(geometry.full_device_z_edges_m, dtype=float)
    full_centers = np.asarray(geometry.full_device_z_centers_m, dtype=float)
    full_B_tilde = np.asarray(geometry.full_device_B_tilde_centers, dtype=float)
    full_volumes = np.asarray(geometry.full_device_cell_volumes_m3, dtype=float)
    # Transfer the confined local distribution before adding directed lost branches
    confined_full = conservative_confined_distribution_on_grid( source_edges_m=np.asarray(geometry.z_edges_m, dtype=float), source_cell_volumes_m3=np.asarray(geometry.cell_volumes_m3, dtype=float), source_distribution_z_v_pitch=np.asarray(local.local_distribution_z_v_pitch, dtype=float), target_edges_m=full_edges, target_cell_volumes_m3=full_volumes,)
    local_speed_grid = local.local_speed_grid or modal.speed_grid
    local_diagnostics = local.local_speed_grid_diagnostics or {}
    metadata.update({"full_device_confined_population_status": "represented", "full_device_fast_ion_distribution_includes_confined": True, "full_device_local_speed_grid_upper_m_s": float(local_speed_grid.faces_m_s[-1]), "full_device_local_speed_required_upper_m_s": local_diagnostics.get( "required_local_physical_speed_max_m_s", float(local_speed_grid.faces_m_s[-1]),), "full_device_local_speed_domain_overflow_possible": not bool(local_diagnostics.get("local_speed_domain_sufficient", True)),})
    if not fixed_boundary_selected:
        metadata.update({"full_device_lost_ion_reconstruction_status": "unavailable_for_selected_loss_model", "full_device_lost_ion_reconstruction_failure_reason": ("selected_loss_model_does_not_supply_the_active_fixed_boundary_Eq63_population"),})
        return confined_full, metadata
    selection = geometry.plasma_boundaries
    profile = modal.electrostatic_profile
    left_H = modal.metadata.get("modal_eq63_left_H_U")
    right_H = modal.metadata.get("modal_eq63_right_H_U")
    left_T_L = modal.metadata.get("lost_ion_parallel_temperature_left_J")
    right_T_L = modal.metadata.get("lost_ion_parallel_temperature_right_J")
    if selection is None or not selection.valid:
        metadata.update({"full_device_lost_ion_reconstruction_status": "unavailable_invalid_plasma_boundary_selection", "full_device_lost_ion_reconstruction_failure_reason": "left_and_right_material_loss_surfaces_are_required",})
        return confined_full, metadata
    left_wall_z = selection.left_selected_loss_surface.z_m
    right_wall_z = selection.right_selected_loss_surface.z_m
    if left_wall_z is None or right_wall_z is None:
        metadata.update({"full_device_lost_ion_reconstruction_status": "unavailable_nonaxial_material_boundary", "full_device_lost_ion_reconstruction_failure_reason": "the_current_one_dimensional_expander_path_requires_axial_material_boundaries",})
        return confined_full, metadata
    if profile is None:
        metadata.update({"full_device_lost_ion_reconstruction_status": "unavailable_no_electrostatic_profile", "full_device_lost_ion_reconstruction_failure_reason": "Eq70_profile_is_required_for_local_energy_mapping",})
        return confined_full, metadata
    if left_H is None or right_H is None or left_T_L is None or right_T_L is None:
        metadata.update({"full_device_lost_ion_reconstruction_status": "unavailable_incomplete_Eq61_Eq64_throat_state", "full_device_lost_ion_reconstruction_failure_reason": "separate_left_right_H_U_and_T_L_are_required",})
        return confined_full, metadata
    left_throat_z = float(geometry.z_edges_m[0])
    right_throat_z = float(geometry.z_edges_m[-1])
    confined_z = np.asarray(geometry.z_centers_m, dtype=float)
    confined_drop = np.maximum(np.asarray(profile.potential_energy_J, dtype=float), 0.0)
    potential_full = np.interp(full_centers, confined_z, confined_drop, left=float(profile.throat_potential_left_energy_J), right=float(profile.throat_potential_right_energy_J))
    potential_full = np.maximum(potential_full, 0.0)
    # Reconstruct left and right Eq 63 branches with the solved local energy mapping
    lost_total, left_lost, right_lost, lost_diagnostics = (reconstruct_full_device_fixed_boundary_eq63(
            invariant_speed_grid=modal.speed_grid,
            local_speed_grid=local_speed_grid,
            pitch_grid=modal.pitch_grid,
            left_H_U=np.asarray(left_H, dtype=float),
            right_H_U=np.asarray(right_H, dtype=float),
            left_parallel_temperature_J=float(left_T_L),
            right_parallel_temperature_J=float(right_T_L),
            z_edges_m=full_edges,
            B_tilde_centers=full_B_tilde,
            potential_drop_magnitude_centers_J=potential_full,
            left_throat_z_m=left_throat_z,
            right_throat_z_m=right_throat_z,
            left_wall_z_m=float(left_wall_z),
            right_wall_z_m=float(right_wall_z),
            mirror_ratio=float(geometry.mirror_ratio),
            half_length_m=float(geometry.half_length_m),
            geometry_factor_G=float(modal.basis.geometry_factor_G),
            quadrature_order=int(quadrature_order),
            particle_mass_kg=modal.species.mass_kg,
        ))
    local_speed_coverage_passed = bool(lost_diagnostics.get("local_speed_coverage_passed", False))
    full_distribution = confined_full + lost_total
    metadata.update({
        "full_device_lost_ion_reconstruction_status": ("central_Eq63_represented_expander_deferred" if local_speed_coverage_passed else "central_Eq63_represented_but_local_speed_domain_insufficient"),
        "full_device_lost_ion_reconstruction_available": True,
        "full_device_lost_ion_reconstruction_failure_reason": (None if local_speed_coverage_passed else "local_physical_speed_grid_does_not_cover_the_central_Eq63_population"),
        "full_device_lost_ion_expander_role": "deferred_to_source_connected_expander_Eq70_stage",
        "full_device_lost_ion_pre_expander_potential_outside_throats_is_authoritative": False,
        "left_lost_ion_distribution": left_lost,
        "right_lost_ion_distribution": right_lost,
        "full_device_lost_ion_distribution": lost_total,
        "full_device_lost_ion_potential_drop_magnitude_J": potential_full,
        "full_device_fast_ion_distribution_includes_directed_lost": True,
        "full_device_lost_ion_reconstruction_matches_active_hot_or_electrostatic_sink": True,
        "full_device_eq63_reference_is_separate_diagnostic": False,
        "full_device_lost_ion_forced_normalization_applied": False,
        "full_device_lost_ion_reconstruction_diagnostics": lost_diagnostics,
        "fixed_boundary_local_speed_coverage_passed": local_speed_coverage_passed,
        "fixed_boundary_required_local_speed_max_m_s": lost_diagnostics.get("required_local_physical_speed_max_m_s"),
        "fixed_boundary_represented_local_speed_max_m_s": lost_diagnostics.get("represented_local_physical_speed_max_m_s"),
        "fixed_boundary_maximum_total_energy_invariant_residual_J": (lost_diagnostics.get("maximum_total_energy_invariant_residual_J")),
        "fixed_boundary_maximum_mu_B0_invariant_residual_J": (lost_diagnostics.get("maximum_mu_B0_invariant_residual_J")),
        "fixed_boundary_left_wrong_sign_population_integral": (lost_diagnostics.get("left_wrong_sign_population_integral")),
        "fixed_boundary_right_wrong_sign_population_integral": (lost_diagnostics.get("right_wrong_sign_population_integral")),
    })

    return full_distribution, metadata