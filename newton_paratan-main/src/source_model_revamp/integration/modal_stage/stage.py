"""Top level modal FBIS integration stage for single species and coupled D plus T operation"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from source_model_revamp.fbis.collision_parameters import collision_parameter_metadata, energy_J_from_keV, rebuild_fbis_collision_parameter_state_with_fast_density
from source_model_revamp.fbis.species import DEUTERON, TRITON
from source_model_revamp.fbis.pairwise_collisions import pairwise_collision_metadata
from source_model_revamp.fbis.modal.lost_ion_distribution import EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY
from source_model_revamp.fbis.modal.source_projection import NO_CONFINED_FBIS_SOURCE, audit_attenuated_modal_source_projection, modal_source_projection_audit_metadata
from source_model_revamp.fbis.modal import ModalFBISBasis
from source_model_revamp.fbis.modal.species_state import FastIonSpeciesState, FastIonSystemState
from source_model_revamp.integration.pipeline_common import distribution_roundoff
from source_model_revamp.integration.pipeline_types import BeamEnsembleResult, GeometryStageResult, KineticStageResult
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.modal_stage.types import _ReusableModalBasisBundle, electron_density_input_is_initial_guess_only
from source_model_revamp.integration.modal_stage.basis_density import _background_positive_charge_on_confined_grid, _modal_numerics, build_reusable_modal_basis
from source_model_revamp.integration.modal_stage.density_closure import _eq42_profile_from_modal_state, _solve_modal_density_closure, _solve_modal_with_collision_state
from source_model_revamp.fbis.modal.density_weighting import EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, compare_eq42_density_profiles, eq42_density_weighting_model
from source_model_revamp.integration.modal_stage.full_device import _expander_trajectory_support_metadata, _full_device_fast_ion_population
from source_model_revamp.integration.modal_stage.end_loss_power import end_loss_power_availability_metadata
from source_model_revamp.integration.modal_stage.metadata_adapter import KineticMetadataContract
from source_model_revamp.integration.modal_stage.prompt_source import _basis_and_assessment, _prompt_only_kinetic_stage, _source_audit_power_metadata
from source_model_revamp.integration.evaluation import FINAL_ONLY_QUALIFICATION_CHECK_NAMES, OperatingPointEvaluationMode
from source_model_revamp.integration.modal_stage.system_stage import build_multispecies_modal_kinetic_stage

def build_modal_kinetic_stage(config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, *, electron_temperature_J: float | None = None, collision_state_recompute_policy: str | None = None, collision_state_is_self_consistent_electron_temperature: bool = False, modal_basis: _ReusableModalBasisBundle | ModalFBISBasis | None = None, eq59_warm_start_state: object | None = None, eq59_warm_start_state_by_species: dict[str, object] | None = None, eq70_initial_potential_energy_J: np.ndarray | None = None, initial_representative_fast_density_m3: float | None = None, initial_fast_D_midpoint_density_m3: float | None = None, initial_fast_D_confined_average_density_m3: float | None = None, initial_fast_T_midpoint_density_m3: float | None = None, initial_fast_T_confined_average_density_m3: float | None = None, initial_electron_midpoint_density_m3: float | None = None, initial_electron_confined_average_density_m3: float | None = None, initial_electron_parent_n0_m3: float | None = None, initial_electron_collision_density_m3: float | None = None, initial_eq42_density_profile: object | None = None, scalar_collision_operator_warm_start_by_species: dict[str, tuple[object, object, object, object]] | None = None, scalar_external_fast_field_warm_start_by_test_species: dict[str, tuple[object, ...]] | None = None, scalar_wall_barrier_energy_warm_start_J: float | None = None, evaluation_mode: OperatingPointEvaluationMode = OperatingPointEvaluationMode.FIXED_TEMPERATURE_FINAL) -> KineticStageResult:
    """
    Build the modal kinetic stage result for the current beam and plasma configuration
    
    NBI supported stationary closure or an active tritium beam is routed to the coupled D plus T system path
    The deuterium benchmark path performs source projection, scalar density closure, Eq 42 basis coupling, Eq 59 and Eq 70 solves, fixed magnetic loss reconstruction, full device mapping, convergence checks, and metadata assembly
    """
    active_beam_species = tuple(sorted(beam.attenuated_source_by_species))
    # Route stationary NBI closure and any active tritium source through the shared D plus T system solve
    if config.plasma_closure.model == "nbi_supported_stationary" or TRITON.species_id in active_beam_species:
        active_bundle = modal_basis if isinstance(modal_basis, _ReusableModalBasisBundle) else None
        initial_midpoint = {key: value for key, value in ((DEUTERON.species_id, initial_fast_D_midpoint_density_m3), (TRITON.species_id, initial_fast_T_midpoint_density_m3)) if value is not None}
        initial_average = {key: value for key, value in ((DEUTERON.species_id, initial_fast_D_confined_average_density_m3), (TRITON.species_id, initial_fast_T_confined_average_density_m3)) if value is not None}
        warm_by_species = dict(eq59_warm_start_state_by_species or {})
        if eq59_warm_start_state is not None and DEUTERON.species_id in active_beam_species and DEUTERON.species_id not in warm_by_species:
            warm_by_species[DEUTERON.species_id] = eq59_warm_start_state
        return build_multispecies_modal_kinetic_stage(config, geometry, beam, electron_temperature_J=electron_temperature_J, collision_state_recompute_policy=collision_state_recompute_policy, collision_state_is_self_consistent_electron_temperature=collision_state_is_self_consistent_electron_temperature, modal_basis=active_bundle, eq59_warm_start_state_by_species=warm_by_species, eq70_initial_potential_energy_J=eq70_initial_potential_energy_J, initial_fast_midpoint_density_m3_by_species=initial_midpoint, initial_fast_average_density_m3_by_species=initial_average, initial_electron_midpoint_density_m3=initial_electron_midpoint_density_m3, initial_electron_confined_average_density_m3=initial_electron_confined_average_density_m3, initial_electron_collision_density_m3=initial_electron_collision_density_m3, initial_eq42_density_profile=initial_eq42_density_profile, scalar_collision_operator_warm_start_by_species=scalar_collision_operator_warm_start_by_species, scalar_external_fast_field_warm_start_by_test_species=scalar_external_fast_field_warm_start_by_test_species, scalar_wall_barrier_energy_warm_start_J=scalar_wall_barrier_energy_warm_start_J, evaluation_mode=evaluation_mode)
    if active_beam_species != (DEUTERON.species_id,):
        raise ValueError("the modal operating point supports deuterium and tritium beam species only")
    evaluation_mode = OperatingPointEvaluationMode.parse(evaluation_mode)
    k = config.kinetic_electrostatic
    T_e_J = float(energy_J_from_keV(config.plasma_closure.electron_temperature_initial_guess_keV) if electron_temperature_J is None else electron_temperature_J)
    T_i_J = float(energy_J_from_keV(config.plasma_closure.background_ion_temperature_keV))
    numerics = _modal_numerics(config)
    if modal_basis is None:
        modal_basis = build_reusable_modal_basis(config, geometry)
    active_basis, _ = _basis_and_assessment(modal_basis)
    source_audit = audit_attenuated_modal_source_projection(attenuated_source=beam.attenuated_source, basis=active_basis, volume_m3=geometry.volume_m3)
    # A fully prompt source has an exact zero confined modal state
    if source_audit.source_state == NO_CONFINED_FBIS_SOURCE:
        prompt_result = _prompt_only_kinetic_stage(config, geometry, beam, electron_temperature_J=T_e_J, ion_temperature_J=T_i_J, numerics=numerics, modal_basis=modal_basis, source_audit=source_audit, collision_state_recompute_policy=collision_state_recompute_policy, collision_state_is_self_consistent_electron_temperature=collision_state_is_self_consistent_electron_temperature)
        return replace(prompt_result, metadata={**prompt_result.metadata, "operating_point_evaluation_mode": evaluation_mode.value, "final_only_qualification_status": "not_applicable_no_confined_source", "full_qualification_count": 0},)
    density_result = _solve_modal_density_closure(config=config, geometry=geometry, beam=beam, electron_temperature_J=T_e_J, ion_temperature_J=T_i_J, numerics=numerics, modal_basis=modal_basis, initial_eq59_warm_start_state=eq59_warm_start_state, initial_representative_fast_density_m3=initial_representative_fast_density_m3, initial_fast_D_midpoint_density_m3=initial_fast_D_midpoint_density_m3, initial_fast_D_confined_average_density_m3=initial_fast_D_confined_average_density_m3, initial_electron_midpoint_density_m3=initial_electron_midpoint_density_m3, initial_electron_confined_average_density_m3=initial_electron_confined_average_density_m3, initial_electron_parent_n0_m3=initial_electron_parent_n0_m3, initial_electron_collision_density_m3=initial_electron_collision_density_m3, initial_eq42_density_profile=initial_eq42_density_profile, initial_eq70_potential_energy_J=eq70_initial_potential_energy_J, evaluation_mode=evaluation_mode)
    collision_state = density_result.collision_state
    if density_result.modal_basis is None:
        raise RuntimeError("Eq 42 density closure did not return its active modal basis")
    modal_basis = density_result.modal_basis
    closure_modal = density_result.modal
    closure_eq42_metadata = {key: value for key, value in closure_modal.metadata.items() if key.startswith("modal_eq42_")}
    closure_potential = (None if closure_modal.electrostatic_profile is None else np.asarray(closure_modal.electrostatic_profile.potential_energy_J, dtype=float,))
    benchmark_density_closure = k.modal_density_closure_model == "egedal_beam_plasma_quasineutral"
    if evaluation_mode.runs_final_qualification:
        assessed_modal = _solve_modal_with_collision_state(config=config, geometry=geometry, beam=beam, electron_temperature_J=T_e_J, ion_temperature_J=T_i_J, electron_midplane_density_m3=density_result.electron_midplane_density_m3, electron_collision_density_m3=density_result.electron_collision_density_m3, eq70_prescribed_electron_midplane_density_m3=(density_result.electron_midplane_density_m3 if benchmark_density_closure else None), collision_state=collision_state, numerics=numerics, modal_basis=modal_basis, eq59_warm_start_state=closure_modal.eq59_warm_start_state, eq59_warm_start_compatibility_policy="closure_continuation", eq70_initial_potential_energy_J=closure_potential, run_numerical_convergence_assessments=True, fixed_boundary_reconstruction_policy="required", evaluation_mode=evaluation_mode)
    else:
        assessed_modal = replace(closure_modal, metadata={**closure_modal.metadata, "operating_point_evaluation_mode": evaluation_mode.value, "final_only_qualification_status": "deferred_to_final_qualification", "final_only_qualification_check_names": list(FINAL_ONLY_QUALIFICATION_CHECK_NAMES)},)
    modal_metadata = dict(assessed_modal.metadata)
    modal_metadata.update(closure_eq42_metadata)
    modal = replace(assessed_modal, metadata=modal_metadata)
    closure_density_state = density_result.density_state
    assessed_density_state = density_result.density_state
    assessed_fixed_point_change = 0.0
    assessed_fixed_point_converged = bool(density_result.converged)
    assessed_history = list(density_result.history)
    assessed_eq42_converged = bool(density_result.eq42_density_basis_converged)
    assessed_eq42_profile_change: float | None = None
    assessed_eq42_history = list(density_result.eq42_density_basis_history or [])
    if eq42_density_weighting_model(k.modal_eq42_density_weighting_model) == EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY:
        active_profile = density_result.eq42_density_profile
        if active_profile is None or assessed_density_state is None:
            assessed_eq42_converged = False
        else:
            assessed_target, assessed_source_converged, assessed_source_status = _eq42_profile_from_modal_state(config=config, geometry=geometry, modal=modal, numerics=numerics, density_state=assessed_density_state)
            assessed_profile_comparison = compare_eq42_density_profiles(active_profile, assessed_target)
            assessed_eq42_profile_change = float(assessed_profile_comparison.volume_weighted_l2_relative_change)
            assessed_eq42_converged = bool(density_result.eq42_density_basis_converged and assessed_source_converged and assessed_target.symmetric_within_tolerance and assessed_eq42_profile_change <= float(numerics.eq42_density_shape_relative_tolerance))
            assessed_eq42_history.append({"iteration": "post_solver_eq42_assessment", "status": "converged" if assessed_eq42_converged else "not_converged", "profile_source_status": assessed_source_status, "density_shape_volume_weighted_l2_relative_change": assessed_eq42_profile_change, "density_shape_relative_tolerance": float(numerics.eq42_density_shape_relative_tolerance), "profile_symmetric": assessed_target.symmetric_within_tolerance, "overall_eq42_outer_iteration_converged": assessed_eq42_converged})
            updated_modal_metadata = dict(modal.metadata)
            updated_modal_metadata["modal_eq42_outer_iteration_converged"] = assessed_eq42_converged
            updated_modal_metadata["modal_eq42_density_shape_final_relative_change"] = assessed_eq42_profile_change
            updated_modal_metadata["modal_eq42_final_assessment_uses_final_density_state"] = True
            updated_modal_metadata["modal_eq42_final_assessment_profile_source"] = assessed_target.source
            updated_modal_metadata["modal_eq42_outer_iteration_history"] = assessed_eq42_history
            if not assessed_eq42_converged:
                updated_modal_metadata["modal_eq42_outer_iteration_failure_reason"] = "final_assessed_density_shape_changed_above_tolerance"
            modal = replace(modal, metadata=updated_modal_metadata)
    density_result = replace(density_result, modal=modal, density_state=assessed_density_state, electron_midplane_density_m3=(density_result.electron_midplane_density_m3 if assessed_density_state is None else assessed_density_state.electron_midplane_density_m3), converged=assessed_fixed_point_converged, relative_error=max(float(density_result.relative_error), assessed_fixed_point_change), history=assessed_history, eq42_density_basis_converged=assessed_eq42_converged, eq42_density_basis_history=assessed_eq42_history)
    active_basis, _ = _basis_and_assessment(modal_basis)
    source_audit = audit_attenuated_modal_source_projection(attenuated_source=beam.attenuated_source, basis=active_basis, volume_m3=geometry.volume_m3)
    background_positive_profile, background_left_throat, background_right_throat, background_profile_source = (_background_positive_charge_on_confined_grid(geometry))
    raw_final_distribution = np.asarray(modal.distribution_v_lambda, dtype=float)
    final_distribution, roundoff_adjusted_count, raw_min_distribution = distribution_roundoff(raw_final_distribution, relative_tolerance=k.negative_distribution_roundoff_tolerance, name="modal fast-D distribution")
    final_inventory_sum = float(np.sum(final_distribution))
    effective_injection_angle_deg = float(beam.metadata.get("beam_effective_injection_angle_deg", 0.0))
    solved_collision_state = rebuild_fbis_collision_parameter_state_with_fast_density(collision_state, modal.density_m3)
    collision_metadata = collision_parameter_metadata(collision_state)
    collision_metadata.update(pairwise_collision_metadata(solved_collision_state.pairwise_collision_state))
    if collision_state_recompute_policy is not None:
        collision_metadata["collision_state_recompute_policy"] = str(collision_state_recompute_policy)
    if collision_state_is_self_consistent_electron_temperature:
        collision_metadata.update({"collision_state_temperature_closure": "self_consistent_electron_temperature_power_balance", "collision_state_is_self_consistent_electron_temperature": True, "plasma_temperature_closure_status": "self_consistent_electron_temperature_power_balance",})
    feedback_model = str(getattr(k, "electrostatic_feedback_model", "magnetic_only_fbis")).strip().lower()
    fixed_boundary_selected = (str(k.ion_loss_closure_model).strip().lower() == EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY)
    zero_potential_reference = feedback_model == "zero_potential_magnetic_reference"
    phi_profile_applicable = not zero_potential_reference
    current_balance_applicable = not zero_potential_reference
    phi_profile_present = modal.electrostatic_profile is not None
    phi_profile_converged = bool(phi_profile_present and modal.electrostatic_profile.converged)
    phi_profile_numerically_valid = bool(phi_profile_converged)
    eq71_low_energy_applicability_valid = bool(modal.metadata.get("modal_phi_z_low_energy_approximation_valid", False) and not modal.metadata.get("modal_phi_z_eq71_closed_interval_intersects_distribution_support", False,))
    eq71_throat_boundary_applicability_valid = bool(modal.metadata.get("modal_phi_z_effective_potential_throat_boundary_valid", False,))
    modal_total_birth_rate_s = float(modal.metadata.get("modal_total_beam_birth_rate_s", 0.0))
    modal_source_birth_rate_s = float(modal.metadata.get("modal_source_birth_rate_s", 0.0))
    modal_projected_source_rate_s = float(modal.metadata.get("modal_projected_source_particle_rate_from_coefficients_s", 0.0))
    source_projection_nonzero = modal_source_birth_rate_s <= 0.0 or modal_projected_source_rate_s > 0.0
    source_projection_relative_error = 0.0 if modal_source_birth_rate_s <= 0.0 else (modal_projected_source_rate_s - modal_source_birth_rate_s) / modal_source_birth_rate_s
    source_projection_tolerance = float(getattr(k, "modal_source_projection_relative_tolerance", 1.0e-2))
    source_projection_converged = bool(source_projection_nonzero and abs(source_projection_relative_error) <= source_projection_tolerance)
    source_audit_metadata = modal_source_projection_audit_metadata(source_audit)
    source_power_metadata = _source_audit_power_metadata(beam, source_audit)
    audit_rate_scale = max(abs(float(source_audit.total_deposited_birth_rate_s)), abs(modal_total_birth_rate_s), 1.0)
    if not np.isclose(float(source_audit.total_deposited_birth_rate_s), modal_total_birth_rate_s, rtol=1.0e-12, atol=1.0e-12 * audit_rate_scale):
        raise RuntimeError("source audit and modal solve disagree on total deposited birth rate")
    if not np.isclose(float(source_audit.confined_physical_birth_rate_s), modal_source_birth_rate_s, rtol=1.0e-12, atol=1.0e-12 * audit_rate_scale):
        raise RuntimeError("source audit and modal solve disagree on confined physical birth rate")
    if not np.isclose(float(source_audit.modal_projected_confined_rate_s), modal_projected_source_rate_s, rtol=1.0e-12, atol=1.0e-12 * audit_rate_scale,):
        raise RuntimeError("source audit and modal solve disagree on projected confined rate")
    mirror_throat_inside_domain = (geometry.metadata.get("fitted_mirror_throat_inside_kinetic_domain") is True)
    beam_injection_geometry_consistent = (beam.metadata.get("beam_injection_geometry_consistent") is True)
    beam_pitch_selection_valid = (beam.metadata.get("beam_pitch_selection_valid") is True)
    boundary_selection_valid = bool(geometry.plasma_boundaries is not None and geometry.plasma_boundaries.valid and geometry.B_T_function is not None and geometry.full_device_z_centers_m is not None)
    left_mirror_ratio = geometry.metadata.get("left_fitted_mirror_ratio")
    right_mirror_ratio = geometry.metadata.get("right_fitted_mirror_ratio")
    if left_mirror_ratio is None or right_mirror_ratio is None:
        end_symmetry_relative_error = 0.0
        symmetric_end_model_applicable = True
    else:
        left_mirror_ratio = float(left_mirror_ratio)
        right_mirror_ratio = float(right_mirror_ratio)
        end_symmetry_relative_error = abs(left_mirror_ratio - right_mirror_ratio) / max(abs(left_mirror_ratio), abs(right_mirror_ratio), np.finfo(float).tiny)
        symmetric_end_model_applicable = bool(end_symmetry_relative_error <= numerics.loss_convention_relative_tolerance)
    geometry_boundary_valid = bool(mirror_throat_inside_domain and beam_pitch_selection_valid and boundary_selection_valid)
    basis_convergence_assessed = (modal.metadata.get("modal_basis_convergence_assessed") is True)
    basis_converged = modal.metadata.get("modal_basis_converged") is True
    eq59_applicable = (modal.metadata.get("modal_rosenbluth_convergence_applicable") is True)
    eq59_convergence_assessed = (modal.metadata.get("modal_rosenbluth_convergence_assessed") is True)
    eq59_converged = modal.metadata.get("modal_rosenbluth_converged") is True
    eq59_base_converged = bool(not eq59_applicable or (modal.metadata.get("modal_rosenbluth_fixed_point_iteration_converged") is True and modal.metadata.get("modal_rosenbluth_nonlinear_residual_converged") is True and modal.metadata.get("modal_rosenbluth_inventory_converged") is True and modal.metadata.get("modal_rosenbluth_effective_energy_converged") is True and modal.metadata.get("modal_rosenbluth_single_grid_tail_indicator_passed") is True))
    retained_mode_applicable = int(modal.metadata.get("modal_n_physical_modes", 0)) >= 2
    retained_mode_assessed = (modal.metadata.get("modal_retained_mode_count_assessed") is True)
    retained_mode_converged = (modal.metadata.get("modal_retained_mode_count_converged") is True)
    local_speed_refinement_applicable = bool(phi_profile_applicable and modal.local_reconstruction is not None)
    local_speed_refinement_assessed = (modal.metadata.get("modal_local_speed_refinement_assessed") is True)
    local_speed_refinement_converged = (modal.metadata.get("modal_local_speed_refinement_converged") is True)
    loss_convention_valid = (modal.metadata.get("modal_loss_convention_validation_passed") is True)
    published_feedback_requested = feedback_model == "egedal_2022_published_electrostatic_fbis_approximation"
    published_population_fraction = modal.metadata.get("published_lambda1_heuristic_population_fraction")
    published_heuristic_valid = bool(fixed_boundary_selected or not published_feedback_requested or (modal.metadata.get("electrostatic_feedback_iteration_converged") is True and published_population_fraction is not None and np.isfinite(float(published_population_fraction)) and float(published_population_fraction) >= 0.0 and modal.metadata.get("published_lambda1_heuristic_higher_modes_modified") is False))
    current_balance_tolerance = float(getattr(k, "modal_current_balance_relative_tolerance", numerics.loss_convention_relative_tolerance))
    current_balance_relative_error = modal.metadata.get("modal_electron_current_balance_relative_error", modal.metadata.get("modal_current_balance_relative_error"),)
    current_balance_valid = bool(np.isfinite(current_balance_tolerance) and current_balance_tolerance >= 0.0 and current_balance_relative_error is not None and np.isfinite(float(current_balance_relative_error)) and abs(float(current_balance_relative_error)) <= current_balance_tolerance)
    expander_metadata = _expander_trajectory_support_metadata(geometry, modal=modal, material_loss_rate_s=float(modal.ion_particle_loss_rate_s),)
    full_device_fast_distribution, full_device_population_metadata = (_full_device_fast_ion_population(geometry=geometry, beam=beam, modal=modal, collision_state=collision_state, quadrature_order=numerics.local_velocity_quadrature_order,))
    expander_population_valid = bool(expander_metadata["expander_unclassified_population_check_passed"])
    full_device_lost_population_applicable = bool(fixed_boundary_selected and not zero_potential_reference and modal.ion_particle_loss_rate_s > 0.0 and boundary_selection_valid)
    full_device_lost_population_required = full_device_lost_population_applicable
    full_device_lost_reference_compatible = bool(fixed_boundary_selected and full_device_population_metadata.get("full_device_lost_ion_reconstruction_matches_active_hot_or_electrostatic_sink") is True and full_device_population_metadata.get("fixed_boundary_local_speed_coverage_passed") is True)
    full_device_lost_population_valid = bool(full_device_lost_population_applicable and full_device_population_metadata["full_device_lost_ion_reconstruction_available"] and full_device_lost_reference_compatible)
    eq42_density_weighting_coupled = bool(modal.metadata.get("modal_eq42_density_weighting_coupled") is True)
    eq42_density_profile_symmetric = bool(modal.metadata.get("modal_eq42_density_profile_symmetric") is True)
    eq42_outer_iteration_applicable = bool(modal.metadata.get("modal_eq42_outer_iteration_applicable") is True)
    eq42_outer_iteration_converged = bool(modal.metadata.get("modal_eq42_outer_iteration_converged") is True)
    eq42_density_basis_valid = bool(eq42_density_weighting_coupled and eq42_density_profile_symmetric and (not eq42_outer_iteration_applicable or eq42_outer_iteration_converged))
    active_global_loss_conservation_applicable = bool(fixed_boundary_selected)
    density_state_required = not benchmark_density_closure
    density_state = density_result.density_state
    density_identity_tolerance = max(128.0 * np.finfo(float).eps, 1.0e-10)
    density_inventory_tolerance = float(numerics.phi_z_relative_tolerance)
    density_state_available = bool(not density_state_required or density_state is not None)
    exact_midplane_quasineutrality_valid = bool(not density_state_required or (density_state is not None and density_state.exact_midplane_quasineutrality_relative_error <= density_identity_tolerance))
    confined_inventory_consistent = bool(not density_state_required or (density_state is not None and density_state.confined_inventory_relative_error <= density_inventory_tolerance))
    collision_composition_consistent = not density_state_required
    midpoint_role_relative_error = None
    collision_role_relative_error = None
    density_roles_separated = not density_state_required
    required_check_status = {
        "geometry_boundary": geometry_boundary_valid,
        "symmetric_end_model_applicability": symmetric_end_model_applicable,
        "source_projection": source_projection_converged,
        "basis": bool(evaluation_mode is OperatingPointEvaluationMode.TEMPERATURE_TRIAL or (basis_convergence_assessed and basis_converged)),
        "retained_mode_count": bool(evaluation_mode is OperatingPointEvaluationMode.TEMPERATURE_TRIAL or not retained_mode_applicable or (retained_mode_assessed and retained_mode_converged)),
        "eq42_density_weighted_basis": eq42_density_basis_valid,
        "eq59": bool(eq59_base_converged if evaluation_mode is OperatingPointEvaluationMode.TEMPERATURE_TRIAL else (not eq59_applicable or (eq59_convergence_assessed and eq59_converged))),
        "eq70_profile": bool(not phi_profile_applicable or phi_profile_numerically_valid),
        "local_speed_grid": bool(evaluation_mode is OperatingPointEvaluationMode.TEMPERATURE_TRIAL or not local_speed_refinement_applicable or (local_speed_refinement_assessed and local_speed_refinement_converged)),
        "eq71_low_energy_applicability": bool(fixed_boundary_selected or not phi_profile_applicable or eq71_low_energy_applicability_valid),
        "eq71_throat_boundary_applicability": bool( fixed_boundary_selected  or not phi_profile_applicable  or eq71_throat_boundary_applicability_valid),
        "published_heuristic": published_heuristic_valid,
        "expander_unclassified_population": bool(zero_potential_reference or expander_population_valid),
        "full_device_lost_population": bool(not full_device_lost_population_applicable or full_device_lost_population_valid),
        "density_closure": bool(density_result.converged),
        "density_state_available": density_state_available,
        "density_roles_separated": density_roles_separated,
        "exact_midplane_quasineutrality": exact_midplane_quasineutrality_valid,
        "electron_ion_confined_inventory": confined_inventory_consistent,
        "collision_composition": collision_composition_consistent,
        "loss_and_global_conservation": bool(not active_global_loss_conservation_applicable or loss_convention_valid),
        "current_balance": bool(not current_balance_applicable or current_balance_valid),
    }
    numerical_check_status = {
        "source_projection": source_projection_converged,
        "basis": bool(evaluation_mode is OperatingPointEvaluationMode.TEMPERATURE_TRIAL or (basis_convergence_assessed and basis_converged)),
        "retained_mode_count": bool(evaluation_mode is OperatingPointEvaluationMode.TEMPERATURE_TRIAL or not retained_mode_applicable or (retained_mode_assessed and retained_mode_converged)),
        "eq42_density_weighted_basis": eq42_density_basis_valid,
        "eq59": bool(eq59_base_converged if evaluation_mode is OperatingPointEvaluationMode.TEMPERATURE_TRIAL else (not eq59_applicable or (eq59_convergence_assessed and eq59_converged))),
        "eq70_profile": bool(not phi_profile_applicable or phi_profile_numerically_valid),
        "local_speed_grid": bool(evaluation_mode is OperatingPointEvaluationMode.TEMPERATURE_TRIAL or not local_speed_refinement_applicable or (local_speed_refinement_assessed and local_speed_refinement_converged)),
        "density_closure": bool(density_result.converged),
        "density_state_available": density_state_available,
        "density_roles_separated": density_roles_separated,
        "exact_midplane_quasineutrality": exact_midplane_quasineutrality_valid,
        "electron_ion_confined_inventory": confined_inventory_consistent,
        "collision_composition": collision_composition_consistent,
    }
    if evaluation_mode is OperatingPointEvaluationMode.TEMPERATURE_TRIAL:
        for check_name in ("basis", "retained_mode_count", "local_speed_grid", "published_heuristic", "expander_unclassified_population"):
            required_check_status.pop(check_name, None)
            numerical_check_status.pop(check_name, None)
        required_check_status.pop("eq59", None)
        numerical_check_status.pop("eq59", None)
        required_check_status["eq59_base_grid"] = eq59_base_converged
        numerical_check_status["eq59_base_grid"] = eq59_base_converged
    numerical_convergence_passed = bool(all(numerical_check_status.values()))
    evaluated_closures_converged = bool(all(required_check_status.values()))
    convergence_failures: list[str] = []
    if not geometry_boundary_valid:
        convergence_failures.append("geometry_or_plasma_boundary_invalid")
    if not symmetric_end_model_applicable:
        convergence_failures.append("asymmetric_branch_kinetics_not_implemented")
    if not source_projection_nonzero:
        convergence_failures.append("nonzero_beam_birth_projected_to_zero_modal_source")
    elif not source_projection_converged:
        convergence_failures.append("modal_source_projection_rate_not_converged")
    if not eq42_density_weighting_coupled:
        convergence_failures.append("modal_eq42_density_weighting_not_coupled")
    elif not eq42_density_profile_symmetric:
        convergence_failures.append("modal_eq42_density_profile_not_symmetric")
    elif eq42_outer_iteration_applicable and not eq42_outer_iteration_converged:
        convergence_failures.append("modal_eq42_density_basis_iteration_not_converged")
    if evaluation_mode.runs_final_qualification and not basis_convergence_assessed:
        convergence_failures.append("modal_basis_grid_convergence_not_assessed")
    elif evaluation_mode.runs_final_qualification and not basis_converged:
        convergence_failures.append("modal_basis_grid_convergence_failed")
    if evaluation_mode.runs_final_qualification and retained_mode_applicable and not retained_mode_assessed:
        convergence_failures.append("modal_retained_mode_count_not_assessed")
    elif evaluation_mode.runs_final_qualification and retained_mode_applicable and not retained_mode_converged:
        convergence_failures.append("modal_retained_mode_count_not_converged")
    if evaluation_mode.runs_final_qualification and eq59_applicable and not eq59_convergence_assessed:
        convergence_failures.append("modal_eq59_nonlinear_convergence_not_assessed")
    elif evaluation_mode.runs_final_qualification and eq59_applicable and not eq59_converged:
        convergence_failures.append("modal_eq59_nonlinear_or_speed_convergence_failed")
    if phi_profile_applicable and not phi_profile_numerically_valid:
        convergence_failures.append("modal_eq70_profile_invalid")
    if evaluation_mode.runs_final_qualification and local_speed_refinement_applicable and not local_speed_refinement_assessed:
        convergence_failures.append("modal_local_speed_refinement_not_assessed")
    elif evaluation_mode.runs_final_qualification and local_speed_refinement_applicable and not local_speed_refinement_converged:
        convergence_failures.append("modal_local_speed_refinement_not_converged")
    if (phi_profile_applicable and not fixed_boundary_selected and not eq71_low_energy_applicability_valid):
        convergence_failures.append("modal_eq71_low_energy_or_closed_interval_applicability_invalid")
    if (phi_profile_applicable and not fixed_boundary_selected and not eq71_throat_boundary_applicability_valid):
        convergence_failures.append("modal_eq71_throat_only_boundary_applicability_invalid")
    if not published_heuristic_valid:
        convergence_failures.append("published_electrostatic_heuristic_status_invalid")
    if not density_result.converged:
        convergence_failures.append("modal_collision_density_closure_not_converged")
    if not density_state_available:
        convergence_failures.append("operating_point_density_state_unavailable")
    if not exact_midplane_quasineutrality_valid:
        convergence_failures.append("exact_midplane_quasineutrality_identity_failed")
    if not confined_inventory_consistent:
        convergence_failures.append("electron_ion_confined_inventory_identity_failed")
    if not collision_composition_consistent:
        convergence_failures.append("collision_composition_not_converged")
    if not density_roles_separated:
        convergence_failures.append("electron_density_roles_not_consistent")
    if current_balance_applicable and not current_balance_valid:
        convergence_failures.append("modal_current_balance_not_converged")
    rosenbluth_history = list(modal.metadata.get("modal_rosenbluth_iteration_history", []))
    rosenbluth_final_change = None
    if rosenbluth_history:
        rosenbluth_final_change = float(rosenbluth_history[-1].get("relative_change", float("nan")))
    birth_lambda_profile = np.asarray(beam.metadata.get("beam_birth_lambda_profile", []), dtype=float)
    if birth_lambda_profile.size:
        if np.any(~np.isfinite(birth_lambda_profile)) or np.any((birth_lambda_profile < 0.0) | (birth_lambda_profile > 1.0)):
            raise ValueError("beam birth Lambda profile must remain finite and inside [0, 1]")
        birth_eta_profile = np.asarray(modal.basis.eta_lambda_map.eta_of_lambda(birth_lambda_profile), dtype=float)
        birth_eta_status = "resolved_from_geometry_dependent_Egedal_eq47_eta_lambda_map"
    else:
        birth_eta_profile = np.empty(0, dtype=float)
        birth_eta_status = "unavailable_no_beam_birth_lambda_profile"
    density_state_metadata = {} if density_state is None else density_state.as_metadata()
    kinetic_state_finite = bool(np.all(np.isfinite(final_distribution)) and (density_state is None or (np.all(np.isfinite(density_state.electron_cell_density_m3)) and np.all(np.isfinite(density_state.fast_deuterium_cell_density_m3)))) and (modal.electrostatic_profile is None or np.all(np.isfinite(modal.electrostatic_profile.potential_energy_J))))
    kinetic_state_physically_evaluable = bool(kinetic_state_finite and (density_state is None or (np.all(density_state.electron_cell_density_m3 >= 0.0) and np.all(density_state.fast_deuterium_cell_density_m3 >= 0.0))))
    kinetic_hard_model_applicability_passed = bool(kinetic_state_physically_evaluable and geometry_boundary_valid and symmetric_end_model_applicable)
    metadata_contract = KineticMetadataContract(exact_midplane_quasineutrality_check_passed=exact_midplane_quasineutrality_valid, modal_source_projection_rate_converged=source_projection_converged, modal_source_projection_relative_error=float(source_projection_relative_error), modal_source_projection_relative_tolerance=source_projection_tolerance, kinetic_min_distribution_value=float(np.min(final_distribution)) if final_distribution.size else 0.0, kinetic_geometry_boundary_check_passed=geometry_boundary_valid, kinetic_symmetric_end_model_applicable=symmetric_end_model_applicable, density_roles_separated=density_roles_separated, electron_ion_confined_inventory_check_passed=confined_inventory_consistent, kinetic_convergence_status="converged" if evaluated_closures_converged else "not_converged", kinetic_convergence_failure_reasons=tuple(convergence_failures))
    metadata = {'kinetic_model': 'egedal_modal_fbis', 'operating_point_evaluation_mode': evaluation_mode.value, 'final_only_qualification_status': 'completed' if evaluation_mode.runs_final_qualification else 'deferred_to_final_qualification', 'full_qualification_count': int(evaluation_mode.runs_final_qualification), 'backend_numerics_settings': k.backend_numerics_metadata(), **collision_metadata, 'modal_density_closure_model': k.modal_density_closure_model, 'modal_density_closure_converged': bool(density_result.converged), 'modal_density_closure_iterations': int(density_result.iterations), 'modal_density_closure_relative_error': float(density_result.relative_error), 'modal_density_closure_relative_tolerance': float(k.modal_density_closure_relative_tolerance), 'modal_density_closure_history': density_result.history, 'modal_density_closure_electron_collision_density_m3': float(density_result.electron_collision_density_m3), 'electron_density_input_is_initial_guess_only': electron_density_input_is_initial_guess_only(k.modal_density_closure_model), 'electron_density_initial_guess_m3': config.plasma_closure.electron_density_initial_guess_m3, 'configured_seed_reference_midpoint_positive_charge_m3': float(config.plasma_closure.background_deuterium_midplane_density_m3 + config.plasma_closure.background_tritium_midplane_density_m3), 'electron_density_initial_guess_valid': config.plasma_closure.electron_density_initial_guess_valid, 'electron_density_initial_guess_source': config.plasma_closure.electron_density_initial_guess_source, 'electron_temperature_initial_guess_keV': config.plasma_closure.electron_temperature_initial_guess_keV, 'background_ion_temperature_keV': config.plasma_closure.background_ion_temperature_keV, 'background_profile_model': geometry.metadata.get('background_profile_model'), 'background_profile_midplane_normalized': geometry.metadata.get('background_profile_midplane_normalized'), 'background_profile_symmetry_relative_error': geometry.metadata.get('background_profile_symmetry_relative_error'), 'background_profile_symmetry_tolerance': geometry.metadata.get('background_profile_symmetry_tolerance'), 'background_profile_valid_and_midplane_normalized': bool(geometry.metadata.get('background_profile_midplane_normalized') is True and geometry.metadata.get('background_profile_finite') is True and (geometry.metadata.get('background_profile_nonnegative') is True) and (geometry.metadata.get('background_profile_symmetric_within_tolerance') is True)), 'background_thermal_ion_end_loss_available': False, 'global_plasma_particle_power_balance_available': False, 'electron_temperature_power_balance_scope': 'fast_D_electron_FBIS_subsystem', 'global_plasma_power_balance_available': False, 'background_positive_charge_profile_m3': None if background_positive_profile is None else background_positive_profile, 'background_left_throat_positive_charge_density_m3': background_left_throat, 'background_right_throat_positive_charge_density_m3': background_right_throat, 'density_roles_separated': density_roles_separated, 'density_identity_relative_tolerance': density_identity_tolerance, 'density_inventory_relative_tolerance': density_inventory_tolerance, 'exact_midplane_quasineutrality_check_passed': exact_midplane_quasineutrality_valid, 'electron_ion_confined_inventory_check_passed': confined_inventory_consistent, 'collision_composition_check_passed': collision_composition_consistent, 'kinetic_convergence_status': 'converged' if evaluated_closures_converged else 'not_converged', 'state_finite': kinetic_state_finite, 'state_physically_evaluable': kinetic_state_physically_evaluable, 'kinetic_eq59_state_finite': bool(np.all(np.isfinite(final_distribution))), 'kinetic_hard_model_applicability_passed': kinetic_hard_model_applicability_passed, 'kinetic_numerical_convergence_passed': numerical_convergence_passed, 'kinetic_iterations': int(density_result.iterations), 'kinetic_all_inner_steps_converged': evaluated_closures_converged, 'modal_rosenbluth_convergence_assessed': eq59_convergence_assessed, 'modal_rosenbluth_converged': eq59_converged, 'modal_source_projection_nonzero_for_nonzero_birth': source_projection_nonzero, 'modal_source_projection_relative_error': float(source_projection_relative_error), 'modal_source_projection_relative_tolerance': source_projection_tolerance, 'kinetic_geometry_contains_fitted_mirror_throat': mirror_throat_inside_domain, 'kinetic_geometry_boundary_check_passed': geometry_boundary_valid, 'kinetic_left_right_mirror_ratio_relative_difference': end_symmetry_relative_error, 'kinetic_symmetric_end_model_applicable': symmetric_end_model_applicable, 'kinetic_boundary_selection_valid': boundary_selection_valid, 'kinetic_beam_pitch_selection_valid': beam_pitch_selection_valid, 'kinetic_beam_injection_geometry_consistent': beam_injection_geometry_consistent, 'kinetic_basis_convergence_assessed': basis_convergence_assessed, 'kinetic_basis_converged': basis_converged, 'kinetic_retained_mode_count_applicable': retained_mode_applicable, 'kinetic_eq59_convergence_applicable': eq59_applicable, 'kinetic_eq59_convergence_assessed': eq59_convergence_assessed, 'kinetic_eq59_converged': eq59_converged, 'kinetic_total_device_loss_convention_valid': loss_convention_valid, 'kinetic_eq70_profile_present': phi_profile_present, 'kinetic_eq70_profile_applicable': phi_profile_applicable, 'kinetic_local_speed_refinement_applicable': local_speed_refinement_applicable, 'kinetic_published_heuristic_requested': bool(published_feedback_requested and (not fixed_boundary_selected)), 'kinetic_eq71_72_loss_boundary_diagnostic_only': bool(fixed_boundary_selected), 'kinetic_published_heuristic_check_passed': published_heuristic_valid, 'kinetic_expander_unclassified_population_check_passed': expander_population_valid, 'kinetic_full_device_lost_population_applicable': full_device_lost_population_applicable, 'kinetic_full_device_lost_population_required': full_device_lost_population_required, 'kinetic_full_device_lost_population_check_passed': full_device_lost_population_valid, 'kinetic_full_device_lost_population_reference_compatible': full_device_lost_reference_compatible, 'closed_interval_population_fraction': None, 'closed_interval_population_fraction_status': 'diagnostic_only_not_part_of_active_fixed_magnetic_loss_boundary' if fixed_boundary_selected else 'unresolved_energy_dependent_loss_operator_required' if bool(modal.metadata.get('modal_phi_z_eq71_closed_interval_intersects_distribution_support', False)) else 'not_activated_for_returned_profile', 'local_well_population_fraction': expander_metadata.get('locally_trapped_expander_fraction'), 'unclassified_population_fraction': expander_metadata.get('unclassified_expander_fraction'), 'kinetic_current_balance_applicable': current_balance_applicable, 'kinetic_current_balance_relative_tolerance': current_balance_tolerance, 'kinetic_effective_injection_angle_deg': effective_injection_angle_deg, 'beam_birth_eta_profile': birth_eta_profile, 'beam_birth_eta_profile_status': birth_eta_status, 'kinetic_raw_final_distribution_min': float(raw_min_distribution), 'kinetic_negative_cell_count': 0, 'kinetic_roundoff_negative_cells_zeroed_for_fusion': int(roundoff_adjusted_count), 'solved_fast_D_inventory_v_lambda_sum': final_inventory_sum, 'solved_fast_D_volume_average_density_m3_from_inventory_sum': modal.density_m3, 'solved_fast_D_inventory_particles': modal.inventory_particles, 'solved_fast_D_ion_particle_loss_rate_s': modal.ion_particle_loss_rate_s, 'solved_fast_D_ion_confinement_time_s': modal.confinement_time_s, 'solved_electrostatic_potential_min_V': modal.electron_wall_potential_relative_to_midplane_V, 'solved_electrostatic_potential_max_V': 0.0, **modal.metadata, **expander_metadata, **full_device_population_metadata, **source_audit_metadata, **source_power_metadata, **density_state_metadata}
    metadata.update({"eq68_root_residual_passed": current_balance_valid, "eq68_root_relative_residual": current_balance_relative_error, **metadata_contract.as_metadata()})
    metadata.update(end_loss_power_availability_metadata(metadata,fixed_boundary_selected=fixed_boundary_selected,current_balance_applicable=current_balance_applicable,current_balance_valid=current_balance_valid,loss_convention_valid=loss_convention_valid,full_device_lost_population_applicable=(full_device_lost_population_applicable),full_device_lost_population_valid=full_device_lost_population_valid))
    if isinstance(modal_basis, _ReusableModalBasisBundle) and modal_basis.basis_cache is not None:
        metadata.update(modal_basis.basis_cache.as_metadata())
    fast_species_state = FastIonSpeciesState(species=modal.species, speed_grid=modal.speed_grid, active=True, source_particle_rate_s=float(beam.total_fast_birth_rate_by_species[modal.species.species_id]), modal_result=modal, collision_state=solved_collision_state, full_device_local_distribution_z_v_pitch=full_device_fast_distribution, status="active_solved")
    fast_ion_system_state = FastIonSystemState(species_states={modal.species.species_id: fast_species_state}, shared_basis=modal.basis, shared_current_balance=modal.current_balance, shared_electrostatic_profile=modal.electrostatic_profile, electron_wall_power_loss_W=modal.electron_wall_power_loss_W, status="active_single_species", metadata={"fast_ion_operating_point_model": "single_species_through_fast_ion_system", "active_fast_species": (modal.species.species_id,), "shared_eq42_basis": True, "shared_eq68_wall_barrier": modal.current_balance is not None, "shared_eq70_potential": modal.electrostatic_profile is not None})

    return KineticStageResult(speed_grid=beam.speed_grid, lambda_grid=beam.lambda_grid, pitch_grid=beam.pitch_grid, final_distribution_v_lambda=final_distribution, metadata=metadata, local_distribution_z_v_lambda=None if modal.local_reconstruction is None else np.asarray(modal.local_reconstruction.local_distribution_z_v_lambda, dtype=float), local_distribution_z_v_pitch=None if modal.local_reconstruction is None else np.asarray(modal.local_reconstruction.local_distribution_z_v_pitch, dtype=float), local_density_m3=None if modal.local_reconstruction is None else np.asarray(modal.local_reconstruction.local_density_m3, dtype=float), full_device_local_distribution_z_v_pitch=full_device_fast_distribution, local_speed_grid=None if modal.local_reconstruction is None else modal.local_reconstruction.local_speed_grid, eq59_warm_start_state=modal.eq59_warm_start_state, operating_point_density_state=density_state, modal_basis_warm_start=modal_basis, eq70_potential_energy_warm_start_J=None if modal.electrostatic_profile is None else np.asarray(modal.electrostatic_profile.potential_energy_J, dtype=float), eq42_density_profile_warm_start=density_result.eq42_density_profile, fast_ion_system_state=fast_ion_system_state, speed_grid_by_species={modal.species.species_id: modal.speed_grid}, final_distribution_v_lambda_by_species={modal.species.species_id: final_distribution}, local_speed_grid_by_species=None if modal.local_reconstruction is None else {modal.species.species_id: modal.local_reconstruction.local_speed_grid}, local_distribution_z_v_lambda_by_species=None if modal.local_reconstruction is None else {modal.species.species_id: np.asarray(modal.local_reconstruction.local_distribution_z_v_lambda, dtype=float)}, local_distribution_z_v_pitch_by_species=None if modal.local_reconstruction is None else {modal.species.species_id: np.asarray(modal.local_reconstruction.local_distribution_z_v_pitch, dtype=float)}, local_density_m3_by_species=None if modal.local_reconstruction is None else {modal.species.species_id: np.asarray(modal.local_reconstruction.local_density_m3, dtype=float)}, full_device_local_distribution_z_v_pitch_by_species=None if full_device_fast_distribution is None else {modal.species.species_id: full_device_fast_distribution}, eq59_warm_start_state_by_species=None if modal.eq59_warm_start_state is None else {modal.species.species_id: modal.eq59_warm_start_state})
