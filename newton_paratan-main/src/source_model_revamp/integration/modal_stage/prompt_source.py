"""Prompt magnetic loss handling when no beam births enter the confined modal domain"""
from __future__ import annotations
from typing import Any
import numpy as np
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState, collision_parameter_metadata
from source_model_revamp.fbis.species import DEUTERON, TRITON
from source_model_revamp.fbis.modal.basis import ModalBasisConvergenceAssessment
from source_model_revamp.fbis.modal.lost_ion_distribution import EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY
from source_model_revamp.fbis.modal.source_projection import ModalSourceProjectionAudit, modal_source_projection_audit_metadata
from source_model_revamp.fbis.modal import ModalFBISBasis, ModalFBISNumerics
from source_model_revamp.fbis.modal.species_state import FastIonSpeciesState, FastIonSystemState
from source_model_revamp.integration.pipeline_types import BeamEnsembleResult, GeometryStageResult, KineticStageResult
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.modal_stage.types import _ReusableModalBasisBundle, electron_density_input_is_initial_guess_only
from source_model_revamp.integration.modal_stage.basis_density import _collision_state_for_density

def _basis_and_assessment(modal_basis: _ReusableModalBasisBundle | ModalFBISBasis) -> tuple[ModalFBISBasis, ModalBasisConvergenceAssessment | None]:
    """Return the physical modal basis and optional basis convergence assessment from a reusable bundle or bare basis"""
    if isinstance(modal_basis, _ReusableModalBasisBundle):
        return modal_basis.basis, modal_basis.convergence_assessment
   
    return modal_basis, None

def _prompt_only_collision_state(config: SourceModelRunConfig, *, electron_temperature_J: float, ion_temperature_J: float) -> FBISCollisionParameterState:
    """
    Build a collision reference for a prompt only source without creating a confined fast ion population
    
    The benchmark closure uses its electron density reference while other closures use the configured D plus T reference densities
    """
    closure = str(config.kinetic_electrostatic.modal_density_closure_model).strip().lower()
    if closure == "egedal_beam_plasma_quasineutral":
        electron_density = float(config.plasma_closure.electron_density_initial_guess_m3)
        ion_densities = np.asarray([electron_density], dtype=float)
        ion_charges = np.asarray([1.0], dtype=float)
        ion_masses = np.asarray([DEUTERON.mass_kg], dtype=float)
    else:
        ion_densities = np.asarray([config.plasma_closure.background_deuterium_midplane_density_m3, config.plasma_closure.background_tritium_midplane_density_m3,],dtype=float,)
        ion_charges = np.asarray([1.0, 1.0], dtype=float)
        ion_masses = np.asarray([DEUTERON.mass_kg, TRITON.mass_kg], dtype=float)
        electron_density = float(np.sum(ion_densities))
        if float(np.sum(ion_densities)) <= 0.0:
            raise ValueError("prompt only collision reference is unavailable when configured seed and reference D and T densities are zero")
    
    return _collision_state_for_density(fast_ion_species=DEUTERON, electron_density_m3=electron_density, ion_densities_m3=ion_densities, ion_charge_numbers=ion_charges, ion_masses_kg=ion_masses, electron_temperature_J=electron_temperature_J, ion_temperature_J=ion_temperature_J,)

def _source_audit_power_metadata(beam: BeamEnsembleResult, audit: ModalSourceProjectionAudit) -> dict[str, Any]:
    """Convert prompt and confined source audit rates into component and total deposited birth powers"""
    component_energies = [ float(component.spatial_source.energy_J) for component in beam.attenuated_source.component_sources]
    prompt_powers = [ item.prompt_magnetic_loss_birth_rate_s * energy for item, energy in zip(audit.component_audits, component_energies, strict=True)]
    confined_powers = [ item.confined_physical_birth_rate_s * energy for item, energy in zip(audit.component_audits, component_energies, strict=True)]
    
    return {
        "modal_component_prompt_loss_birth_power_W": prompt_powers,
        "modal_component_confined_birth_power_W": confined_powers,
        "modal_prompt_loss_birth_power_W": float(sum(prompt_powers)),
        "modal_confined_physical_birth_power_W": float(sum(confined_powers)),
    }

def _prompt_only_kinetic_stage( config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, *, electron_temperature_J: float, ion_temperature_J: float, numerics: ModalFBISNumerics, modal_basis: _ReusableModalBasisBundle | ModalFBISBasis, source_audit: ModalSourceProjectionAudit, collision_state_recompute_policy: str | None, collision_state_is_self_consistent_electron_temperature: bool) -> KineticStageResult:
    """
    Return the exact zero confined kinetic state when all beam births are prompt magnetic losses
    
    The result preserves source partition, geometry, basis, collision reference, and applicability metadata while leaving Eq 68, Eq 70, and confined loss reconstruction unavailable
    """
    basis, assessment = _basis_and_assessment(modal_basis)
    fixed_boundary_selected = (str(config.kinetic_electrostatic.ion_loss_closure_model).strip().lower() == EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY)
    collision_state = _prompt_only_collision_state(config, electron_temperature_J=electron_temperature_J, ion_temperature_J=ion_temperature_J)
    collision_metadata = collision_parameter_metadata(collision_state)
    collision_metadata.update({ "collision_state_scope": "background_reference_no_confined_fast_ion_collision_solve", "collision_state_recompute_policy": str( collision_state_recompute_policy or "computed_once_for_prompt_only_background_reference"), "collision_state_is_self_consistent_electron_temperature": bool(collision_state_is_self_consistent_electron_temperature),})
    if collision_state_is_self_consistent_electron_temperature:
        collision_metadata["collision_state_temperature_closure"] = "self_consistent_electron_temperature_trial_background_reference"
        collision_metadata["plasma_temperature_closure_status"] = "self_consistent_electron_temperature_trial_background_reference"
    n_speed = beam.speed_grid.centers_m_s.size
    n_lambda = beam.lambda_grid.centers.size
    n_pitch = beam.pitch_grid.centers.size
    n_z = geometry.z_centers_m.size
    zero_v_lambda = np.zeros((n_speed, n_lambda), dtype=float)
    zero_z_v_lambda = np.zeros((n_z, n_speed, n_lambda), dtype=float)
    zero_z_v_pitch = np.zeros((n_z, n_speed, n_pitch), dtype=float)
    zero_density = np.zeros(n_z, dtype=float)
    source_metadata = modal_source_projection_audit_metadata(source_audit)
    source_power_metadata = _source_audit_power_metadata(beam, source_audit)
    prompt_rate = float(source_audit.prompt_magnetic_loss_birth_rate_s)
    prompt_power = float(source_power_metadata["modal_prompt_loss_birth_power_W"])
    total_rate = float(source_audit.total_deposited_birth_rate_s)
    partition_scale = max(abs(total_rate), 1.0)
    partition_conserved = bool( abs(float(source_audit.total_partition_residual_s)) <= 1.0e-12 * partition_scale)
    assessment_present = assessment is not None
    basis_converged = bool(assessment is not None and assessment.converged)
    published_feedback_requested = str(config.kinetic_electrostatic.electrostatic_feedback_model).strip().lower() not in {"magnetic_only_fbis"}
    mirror_throat_inside_domain = (geometry.metadata.get("fitted_mirror_throat_inside_kinetic_domain") is True)
    beam_pitch_selection_valid = (beam.metadata.get("beam_pitch_selection_valid") is True)
    boundary_selection_valid = bool(geometry.plasma_boundaries is not None and geometry.plasma_boundaries.valid and geometry.B_T_function is not None and geometry.full_device_z_centers_m is not None)
    geometry_boundary_valid = bool(mirror_throat_inside_domain and beam_pitch_selection_valid and boundary_selection_valid)
    metadata: dict[str, Any] = {'kinetic_model': 'egedal_modal_fbis', 'backend_numerics_settings': config.kinetic_electrostatic.backend_numerics_metadata(), 'ion_loss_closure_model': config.kinetic_electrostatic.ion_loss_closure_model, 'active_ion_loss_boundary_model': 'fixed_magnetic_boundary_Lambda_M_equals_1_over_RM' if fixed_boundary_selected else None, 'active_ion_loss_uses_Eq71_moving_boundary': False, 'Eq71_72_loss_boundary_role': 'diagnostic_local_accessibility_and_quasineutral_reconstruction_only' if fixed_boundary_selected else None, 'lost_ion_parallel_temperature_model': config.kinetic_electrostatic.lost_ion_parallel_temperature_model, 'lost_ion_parallel_temperature_left_keV': None, 'lost_ion_parallel_temperature_right_keV': None, **collision_metadata, 'modal_density_closure_model': config.kinetic_electrostatic.modal_density_closure_model, 'modal_density_closure_applicable': False, 'modal_density_closure_converged': False, 'modal_density_closure_iterations': 0, 'modal_density_closure_relative_error': None, 'modal_density_closure_history': [], 'modal_density_closure_status': 'not_applicable_no_confined_fbis_source', 'electron_density_input_is_initial_guess_only': electron_density_input_is_initial_guess_only(config.kinetic_electrostatic.modal_density_closure_model), 'electron_density_initial_guess_m3': config.plasma_closure.electron_density_initial_guess_m3, 'configured_seed_reference_midpoint_positive_charge_m3': float(config.plasma_closure.background_deuterium_midplane_density_m3 + config.plasma_closure.background_tritium_midplane_density_m3), 'electron_density_initial_guess_valid': config.plasma_closure.electron_density_initial_guess_valid, 'electron_density_initial_guess_source': config.plasma_closure.electron_density_initial_guess_source, 'electron_midplane_density_m3': None, 'electron_confined_volume_average_density_m3': None, 'electron_parent_maxwellian_n0_m3': None, 'electron_collision_density_m3': float(collision_state.coulomb_logs.electron_density_m3), 'electron_density_profile_m3': None, 'electron_profile_scope': 'unavailable_no_confined_fbis_source', 'expander_electron_profile_available': False, 'electron_temperature_initial_guess_keV': config.plasma_closure.electron_temperature_initial_guess_keV, 'background_ion_temperature_keV': config.plasma_closure.background_ion_temperature_keV, 'background_thermal_ion_end_loss_available': False, 'global_plasma_particle_power_balance_available': False, 'electron_temperature_power_balance_scope': 'fast_D_electron_FBIS_subsystem', 'global_plasma_power_balance_available': False, 'kinetic_convergence_status': 'not_applicable_no_confined_fbis_source', 'kinetic_convergence_failure_reason': 'no_confined_fbis_source_kinetic_and_electrostatic_blocks_not_applicable_or_unavailable', 'kinetic_iterations': 0, 'kinetic_all_inner_steps_converged': False, 'kinetic_geometry_boundary_check_passed': geometry_boundary_valid, 'kinetic_geometry_contains_fitted_mirror_throat': mirror_throat_inside_domain, 'kinetic_boundary_selection_valid': boundary_selection_valid, 'kinetic_beam_pitch_selection_valid': beam_pitch_selection_valid, 'kinetic_beam_injection_geometry_consistent': beam.metadata.get('beam_injection_geometry_consistent') is True, 'kinetic_symmetric_end_model_applicable': True, 'kinetic_left_right_mirror_ratio_relative_difference': 0.0, 'kinetic_basis_convergence_assessed': assessment_present, 'kinetic_basis_converged': basis_converged, 'modal_basis_convergence_assessed': assessment_present, 'modal_basis_converged': basis_converged, 'modal_basis_convergence_failure_reason': None if basis_converged else None if assessment is None else assessment.failure_reason, 'modal_basis_eigenvalues': np.asarray(basis.physical_basis.eigenvalues, dtype=float), 'modal_basis_gram_condition_number': float(basis.physical_basis.gram_condition_number), 'kinetic_eq59_convergence_applicable': False, 'kinetic_eq59_convergence_assessed': False, 'kinetic_eq59_converged': True, 'modal_rosenbluth_convergence_applicable': False, 'modal_rosenbluth_convergence_assessed': False, 'modal_rosenbluth_converged': True, 'modal_rosenbluth_convergence_status': 'not_applicable_no_confined_fbis_source', 'modal_rosenbluth_convergence_failure_reason': None, 'modal_rosenbluth_final_relative_residual': None, 'modal_rosenbluth_high_speed_tail_population_fraction': 0.0, 'modal_rosenbluth_high_speed_tail_population_tolerance': numerics.rosenbluth_high_speed_tail_population_tolerance, 'modal_source_projection_convergence_applicable': False, 'modal_source_projection_nonzero_for_nonzero_birth': True, 'modal_source_projection_rate_converged': True, 'modal_source_projection_relative_error': 0.0, 'modal_source_projection_relative_tolerance': config.kinetic_electrostatic.modal_source_projection_relative_tolerance, 'modal_source_projection_status': 'not_applicable_zero_confined_physical_birth_rate', 'kinetic_eq70_profile_present': False, 'kinetic_eq70_profile_applicable': False, 'modal_phi_z_converged': False, 'modal_phi_z_profile_V': None, 'modal_phi_z_failure_reason': 'unavailable_no_confined_fbis_source_electrostatic_closure_not_solved', 'kinetic_published_heuristic_requested': bool(published_feedback_requested and (not fixed_boundary_selected)), 'kinetic_eq71_72_loss_boundary_diagnostic_only': bool(fixed_boundary_selected), 'kinetic_published_heuristic_check_passed': bool(fixed_boundary_selected or not published_feedback_requested), 'kinetic_expander_unclassified_population_check_passed': False, 'kinetic_full_device_lost_population_required': prompt_rate > 0.0, 'kinetic_full_device_lost_population_check_passed': False, 'kinetic_full_device_lost_population_reference_compatible': False, 'kinetic_current_balance_relative_tolerance': float(getattr(config.kinetic_electrostatic, 'modal_current_balance_relative_tolerance', numerics.loss_convention_relative_tolerance)), 'kinetic_total_device_loss_convention_valid': False, 'modal_loss_convention_validation_passed': False, 'modal_loss_convention_failure_reason': 'electrostatic_and_electron_end_loss_closure_unavailable_for_prompt_only_source', 'modal_prompt_birth_partition_conservation_check_passed': partition_conserved, 'modal_prompt_birth_partition_relative_error': float(source_audit.total_partition_residual_s / partition_scale), 'modal_total_beam_birth_rate_s': total_rate, 'modal_total_beam_deposited_birth_power_W': float(beam.total_deposited_birth_power_W), 'modal_source_birth_rate_s': 0.0, 'modal_source_deposited_birth_power_W': 0.0, 'modal_prompt_loss_birth_rate_s': prompt_rate, 'modal_prompt_loss_birth_power_W': prompt_power, 'modal_prompt_loss_birth_fraction': 1.0 if total_rate > 0.0 else 0.0, 'modal_confined_birth_fraction': 0.0, 'modal_projected_source_particle_rate_from_coefficients_s': 0.0, 'modal_projected_source_component_rates_s': [0.0 for _ in source_audit.component_audits], 'modal_ion_particle_loss_rate_s': prompt_rate, 'modal_ion_midplane_kinetic_power_loss_W': prompt_power, 'modal_ion_wall_power_loss_W': None, 'modal_electron_wall_power_loss_W': None, 'modal_end_loss_power_terms_available': False, 'modal_end_loss_power_unavailability_reason': 'electrostatic_barrier_and_electron_current_balance_not_solved_for_prompt_only_source', 'electrostatic_model': 'unavailable_no_confined_fbis_source', 'electron_loss_model': 'unavailable_no_confined_fbis_source', 'closed_interval_population_fraction': None, 'closed_interval_population_fraction_status': 'not_applicable_no_confined_fbis_source', 'local_well_population_fraction': None, 'unclassified_population_fraction': None, 'kinetic_min_distribution_value': 0.0, 'kinetic_raw_final_distribution_min': 0.0, 'kinetic_negative_cell_count': 0, 'kinetic_roundoff_negative_cells_zeroed_for_fusion': 0, 'solved_fast_D_inventory_v_lambda_sum': 0.0, 'solved_fast_D_volume_average_density_m3_from_inventory_sum': 0.0, 'solved_fast_D_inventory_particles': 0.0, 'solved_fast_D_ion_particle_loss_rate_s': prompt_rate, 'solved_fast_D_ion_confinement_time_s': None, 'solved_electrostatic_potential_min_V': None, 'solved_electrostatic_potential_max_V': None, 'kinetic_confined_arrays_are_exact_zero': True, 'kinetic_prompt_only_source_is_hard_invalid': False, **source_metadata, **source_power_metadata}
    inactive_species_state = FastIonSpeciesState(species=DEUTERON, speed_grid=beam.speed_grid_by_species[DEUTERON.species_id], active=False, source_particle_rate_s=0.0, modal_result=None, status="inactive_zero_confined_source")
    fast_ion_system_state = FastIonSystemState(species_states={DEUTERON.species_id: inactive_species_state}, shared_basis=basis)

    return KineticStageResult(speed_grid=beam.speed_grid, lambda_grid=beam.lambda_grid, pitch_grid=beam.pitch_grid, final_distribution_v_lambda=zero_v_lambda, metadata=metadata, local_distribution_z_v_lambda=zero_z_v_lambda, local_distribution_z_v_pitch=zero_z_v_pitch, local_density_m3=zero_density, full_device_local_distribution_z_v_pitch=None, local_speed_grid=beam.speed_grid, eq59_warm_start_state=None, fast_ion_system_state=fast_ion_system_state)
