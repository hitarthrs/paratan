"""
Final source model metadata assembly and plotting contract

The stage filters stage metadata into the canonical run record, adds missing physics availability state, and declares which grids own each population and plotting quantity
"""
from __future__ import annotations
from collections.abc import Mapping
from typing import Any
import numpy as np
from source_model_revamp.integration.metadata import REQUIRED_RUN_METADATA_KEYS, require_metadata_keys
from source_model_revamp.integration.modal_stage.metadata_adapter import KINETIC_STATUS_METADATA_KEYS
from source_model_revamp.integration.pipeline_types import SourceModelPipelineResult

_HISTORY_MARKERS = ("_history", "_histories", "trial_records", "iteration_records")
_CANONICAL_EXACT_KEYS = set(REQUIRED_RUN_METADATA_KEYS) | set(KINETIC_STATUS_METADATA_KEYS) | {
    "kinetic_convergence_failure_reasons",
    "kinetic_convergence_additional_failure_count",
    "total_current_numerical_fixed_point_converged",
    "equivalent_end_model_applicable",
    "represented_populations",
    "excluded_populations",
    "limitations",
    "beam_collection_model",
    "beam_count",
    "enabled_beam_count",
    "beam_order",
    "beam_ids",
    "enabled_beam_ids",
    "beam_species_by_id",
    "beam_power_W_by_id",
    "beam_energy_keV_by_id",
    "beam_component_ids_by_id",
    "beam_fast_birth_rate_s_by_id",
    "beam_deposited_birth_power_W_by_id",
    "beam_shine_through_power_W_by_id",
    "beam_deposition_by_id",
    "beam_source_aggregation",
    "beam_source_aggregation_order",
    "beam_active_projectile_species",
    "beam_effective_birth_energy_keV",
    "beam_representative_injection_angle_deg",
    "beam_atomic_data_model",
    "beam_atomic_data_reference",
    "beam_atomic_data_collisiondb_qid",
    "beam_atomic_data_collisiondb_doi",
    "beam_atomic_data_collisiondb_method",
    "beam_atomic_data_collisiondb_energy_frame",
    "beam_atomic_data_collisiondb_energy_units",
    "beam_atomic_data_collisiondb_cross_section_units",
    "beam_atomic_data_collisiondb_data_sha256",
    "beam_atomic_data_collisiondb_source_pdf_sha256",
    "beam_atomic_data_collisiondb_min_energy_keV_per_u",
    "beam_atomic_data_collisiondb_max_energy_keV_per_u",
    "beam_atomic_data_ionization_low_energy_truncation_max_relative_rate",
    "beam_atomic_data_ionization_high_energy_truncation_max_relative_rate",
    "beam_stopping_density_model",
    "beam_stopping_species_resolved",
    "beam_stopping_target_species",
    "beam_stopping_model_limitation",
    "beam_electron_temperature_keV",
    "beam_ion_temperature_keV",
    "beam_electron_ionization_birth_rate_s",
    "beam_deuterium_ionization_birth_rate_s",
    "beam_tritium_ionization_birth_rate_s",
    "beam_ionization_birth_rate_s",
    "beam_deuterium_charge_exchange_fast_birth_rate_s",
    "beam_tritium_charge_exchange_fast_birth_rate_s",
    "beam_charge_exchange_fast_birth_rate_s",
    "beam_charge_exchange_target_ion_pumpout_rate_s",
    "beam_deuterium_target_pumpout_rate_s",
    "beam_tritium_target_pumpout_rate_s",
    "beam_deuterium_target_pumpout_energy_W",
    "beam_tritium_target_pumpout_energy_W",
    "beam_charge_exchange_target_pumpout_energy_W",
    "beam_net_plasma_fueling_rate_s",
    "beam_particle_bookkeeping_scope",
    "beam_charge_exchange_target_pumpout_kinetic_model",
    "beam_primary_projectile_species",
    "beam_total_fast_birth_rate_s_by_species",
    "beam_total_deposited_birth_power_W_by_species",
    "beam_source_construction_status",
    "species_density_relative_change_by_species",
    "species_density_change_metrics_by_species",
    "nbi_particle_authority_chain_by_species",
}
_PLOTTING_METADATA_KEYS = {
    "B0_T",
    "B_tilde_centers",
    "axial_bin_rates_s",
    "active_fast_species",
    "background_positive_charge_profile_m3",
    "background_positive_charge_density_profile_m3",
    "beam_attenuation_beam_end_m",
    "beam_attenuation_beam_radius_m",
    "beam_attenuation_beam_start_m",
    "beam_deposited_birth_power_W",
    "beam_effective_injection_angle_deg",
    "beam_energy_keV",
    "beam_injection_angle_deg",
    "beam_particle_mass_kg",
    "beam_path_beam_radius_profile_m",
    "beam_path_center_points_m",
    "beam_axial_birth_rate_density_m3_s_by_species",
    "beam_axial_birth_rate_s_by_species",
    "beam_axial_birth_power_density_W_m3_by_species",
    "beam_axial_birth_power_W_by_species",
    "beam_canonical_order",
    "beam_species_birth_rate_profile_identity_relative_error",
    "beam_species_birth_power_profile_identity_relative_error",
    "beam_species_profile_identity_relative_tolerance",
    "beam_species_birth_rate_profile_identity_check_passed",
    "beam_species_birth_power_profile_identity_check_passed",
    "beam_power_W",
    "beam_radius_m",
    "beam_source_axial_physical_birth_rate_density_m3_s",
    "beam_source_axial_physical_birth_rate_s",
    "beam_source_axial_physical_birth_power_density_W_m3",
    "beam_source_axial_physical_birth_power_W",
    "beam_birth_profile_support_mask",
    "bottleneck_vacuum_radius_m",
    "candidate_plasma_facing_surfaces",
    "cell_volumes_m3",
    "central_cell_vacuum_radius_m",
    "coil_group",
    "coil_radius_m",
    "coil_z_center_m",
    "collision_state_electron_temperature_keV",
    "coordinate_centers_m",
    "coordinate_edges_m",
    "device_domains",
    "electron_temperature_solved_keV",
    "electron_density_profile_m3",
    "end_cell_vacuum_radius_m",
    "energy_bin_rates_s",
    "energy_edges_J",
    "energy_edges_MeV",
    "energy_edges_eV",
    "fast_ion_species_mass_kg",
    "fast_ion_species_mass_kg_by_species",
    "fast_ion_species_charge_number_by_species",
    "fast_ion_speed_grid_m_s_by_species",
    "fast_ion_speed_grid_faces_m_s_by_species",
    "fast_ion_distribution_v_lambda_by_species",
    "fast_ion_local_speed_grid_m_s_by_species",
    "fast_ion_local_speed_grid_faces_m_s_by_species",
    "fast_ion_local_distribution_z_v_lambda_by_species",
    "fast_ion_local_distribution_z_v_pitch_by_species",
    "fast_ion_local_density_m3_by_species",
    "fast_ion_full_device_distribution_z_v_pitch_by_species",
    "fast_ion_full_device_density_m3_by_species",
    "final_fast_D_distribution_v_lambda",
    "fitted_mirror_ratio",
    "fitted_mirror_throat_z_m",
    "full_device_z_centers_m",
    "full_device_z_edges_m",
    "full_device_cell_volumes_m3",
    "deuterium_background_density_profile_m3",
    "tritium_background_density_profile_m3",
    "background_profile_shape",
    "background_profile_model",
    "background_population_scope",
    "fusion_channel_axial_neutron_rate_s",
    "fusion_total_axial_neutron_rate_s",
    "fusion_total_neutron_rate_density_m3_s",
    "fusion_total_fusion_power_density_W_m3",
    "half_length_m",
    "lambda_centers",
    "lambda_faces",
    "left_fitted_mirror_throat_z_m",
    "left_wall_power_W",
    "magnetic_field_derived_Bmax_T",
    "magnetic_field_derived_throat_z_m",
    "magnetic_field_effective_Bmax_T",
    "magnetic_field_visual_B_T",
    "magnetic_field_visual_B_tilde",
    "magnetic_field_visual_coordinate_m",
    "magnetic_field_visual_flux_tube_radius_m",
    "mirror_ratio",
    "modal_component_source_coefficients",
    "modal_confined_physical_birth_power_W",
    "modal_effective_ion_temperature_keV",
    "modal_effective_ion_temperature_over_beam_energy",
    "modal_eigenvalues",
    "modal_electron_wall_power_loss_W",
    "modal_eta_lambda_grid",
    "modal_eta_of_lambda_grid",
    "modal_eta_trapped_passing_boundary",
    "modal_ion_confinement_time_s",
    "modal_ion_wall_power_loss_W",
    "modal_lambda_boundary",
    "modal_local_density_profile_m3",
    "modal_local_distribution_z_v_pitch",
    "modal_local_speed_grid_faces_m_s",
    "modal_phi_z_electron_density_profile_m3",
    "modal_phi_z_ion_density_profile_m3",
    "modal_phi_z_profile_V",
    "modal_phi_z_converged",
    "modal_physical_eigenfunctions_eta",
    "modal_physical_eta_grid",
    "modal_source_mean_birth_energy_keV",
    "modal_source_mean_birth_speed_m_s",
    "modal_total_wall_power_loss_W",
    "modal_wall_barrier_energy_J",
    "modal_wall_barrier_over_Te",
    "neutron_fusion_profile_total_rate_s",
    "paratan_layout",
    "pitch_centers",
    "pitch_edges",
    "pitch_faces",
    "plasma_radius_m",
    "radius_at_z_max_m",
    "radius_at_z_min_m",
    "right_fitted_mirror_throat_z_m",
    "right_wall_power_W",
    "solved_fast_D_ion_confinement_time_s",
    "source_geometry_linkage_source_z_max_m",
    "source_geometry_linkage_source_z_min_m",
    "source_rate_matrix_s",
    "speed_centers_m_s",
    "speed_faces_m_s",
    "total_neutron_rate_s",
    "total_positive_charge_profile_m3",
    "z_edges_cm",
    "z_edges_m",
}
_CANONICAL_STAGE_KEYS = {
    "fast_ion_species_id",
    "fast_ion_species_name",
    "fast_ion_species_symbol",
    "fast_ion_species_mass_kg",
    "fast_ion_species_charge_number",
    "ion_loss_closure_model",
    "modal_eq61_left_one_end_particle_loss_rate_s",
    "modal_eq61_right_one_end_particle_loss_rate_s",
    "modal_two_end_eq61_loss_rate_s",
    "modal_eq63_left_throat_crossing_rate_s",
    "modal_eq63_right_throat_crossing_rate_s",
    "modal_eq63_left_rate_relative_error",
    "modal_eq63_right_rate_relative_error",
    "modal_eq64_matching_model",
    "modal_particle_balance_relative_error_loss_minus_source_over_source",
    "modal_power_balance_relative_error",
    "modal_current_balance_relative_error",
    "modal_electron_current_balance_relative_error",
    "modal_confined_ion_collisional_confinement_time_s",
    "solved_fast_D_inventory_particles",
    "kinetic_convergence_failure_reason",
    "modal_density_closure_model",
    "modal_density_closure_applicable",
    "modal_density_closure_converged",
    "modal_density_closure_relative_error",
    "modal_density_closure_relative_tolerance",
    "modal_density_iteration_relative_tolerance",
    "kinetic_local_speed_refinement_applicable",
    "species_eq59_converged",
    "species_local_reconstruction_convergence",
    "all_species_local_reconstruction_converged",
    "fast_D_fast_T_cross_collisions_available",
    "shared_eq70_production_state",
    "operating_point_converged",
    "operating_point_failure_reason",
    "electron_temperature_solve_failure_reason",
    "electron_temperature_valid_bracket_keV",
    "total_charged_product_power_W",
    "openmc_file_source_path",
    "openmc_file_source_validity_path",
    "full_device_expander_population_model",
    "full_device_lost_fast_ion_population_model",
    "full_device_fast_ion_distribution_includes_source_connected_expander",
    "full_device_fast_ion_distribution_includes_finite_return",
    "full_device_fast_ion_distribution_excludes_source_disconnected_local_wells",
    "full_device_linear_expander_potential_is_authoritative",
    "fusion_source_disconnected_expander_population_included",
    "backend_numerics_settings",
    "kinetic_basis_convergence_assessed",
    "kinetic_basis_converged",
    "modal_n_lambda_grid",
    "modal_n_eta_grid",
}
_CANONICAL_PREFIXES = ('expander_', 'terminal_current_', 'total_current_balance_', 'electron_temperature_', 'electron_energy_balance_', 'represented_electron_energy_', 'fast_ion_electron_transfer_', 'physical_electron_energy_', 'global_plasma_power_balance_', 'global_device_power_balance_', 'electron_tail_integral_', 'electron_mean_wall_kinetic_energy_', 'eq68_', 'eq69_', 'fast_d_', 'fast_t_', 'beam_density_', 'beam_current_outer_iteration_', 'fusion_alpha_heating_', 'fusion_expander_', 'radiative_electron_losses_', 'beam_ionization_electron_energy_loss_', 'radial_electron_transport_', 'fusion_channel_', 'fusion_component_', 'fusion_population_', 'fusion_reaction_', 'fusion_dd_', 'fusion_dt_', 'fusion_fast_t_', 'fusion_tt_', 'modal_basis_', 'modal_rosenbluth_', 'neutron_component_', 'neutron_dt_', 'runtime_')
_PLOTTING_CONTRACT_VERSION = 1
_PLOTTING_SPECIES_ORDER = ("deuterium", "tritium")

def _nonempty_mapping(value: object) -> bool:
    """Return whether a value is a nonempty mapping"""
    return isinstance(value, Mapping) and bool(value)

def _plotting_metadata_contract(metadata: Mapping[str, Any]) -> dict[str, Any]:
    """
    Build the stable plotting and population ownership contract from assembled stage metadata
    
    The contract identifies active fast species, required full device distribution and density keys, canonical beam ordering, grid ownership, and population roles
    """
    canonical_beams = metadata.get("beam_canonical_order")
   
    if not isinstance(canonical_beams, (list, tuple)):
        beam_records = metadata.get("beam_deposition_by_id")
        canonical_beams = sorted(str(key) for key in beam_records) if isinstance(beam_records, Mapping) else []
   
    active_fast_species = metadata.get("active_fast_species")
  
    if not isinstance(active_fast_species, (list, tuple)):
        active_fast_species = []
   
    active_fast_set = {str(value) for value in active_fast_species}
    full_device_distribution = metadata.get("fast_ion_full_device_distribution_z_v_pitch_by_species")
    full_device_density = metadata.get("fast_ion_full_device_density_m3_by_species")
    full_device_fast_species_available = bool(active_fast_set and isinstance(full_device_distribution, Mapping) and isinstance(full_device_density, Mapping) and set(str(key) for key in full_device_distribution) == active_fast_set and set(str(key) for key in full_device_density) == active_fast_set)
    nbi_supported = str(metadata.get("plasma_closure_model", "")).strip().lower() == "nbi_supported_stationary"
    grid_ownership = {"beam_axial_profiles": "confined_kinetic", "fast_confined_profiles": "confined_kinetic", "fast_full_device_profiles": "full_device_population", "source_connected_expander_profiles": "full_device_population", "electrostatic_profiles": "confined_kinetic", "fusion_component_profiles": "full_device_population"}
  
    if nbi_supported:
        grid_ownership['startup_seed_profiles'] = 'full_device_population'
        population_roles = {'startup_seed_deuterium': {'species_id': 'deuterium', 'population_kind': 'startup_seed', 'role': 'initialization_only', 'authoritative_for': []}, 'startup_seed_tritium': {'species_id': 'tritium', 'population_kind': 'startup_seed', 'role': 'initialization_only', 'authoritative_for': []}, 'fast_deuterium': {'species_id': 'deuterium', 'population_kind': 'fast', 'role': 'species_resolved_fbis_population', 'authoritative_for': ['beam_stopping_after_first_iteration', 'Eq42_confined_collision_weighting', 'Eq70', 'fusion']}, 'fast_tritium': {'species_id': 'tritium', 'population_kind': 'fast', 'role': 'species_resolved_fbis_population', 'authoritative_for': ['beam_stopping_after_first_iteration', 'Eq42_confined_collision_weighting', 'Eq70', 'fusion']}, 'source_connected_expander_fast_deuterium': {'species_id': 'deuterium', 'population_kind': 'fast', 'role': 'source_connected_U_mu_residence_mapping'}, 'source_connected_expander_fast_tritium': {'species_id': 'tritium', 'population_kind': 'fast', 'role': 'source_connected_U_mu_residence_mapping'}}
    else:
        population_roles = {'fast_deuterium': {'species_id': 'deuterium', 'population_kind': 'fast', 'role': 'species_resolved_fbis_population'}, 'fast_tritium': {'species_id': 'tritium', 'population_kind': 'fast', 'role': 'species_resolved_fbis_population'}, 'source_connected_expander_fast_deuterium': {'species_id': 'deuterium', 'population_kind': 'fast', 'role': 'source_connected_U_mu_residence_mapping'}, 'source_connected_expander_fast_tritium': {'species_id': 'tritium', 'population_kind': 'fast', 'role': 'source_connected_U_mu_residence_mapping'}}
  
    return {
        "version": _PLOTTING_CONTRACT_VERSION,
        "canonical_species_order": list(_PLOTTING_SPECIES_ORDER),
        "canonical_beam_order": [str(value) for value in canonical_beams],
        "active_fast_species": [str(value) for value in active_fast_species],
        "grid_definitions": {
            "confined_kinetic": {"edges_key": "coordinate_edges_m", "centers_key": "coordinate_centers_m", "cell_volumes_key": "cell_volumes_m3", "scope": "confined_throat_to_throat"},
            "full_device_population": {"edges_key": "full_device_z_edges_m", "centers_key": "full_device_z_centers_m", "cell_volumes_key": "full_device_cell_volumes_m3", "scope": "full_device"},
            "magnetic_visual": {"centers_key": "magnetic_field_visual_coordinate_m", "scope": "full_device_visual"},
        },
        "grid_ownership": grid_ownership,
        "population_roles": population_roles,
        "source_keys": {
            "beam_records": "beam_deposition_by_id",
            "beam_profiles_by_species": "beam_axial_birth_rate_density_m3_s_by_species",
            "beam_profile_total": "beam_source_axial_physical_birth_rate_density_m3_s",
            "startup_seed_positive_charge": "background_positive_charge_density_profile_m3",
            "fast_density_by_species": "fast_ion_local_density_m3_by_species",
            "fast_full_device_density_by_species": "fast_ion_full_device_density_m3_by_species",
            "fusion_component_profiles": "fusion_channel_axial_neutron_rate_s",
            "fusion_profile_total": "fusion_total_axial_neutron_rate_s",
            "current_components": "total_current_balance_components",
        },
        "available_blocks": {'multi_beam': _nonempty_mapping(metadata.get('beam_deposition_by_id')), 'beam_profiles_by_species': _nonempty_mapping(metadata.get('beam_axial_birth_rate_density_m3_s_by_species')), 'fast_species': _nonempty_mapping(metadata.get('fast_ion_local_density_m3_by_species')), 'fast_full_device_species': full_device_fast_species_available, 'source_connected_expander': bool(metadata.get('expander_population_included_in_fusion_input')), 'fusion_registry': _nonempty_mapping(metadata.get('fusion_component_population_pairs')), 'total_current_balance': _nonempty_mapping(metadata.get('total_current_balance_components'))},
    }

def _include_metadata_key(key: str, *, write_iteration_histories: bool) -> bool:
    """Return whether one metadata key belongs in the canonical run record for the selected iteration history setting"""
    lowered = key.lower()
    if any(marker in lowered for marker in _HISTORY_MARKERS):
        return write_iteration_histories
 
    return bool(key in _CANONICAL_EXACT_KEYS or key in _PLOTTING_METADATA_KEYS or key in _CANONICAL_STAGE_KEYS or lowered.startswith(_CANONICAL_PREFIXES))

def _missing_physics_metadata(metadata: Mapping[str, Any]) -> dict[str, object]:
    """Return missing physics availability metadata after considering whether the expander closure is present and closed"""
    remaining = [str(value) for value in metadata.get("electron_energy_balance_missing_physics_terms", ())]
    background_available = metadata.get("represented_background_particle_maintenance_available") is True and metadata.get("represented_background_energy_maintenance_available") is True
    expander_available = metadata.get("expander_population_included_in_fusion_input") is True and metadata.get("terminal_current_balance_authoritative") is True and metadata.get("electron_energy_balance_uses_terminal_ion_power") is True
    expander_closed = expander_available and metadata.get("expander_population_classification_passed") is True and metadata.get("expander_nonempirical_electron_trapping_closure_unavailable") is not True
  
    if expander_closed and "expander_population_maintenance" in remaining:
        remaining.remove("expander_population_maintenance")
  
    elif expander_available and "expander_population_maintenance" in remaining:
        index = remaining.index("expander_population_maintenance")
        remaining[index] = "source_disconnected_or_unclosed_expander_population_maintenance"
   
    return {"electron_energy_balance_missing_physics_terms": tuple(dict.fromkeys(remaining))}

def build_pipeline_metadata(result: SourceModelPipelineResult) -> dict[str, Any]:
    """
    Collect filtered metadata from geometry, beam, operating point, kinetic, expander, fusion, neutron, and OpenMC stages
    
    Required run keys are checked after adding missing physics state and the plotting metadata contract
    """
    metadata: dict[str, Any] = {"source_model_backend": "source_model_revamp", "source_model_workflow_model": result.config.model}
    write_histories = bool(result.config.output.write_iteration_histories)
  
    for stage in (result.geometry, result.beam, result.operating_point, result.kinetic, result.expander, result.fusion, result.neutrons, result.openmc_export):
        if stage is None:
            continue
        for key, value in getattr(stage, "metadata").items():
            if _include_metadata_key(str(key), write_iteration_histories=write_histories):
                metadata[str(key)] = value
  
    metadata.update(_missing_physics_metadata(metadata))
    metadata["plotting_metadata_contract_version"] = _PLOTTING_CONTRACT_VERSION
    metadata["plotting_metadata_contract"] = _plotting_metadata_contract(metadata)
    metadata["iteration_histories_written"] = write_histories
    require_metadata_keys(metadata, REQUIRED_RUN_METADATA_KEYS)
   
    return metadata
