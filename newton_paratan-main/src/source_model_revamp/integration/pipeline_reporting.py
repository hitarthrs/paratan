"""Represented population, exclusion, and limitation reporting for the final run record"""
from __future__ import annotations
from source_model_revamp.integration.pipeline_types import KineticStageResult
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.plasma_closure import uses_startup_seed_only

def represented_population_scope(config: SourceModelRunConfig, kinetic: KineticStageResult) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    """
    Return explicit represented populations, excluded populations, and model limitations for the solved closure state
    
    The report distinguishes confined, fixed boundary lost, and source connected expander fast ion populations and records additional NBI supported stationary closure limitations when that closure is active
    """
    included: list[str] = []
    startup_seed_only = uses_startup_seed_only(config.plasma_closure)
    system = kinetic.fast_ion_system_state
    source_species = () if system is None else system.source_species_ids
    active_species = () if system is None else system.active_species_ids
    prompt_only = str(kinetic.metadata.get("kinetic_source_state", "")).strip().lower() == "no_confined_fbis_source"
    if prompt_only:
        for species_id in source_species:
            included.append(f"exact_zero_confined_fast_{species_id}_population_no_confined_fbis_source")
    else:
        for species_id in active_species:
            included.append(f"confined_fast_{species_id}_on_geometry_confined_support")
    expander_population_included = bool(kinetic.metadata.get("expander_population_included_in_fusion_input", False))
    directed_population_included = bool(kinetic.metadata.get("full_device_fast_ion_distribution_includes_directed_lost", False))
    if expander_population_included:
        for species_id in active_species:
            included.append(f"source_connected_expander_fast_{species_id}_U_mu_residence_mapping")
    elif directed_population_included:
        for species_id in active_species:
            included.append(f"egedal_hot_fixed_boundary_directed_lost_fast_{species_id}")
    excluded = [
        "cold_magnetic_eq63_reference_as_active_population",
        "source_disconnected_local_well_population_without_collision_closure",
        "unclassified_expander_population",
        "confined_fast_distribution_extrapolated_into_expanders",
        "material_reflected_fast_population_without_Eq42_Eq59_return_coupling",
        "prompt_expander_population_without_side_resolved_throat_distribution",
    ]
    limitations = [
        "eq71_moving_boundary_loss_population_not_derived_or_validated",
        "results_include_only_finite_explicitly_represented_populations",
        "implemented_cold_zero_potential_eq63_diagnostic_is_never_an_active_downstream_population",
    ]
    if startup_seed_only:
        limitations.append("startup_D_T_profiles_initialize_the_nonzero_branch_but_are_excluded_from_final_density_and_fusion")
        limitations.append("heavy_particle_beam_targets_use_the_previous_converged_kinetic_D_T_state_after_the_first_iteration")
        limitations.append("primary_beam_charge_exchange_removes_target_ions_with_a_reaction_weighted_speed_resolved_pitch_averaged_sink")
        limitations.append("charge_exchange_target_sink_does_not_retain_axial_or_pitch_localization_inside_the_reduced_modal_solver")
        limitations.append("fast_D_T_cross_collisions_use_the_other_species_isotropic_first_pitch_mode_and_the_existing_representative_ion_Coulomb_log")
        limitations.append("stationary_species_particle_balance_requires_NBI_ionization_and_cross_isotope_CX_transfer_to_balance_terminal_loss")
        limitations.append("fusion_uses_only_converged_beam_born_kinetic_D_T_populations_in_the_NBI_supported_closure")
        limitations.append("fast_D_and_fast_T_are_origin_labels_and_include_slowed_and_warm_ions")
        limitations.append("fusion_burnup_is_not_an_active_kinetic_sink_and_must_pass_the_configured_negligibility_gate")
    if expander_population_included:
        limitations.append("source_connected_expander_ion_populations_use_collisionless_U_mu_residence_mapping")
        limitations.append("Egedal_fixed_fast_boundary_requires_material_fast_expander_return_to_be_negligible")
        limitations.append("source_disconnected_expander_wells_are_detected_but_not_populated")
        limitations.append("complete_expander_electron_kinetics_beyond_Eq68_current_and_Eq70_density_closures_are_not_modeled")
        limitations.append("collisional_expander_transport_is_not_modeled")
        limitations.append("expander_terminal_material_roles_are_restricted_to_absorbing_collector_and_end_ring")
    else:
        limitations.append("source_connected_expander_population_unavailable_for_this_run")
    if not directed_population_included:
        for species_id in active_species:
            excluded.append(f"active_hot_electrostatic_lost_fast_{species_id}")
        limitations.append("active_hot_fixed_boundary_lost_population_unavailable_for_this_run")
   
    return tuple(included), tuple(excluded), tuple(limitations)


__all__ = ["represented_population_scope"]
