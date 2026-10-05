"""
Shared ion current rebuild and postclosure consistency helpers

These helpers rebuild the coupled kinetic state when explicit ion current components are supplied and assess whether that closure changed the beam attenuation density state
"""
from __future__ import annotations
from dataclasses import replace
from collections.abc import Mapping
import numpy as np
from source_model_revamp.constants import J_TO_KEV
from source_model_revamp.fbis.species import DEUTERON, TRITON
from source_model_revamp.integration.evaluation import OperatingPointEvaluationMode
from source_model_revamp.integration.modal_stage.system_stage import build_multispecies_modal_kinetic_stage
from source_model_revamp.integration.pipeline_types import BeamEnsembleResult, GeometryStageResult, KineticStageResult
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.plasma_closure import uses_startup_seed_only

def _temperature_metadata_for_rebuild(kinetic: KineticStageResult, electron_temperature_J: float) -> dict[str, object]:
    """Return electron temperature metadata retained when the kinetic state is rebuilt for current closure"""
    metadata = {key: value for key, value in kinetic.metadata.items() if "electron_temperature" in key or key in {"collision_state_electron_temperature_keV", "collision_state_ion_temperature_keV", "collision_state_temperature_closure"}}
    temperature_keV = float(electron_temperature_J) * J_TO_KEV
    metadata["collision_state_electron_temperature_keV"] = temperature_keV
    metadata["electron_temperature_solved_keV"] = temperature_keV
    metadata["operating_point_electron_temperature_J"] = float(electron_temperature_J)
    return metadata

def _initial_density_inputs(kinetic: KineticStageResult) -> dict[str, object]:
    """Extract solved fast D, fast T, and electron scalar density values used to warm start a current rebuild"""
    density = kinetic.operating_point_density_state
    if density is None:
        return {}
   
    return {
        "initial_fast_midpoint_density_m3_by_species": {DEUTERON.species_id: float(density.fast_deuterium_midplane_density_m3), TRITON.species_id: float(density.fast_tritium_midplane_density_m3)},
        "initial_fast_average_density_m3_by_species": {DEUTERON.species_id: float(density.fast_deuterium_confined_volume_average_density_m3), TRITON.species_id: float(density.fast_tritium_confined_volume_average_density_m3)},
        "initial_electron_midpoint_density_m3": float(density.electron_midplane_density_m3),
        "initial_electron_confined_average_density_m3": float(density.electron_confined_volume_average_density_m3),
        "initial_electron_collision_density_m3": float(density.electron_collision_density_m3),
    }

def _primary_species_id(kinetic: KineticStageResult) -> str | None:
    """Return the preferred active fast species identifier for flattened compatibility metadata"""
    system = kinetic.fast_ion_system_state
    if system is None or not system.active_species_ids:
        return None
 
    return DEUTERON.species_id if DEUTERON.species_id in system.active_species_ids else system.active_species_ids[0]

def _flatten_primary_species_metadata(kinetic: KineticStageResult) -> dict[str, object]:
    """Return selected primary species kinetic metadata fields at the shared kinetic metadata level"""
    primary_id = _primary_species_id(kinetic)
    if primary_id is None:
        return {}
  
    metadata = kinetic.metadata
    flattened: dict[str, object] = {}
   
    for key in ("species_collision_metadata", "species_source_projection_audits", "species_full_device_population_metadata"):
        values = metadata.get(key)
        if isinstance(values, Mapping):
            item = values.get(primary_id)
            if isinstance(item, Mapping):
                flattened.update(dict(item))
  
    system = kinetic.fast_ion_system_state
    state = system.species_states.get(primary_id)
    modal = None if state is None else state.modal_result
   
    if modal is not None:
        modal_metadata = modal.metadata
        flattened.update({
            "kinetic_basis_convergence_assessed": modal_metadata.get("modal_basis_convergence_assessed") is True,
            "kinetic_basis_converged": modal_metadata.get("modal_basis_converged") is True,
            "kinetic_retained_mode_count_applicable": int(modal_metadata.get("modal_n_physical_modes", 0)) > 1,
            "kinetic_eq59_convergence_applicable": True,
            "kinetic_eq59_convergence_assessed": modal_metadata.get("modal_rosenbluth_convergence_assessed") is True,
            "kinetic_eq59_converged": modal_metadata.get("modal_rosenbluth_converged") is True,
            "kinetic_local_speed_refinement_applicable": modal_metadata.get("modal_local_speed_refinement_assessed") is not None,
        })
  
    return flattened

def _merge_rebuilt_metadata(previous: KineticStageResult, rebuilt: KineticStageResult, electron_temperature_J: float) -> dict[str, object]:
    """Merge prior kinetic metadata, rebuilt system metadata, primary species compatibility fields, and retained temperature metadata"""
    metadata = dict(previous.metadata)
    metadata.update(rebuilt.metadata)
    metadata.update(_flatten_primary_species_metadata(rebuilt))
    metadata.update(_temperature_metadata_for_rebuild(previous, electron_temperature_J))
    metadata["total_current_balance_rebuild_preserved_prior_diagnostics"] = True
    metadata["total_current_balance_rebuild_current_state_overrides_prior_state"] = True
  
    return metadata

def _postclosure_beam_density_consistency(config: SourceModelRunConfig, initial: KineticStageResult, final: KineticStageResult) -> dict[str, object]:
    """
    Compare the electron density state before and after a current closure rebuild
    
    The checks reuse the beam density pointwise, volume L2, and absolute reference tolerances so a current closure cannot silently invalidate the attenuation target
    """
    initial_density = initial.operating_point_density_state
    final_density = final.operating_point_density_state
   
    if initial_density is None or final_density is None:
        return {
            "total_current_balance_postclosure_beam_density_consistency_assessed": False,
            "total_current_balance_postclosure_beam_density_consistency_passed": False,
            "total_current_balance_postclosure_beam_density_consistency_failure_reason": "operating point density state unavailable",
        }
   
    initial_profile = np.asarray(initial_density.electron_cell_density_m3, dtype=float)
    final_profile = np.asarray(final_density.electron_cell_density_m3, dtype=float)
    volumes = np.asarray(final_density.cell_volumes_m3, dtype=float)
  
    if initial_profile.shape != final_profile.shape or initial_profile.shape != volumes.shape:
        return {
            "total_current_balance_postclosure_beam_density_consistency_assessed": False,
            "total_current_balance_postclosure_beam_density_consistency_passed": False,
            "total_current_balance_postclosure_beam_density_consistency_failure_reason": "electron density grids changed",
        }
  
    controls = config.beam_density_coupling
    reference = max(float(np.max(np.abs(initial_profile))), float(np.max(np.abs(final_profile))), 1.0)
    floor = max(float(controls.density_floor_absolute_m3), float(controls.density_floor_reference_fraction) * reference)
    scale = np.maximum(np.maximum(np.abs(initial_profile), np.abs(final_profile)), floor)
    difference = np.abs(final_profile - initial_profile)
    pointwise = float(np.max(difference / scale))
    numerator = float(np.sum(volumes * difference**2))
    denominator = float(np.sum(volumes * scale**2))
    volume_L2 = float(np.sqrt(numerator / max(denominator, np.finfo(float).tiny)))
    absolute_reference = float(np.max(difference) / reference)
    passed = bool(pointwise <= float(controls.electron_profile_relative_tolerance) and volume_L2 <= float(controls.profile_volume_L2_tolerance) and absolute_reference <= float(controls.profile_absolute_reference_tolerance))
  
    return {
        "total_current_balance_postclosure_beam_density_consistency_assessed": True,
        "total_current_balance_postclosure_beam_density_consistency_passed": passed,
        "total_current_balance_postclosure_beam_density_consistency_failure_reason": None if passed else "Pass 11 electron density change exceeds the beam density fixed point tolerances",
        "total_current_balance_postclosure_electron_density_pointwise_relative_change": pointwise,
        "total_current_balance_postclosure_electron_density_volume_L2_relative_change": volume_L2,
        "total_current_balance_postclosure_electron_density_absolute_reference_change": absolute_reference,
        "total_current_balance_postclosure_electron_density_pointwise_relative_tolerance": float(controls.electron_profile_relative_tolerance),
        "total_current_balance_postclosure_electron_density_volume_L2_relative_tolerance": float(controls.profile_volume_L2_tolerance),
        "total_current_balance_postclosure_electron_density_absolute_reference_tolerance": float(controls.profile_absolute_reference_tolerance),
        "total_current_balance_beam_source_held_fixed": True,
    }

def postclosure_beam_density_consistency_metadata(config: SourceModelRunConfig, initial: KineticStageResult, final: KineticStageResult) -> dict[str, object]:
    """Return beam density consistency metadata for a current closure rebuild"""
    return _postclosure_beam_density_consistency(config, initial, final)

def rebuild_kinetic_for_ion_current_components(*, config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, kinetic: KineticStageResult, electron_temperature_J: float, ion_loss_rates_s_by_component: Mapping[str, float], evaluation_mode: OperatingPointEvaluationMode | str, collision_state_recompute_policy: str) -> KineticStageResult:
    """
    Rebuild the coupled modal kinetic state with explicit species ion current component rates
    
    The current kinetic solution supplies density, Eq 59, Eq 70, Eq 42, collision, and barrier warm states when compatible
    """
    mode = OperatingPointEvaluationMode.parse(evaluation_mode)
    rates = {str(key): float(value) for key, value in dict(ion_loss_rates_s_by_component).items()}
    rebuilt = build_multispecies_modal_kinetic_stage(config, geometry, beam, electron_temperature_J=float(electron_temperature_J), collision_state_recompute_policy=str(collision_state_recompute_policy), collision_state_is_self_consistent_electron_temperature=bool(kinetic.metadata.get('collision_state_is_self_consistent_electron_temperature', False)), modal_basis=kinetic.modal_basis_warm_start, eq59_warm_start_state_by_species=None if kinetic.eq59_warm_start_state_by_species is None else dict(kinetic.eq59_warm_start_state_by_species), eq70_initial_potential_energy_J=kinetic.eq70_potential_energy_warm_start_J, initial_eq42_density_profile=kinetic.eq42_density_profile_warm_start, scalar_collision_operator_warm_start_by_species=None if kinetic.scalar_collision_operator_warm_start_by_species is None else dict(kinetic.scalar_collision_operator_warm_start_by_species), scalar_external_fast_field_warm_start_by_test_species=None if kinetic.scalar_external_fast_field_warm_start_by_test_species is None else {key: tuple(values) for key, values in kinetic.scalar_external_fast_field_warm_start_by_test_species.items()}, scalar_wall_barrier_energy_warm_start_J=kinetic.scalar_wall_barrier_energy_warm_start_J, ion_loss_rates_s_by_component_override=rates, evaluation_mode=mode, **_initial_density_inputs(kinetic))
    metadata = _merge_rebuilt_metadata(kinetic, rebuilt, electron_temperature_J)
    metadata.update({"explicit_ion_current_component_rebuild": True, "explicit_ion_current_component_rates_s": rates, "explicit_ion_current_component_rebuild_policy": str(collision_state_recompute_policy)})
   
    return replace(rebuilt, metadata=metadata)

def mark_pre_expander_current_candidate(config: SourceModelRunConfig, kinetic: KineticStageResult) -> KineticStageResult:
    """Mark throat level Eq 68 current balance as a pre expander candidate rather than a terminal current qualification"""
    if not uses_startup_seed_only(config.plasma_closure):
        return kinetic
    metadata = {**kinetic.metadata, 'total_current_balance_model': 'shared_Eq68_fast_and_prompt_kinetic_ion_current_pre_expander_candidate', 'total_current_balance_scope': 'kinetic_D_kinetic_T_prompt_D_prompt_T_pre_expander', 'total_current_balance_applicable': True, 'total_current_balance_candidate_active': True, 'total_current_balance_numerical_fixed_point_converged': False, 'total_current_balance_end_model_applicable': False, 'total_current_balance_converged': False, 'total_current_numerical_fixed_point_converged': False, 'equivalent_end_model_applicable': False, 'total_current_balance_failure_reason': 'terminal_expander_current_not_yet_closed', 'total_current_balance_failure_reasons': ('terminal_expander_current_not_yet_closed',), 'startup_seed_ion_current_included': False, 'global_current_balance_claimed': False, 'physical_electron_energy_equation_active': True}
  
    return replace(kinetic, metadata=metadata)
