"""
Source connected expander population and terminal current closure

The stage couples Eq 63 throat crossing rates to expander trajectories and iterates the terminal ion current with the shared Eq 68 barrier and Eq 70 potential state
"""
from __future__ import annotations
from time import perf_counter
from collections.abc import Mapping
from dataclasses import replace
import numpy as np
from source_model_revamp.constants import ELECTRON_CHARGE_C
from source_model_revamp.expander.solver import solve_expander_system
from source_model_revamp.expander.types import ExpanderSystemState
from source_model_revamp.fbis.species import DEUTERON, TRITON
from source_model_revamp.fbis.velocity_space_grid import gyrotropic_velocity_cell_volumes
from source_model_revamp.fusion.populations import FAST_D, FAST_T
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.plasma_closure import uses_startup_seed_only
from source_model_revamp.integration.evaluation import OperatingPointEvaluationMode
from source_model_revamp.integration.full_device_populations import conservative_confined_distribution_on_grid
from source_model_revamp.integration.pipeline_stages.ion_current_helpers import postclosure_beam_density_consistency_metadata, rebuild_kinetic_for_ion_current_components
from source_model_revamp.integration.pipeline_types import BeamEnsembleResult, ExpanderStageResult, GeometryStageResult, KineticStageResult
from source_model_revamp.runtime_profile import finalize_runtime_profile, new_runtime_profile, record_runtime

def _relative_scalar_change(value_a: float, value_b: float) -> float:
    """Return the symmetric relative change between two finite scalar values"""
    a = float(value_a)
    b = float(value_b)

    return abs(a - b) / max(abs(a), abs(b), np.finfo(float).tiny)

def _relative_array_change(value_a: np.ndarray, value_b: np.ndarray) -> float:
    """Return the largest normalized absolute change between two finite arrays"""
    a = np.asarray(value_a, dtype=float)
    b = np.asarray(value_b, dtype=float)
  
    if a.shape != b.shape:
        return float("inf")
    numerator = float(np.max(np.abs(a - b))) if a.size else 0.0
    denominator = max(float(np.max(np.abs(a))) if a.size else 0.0, float(np.max(np.abs(b))) if b.size else 0.0, np.finfo(float).tiny)

    return numerator / denominator

def _prompt_component_rates(kinetic: KineticStageResult) -> dict[str, float]:
    """Return prompt magnetic loss rates by fast ion species in s⁻¹"""
    system = kinetic.fast_ion_system_state
  
    if system is None:
        return {}
    rates: dict[str, float] = {}
  
    for species_id, state in system.species_states.items():
        if state.modal_result is not None:
            rate = float(state.modal_result.metadata.get("modal_prompt_loss_birth_rate_s", 0.0))
        else:
            rate = float(state.prompt_only_particle_loss_rate_s)
        if rate > 0.0:
            rates[f"prompt_{species_id}"] = rate
  
    return rates

def _throat_component_rates(kinetic: KineticStageResult) -> dict[str, float]:
    """Return nonprompt fixed magnetic boundary throat crossing rates by fast ion species in s⁻¹"""
    rates: dict[str, float] = {}
    system = kinetic.fast_ion_system_state
   
    if system is not None:
        for species_id, state in system.species_states.items():
            if state.modal_result is not None:
                prompt = float(state.modal_result.metadata.get("modal_prompt_loss_birth_rate_s", 0.0))
                rates[f"fast_{species_id}"] = max(0.0, float(state.modal_result.ion_particle_loss_rate_s) - prompt)
                if prompt > 0.0:
                    rates[f"prompt_{species_id}"] = prompt
            elif float(state.prompt_only_particle_loss_rate_s) > 0.0:
                rates[f"prompt_{species_id}"] = float(state.prompt_only_particle_loss_rate_s)
  
    return rates

def _terminal_current_target_rates(expander: ExpanderSystemState, kinetic: KineticStageResult) -> dict[str, float]:
    """Return species terminal current target rates including source connected expander losses and unresolved prompt loss"""
    rates = {str(key): float(value) for key, value in expander.terminal_particle_rate_s_by_component.items()}
    rates.update(_prompt_component_rates(kinetic))
  
    return {key: max(0.0, value) for key, value in rates.items()}

def _relax_rates(current: Mapping[str, float], target: Mapping[str, float], relaxation: float) -> dict[str, float]:
    """Logarithmically relax positive species rates toward terminal target rates"""
    keys = sorted(set(current) | set(target))
   
    return {key: max(0.0, (1.0 - relaxation) * float(current.get(key, 0.0)) + relaxation * float(target.get(key, 0.0))) for key in keys}

def _component_changes(current: Mapping[str, float], target: Mapping[str, float]) -> dict[str, float]:
    """Return species resolved relative changes between current and target rate mappings"""
    keys = sorted(set(current) | set(target))
  
    return {key: _relative_scalar_change(float(current.get(key, 0.0)), float(target.get(key, 0.0))) for key in keys}

def _barrier_energy_J(kinetic: KineticStageResult) -> float:
    """Return the shared Eq 68 electron barrier energy in J"""
    system = kinetic.fast_ion_system_state
    if system is None or system.shared_current_balance is None:
        return float("nan")
  
    return float(system.shared_current_balance.electron_balance.barrier_energy_J)

def _eq70_profile_J(kinetic: KineticStageResult) -> np.ndarray:
    """Return the shared Eq 70 electron potential energy profile on electrostatic nodes in J"""
    system = kinetic.fast_ion_system_state
    if system is None or system.shared_electrostatic_profile is None:
        return np.asarray([], dtype=float)
  
    return np.asarray(system.shared_electrostatic_profile.potential_energy_J, dtype=float)

def _expander_potential_change(current: ExpanderSystemState, previous: ExpanderSystemState) -> float:
    """Return the relative change between two expander potential energy profiles"""
    changes = []
    for side in ("left", "right"):
        current_state = current.potential_by_side[side]
        previous_state = previous.potential_by_side[side]
        changes.append(_relative_array_change(np.asarray(current_state.potential_drop_centers_outward_J, dtype=float), np.asarray(previous_state.potential_drop_centers_outward_J, dtype=float)))
  
    return max(changes, default=0.0)

def _active_total_current_values(*, terminal_state_unavailable: bool, pass11_metadata: Mapping[str, object], pass11_target_current: object, pass11_electron_current: object, terminal_target_current: float, Eq68_electron_current: float, terminal_component_current: float, terminal_component_identity_error: float, terminal_component_identity_passed: bool) -> dict[str, object]:
    """
    Select the total current diagnostics that are authoritative for the active closure state
    
    Terminal values are used when the expander state is available and throat level Eq 68 values remain a fallback when terminal closure is unavailable
    """
    if not terminal_state_unavailable:
        target = float(terminal_target_current)
        electron = float(Eq68_electron_current)
        component = float(terminal_component_current)
        component_error = float(terminal_component_identity_error)
        component_passed = bool(terminal_component_identity_passed)
    else:
        target = float(pass11_target_current if pass11_target_current is not None else terminal_target_current)
        electron = float(pass11_electron_current if pass11_electron_current is not None else Eq68_electron_current)
        component = float(pass11_metadata.get("total_current_balance_component_current_sum_A", target))
        component_error = float(pass11_metadata.get("total_current_balance_component_current_identity_relative_error", _relative_scalar_change(component, target)))
        component_passed = pass11_metadata.get("total_current_balance_component_current_identity_check_passed") is True
    
    residual = electron - target
    relative_residual = abs(residual) / max(abs(electron), abs(target), np.finfo(float).tiny)
   
    return {
        "target_ion_current_A": target,
        "electron_current_A": electron,
        "current_residual_A": residual,
        "current_relative_residual": relative_residual,
        "component_current_sum_A": component,
        "component_current_identity_relative_error": component_error,
        "component_current_identity_check_passed": component_passed,
        "Pass11_throat_diagnostic": bool(terminal_state_unavailable),
    }

def _terminal_end_metadata(expander: ExpanderSystemState, *, asymmetry_tolerance: float, unresolved_tolerance: float) -> dict[str, object]:
    """Return left and right terminal current, asymmetry, unresolved prompt fraction, and associated qualification metadata"""
    left_current = 0.0
    right_current = 0.0
  
    for branch in expander.branches_by_population_side.values():
        current = ELECTRON_CHARGE_C * float(branch.terminal_particle_rate_s)
        if branch.side == "left":
            left_current += current
        else:
            right_current += current
   
    total = left_current + right_current
    asymmetry = abs(left_current - right_current) / max(abs(left_current), abs(right_current), np.finfo(float).tiny)
    unresolved_rate = float(expander.prompt_unresolved_particle_rate_s)
    unresolved_fraction = unresolved_rate / max(float(expander.total_throat_particle_rate_s), np.finfo(float).tiny)
    prompt_mapping_available = bool(unresolved_rate <= np.finfo(float).tiny)
    applicable = prompt_mapping_available
    passed = bool(applicable and asymmetry <= float(asymmetry_tolerance))
  
    return {
        "terminal_current_left_ion_current_A": left_current,
        "terminal_current_right_ion_current_A": right_current,
        "terminal_current_side_sum_A": total,
        "terminal_current_end_asymmetry_relative_error": asymmetry,
        "terminal_current_end_asymmetry_relative_tolerance": float(asymmetry_tolerance),
        "terminal_current_end_asymmetry_applicable": applicable,
        "terminal_current_end_asymmetry_check_passed": passed,
        "terminal_current_prompt_distribution_available": prompt_mapping_available,
        "terminal_current_prompt_unresolved_particle_rate_s": unresolved_rate,
        "terminal_current_prompt_unresolved_fraction": unresolved_fraction,
        "terminal_current_prompt_unresolved_fraction_tolerance": 0.0,
        "terminal_current_prompt_unresolved_materiality_tolerance": float(unresolved_tolerance),
    }

def _mark_terminal_current_state(kinetic: KineticStageResult, expander: ExpanderSystemState, metadata: Mapping[str, object]) -> KineticStageResult:
    """Attach terminal current closure metadata to the kinetic result, shared fast ion system, and modal species results"""
    combined = dict(metadata)
    system = kinetic.fast_ion_system_state
  
    if system is not None:
        states = {}
        for species_id, state in system.species_states.items():
            modal = state.modal_result
            if modal is None:
                states[species_id] = state
            else:
                states[species_id] = replace(state, modal_result=replace(modal, metadata={**modal.metadata, **combined}))
        system = replace(system, species_states=states, metadata={**system.metadata, **combined})
    kinetic = replace(kinetic, fast_ion_system_state=system, metadata={**kinetic.metadata, **combined})
 
    return kinetic

def _full_fast_distributions(geometry: GeometryStageResult, kinetic: KineticStageResult, expander: ExpanderSystemState) -> dict[str, np.ndarray]:
    """Return full device fast D and fast T local distributions formed from confined, central Eq 63, and directed expander populations"""
    if geometry.full_device_z_edges_m is None or geometry.full_device_cell_volumes_m3 is None or kinetic.local_distribution_z_v_pitch_by_species is None:
        raise ValueError("full device grid and local fast distributions are required")
   
    full_edges = np.asarray(geometry.full_device_z_edges_m, dtype=float)
    full_volumes = np.asarray(geometry.full_device_cell_volumes_m3, dtype=float)
    result: dict[str, np.ndarray] = {}
    system = kinetic.fast_ion_system_state
  
    if system is None:
        raise ValueError("full device fast distribution requires the solved fast ion system state")
    full_centers = np.asarray(geometry.full_device_z_centers_m, dtype=float)
    central_mask = (full_centers >= float(geometry.z_edges_m[0])) & (full_centers <= float(geometry.z_edges_m[-1]))
  
    for species, population_id in ((DEUTERON, FAST_D), (TRITON, FAST_T)):
        local = kinetic.local_distribution_z_v_pitch_by_species.get(species.species_id)
        if local is None:
            continue
        confined = conservative_confined_distribution_on_grid(source_edges_m=np.asarray(geometry.z_edges_m, dtype=float), source_cell_volumes_m3=np.asarray(geometry.cell_volumes_m3, dtype=float), source_distribution_z_v_pitch=np.asarray(local, dtype=float), target_edges_m=full_edges, target_cell_volumes_m3=full_volumes)
        central_lost = np.zeros_like(confined)
        population_metadata = dict(kinetic.metadata.get("species_full_device_population_metadata", {})).get(species.species_id, {})
        represented_lost = population_metadata.get("full_device_lost_ion_distribution") if isinstance(population_metadata, Mapping) else None
      
        if represented_lost is not None:
            lost_array = np.asarray(represented_lost, dtype=float)
            if lost_array.shape != confined.shape:
                raise ValueError(f"{population_id} central Eq 63 distribution does not match the full device fast grid")
            central_lost[central_mask] = lost_array[central_mask]
        directed = expander.full_device_distribution_by_population.get(population_id)
      
        if directed is None:
            result[species.species_id] = confined + central_lost
            continue
        directed_array = np.asarray(directed, dtype=float)
      
        if directed_array.shape != confined.shape:
            raise ValueError(f"{population_id} expander distribution does not match the full device fast grid")
        result[species.species_id] = confined + central_lost + directed_array
   
    return result

def _attach_full_device_populations(geometry: GeometryStageResult, kinetic: KineticStageResult, expander: ExpanderSystemState) -> KineticStageResult:
    """
    Attach full device fast ion distributions and density profiles to the kinetic result
    
    The metadata distinguishes central Eq 63 throat populations from source connected expander populations and records finite return and material connection state
    """
    full_by_species = _full_fast_distributions(geometry, kinetic, expander)
  
    system = kinetic.fast_ion_system_state
    if system is None:
        raise ValueError("expander integration requires a fast ion system state")
   
    species_states = dict(system.species_states)
    population_by_species = {DEUTERON.species_id: FAST_D, TRITON.species_id: FAST_T}
    species_population_metadata = dict(kinetic.metadata.get("species_full_device_population_metadata", {}))
    full_device_density_by_species: dict[str, np.ndarray] = {}
    directed_by_species: dict[str, bool] = {}
    finite_return_by_species: dict[str, bool] = {}
    for species_id, state in tuple(species_states.items()):
        distribution = full_by_species.get(species_id)
        population_id = population_by_species.get(species_id)
        directed = None if population_id is None else expander.full_device_distribution_by_population.get(population_id)
        directed_included = bool(directed is not None and np.any(np.asarray(directed, dtype=float) > 0.0))
        finite_return_included = bool(population_id is not None and any(branch.population_id == population_id and branch.reflected_particle_rate_s > 0.0 for branch in expander.branches_by_population_side.values()))
        directed_by_species[species_id] = directed_included
        finite_return_by_species[species_id] = finite_return_included
      
        if distribution is not None and state.modal_result is not None:
            local_speed_grid = state.modal_result.local_reconstruction.local_speed_grid
            if local_speed_grid is None:
                local_speed_grid = state.modal_result.speed_grid
            velocity_measure = gyrotropic_velocity_cell_volumes(local_speed_grid, kinetic.pitch_grid)
            full_device_density_by_species[species_id] = np.sum(np.maximum(np.asarray(distribution, dtype=float), 0.0) * velocity_measure[None, :, :], axis=(1, 2))
      
        existing_population_metadata = dict(species_population_metadata.get(species_id, {}))
        central_eq63_included = existing_population_metadata.get("full_device_lost_ion_distribution") is not None
        species_population_metadata[species_id] = {
            **existing_population_metadata,
            **expander.metadata,
            "full_device_lost_ion_reconstruction_status": "Egedal_Eq63_central_and_source_connected_expander_with_Eq70_center_nodal_potential",
            "full_device_lost_ion_reconstruction_available": distribution is not None,
            "full_device_lost_ion_reconstruction_failure_reason": None if distribution is not None else "Pass13_full_device_fast_distribution_unavailable",
            "full_device_fast_ion_distribution_includes_confined": distribution is not None,
            "full_device_fast_ion_distribution_includes_central_Eq63_lost": central_eq63_included,
            "full_device_fast_ion_distribution_includes_directed_lost": directed_included or central_eq63_included,
            "full_device_fast_ion_distribution_includes_source_connected_expander": directed_included,
            "full_device_fast_ion_distribution_includes_finite_return": finite_return_included,
            "full_device_fast_ion_distribution_excludes_source_disconnected_local_wells": True,
            "full_device_lost_ion_potential_model": "Egedal_Eq70_center_nodal_expander",
            "full_device_linear_expander_potential_is_authoritative": False,
        }
        modal = state.modal_result

        if modal is not None:
            modal_metadata = {**modal.metadata, **species_population_metadata[species_id]}
            modal = replace(modal, metadata=modal_metadata)

        species_states[species_id] = replace(state, modal_result=modal, full_device_local_distribution_z_v_pitch=distribution)
    active_species = tuple(species_id for species_id, state in species_states.items() if state.active)
    full_device_includes_confined = bool(active_species and all(full_by_species.get(species_id) is not None for species_id in active_species))
    central_eq63_by_species = {species_id: bool(species_population_metadata.get(species_id, {}).get("full_device_fast_ion_distribution_includes_central_Eq63_lost", False)) for species_id in active_species}
    full_device_includes_central_eq63 = bool(active_species and all(central_eq63_by_species.values()))
    combined_directed_by_species = {species_id: bool(directed_by_species.get(species_id, False) or central_eq63_by_species.get(species_id, False)) for species_id in active_species}
    full_device_includes_directed = bool(active_species and all(combined_directed_by_species.values()))
    full_device_includes_any_directed = bool(any(combined_directed_by_species.values()))
    full_device_includes_source_connected_expander = bool(any(directed_by_species.get(species_id, False) for species_id in active_species))
    full_device_includes_finite_return = bool(any(finite_return_by_species.get(species_id, False) for species_id in active_species))
    system_metadata = {
        **system.metadata,
        **expander.metadata,
        "full_device_expander_population_model": "source_connected_U_mu_conserving_residence_mapping",
        "species_full_device_population_metadata": species_population_metadata,
        "full_device_fast_ion_distribution_includes_confined": full_device_includes_confined,
        "full_device_fast_ion_distribution_includes_central_Eq63_lost": full_device_includes_central_eq63,
        "full_device_fast_ion_distribution_includes_central_Eq63_lost_by_species": central_eq63_by_species,
        "full_device_fast_ion_distribution_includes_directed_lost": full_device_includes_directed,
        "full_device_fast_ion_distribution_includes_any_directed_lost": full_device_includes_any_directed,
        "full_device_fast_ion_distribution_includes_directed_lost_by_species": combined_directed_by_species,
        "full_device_fast_ion_distribution_includes_source_connected_expander": full_device_includes_source_connected_expander,
        "full_device_fast_ion_distribution_includes_source_connected_expander_by_species": directed_by_species,
        "full_device_fast_ion_distribution_includes_finite_return": full_device_includes_finite_return,
        "full_device_fast_ion_distribution_includes_finite_return_by_species": finite_return_by_species,
        "full_device_fast_ion_distribution_excludes_source_disconnected_local_wells": True,
        "fast_ion_full_device_distribution_z_v_pitch_by_species": full_by_species,
        "fast_ion_full_device_density_m3_by_species": full_device_density_by_species,
    }
    system = replace(system, species_states=species_states, metadata=system_metadata)
    metadata = {
        **kinetic.metadata,
        **expander.metadata,
        "full_device_expander_population_model": "source_connected_U_mu_conserving_residence_mapping",
        "full_device_lost_fast_ion_population_model": "Egedal_Eq63_source_connected_U_mu_mapping",
        "species_full_device_population_metadata": species_population_metadata,
        "full_device_fast_ion_distribution_includes_confined": full_device_includes_confined,
        "full_device_fast_ion_distribution_includes_central_Eq63_lost": full_device_includes_central_eq63,
        "full_device_fast_ion_distribution_includes_central_Eq63_lost_by_species": central_eq63_by_species,
        "full_device_fast_ion_distribution_includes_directed_lost": full_device_includes_directed,
        "full_device_fast_ion_distribution_includes_any_directed_lost": full_device_includes_any_directed,
        "full_device_fast_ion_distribution_includes_directed_lost_by_species": combined_directed_by_species,
        "full_device_fast_ion_distribution_includes_source_connected_expander": full_device_includes_source_connected_expander,
        "full_device_fast_ion_distribution_includes_source_connected_expander_by_species": directed_by_species,
        "full_device_fast_ion_distribution_includes_finite_return": full_device_includes_finite_return,
        "full_device_fast_ion_distribution_includes_finite_return_by_species": finite_return_by_species,
        "full_device_fast_ion_distribution_excludes_source_disconnected_local_wells": True,
        "full_device_linear_expander_potential_is_authoritative": False,
        "expander_population_included_in_fusion_input": True,
        "fast_ion_full_device_distribution_z_v_pitch_by_species": full_by_species,
        "fast_ion_full_device_density_m3_by_species": full_device_density_by_species,
        "global_plasma_power_balance_available": False,
        "global_device_power_balance_claimed": False,
    }
    primary_species = DEUTERON.species_id if DEUTERON.species_id in full_by_species else TRITON.species_id if TRITON.species_id in full_by_species else None
   
    return replace(kinetic, metadata=metadata, fast_ion_system_state=system, full_device_local_distribution_z_v_pitch=None if primary_species is None else full_by_species[primary_species], full_device_local_distribution_z_v_pitch_by_species=full_by_species or None)

def build_expander_stage(config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, kinetic: KineticStageResult, *, electron_temperature_J: float, evaluation_mode: OperatingPointEvaluationMode | str) -> tuple[KineticStageResult, ExpanderStageResult]:
    """
    Iterate the source connected expander and terminal current closure
    
    Each iteration rebuilds the shared kinetic current state with species terminal target rates, solves the expander populations, relaxes the next current targets, and checks component rates, Eq 68 barrier, Eq 70 profile, expander potential, and current residuals before a final consistency reconstruction
    """
    expander_runtime_started = perf_counter()
    runtime_profile = new_runtime_profile()
    runtime_kinetic_evaluations: list[dict[str, object]] = []
    mode = OperatingPointEvaluationMode.parse(evaluation_mode)
    startup_seed_only = uses_startup_seed_only(config.plasma_closure)
    controls = config.kinetic_electrostatic
    tolerance = float(controls.total_current_balance_relative_tolerance)
    relaxation = float(controls.total_current_balance_relaxation)
    asymmetry_tolerance = float(controls.total_current_balance_end_asymmetry_relative_tolerance)
    unresolved_tolerance = float(controls.expander_fast_nonterminal_fraction_tolerance)
    pass11_scope = str(kinetic.metadata.get("current_balance_scope") or "fast_D_fast_T_prompt_D_prompt_T_throat_crossing_current")
    pass11_numerical = kinetic.metadata.get("total_current_numerical_fixed_point_converged", kinetic.metadata.get("total_current_balance_numerical_fixed_point_converged")) is True
    pass11_qualified = kinetic.metadata.get("total_current_balance_converged") is True
    pass11_target_current = kinetic.metadata.get("total_current_balance_target_ion_current_A")
    pass11_electron_current = kinetic.metadata.get("total_current_balance_electron_current_A")
    initial_expander_started = perf_counter()
    initial_expander = solve_expander_system(config, geometry, kinetic)
   
    record_runtime(runtime_profile, "initial_expander_solve", perf_counter() - initial_expander_started)
  
    current_kinetic = kinetic
    current_expander = initial_expander
    current_target = _throat_component_rates(current_kinetic)
    previous_barrier = _barrier_energy_J(current_kinetic)
    previous_eq70 = _eq70_profile_J(current_kinetic)
    history: list[dict[str, object]] = []
    converged = False
    terminal_iteration_attempted = False
    closure_unavailable_initial = False
    # Iterate species terminal rates with the barrier, Eq 70 profile, and expander potential as one closure
    for iteration in range(1, max(1, int(controls.total_current_balance_iterations)) + 1):
        terminal_iteration_attempted = True
        observed_target = _terminal_current_target_rates(current_expander, current_kinetic)
        relaxed_target = _relax_rates(current_target, observed_target, relaxation)
        candidate_kinetic_started = perf_counter()
        candidate_kinetic = rebuild_kinetic_for_ion_current_components(config=config, geometry=geometry, beam=beam, kinetic=current_kinetic, electron_temperature_J=float(electron_temperature_J), ion_loss_rates_s_by_component=relaxed_target, evaluation_mode=OperatingPointEvaluationMode.TEMPERATURE_TRIAL, collision_state_recompute_policy='rebuilt_for_Pass13_terminal_current_fixed_point')
        candidate_kinetic_runtime_s = perf_counter() - candidate_kinetic_started
       
        record_runtime(runtime_profile, "terminal_current_kinetic_rebuild", candidate_kinetic_runtime_s)
      
        runtime_kinetic_evaluations.append({"kind": "terminal_current_trial", "iteration": iteration, "runtime_s": candidate_kinetic_runtime_s, "profile": candidate_kinetic.metadata.get("runtime_kinetic_profile")})
        candidate_expander_started = perf_counter()
        candidate_expander = solve_expander_system(config, geometry, candidate_kinetic)
       
        record_runtime(runtime_profile, "terminal_current_expander_solve", perf_counter() - candidate_expander_started)
       
        candidate_observed = _terminal_current_target_rates(candidate_expander, candidate_kinetic)
        component_changes = _component_changes(relaxed_target, candidate_observed)
        component_change = max(component_changes.values(), default=0.0)
        barrier = _barrier_energy_J(candidate_kinetic)
        eq70 = _eq70_profile_J(candidate_kinetic)
        barrier_change = _relative_scalar_change(barrier, previous_barrier)
        eq70_change = _relative_array_change(eq70, previous_eq70)
        potential_change = _expander_potential_change(candidate_expander, current_expander)
        Eq68_internal_current_residual = float("inf")
        terminal_ion_to_Eq68_target_residual = float("inf")
        candidate_system = candidate_kinetic.fast_ion_system_state

        if candidate_system is not None and candidate_system.shared_current_balance is not None:
            balance = candidate_system.shared_current_balance
            Eq68_target_current = float(balance.ion_current_loss.total_ion_current_A)
            Eq68_internal_current_residual = abs(float(balance.current_residual_A)) / max(abs(Eq68_target_current), np.finfo(float).tiny)
            terminal_ion_to_Eq68_target_residual = _relative_scalar_change(float(candidate_expander.terminal_current_A), Eq68_target_current)

        iteration_converged = bool(component_change <= tolerance and barrier_change <= tolerance and (eq70_change <= tolerance) and (potential_change <= tolerance) and (Eq68_internal_current_residual <= tolerance) and (terminal_ion_to_Eq68_target_residual <= tolerance) and candidate_expander.potential_converged)
        history.append({
            "iteration": iteration,
            "component_relative_changes": component_changes,
            "maximum_component_relative_change": component_change,
            "barrier_relative_change": barrier_change,
            "Eq70_profile_relative_change": eq70_change,
            "expander_potential_relative_change": potential_change,
            "Eq68_internal_current_relative_residual": Eq68_internal_current_residual,
            "terminal_ion_to_Eq68_target_current_relative_residual": terminal_ion_to_Eq68_target_residual,
            "expander_potential_converged": candidate_expander.potential_converged,
            "iteration_converged": iteration_converged,
        })
        current_kinetic = candidate_kinetic
        current_expander = candidate_expander
        current_target = relaxed_target
        previous_barrier = barrier
        previous_eq70 = eq70

        if iteration_converged:
            converged = True
            break

    # Perform one unrelaxed rebuild at the final observed terminal rates
    final_target = _terminal_current_target_rates(current_expander, current_kinetic)
    final_kinetic_started = perf_counter()
    final_kinetic = rebuild_kinetic_for_ion_current_components(config=config, geometry=geometry, beam=beam, kinetic=current_kinetic, electron_temperature_J=float(electron_temperature_J), ion_loss_rates_s_by_component=final_target, evaluation_mode=mode, collision_state_recompute_policy='rebuilt_for_final_Pass13_terminal_current_consistency')
    final_kinetic_runtime_s = perf_counter() - final_kinetic_started

    record_runtime(runtime_profile, "final_terminal_current_kinetic_rebuild", final_kinetic_runtime_s)

    runtime_kinetic_evaluations.append({"kind": "terminal_current_final", "iteration": len(history) + 1, "runtime_s": final_kinetic_runtime_s, "profile": final_kinetic.metadata.get("runtime_kinetic_profile")})
    final_expander_started = perf_counter()
    final_expander = solve_expander_system(config, geometry, final_kinetic)

    record_runtime(runtime_profile, "final_expander_solve", perf_counter() - final_expander_started)

    beam_density_consistency = postclosure_beam_density_consistency_metadata(config, kinetic, final_kinetic)
    final_observed = _terminal_current_target_rates(final_expander, final_kinetic)
    final_component_changes = _component_changes(final_target, final_observed)
    final_component_consistency = max(final_component_changes.values(), default=0.0)
    final_consistency = final_component_consistency
    final_system = final_kinetic.fast_ion_system_state
    Eq68_internal_current_residual = float("inf")
    active_Eq68_target_current = float("nan")
    Eq68_electron_current = float("nan")

    if final_system is not None and final_system.shared_current_balance is not None:
        balance = final_system.shared_current_balance
        active_Eq68_target_current = float(balance.ion_current_loss.total_ion_current_A)
        Eq68_electron_current = float(balance.electron_balance.electron_current_A)
        Eq68_internal_current_residual = abs(float(balance.current_residual_A)) / max(abs(active_Eq68_target_current), np.finfo(float).tiny)

    terminal_ion_to_Eq68_target_residual = _relative_scalar_change(float(final_expander.terminal_current_A), active_Eq68_target_current)
    end_metadata = _terminal_end_metadata(final_expander, asymmetry_tolerance=asymmetry_tolerance, unresolved_tolerance=unresolved_tolerance)
    prompt_unresolved = end_metadata["terminal_current_prompt_distribution_available"] is not True
    fast_state_numerical = bool(final_kinetic.metadata.get("kinetic_numerical_convergence_passed"))
    expander_intrinsic = bool(final_expander.potential_converged and final_expander.fixed_fast_boundary_applicable and final_expander.metadata.get("expander_branch_conservation_passed") is True and not prompt_unresolved)
    final_iteration = history[-1] if history else {}
    outer_fixed_point_converged = bool(history and float(final_iteration.get('maximum_component_relative_change', float('inf'))) <= tolerance and (float(final_iteration.get('barrier_relative_change', float('inf'))) <= tolerance) and (float(final_iteration.get('Eq70_profile_relative_change', float('inf'))) <= tolerance) and (float(final_iteration.get('expander_potential_relative_change', float('inf'))) <= tolerance) and (float(final_iteration.get('terminal_ion_to_Eq68_target_current_relative_residual', float('inf'))) <= tolerance) and (final_consistency <= tolerance))
    numerical_fixed_point = bool(converged and final_consistency <= tolerance and (Eq68_internal_current_residual <= tolerance) and (terminal_ion_to_Eq68_target_residual <= tolerance) and final_expander.potential_converged)
    end_model_applicable = bool(expander_intrinsic and end_metadata["terminal_current_end_asymmetry_check_passed"] is True)
    qualification_base = bool(numerical_fixed_point and end_model_applicable and fast_state_numerical and (final_expander.status == 'qualified'))
    failure_reasons: list[str] = []

    if not converged:
        failure_reasons.append("terminal_current_outer_fixed_point_not_converged")
    if final_component_consistency > tolerance:
        failure_reasons.append("final_terminal_component_rates_inconsistent")
    if not final_expander.potential_converged:
        failure_reasons.append("expander_quasineutral_potential_not_converged")
    if not final_expander.fixed_fast_boundary_applicable:
        failure_reasons.append("fast_ion_fixed_boundary_inapplicable_due_to_material_expander_return")
    if prompt_unresolved:
        failure_reasons.append("prompt_ion_terminal_distribution_unavailable")
    if end_metadata["terminal_current_end_asymmetry_check_passed"] is not True:
        failure_reasons.append("terminal_side_resolved_ion_current_asymmetry_outside_tolerance")
    if mode.runs_final_qualification and not fast_state_numerical:
        failure_reasons.append("final_fast_ion_state_not_converged")

    failure_reasons.extend(str(value) for value in final_expander.metadata.get("expander_failure_reasons", ()))
    failure_reasons = list(dict.fromkeys(failure_reasons))
    terminal_scope = 'fast_D_fast_T_prompt_D_prompt_T_terminal_surface_current'
    prompt_fallback_scope = 'fast_D_fast_T_terminal_surface_plus_prompt_D_prompt_T_throat_fallback_diagnostic'
    terminal_state_unavailable = bool(not final_expander.potential_converged)
    active_scope = pass11_scope if terminal_state_unavailable else prompt_fallback_scope if prompt_unresolved else terminal_scope
    terminal_component_rates = final_observed
    terminal_component_current = ELECTRON_CHARGE_C * float(sum(terminal_component_rates.values()))
    terminal_current_residual_A = Eq68_electron_current - terminal_component_current
    terminal_current_relative_residual = abs(terminal_current_residual_A) / max(abs(Eq68_electron_current), abs(terminal_component_current), np.finfo(float).tiny)

    if terminal_current_relative_residual > tolerance:
        failure_reasons.append("Eq68_terminal_electron_current_not_equal_to_terminal_ion_current")

    terminal_component_identity_error = _relative_scalar_change(terminal_component_current, active_Eq68_target_current)

    if terminal_component_identity_error > tolerance:
        failure_reasons.append("terminal_component_current_not_equal_to_active_Eq68_target")

    terminal_component_identity_passed = bool(not terminal_state_unavailable and terminal_component_identity_error <= tolerance)
    terminal_side_sum_error = _relative_scalar_change(float(end_metadata["terminal_current_side_sum_A"]), terminal_component_current)
    terminal_side_sum_passed = bool(not terminal_state_unavailable and end_metadata["terminal_current_end_asymmetry_applicable"] and terminal_side_sum_error <= tolerance)

    if end_metadata["terminal_current_end_asymmetry_applicable"] and not terminal_side_sum_passed:
        failure_reasons.append("terminal_side_resolved_current_sum_does_not_match_component_sum")

    qualified = bool(qualification_base and terminal_current_relative_residual <= tolerance and terminal_component_identity_passed and terminal_side_sum_passed)
    terminal_current_authoritative = bool(numerical_fixed_point and not prompt_unresolved and terminal_component_identity_passed and terminal_side_sum_passed and terminal_current_relative_residual <= tolerance)
    failure_reasons = list(dict.fromkeys(failure_reasons))
    active_total_current = _active_total_current_values(terminal_state_unavailable=terminal_state_unavailable, pass11_metadata=kinetic.metadata, pass11_target_current=pass11_target_current, pass11_electron_current=pass11_electron_current, terminal_target_current=terminal_component_current, Eq68_electron_current=Eq68_electron_current, terminal_component_current=terminal_component_current, terminal_component_identity_error=terminal_component_identity_error, terminal_component_identity_passed=terminal_component_identity_passed)
    unclassified_fraction = float(final_expander.total_unclassified_particle_rate_s) / max(float(final_expander.total_throat_particle_rate_s), np.finfo(float).tiny)
    population_classification_fraction = max(unclassified_fraction, float(final_expander.fast_nonterminal_fraction))
    population_classification_passed = final_expander.metadata.get("expander_population_classification_passed") is True
    lost_population_status = "Egedal_Eq63_source_connected_expander_mapping"
    metadata = {**final_expander.metadata, 'terminal_current_balance_model': 'Eq68_electron_current_with_Eq70_center_nodal_expander_density', 'startup_seed_excluded_from_terminal_current': startup_seed_only, 'terminal_current_balance_scope': terminal_scope, 'terminal_current_balance_active_scope': active_scope, 'terminal_current_balance_evaluation_mode': mode.value, 'terminal_current_balance_attempted': terminal_iteration_attempted, 'terminal_current_balance_initial_expander_closure_unavailable': closure_unavailable_initial, 'terminal_current_balance_skipped_due_to_unavailable_expander_closure': False, 'terminal_current_balance_aborted_due_to_unavailable_expander_closure': False, 'terminal_current_balance_authoritative': terminal_current_authoritative, 'terminal_current_balance_iterations': len(history), 'terminal_current_balance_iteration_history': tuple(history), 'terminal_current_balance_relative_tolerance': tolerance, 'terminal_current_balance_relaxation': relaxation, 'terminal_current_balance_final_component_relative_changes': final_component_changes, 'terminal_current_balance_final_component_consistency_relative_change': final_component_consistency, 'terminal_current_balance_final_consistency_relative_change': final_consistency, 'terminal_current_balance_target_component_rates_s': terminal_component_rates, 'terminal_current_balance_observed_component_rates_s': final_observed, 'terminal_current_balance_target_ion_current_A': terminal_component_current, 'terminal_current_balance_active_Eq68_target_ion_current_A': active_Eq68_target_current, 'terminal_current_balance_electron_current_A': Eq68_electron_current, 'terminal_current_balance_Eq68_electron_current_A': Eq68_electron_current, 'terminal_current_balance_terminal_ion_to_Eq68_target_current_relative_residual': terminal_ion_to_Eq68_target_residual, 'terminal_current_balance_current_residual_A': terminal_current_residual_A, 'terminal_current_balance_current_relative_residual': terminal_current_relative_residual, 'terminal_current_balance_active_Eq68_internal_current_relative_residual': Eq68_internal_current_residual, 'terminal_current_balance_component_current_sum_A': terminal_component_current, 'terminal_current_balance_component_current_identity_relative_error': terminal_component_identity_error, 'terminal_current_balance_component_current_identity_check_passed': terminal_component_identity_passed, 'terminal_current_balance_side_sum_identity_relative_error': terminal_side_sum_error, 'terminal_current_balance_side_sum_identity_relative_tolerance': tolerance, 'terminal_current_balance_side_sum_identity_applicable': bool(not terminal_state_unavailable and end_metadata['terminal_current_end_asymmetry_applicable']), 'terminal_current_balance_side_sum_identity_check_passed': terminal_side_sum_passed, 'terminal_current_balance_numerical_fixed_point_converged': numerical_fixed_point, 'terminal_current_balance_end_model_applicable': end_model_applicable, 'terminal_current_balance_qualified': qualified, 'eq68_root_residual_passed': bool(Eq68_internal_current_residual <= tolerance), 'eq68_root_relative_residual': Eq68_internal_current_residual, 'total_current_numerical_fixed_point_converged': outer_fixed_point_converged, 'equivalent_end_model_applicable': end_model_applicable, 'terminal_current_balance_failure_reason': None if qualified else ';'.join(failure_reasons or ['terminal_current_balance_not_qualified']), 'terminal_current_balance_failure_reasons': tuple(failure_reasons), 'terminal_current_balance_prompt_rate_uses_throat_fallback_for_diagnostic_iteration': bool(final_expander.prompt_unresolved_particle_rate_s > 0.0), 'terminal_current_balance_prompt_rate_fallback_is_qualified': False, 'throat_crossing_and_terminal_loss_are_distinct': True, 'Pass11_throat_current_balance_scope': pass11_scope, 'Pass11_throat_current_balance_numerical_fixed_point_converged': pass11_numerical, 'Pass11_throat_current_balance_qualified': pass11_qualified, 'Pass11_throat_current_balance_target_ion_current_A': pass11_target_current, 'Pass11_throat_current_balance_electron_current_A': pass11_electron_current, 'current_balance_scope': active_scope, 'active_current_boundary_definition': 'Pass11_throat_crossing_diagnostic_pending_expander_terminal_state' if terminal_state_unavailable else 'terminal_surface_plus_prompt_throat_fallback_diagnostic' if prompt_unresolved else 'terminal_surface', 'total_current_balance_active_state_is_Pass11_throat_diagnostic': active_total_current['Pass11_throat_diagnostic'], 'total_current_balance_numerical_fixed_point_converged': numerical_fixed_point, 'total_current_balance_end_model_applicable': end_model_applicable, 'total_current_balance_current_relative_residual': active_total_current['current_relative_residual'], 'total_current_balance_current_residual_A': active_total_current['current_residual_A'], 'total_current_balance_target_ion_current_A': active_total_current['target_ion_current_A'], 'total_current_balance_electron_current_A': active_total_current['electron_current_A'], 'total_current_balance_component_current_sum_A': active_total_current['component_current_sum_A'], 'total_current_balance_component_current_identity_relative_error': active_total_current['component_current_identity_relative_error'], 'total_current_balance_component_current_identity_check_passed': active_total_current['component_current_identity_check_passed'], 'total_current_balance_end_asymmetry_applicable': bool(not terminal_state_unavailable and end_metadata['terminal_current_end_asymmetry_applicable']), 'total_current_balance_end_asymmetry_check_passed': bool(not terminal_state_unavailable and end_metadata['terminal_current_end_asymmetry_check_passed']), 'total_current_balance_end_asymmetry_relative_error': float(end_metadata['terminal_current_end_asymmetry_relative_error']), 'total_current_balance_end_asymmetry_relative_tolerance': float(end_metadata['terminal_current_end_asymmetry_relative_tolerance']), 'total_current_balance_side_resolved_left_ion_current_A': float(end_metadata['terminal_current_left_ion_current_A']), 'total_current_balance_side_resolved_right_ion_current_A': float(end_metadata['terminal_current_right_ion_current_A']), 'total_current_balance_side_sum_identity_applicable': bool(not terminal_state_unavailable and end_metadata['terminal_current_end_asymmetry_applicable']), 'total_current_balance_side_sum_identity_check_passed': terminal_side_sum_passed, 'total_current_balance_side_sum_identity_relative_error': terminal_side_sum_error, 'total_current_balance_side_sum_identity_relative_tolerance': tolerance, 'total_current_balance_converged': qualified, 'total_current_balance_failure_reason': None if qualified else ';'.join(failure_reasons or ['terminal_current_balance_not_qualified']), 'total_current_balance_failure_reasons': tuple(failure_reasons), 'total_confined_plasma_current_balance_claimed': qualified, 'full_device_wall_current_balance_claimed': False, 'global_current_balance_claimed': False, 'terminal_current_used_by_Eq68': terminal_current_authoritative, 'terminal_plus_prompt_throat_fallback_used_by_Eq68': bool(not terminal_state_unavailable and prompt_unresolved), 'Pass11_throat_current_retained_as_diagnostic_only': bool(not terminal_state_unavailable), 'Pass11_throat_current_retained_as_active_diagnostic_state': terminal_state_unavailable, 'modal_current_balance_relative_error': Eq68_internal_current_residual, 'kinetic_expander_unclassified_population_check_passed': population_classification_passed, 'unclassified_population_fraction': population_classification_fraction, 'expander_electron_profile_available': True, 'expander_electron_profile_complete': bool(not terminal_state_unavailable), 'kinetic_full_device_lost_population_check_passed': population_classification_passed, 'full_device_lost_ion_reconstruction_status': lost_population_status, 'full_device_lost_ion_reconstruction_available': True, 'global_plasma_power_balance_available': False, 'global_device_power_balance_claimed': False, **beam_density_consistency, **end_metadata}
    final_population_attachment_started = perf_counter()
    final_kinetic = _mark_terminal_current_state(final_kinetic, final_expander, metadata)
    final_kinetic = _attach_full_device_populations(geometry, final_kinetic, final_expander)

    record_runtime(runtime_profile, "final_population_attachment", perf_counter() - final_population_attachment_started)

    runtime_metadata = {"runtime_expander_profile": finalize_runtime_profile(runtime_profile, total_s=perf_counter() - expander_runtime_started), "runtime_expander_kinetic_evaluations": tuple(runtime_kinetic_evaluations)}
    metadata = {**metadata, **runtime_metadata}
    final_kinetic = replace(final_kinetic, metadata={**final_kinetic.metadata, **runtime_metadata})
    expander_metadata = {**final_expander.metadata, **metadata}
    expander_result = ExpanderStageResult(system_state=replace(final_expander, metadata=expander_metadata), metadata=expander_metadata)
  
    return final_kinetic, expander_result

__all__ = ["build_expander_stage"]
