"""Shared multi species modal FBIS orchestration for fast deuterium and tritium"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import replace
from time import perf_counter
from typing import Any
import numpy as np
from source_model_revamp.electrostatic.current_balance import AmbipolarCurrentBalance, build_ion_current_loss, egedal_tail_refilling_current_balance_for_actual_midplane_density, egedal_tail_refilling_current_balance_for_parent_maxwellian_density
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState, rebuild_fbis_collision_parameter_state_with_fast_density
from source_model_revamp.fbis.modal.basis import build_modal_fbis_basis
from source_model_revamp.fbis.modal.electrostatic_feedback import rescaled_magnetic_profile_for_reference_ratio
from source_model_revamp.fbis.modal.electrostatic.system_iteration import ModalElectrostaticSpeciesInput, ModalElectrostaticSystemResult, _solve_system_phi_profile_quasineutrality
from source_model_revamp.fbis.modal.local.speed import assess_local_physical_speed_grid_convergence, derive_local_physical_speed_grid
from source_model_revamp.fbis.modal.lost_ion_distribution import BALDWIN_1972_THROAT_DENSITY_CLOSURE, EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY, build_fixed_boundary_eq63_branch, calculate_baldwin_throat_parallel_temperature, ion_loss_model
from source_model_revamp.fbis.modal.models import EGEDAL_2022_PUBLISHED_ELECTROSTATIC_FBIS_APPROXIMATION, MAGNETIC_ONLY_FBIS, ZERO_POTENTIAL_MAGNETIC_REFERENCE, electrostatic_feedback_model as canonical_electrostatic_feedback_model
from source_model_revamp.fbis.modal.solver.solver_core import _solve_modal_fbis_state
from source_model_revamp.fbis.modal.solver.solver_loss_rates import _conservative_eq61_mode_slopes_dI_dlambda
from source_model_revamp.fbis.modal.species_state import FastIonSpeciesRequest, FastIonSpeciesState, FastIonSystemState
from source_model_revamp.fbis.modal.types import ModalFBISBasis, ModalFBISResult, ModalLocalReconstruction
from source_model_revamp.runtime_profile import record_eq70_confined_interpolation_audit_call, record_runtime

def _validate_collision_state(request: FastIonSpeciesRequest, collision_state: FBISCollisionParameterState | None) -> None:
    """Validate that the collision state and reduced Eq 59 projection belong to the requested species"""
    if collision_state is None:
        raise ValueError("an active fast ion species solve requires a collision state")
    if collision_state.fast_ion_species != request.species:
        raise ValueError("the collision state fast ion species must match the requested species")
    if collision_state.pairwise_collision_state.test_species != request.species:
        raise ValueError("the pairwise collision test species must match the requested species")
    projection = collision_state.pairwise_collision_state.reduced_eq59_projection
    if projection.spitzer_slowing_down_time_s != collision_state.spitzer_slowing_down_time_s or projection.critical_velocity_m_s != collision_state.critical_velocity_m_s or projection.beta_m != collision_state.beta_m:
        raise ValueError("the legacy scalar collision diagnostic must exactly reproduce the stored scalar collision state")

def _solve_fast_ion_species_direct(*, request: FastIonSpeciesRequest, collision_state: FBISCollisionParameterState | None, solver_arguments: Mapping[str, Any]) -> FastIonSpeciesState:
    """Solve one species for a fixed collision state or return an inactive zero source state"""
    source_rate = float(request.attenuated_source.total_birth_rate_s)
    if source_rate <= 0.0:
        return FastIonSpeciesState(species=request.species, speed_grid=request.speed_grid, active=False, source_particle_rate_s=source_rate, modal_result=None, status="inactive_zero_source")
    _validate_collision_state(request, collision_state)
    arguments = dict(solver_arguments)
    external_by_test_species = dict(arguments.pop("_external_fast_field_collision_states_by_test_species", {}) or {})
    arguments.pop("_runtime_eq70_confined_interpolation_audit_enabled", None)
    arguments.pop("_runtime_eq70_solve_context", None)
    arguments.pop("_runtime_eq70_fast_system_solve_index", None)
    if "_external_fast_field_collision_states" not in arguments:
        arguments["_external_fast_field_collision_states"] = tuple(external_by_test_species.get(request.species.species_id, ()))
    modal_result = _solve_modal_fbis_state(fast_ion_species=request.species, speed_grid=request.speed_grid, attenuated_source=request.attenuated_source, collision_state=collision_state, _charge_exchange_sink_state=request.charge_exchange_sink_state, **arguments)
    operator_state = modal_result.eq59_collision_operator_state
    if operator_state is not None:
        if operator_state.test_species != request.species:
            raise ValueError("the final Eq 59 collision operator species must match the requested species")
        density_scale = max(abs(float(operator_state.fast_self_physical_density_m3)), abs(float(modal_result.density_m3)), 1.0)
        density_tolerance = 256.0 * np.finfo(float).eps * density_scale
        if abs(float(operator_state.fast_self_physical_density_m3) - float(modal_result.density_m3)) > density_tolerance:
            raise ValueError("the final Eq 59 fast self density must match the solved modal density")
    solved_collision_state = rebuild_fbis_collision_parameter_state_with_fast_density(collision_state, modal_result.density_m3)
   
    return FastIonSpeciesState(species=request.species, speed_grid=request.speed_grid, active=True, source_particle_rate_s=source_rate, modal_result=modal_result, collision_state=solved_collision_state, status="active_solved")

def _shared_current_balance(*, states: Mapping[str, FastIonSpeciesState], electron_midplane_density_m3: float, electron_collision_density_m3: float, electron_temperature_J: float, electron_collision_frequency_s: float, midplane_area_m2: float, half_length_m: float, mirror_ratio: float, basis: ModalFBISBasis, ion_loss_rates_s_by_component_override: Mapping[str, float] | None=None) -> AmbipolarCurrentBalance:
    """Build one Eq 68 electron barrier from the summed fast and prompt ion loss currents"""
    active = tuple(states[key] for key in sorted(states) if states[key].modal_result is not None)
    if not active:
        raise ValueError("shared current balance requires at least one confined fast species")
    electron_frequency = float(electron_collision_frequency_s)
    if not np.isfinite(electron_frequency) or electron_frequency <= 0.0:
        raise ValueError("shared electron collision frequency must be finite and positive")
    override = None if ion_loss_rates_s_by_component_override is None else {str(key): float(value) for key, value in dict(ion_loss_rates_s_by_component_override).items()}
    allowed_components = {'fast_deuterium', 'fast_tritium', 'prompt_deuterium', 'prompt_tritium'}
    if override is not None:
        if not set(override).issubset(allowed_components):
            raise ValueError("ion loss current override contains an unsupported component")
        if any(not np.isfinite(value) or value < 0.0 for value in override.values()):
            raise ValueError("ion loss current override rates must be finite and nonnegative")
        if not override:
            raise ValueError("ion loss current override must contain at least one component")
        rates = [override[key] for key in sorted(override)]
        charges = [1.0 for _ in rates]
    else:
        source_states = tuple(states[key] for key in sorted(states) if states[key].modal_result is not None or float(states[key].prompt_only_particle_loss_rate_s) > 0.0)
        rates = [float(state.modal_result.ion_particle_loss_rate_s) if state.modal_result is not None else float(state.prompt_only_particle_loss_rate_s) for state in source_states]
        charges = [float(state.species.charge_number) for state in source_states]
    ion_current_loss = build_ion_current_loss(ion_particle_loss_rates_s=np.asarray(rates, dtype=float), ion_charge_numbers=np.asarray(charges, dtype=float))
    first_slope = float(_conservative_eq61_mode_slopes_dI_dlambda(basis)[0])
    electron_balance = egedal_tail_refilling_current_balance_for_actual_midplane_density(target_current_A=float(ion_current_loss.total_ion_current_A), electron_midplane_density_m3=float(electron_midplane_density_m3), electron_collision_density_m3=float(electron_collision_density_m3), electron_temperature_J=float(electron_temperature_J), midplane_area_m2=float(midplane_area_m2), half_length_m=float(half_length_m), mirror_ratio=float(mirror_ratio), loss_geometry_factor=float(basis.geometry_factor_G), loss_cone_slope_dI_dlambda=first_slope, electron_collision_frequency_s=electron_frequency, number_of_ends=2.0)
  
    return AmbipolarCurrentBalance(ion_current_loss=ion_current_loss, electron_balance=electron_balance, current_residual_A=electron_balance.current_residual_A, electron_wall_potential_relative_to_midplane_V=electron_balance.wall_potential_relative_to_midplane_V)

def _attach_baldwin_eq63_branches(*, states: Mapping[str, FastIonSpeciesState], solver_arguments: Mapping[str, Any]) -> dict[str, FastIonSpeciesState]:
    """Attach left and right fixed magnetic boundary Eq 63 branches before the shared Eq 70 solve"""
    model = str(solver_arguments.get("lost_ion_parallel_temperature_model", BALDWIN_1972_THROAT_DENSITY_CLOSURE)).strip().lower()
    if model != BALDWIN_1972_THROAT_DENSITY_CLOSURE:
        raise ValueError("production lost ion parallel temperature model must be baldwin_1972_throat_density_closure")
    result: dict[str, FastIonSpeciesState] = {}
    for species_id, state in states.items():
        modal = state.modal_result
        if modal is None:
            result[species_id] = state
            continue
        if state.collision_state is None:
            raise ValueError("Baldwin Eq 63 reconstruction requires the solved collision state")
        one_end_rate_v_s = 0.5 * np.asarray(modal.boundary_flux_density_v, dtype=float) * np.asarray(modal.speed_grid.shell_volumes_m3_s3, dtype=float) * float(solver_arguments["volume_m3"])
        common = dict(
            speed_grid=modal.speed_grid,
            one_end_eq61_rate_v_s=one_end_rate_v_s,
            dfdlambda_boundary_v=np.asarray(modal.dF_dlambda_boundary_v, dtype=float),
            mirror_ratio=float(solver_arguments["mirror_ratio"]),
            midplane_area_m2=float(solver_arguments["midplane_area_m2"]),
            half_length_m=float(solver_arguments["half_length_m"]),
            geometry_factor_G=float(modal.basis.geometry_factor_G),
            particle_mass_kg=float(state.species.mass_kg),
            collision_state=state.collision_state,
            collision_operator_state=modal.eq59_collision_operator_state,
            volume_geometry_integral=float(modal.basis.volume_geometry_integral),
            volume_m3=float(solver_arguments["volume_m3"]),
        )
        left_temperature = calculate_baldwin_throat_parallel_temperature(side="left", **common)
        right_temperature = calculate_baldwin_throat_parallel_temperature(side="right", **common)
        left_branch = build_fixed_boundary_eq63_branch(
            side="left",
            speed_grid=modal.speed_grid,
            lambda_grid=modal.lambda_grid,
            one_end_eq61_rate_v_s=one_end_rate_v_s,
            mirror_ratio=float(solver_arguments["mirror_ratio"]),
            midplane_area_m2=float(solver_arguments["midplane_area_m2"]),
            geometry_factor_G=float(modal.basis.geometry_factor_G),
            parallel_temperature_J=float(left_temperature.eq63_parallel_temperature_J),
            particle_mass_kg=float(state.species.mass_kg),
        )
        right_branch = build_fixed_boundary_eq63_branch(
            side="right",
            speed_grid=modal.speed_grid,
            lambda_grid=modal.lambda_grid,
            one_end_eq61_rate_v_s=one_end_rate_v_s,
            mirror_ratio=float(solver_arguments["mirror_ratio"]),
            midplane_area_m2=float(solver_arguments["midplane_area_m2"]),
            geometry_factor_G=float(modal.basis.geometry_factor_G),
            parallel_temperature_J=float(right_temperature.eq63_parallel_temperature_J),
            particle_mass_kg=float(state.species.mass_kg),
        )
        metadata = dict(modal.metadata)
        metadata.update({
            "lost_ion_parallel_temperature_model": BALDWIN_1972_THROAT_DENSITY_CLOSURE,
            "modal_lost_ion_parallel_temperature_model": BALDWIN_1972_THROAT_DENSITY_CLOSURE,
            "fixed_boundary_lost_reconstruction_status": "represented",
            "fixed_boundary_lost_reconstruction_failure_reason": None,
            "fixed_boundary_lost_reconstruction_complete": True,
            "modal_lost_ion_distribution_model": "egedal_eq63_directed_left_and_right_fixed_magnetic_boundary_with_baldwin_1972_throat_density_closure",
            "modal_eq63_eq64_active_implementation_status": "implemented_with_analytic_continuous_loss_cone_flux_matching",
            "modal_eq63_pitch_integration_model": "analytic_exact_magnetic_loss_cone_boundary",
            "modal_eq61_left_one_end_rate_v_s": one_end_rate_v_s,
            "modal_eq61_right_one_end_rate_v_s": one_end_rate_v_s,
            "lost_ion_parallel_temperature_left_J": float(left_temperature.eq63_parallel_temperature_J),
            "lost_ion_parallel_temperature_right_J": float(right_temperature.eq63_parallel_temperature_J),
            "lost_ion_parallel_temperature_left_diagnostics": dict(left_temperature.diagnostics or {}),
            "lost_ion_parallel_temperature_right_diagnostics": dict(right_temperature.diagnostics or {}),
            "modal_eq63_left_throat_crossing_rate_s": float(left_branch.reconstructed_total_rate_s),
            "modal_eq63_right_throat_crossing_rate_s": float(right_branch.reconstructed_total_rate_s),
            "modal_eq63_left_rate_relative_error": float(left_branch.relative_rate_error),
            "modal_eq63_right_rate_relative_error": float(right_branch.relative_rate_error),
            "modal_eq63_left_H_U": left_branch.H_U,
            "modal_eq63_right_H_U": right_branch.H_U,
            "modal_eq63_left_throat_distribution_v_lambda": left_branch.throat_distribution_v_lambda,
            "modal_eq63_right_throat_distribution_v_lambda": right_branch.throat_distribution_v_lambda,
            "modal_eq63_left_throat_rate_v_lambda_s": left_branch.throat_rate_v_lambda_s,
            "modal_eq63_right_throat_rate_v_lambda_s": right_branch.throat_rate_v_lambda_s,
            "modal_directed_eq63_throat_rate_v_lambda_s": right_branch.throat_rate_v_lambda_s,
            "modal_directed_eq63_throat_rate_per_end_s": float(right_branch.reconstructed_total_rate_s),
            "modal_eq63_reference_compatible": True,
            "modal_active_hot_electrostatic_eq63_population_available": True,
            "modal_eq64_matching_model": right_branch.matching_model,
            "baldwin_1972_fixed_magnetic_boundary_specialization": True,
            "baldwin_1972_full_anisotropic_collision_operator_claimed": False,
        })
        result[species_id] = replace(state, modal_result=replace(modal, metadata=metadata))
   
    return result

def _current_balance_from_floating_profile(*, current_balance: AmbipolarCurrentBalance, electrostatic: ModalElectrostaticSystemResult, solver_arguments: Mapping[str, Any]) -> AmbipolarCurrentBalance:
    """Recompute Eq 68 from the Eq 70 parent Maxwellian n0 and exact midplane reference"""
    electron = current_balance.electron_balance
    if electron.electron_collision_frequency_s is None or electron.loss_geometry_factor is None or electron.loss_cone_slope_dI_dlambda is None:
        raise ValueError("Eq 68 state is missing the collision or geometry quantities needed for floating closure")
    profile = electrostatic.profile
    updated_electron = egedal_tail_refilling_current_balance_for_parent_maxwellian_density(
        target_current_A=float(current_balance.ion_current_loss.total_ion_current_A),
        electron_parent_maxwellian_n0_m3=float(profile.electron_parent_maxwellian_n0_m3),
        electron_midplane_density_m3=float(profile.electron_midplane_density_m3),
        electron_collision_density_m3=float(profile.electron_collision_density_m3),
        electron_temperature_J=float(solver_arguments["electron_temperature_J"]),
        midplane_area_m2=float(solver_arguments["midplane_area_m2"]),
        half_length_m=float(solver_arguments["half_length_m"]),
        mirror_ratio=float(solver_arguments["mirror_ratio"]),
        loss_geometry_factor=float(electron.loss_geometry_factor),
        loss_cone_slope_dI_dlambda=float(electron.loss_cone_slope_dI_dlambda),
        electron_collision_frequency_s=float(electron.electron_collision_frequency_s),
        exact_midplane_potential_drop_J=float(profile.exact_midplane_potential_energy_J),
        number_of_ends=float(electron.number_of_ends),
    )
 
    return AmbipolarCurrentBalance(
        ion_current_loss=current_balance.ion_current_loss,
        electron_balance=updated_electron,
        current_residual_A=updated_electron.current_residual_A,
        electron_wall_potential_relative_to_midplane_V=updated_electron.wall_potential_relative_to_midplane_V,
    )

def _active_fast_cross_pair_ids(state: FastIonSpeciesState) -> tuple[str, ...]:
    """Return identifiers for active reduced fast cross species collision pairs"""
    modal = state.modal_result
    operator = None if modal is None else modal.eq59_collision_operator_state
    if operator is None:
        return ()
   
    return tuple(sorted(
        contribution.pair_id
        for contribution in operator.pair_contributions_by_id.values()
        if contribution.active
        and contribution.field_population_id.startswith("fast_")
        and contribution.field_population_id != "fast_self"
    ))

def _invariant_population_weights(modal: ModalFBISResult) -> np.ndarray:
    """Return nonnegative invariant particle weights by speed cell from f(v, η)"""
    eta_integral = np.trapezoid(np.maximum(np.asarray(modal.distribution_v_eta, dtype=float), 0.0), np.asarray(modal.basis.physical_basis.eta_grid, dtype=float), axis=1)
   
    return np.asarray(modal.speed_grid.shell_volumes_m3_s3, dtype=float) * eta_integral

def _shared_electrostatic_state(*, states: Mapping[str, FastIonSpeciesState], current_balance: AmbipolarCurrentBalance, solver_arguments: Mapping[str, Any], initial_potential_energy_J: np.ndarray | None, audit_context: Mapping[str, object] | None = None) -> ModalElectrostaticSystemResult | None:
    """Build species inputs and solve the shared Eq 70 to Eq 72 electrostatic state"""
    pitch_grid = solver_arguments.get("pitch_grid")
    cell_volumes_m3 = solver_arguments.get("cell_volumes_m3")
   
    if pitch_grid is None or cell_volumes_m3 is None:
        return None
  
    basis = next(state.modal_result.basis for state in states.values() if state.modal_result is not None)
    species_inputs: dict[str, ModalElectrostaticSpeciesInput] = {}
    for species_id in sorted(states):
        state = states[species_id]
        modal = state.modal_result
        if modal is None:
            continue
        metadata = dict(modal.metadata)
        species_inputs[species_id] = ModalElectrostaticSpeciesInput(
            species_id=species_id,
            speed_grid=modal.speed_grid,
            lambda_grid=modal.lambda_grid,
            pitch_grid=pitch_grid,
            base_distribution_v_lambda=np.asarray(modal.distribution_v_lambda, dtype=float),
            particle_mass_kg=state.species.mass_kg,
            eta_to_local_phase_space_normalization=float(basis.volume_geometry_integral / basis.eta_lambda_map.lambda_normalization),
            invariant_energy_cell_population_weights=_invariant_population_weights(modal),
            charge_number=float(state.species.charge_number),
            eq63_left_H_U=metadata.get("modal_eq63_left_H_U"),
            eq63_right_H_U=metadata.get("modal_eq63_right_H_U"),
            eq63_left_throat_rate_v_lambda_s=metadata.get("modal_eq63_left_throat_rate_v_lambda_s"),
            eq63_right_throat_rate_v_lambda_s=metadata.get("modal_eq63_right_throat_rate_v_lambda_s"),
            eq63_left_parallel_temperature_J=metadata.get("lost_ion_parallel_temperature_left_J"),
            eq63_right_parallel_temperature_J=metadata.get("lost_ion_parallel_temperature_right_J"),
            eq63_geometry_factor_G=float(modal.basis.geometry_factor_G) if metadata.get("modal_eq63_left_H_U") is not None else None,
        )
    zeta_faces = np.asarray(solver_arguments["zeta_faces"], dtype=float)
    zeta = 0.5 * (zeta_faces[:-1] + zeta_faces[1:])
    runtime_profile = solver_arguments.get("_runtime_profile_accumulator")
    audit_enabled = bool(solver_arguments.get("_runtime_eq70_confined_interpolation_audit_enabled", False)) and isinstance(runtime_profile, dict)
    confined_interpolation_audit: dict[str, object] | None = {} if audit_enabled else None
    result = _solve_system_phi_profile_quasineutrality(species_inputs=species_inputs, zeta=zeta, B_tilde=np.asarray(solver_arguments['B_tilde_midpoints'], dtype=float), cell_volumes_m3=np.asarray(cell_volumes_m3, dtype=float), mirror_ratio=float(solver_arguments['mirror_ratio']), wall_barrier_energy_J=float(current_balance.electron_balance.barrier_energy_J), electron_temperature_J=float(solver_arguments['electron_temperature_J']), iterations=int(solver_arguments['numerics'].phi_z_iterations), root_scan_points=int(solver_arguments['numerics'].phi_z_root_scan_points), relaxation=float(solver_arguments['numerics'].phi_z_relaxation), relative_tolerance=float(solver_arguments['numerics'].phi_z_relative_tolerance), electron_midplane_density_m3=solver_arguments.get('eq70_prescribed_electron_midplane_density_m3'), electron_collision_density_m3=float(solver_arguments['electron_collision_density_m3']), background_positive_charge_density_m3=solver_arguments.get('background_positive_charge_density_m3'), background_midplane_positive_charge_density_m3=solver_arguments.get('background_midplane_positive_charge_density_m3'), background_left_throat_positive_charge_density_m3=float(solver_arguments.get('background_left_throat_positive_charge_density_m3', 0.0)), background_right_throat_positive_charge_density_m3=float(solver_arguments.get('background_right_throat_positive_charge_density_m3', 0.0)), target_volume_averaged_ion_density_m3=None, local_velocity_quadrature_order=int(solver_arguments['numerics'].local_velocity_quadrature_order), low_energy_weight_fraction_tolerance=float(solver_arguments['numerics'].phi_z_low_energy_weight_fraction_tolerance), initial_potential_energy_J=initial_potential_energy_J, zeta_faces=zeta_faces, _confined_interpolation_audit=confined_interpolation_audit)
    if confined_interpolation_audit is not None and isinstance(runtime_profile, dict):
        record_eq70_confined_interpolation_audit_call(runtime_profile, {
            **confined_interpolation_audit,
            "fast_ion_system_solve_index": None if audit_context is None else audit_context.get("fast_ion_system_solve_index"),
            "solve_context": None if audit_context is None else audit_context.get("solve_context"),
            "mixed_pass_kind": None if audit_context is None else audit_context.get("mixed_pass_kind"),
            "coupling_phase": None if audit_context is None else audit_context.get("coupling_phase"),
            "output_converged": bool(result.profile.converged),
            "output_failure_reason": result.profile.failure_reason,
        })
  
    return result

def _finalize_species_result(*, state: FastIonSpeciesState, current_balance: AmbipolarCurrentBalance, electrostatic: ModalElectrostaticSystemResult | None, solver_arguments: Mapping[str, Any]) -> FastIonSpeciesState:
    """Attach the shared barrier, local reconstruction, directed loss state, and final metadata to one species"""
    modal = state.modal_result
    barrier = float(current_balance.electron_balance.barrier_energy_J)
    if modal is None:
        if float(state.prompt_only_particle_loss_rate_s) <= 0.0:
            return state
        prompt_wall_power = float(state.prompt_only_midplane_kinetic_power_loss_W) + float(state.prompt_only_particle_loss_rate_s) * barrier
        return replace(state, prompt_only_ion_wall_power_loss_W=prompt_wall_power, status="prompt_only_with_shared_barrier")
    ion_barrier_power_W = float(modal.ion_particle_loss_rate_s) * barrier
    ion_wall_power_W = float(modal.ion_midplane_kinetic_power_loss_W) + ion_barrier_power_W
    local_reconstruction = None
    metadata = dict(modal.metadata)
    fixed_boundary_status = "not_evaluated_no_pitch_grid"
    if electrostatic is not None:
        species_id = state.species.species_id
        local_lambda = np.asarray(electrostatic.local_distribution_z_v_lambda_by_species[species_id], dtype=float)
        local_pitch = np.asarray(electrostatic.local_distribution_z_v_pitch_by_species[species_id], dtype=float)
        local_density = np.asarray(electrostatic.local_density_m3_by_species[species_id], dtype=float)
        local_speed_grid = electrostatic.local_speed_grid_by_species[species_id]
        _, local_grid_diagnostics = derive_local_physical_speed_grid(invariant_speed_grid=modal.speed_grid, maximum_potential_drop_magnitude_J=barrier, particle_mass_kg=state.species.mass_kg)
        if bool(solver_arguments.get("_run_numerical_convergence_assessments", True)):
            local_refinement = assess_local_physical_speed_grid_convergence(invariant_speed_grid=modal.speed_grid, local_speed_grid=local_speed_grid, lambda_grid=modal.lambda_grid, pitch_grid=solver_arguments["pitch_grid"], base_distribution_v_lambda=np.asarray(modal.distribution_v_lambda, dtype=float), mirror_ratio=float(solver_arguments["mirror_ratio"]), zeta=np.asarray(electrostatic.profile.zeta, dtype=float), B_tilde=np.asarray(electrostatic.profile.B_tilde, dtype=float), potential_drop_magnitude_J=np.asarray(electrostatic.profile.potential_energy_J, dtype=float), left_throat_potential_drop_magnitude_J=float(electrostatic.profile.throat_potential_left_energy_J), right_throat_potential_drop_magnitude_J=float(electrostatic.profile.throat_potential_right_energy_J), eta_to_local_phase_space_normalization=float(modal.basis.volume_geometry_integral / modal.basis.eta_lambda_map.lambda_normalization), cell_volumes_m3=np.asarray(solver_arguments["cell_volumes_m3"], dtype=float), quadrature_order=int(solver_arguments["numerics"].local_velocity_quadrature_order), relative_tolerance=float(solver_arguments["numerics"].local_speed_refinement_relative_tolerance), reference_local_distribution_z_v_pitch=local_pitch, particle_mass_kg=state.species.mass_kg)
        else:
            local_refinement = {"local_grid_refinement_assessed": False, "local_grid_refinement_converged": False, "local_grid_refinement_status": "deferred_during_intermediate_closure_iteration", "local_grid_refinement_failure_reason": "intermediate_closure_state", "local_grid_refinement_relative_tolerance": float(solver_arguments["numerics"].local_speed_refinement_relative_tolerance), "local_grid_refinement_history": [], "eq71_72_population_change_included_in_grid_error": False}
        represented_inventory = float(np.sum(local_density * np.asarray(solver_arguments["cell_volumes_m3"], dtype=float)))
        local_grid_diagnostics = {**local_grid_diagnostics, "remap_inventory_relative_error": ((represented_inventory - float(modal.inventory_particles)) / float(modal.inventory_particles) if float(modal.inventory_particles) > 0.0 else 0.0), "remap_energy_relative_error": None, "remap_energy_status": "not_a_conservative_energy_remap_Eq71_72_is_an_energy_mapping", **local_refinement}
        representative_T_L = None
        loss_model = ion_loss_model(solver_arguments.get("ion_loss_closure_model", EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY))
        policy = str(solver_arguments.get("_fixed_boundary_reconstruction_policy", "required")).strip().lower()
        if loss_model == EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY and policy != "deferred_eq42_iteration":
            left_T = metadata.get("lost_ion_parallel_temperature_left_J")
            right_T = metadata.get("lost_ion_parallel_temperature_right_J")
            left_H = metadata.get("modal_eq63_left_H_U")
            right_H = metadata.get("modal_eq63_right_H_U")
            if left_T is None or right_T is None or left_H is None or right_H is None:
                if policy == "required":
                    raise ValueError("production fixed boundary Eq 63 reconstruction requires the Baldwin throat density closure before Eq 70")
                fixed_boundary_status = "unavailable_baldwin_1972_throat_density_closure"
            else:
                representative_T_L = 0.5 * (float(left_T) + float(right_T))
                fixed_boundary_status = "represented"
        elif policy == "deferred_eq42_iteration":
            fixed_boundary_status = "deferred_during_eq42_density_basis_iteration"
        local_reconstruction = ModalLocalReconstruction(local_distribution_z_v_lambda=local_lambda, local_distribution_z_v_pitch=local_pitch, lost_ion_distribution_z_v_lambda=None, lost_ion_parallel_temperature_J=representative_T_L, lost_ion_geometry_profile_G_z=None, local_density_m3=local_density, local_speed_grid=local_speed_grid, local_speed_grid_diagnostics=local_grid_diagnostics)
        profile = electrostatic.profile
        metadata.update({"modal_phi_z_converged": bool(profile.converged), "modal_phi_z_profile_V": np.asarray(profile.potential_relative_to_midplane_V, dtype=float), "modal_phi_z_failure_reason": profile.failure_reason, "modal_electron_density_profile_m3": np.asarray(profile.electron_density_m3, dtype=float), "modal_fast_ion_density_profile_m3": np.asarray(profile.fast_ion_density_by_species_m3[species_id], dtype=float), "modal_total_ion_density_profile_m3": np.asarray(profile.ion_density_m3, dtype=float), "modal_electron_midplane_density_m3": float(profile.electron_midplane_density_m3), "modal_electron_collision_density_m3": float(profile.electron_collision_density_m3), "modal_electron_volume_average_density_m3": float(profile.electron_volume_average_density_m3), "modal_electron_parent_maxwellian_n0_m3": float(profile.electron_parent_maxwellian_n0_m3), "modal_phi_z_exact_midplane_quasineutrality_valid": bool(profile.exact_midplane_quasineutrality_valid), "modal_phi_z_exact_throat_quasineutrality_valid": bool(profile.exact_throat_quasineutrality_valid), "modal_phi_z_direct_throat_roots_valid": bool(profile.direct_throat_roots_valid), "modal_phi_z_direct_central_roots_valid": bool(profile.direct_central_roots_valid), "modal_phi_z_direct_central_no_root_cell_count": int(profile.direct_central_no_root_cell_count), "modal_phi_z_direct_central_interior_scan_cell_count": int(profile.direct_central_interior_scan_cell_count), "modal_phi_z_direct_central_multiple_root_cell_count": int(profile.direct_central_multiple_root_cell_count), "modal_phi_z_direct_central_maximum_root_count": int(profile.direct_central_maximum_root_count), "modal_phi_z_direct_central_root_scan_points": int(profile.direct_central_root_scan_points), "modal_phi_z_eq71_low_energy_approximation_valid": bool(profile.low_energy_approximation_valid), "modal_phi_z_effective_potential_throat_boundary_valid": bool(profile.effective_potential_throat_boundary_valid), "modal_phi_z_floating_nonnegative_reference_active": bool(profile.floating_nonnegative_reference_active), "modal_phi_z_floating_zero_reference_kind": profile.floating_zero_reference_kind, "modal_phi_z_floating_zero_reference_zeta": profile.floating_zero_reference_zeta, "modal_phi_z_floating_zero_reference_B_tilde": profile.floating_zero_reference_B_tilde, "modal_phi_z_exact_midplane_potential_energy_J": float(profile.exact_midplane_potential_energy_J), "modal_phi_z_floating_scalar_closure_converged": bool(profile.floating_scalar_closure_converged), "modal_phi_z_direct_exact_midplane_root_valid": bool(profile.direct_exact_midplane_root_valid), "modal_phi_z_exact_throat_population_qualification_status": "confined_plus_directed_Eq63_fixed_boundary_with_expander_return_fail_closed", "modal_phi_z_exact_throat_population_qualification_available": True, "modal_phi_z_exact_throat_population_classes_included": ["confined_fast", "directed_outbound_Eq63", "configured_background_positive_charge"], "modal_phi_z_exact_throat_population_classes_excluded": ["material_returning_fast_population_requires_fixed_boundary_failure"], "modal_local_speed_grid_diagnostics": local_grid_diagnostics})
        metadata.update({
            "electron_midplane_density_m3": float(profile.electron_midplane_density_m3),
            "electron_midplane_density_used_by_Eq68_m3": float(current_balance.electron_balance.electron_midplane_density_m3),
            "electron_parent_maxwellian_n0_m3": float(profile.electron_parent_maxwellian_n0_m3),
            "electron_collision_density_m3": float(profile.electron_collision_density_m3),
            "electron_confined_volume_average_density_m3": float(profile.electron_volume_average_density_m3),
            "electron_density_profile_m3": [float(x) for x in np.asarray(profile.electron_density_m3, dtype=float)],
            "electron_profile_scope": "confined_throat_to_throat",
        })
    global_current_claimed = bool(solver_arguments.get("_global_current_balance_claimed", False))
    current_scope = str(solver_arguments.get("_current_balance_scope", "fast_D_fast_T_electron_FBIS_subsystem"))
    current_relative_error = float(current_balance.current_residual_A) / max(abs(float(current_balance.ion_current_loss.total_ion_current_A)), 1.0e-300)
    source_mean_energy_J = float(metadata.get("modal_source_mean_birth_energy_J", 0.0))
    barrier_over_birth = barrier / source_mean_energy_J if source_mean_energy_J > 0.0 else None
    fast_cross_pair_ids = _active_fast_cross_pair_ids(state)
    metadata.update({'fast_ion_operating_point_model': 'coupled_fast_D_fast_T_shared_electrostatic_state', 'active_fast_species': tuple(sorted(electrostatic.profile.fast_ion_density_by_species_m3)) if electrostatic is not None else (state.species.species_id,), 'shared_eq42_basis': True, 'shared_eq68_wall_barrier': True, 'shared_eq70_potential': electrostatic is not None, 'fast_D_fast_T_cross_collisions_available': bool(fast_cross_pair_ids), 'fast_D_fast_T_cross_collision_pair_ids': fast_cross_pair_ids, 'fast_D_fast_T_cross_collision_model': None if not fast_cross_pair_ids else 'reduced_isotropic_first_mode_Rosenbluth_cross_species', 'global_current_balance_claimed': global_current_claimed, 'total_current_balance_candidate_active': False, 'physical_electron_energy_equation_active': False, 'fast_T_fusion_channels_available': True, 'current_balance_scope': current_scope, 'modal_ion_barrier_power_loss_W': ion_barrier_power_W, 'modal_ion_wall_power_loss_W': ion_wall_power_W, 'modal_total_wall_power_loss_W': None, 'modal_per_lost_ion_barrier_energy_over_birth_energy': barrier_over_birth, 'modal_per_lost_ion_total_wall_energy_over_birth_energy': None, 'modal_electron_wall_power_loss_W': None, 'modal_shared_electron_wall_power_loss_W': float(current_balance.electron_balance.electron_wall_power_loss_W), 'modal_electron_current_balance_residual_A': float(current_balance.current_residual_A), 'modal_electron_current_balance_relative_error': current_relative_error, 'modal_wall_barrier_energy_J': barrier, 'modal_wall_barrier_over_Te': float(current_balance.electron_balance.normalized_barrier), 'modal_electron_wall_potential_relative_to_midplane_V': float(current_balance.electron_wall_potential_relative_to_midplane_V), 'modal_normalized_wall_barrier': float(current_balance.electron_balance.normalized_barrier), 'modal_current_balance_relative_error': current_relative_error, 'fixed_boundary_lost_reconstruction_status': fixed_boundary_status, 'fixed_boundary_lost_reconstruction_complete': fixed_boundary_status == 'represented'})
    updated_modal = replace(modal, pitch_grid=solver_arguments.get("pitch_grid"), ion_wall_power_loss_W=ion_wall_power_W, electron_wall_power_loss_W=None, wall_barrier_energy_J=barrier, normalized_wall_barrier=float(current_balance.electron_balance.normalized_barrier), electron_wall_potential_relative_to_midplane_V=float(current_balance.electron_wall_potential_relative_to_midplane_V), electrostatic_profile=None if electrostatic is None else electrostatic.profile, local_reconstruction=local_reconstruction, metadata=metadata, current_balance=current_balance)
    solved_collision_state = rebuild_fbis_collision_parameter_state_with_fast_density(state.collision_state, updated_modal.density_m3)
   
    return replace(state, modal_result=updated_modal, collision_state=solved_collision_state)

def _solve_mixed_system_once(*, requests: Mapping[str, FastIonSpeciesRequest], collision_states: Mapping[str, FBISCollisionParameterState], prompt_only_loss_rates_s: Mapping[str, float], prompt_only_midplane_power_W: Mapping[str, float], ion_loss_rates_s_by_component_override: Mapping[str, float] | None, solver_arguments: Mapping[str, Any], throat_drop_override_J: float | None, published_reference_lambda1: float | None, initial_warm_starts: Mapping[str, object] | None, warm_start_compatibility_policies: Mapping[str, str] | None, initial_potential_energy_J: np.ndarray | None, run_numerical_convergence_assessments: bool, runtime_mixed_pass_kind: str) -> tuple[dict[str, FastIonSpeciesState], AmbipolarCurrentBalance | None, ModalElectrostaticSystemResult | None]:
    """Solve one coupled species pass with fixed collision states and reconcile the Eq 68 and Eq 70 references"""
    runtime_profile = solver_arguments.get("_runtime_profile_accumulator")
    preliminary_states: dict[str, FastIonSpeciesState] = {}
    shared_basis = solver_arguments.get("basis")
    for species_id in sorted(requests):
        request = requests[species_id]
        source_rate = float(request.attenuated_source.total_birth_rate_s)
        prompt_rate = float(prompt_only_loss_rates_s.get(species_id, 0.0))
        if prompt_rate > 0.0:
            preliminary_states[species_id] = FastIonSpeciesState(species=request.species, speed_grid=request.speed_grid, active=False, source_particle_rate_s=source_rate, modal_result=None, prompt_only_particle_loss_rate_s=prompt_rate, prompt_only_midplane_kinetic_power_loss_W=float(prompt_only_midplane_power_W.get(species_id, 0.0)), status="prompt_only_no_confined_source")
            continue
        if source_rate <= 0.0:
            preliminary_states[species_id] = _solve_fast_ion_species_direct(request=request, collision_state=None, solver_arguments={})
            continue
        local_arguments = dict(solver_arguments)
        local_arguments["basis"] = shared_basis
        local_arguments["pitch_grid"] = None
        local_arguments["cell_volumes_m3"] = None
        local_arguments["electrostatic_feedback_model"] = MAGNETIC_ONLY_FBIS if throat_drop_override_J is None else EGEDAL_2022_PUBLISHED_ELECTROSTATIC_FBIS_APPROXIMATION
        local_arguments["_published_throat_drop_override_J"] = throat_drop_override_J
        local_arguments["_published_reference_lambda1"] = published_reference_lambda1
        local_arguments["_eq59_warm_start_state"] = None if initial_warm_starts is None else initial_warm_starts.get(species_id)
        default_policy = "exact_physics" if local_arguments["_eq59_warm_start_state"] is None else "closure_continuation"
        local_arguments["_eq59_warm_start_compatibility_policy"] = default_policy if warm_start_compatibility_policies is None else warm_start_compatibility_policies.get(species_id, default_policy)
        local_arguments["_eq70_initial_potential_energy_J"] = None
        local_arguments["_run_numerical_convergence_assessments"] = bool(run_numerical_convergence_assessments)
        species_runtime_started = perf_counter()
        state = _solve_fast_ion_species_direct(request=request, collision_state=collision_states[species_id], solver_arguments=local_arguments)
        record_runtime(runtime_profile, f"fast_species_preliminary_{species_id}", perf_counter() - species_runtime_started, nested=True)
        preliminary_states[species_id] = state
        if state.modal_result is not None:
            shared_basis = state.modal_result.basis
    active = tuple(state for state in preliminary_states.values() if state.modal_result is not None)
    if not active:
        return preliminary_states, None, None
    electron_frequencies = np.asarray([float(collision_states[state.species.species_id].electron_electron_collision_frequency_s) for state in active], dtype=float)
    if not np.allclose(electron_frequencies, electron_frequencies[0], rtol=4.0e-15, atol=0.0):
        raise ValueError("active fast species input collision states must use one electron collision frequency")
    # The first Eq 68 barrier uses the supplied midplane electron density before Eq 70 determines its parent Maxwellian n0
    current_runtime_started = perf_counter()
    current_balance = _shared_current_balance(states=preliminary_states, electron_midplane_density_m3=float(solver_arguments['electron_midplane_density_m3']), electron_collision_density_m3=float(solver_arguments['electron_collision_density_m3']), electron_temperature_J=float(solver_arguments['electron_temperature_J']), electron_collision_frequency_s=float(electron_frequencies[0]), midplane_area_m2=float(solver_arguments['midplane_area_m2']), half_length_m=float(solver_arguments['half_length_m']), mirror_ratio=float(solver_arguments['mirror_ratio']), basis=shared_basis, ion_loss_rates_s_by_component_override=ion_loss_rates_s_by_component_override)
    record_runtime(runtime_profile, "eq68_shared_current_balance", perf_counter() - current_runtime_started, nested=True)
    # Directed Eq 63 branches are attached before Eq 70 so their charge can enter the shared quasineutrality closure
    baldwin_runtime_started = perf_counter()
    preliminary_states = _attach_baldwin_eq63_branches(states=preliminary_states, solver_arguments=solver_arguments)
    record_runtime(runtime_profile, "baldwin_eq63_throat_closure", perf_counter() - baldwin_runtime_started, nested=True)
    numerics = solver_arguments["numerics"]
    coupling_limit = max(int(numerics.phi_z_iterations), 1)
    coupling_tolerance = float(numerics.phi_z_relative_tolerance)
    electrostatic = None
    coupling_history: list[dict[str, object]] = []
    potential_seed = initial_potential_energy_J
    coupling_iteration_converged = False
    coupling_converged = False
    # Eq 70 changes n0 and the exact midplane reference so Eq 68 and Eq 70 are iterated to a common scalar state
    for iteration_index in range(1, coupling_limit + 1):
        eq70_runtime_started = perf_counter()
        electrostatic = _shared_electrostatic_state(states=preliminary_states, current_balance=current_balance, solver_arguments=solver_arguments, initial_potential_energy_J=potential_seed, audit_context={"fast_ion_system_solve_index": solver_arguments.get("_runtime_eq70_fast_system_solve_index"), "solve_context": solver_arguments.get("_runtime_eq70_solve_context"), "mixed_pass_kind": runtime_mixed_pass_kind, "coupling_phase": f"floating_iteration_{iteration_index}"})
        record_runtime(runtime_profile, "eq70_shared_profile", perf_counter() - eq70_runtime_started, nested=True)
        if electrostatic is None:
            coupling_converged = True
            break
        floating_runtime_started = perf_counter()
        updated_balance = _current_balance_from_floating_profile(current_balance=current_balance, electrostatic=electrostatic, solver_arguments=solver_arguments)
        record_runtime(runtime_profile, "eq68_floating_reference_update", perf_counter() - floating_runtime_started, nested=True)
        old_wall = float(current_balance.electron_balance.barrier_energy_J)
        new_wall = float(updated_balance.electron_balance.barrier_energy_J)
        wall_change = abs(new_wall - old_wall) / max(abs(new_wall), abs(old_wall), float(solver_arguments["electron_temperature_J"]), np.finfo(float).tiny)
        current_change = abs(float(updated_balance.current_residual_A)) / max(abs(float(updated_balance.ion_current_loss.total_ion_current_A)), 1.0)
        coupling_iteration_converged = bool(electrostatic.profile.floating_scalar_closure_converged and wall_change <= coupling_tolerance)
        coupling_history.append({
            "iteration": int(iteration_index),
            "wall_barrier_energy_J_used_for_Eq70": old_wall,
            "wall_barrier_energy_J_recomputed_from_parent_n0": new_wall,
            "wall_barrier_relative_change": float(wall_change),
            "electron_parent_maxwellian_n0_m3": float(electrostatic.profile.electron_parent_maxwellian_n0_m3),
            "exact_midplane_potential_energy_J": float(electrostatic.profile.exact_midplane_potential_energy_J),
            "Eq68_current_relative_residual": float(current_change),
            "Eq70_converged": bool(electrostatic.profile.converged),
            "Eq70_floating_scalar_closure_converged": bool(electrostatic.profile.floating_scalar_closure_converged),
            "iteration_fixed_point_converged": bool(coupling_iteration_converged),
        })
        current_balance = updated_balance
        potential_seed = np.asarray(electrostatic.profile.potential_energy_J, dtype=float)
        if coupling_iteration_converged:
            break
    if electrostatic is not None:
        final_eq70_runtime_started = perf_counter()
        final_electrostatic = _shared_electrostatic_state(states=preliminary_states, current_balance=current_balance, solver_arguments=solver_arguments, initial_potential_energy_J=potential_seed, audit_context={"fast_ion_system_solve_index": solver_arguments.get("_runtime_eq70_fast_system_solve_index"), "solve_context": solver_arguments.get("_runtime_eq70_solve_context"), "mixed_pass_kind": runtime_mixed_pass_kind, "coupling_phase": "final_reconciliation"})
        record_runtime(runtime_profile, "eq70_shared_profile", perf_counter() - final_eq70_runtime_started, nested=True)
        if final_electrostatic is not None:
            electrostatic = final_electrostatic
            final_floating_runtime_started = perf_counter()
            final_balance = _current_balance_from_floating_profile(current_balance=current_balance, electrostatic=electrostatic, solver_arguments=solver_arguments)
            record_runtime(runtime_profile, "eq68_floating_reference_update", perf_counter() - final_floating_runtime_started, nested=True)
            used_wall = float(current_balance.electron_balance.barrier_energy_J)
            reconciled_wall = float(final_balance.electron_balance.barrier_energy_J)
            final_wall_change = abs(reconciled_wall - used_wall) / max(abs(reconciled_wall), abs(used_wall), float(solver_arguments["electron_temperature_J"]), np.finfo(float).tiny)
            current_balance = final_balance
            coupling_iteration_converged = bool(electrostatic.profile.floating_scalar_closure_converged and final_wall_change <= coupling_tolerance)
            coupling_converged = bool(electrostatic.profile.converged and coupling_iteration_converged)
            coupling_history.append({
                "iteration": "final_reconciliation",
                "wall_barrier_energy_J_used_for_Eq70": used_wall,
                "wall_barrier_energy_J_recomputed_from_parent_n0": reconciled_wall,
                "wall_barrier_relative_change": float(final_wall_change),
                "electron_parent_maxwellian_n0_m3": float(electrostatic.profile.electron_parent_maxwellian_n0_m3),
                "exact_midplane_potential_energy_J": float(electrostatic.profile.exact_midplane_potential_energy_J),
                "Eq68_current_relative_residual": float(abs(float(final_balance.current_residual_A)) / max(abs(float(final_balance.ion_current_loss.total_ion_current_A)), 1.0)),
                "Eq70_converged": bool(electrostatic.profile.converged),
                "Eq70_floating_scalar_closure_converged": bool(electrostatic.profile.floating_scalar_closure_converged),
                "iteration_fixed_point_converged": bool(coupling_iteration_converged),
            })
            if not coupling_iteration_converged:
                failure = electrostatic.profile.failure_reason
                reason = "eq68_floating_reference_reconciliation_not_converged"
                failure = reason if not failure else f"{failure};{reason}"
                electrostatic = replace(electrostatic, profile=replace(electrostatic.profile, converged=False, failure_reason=failure))
    final_arguments = dict(solver_arguments)
    final_arguments["_run_numerical_convergence_assessments"] = bool(run_numerical_convergence_assessments)
    finalized_states = {}
    for key, state in preliminary_states.items():
        finalize_runtime_started = perf_counter()
        finalized_states[key] = _finalize_species_result(state=state, current_balance=current_balance, electrostatic=electrostatic, solver_arguments=final_arguments)
        record_runtime(runtime_profile, f"species_finalize_{key}", perf_counter() - finalize_runtime_started, nested=True)
    if coupling_history:
        finalized_states = {
            key: replace(
                state,
                modal_result=None if state.modal_result is None else replace(
                    state.modal_result,
                    metadata={
                        **state.modal_result.metadata,
                        "modal_eq68_eq70_floating_reference_coupling_history": tuple(coupling_history),
                        "modal_eq68_eq70_floating_reference_iteration_converged": bool(coupling_iteration_converged),
                        "modal_eq68_eq70_floating_reference_coupling_converged": bool(coupling_converged),
                    },
                ),
            )
            for key, state in finalized_states.items()
        }
  
    return finalized_states, current_balance, electrostatic

def solve_fast_ion_system(*, requests: Mapping[str, FastIonSpeciesRequest], collision_states_by_species: Mapping[str, FBISCollisionParameterState], eq59_warm_start_states_by_species: Mapping[str, object] | None=None, eq59_warm_start_compatibility_policy_by_species: Mapping[str, str] | None=None, prompt_only_loss_rates_s_by_species: Mapping[str, float] | None=None, prompt_only_midplane_power_W_by_species: Mapping[str, float] | None=None, ion_loss_rates_s_by_component_override: Mapping[str, float] | None=None, **solver_arguments: Any) -> FastIonSystemState:
    """Solve active fast species with one Eq 42 basis, Eq 68 barrier, Eq 70 potential, and optional feedback iteration"""
    request_map = dict(requests)
    collision_map = dict(collision_states_by_species)
    warm_start_policies = {str(key): str(value) for key, value in dict(eq59_warm_start_compatibility_policy_by_species or {}).items()}
    prompt_rates = {str(key): float(value) for key, value in dict(prompt_only_loss_rates_s_by_species or {}).items()}
    prompt_powers = {str(key): float(value) for key, value in dict(prompt_only_midplane_power_W_by_species or {}).items()}
    current_override = None if ion_loss_rates_s_by_component_override is None else {str(key): float(value) for key, value in dict(ion_loss_rates_s_by_component_override).items()}
    if not request_map:
        return FastIonSystemState(species_states={}, shared_basis=solver_arguments.get("basis"), status="inactive_no_species", metadata={"fast_ion_operating_point_model": "inactive_no_species"})
    for key, request in request_map.items():
        if str(key) != request.species.species_id:
            raise ValueError("fast ion system request keys must match species identifiers")
        prompt_rate = float(prompt_rates.get(key, 0.0))
        prompt_power = float(prompt_powers.get(key, 0.0))
        source_rate = float(request.attenuated_source.total_birth_rate_s)
        collision_state = collision_map.get(key)
        if collision_state is not None and not isinstance(collision_state, FBISCollisionParameterState):
            raise TypeError("collision states must be FBISCollisionParameterState values")
        policy = warm_start_policies.get(key)
        if policy is not None and policy not in {"exact_physics", "closure_continuation"}:
            raise ValueError("Eq 59 warm start compatibility policy must be exact_physics or closure_continuation")
        if not np.isfinite(prompt_rate) or prompt_rate < 0.0 or prompt_rate > source_rate + 128.0 * np.finfo(float).eps * max(source_rate, 1.0):
            raise ValueError("prompt only loss rate must be finite and no larger than the species source rate")
        if prompt_rate > 0.0 and not np.isclose(prompt_rate, source_rate, rtol=1.0e-12, atol=1.0e-12 * max(source_rate, 1.0)):
            raise ValueError("prompt only loss rate must represent the complete species source")
        if not np.isfinite(prompt_power) or prompt_power < 0.0:
            raise ValueError("prompt only loss power must be finite and nonnegative")
        if prompt_rate <= 0.0 and source_rate > 0.0:
            _validate_collision_state(request, collision_map.get(key))
    prompt_ids = tuple(sorted(key for key, rate in prompt_rates.items() if rate > 0.0))
    active_ids = tuple(sorted(key for key, request in request_map.items() if float(request.attenuated_source.total_birth_rate_s) > 0.0 and key not in prompt_ids))
    if not active_ids:
        states = {}
        for species_id in sorted(request_map):
            request = request_map[species_id]
            source_rate = float(request.attenuated_source.total_birth_rate_s)
            states[species_id] = FastIonSpeciesState(species=request.species, speed_grid=request.speed_grid, active=False, source_particle_rate_s=source_rate, modal_result=None, prompt_only_particle_loss_rate_s=float(prompt_rates.get(species_id, 0.0)), prompt_only_midplane_kinetic_power_loss_W=float(prompt_powers.get(species_id, 0.0)), status="prompt_only_no_confined_source" if float(prompt_rates.get(species_id, 0.0)) > 0.0 else "inactive_zero_source")
        return FastIonSystemState(species_states=states, shared_basis=solver_arguments.get("basis"), status="no_confined_fbis_source", metadata={"fast_ion_operating_point_model": "no_confined_fbis_source", "source_fast_species": tuple(sorted(key for key, state in states.items() if state.source_particle_rate_s > 0.0)), "active_fast_species": ()})
    feedback_model = canonical_electrostatic_feedback_model(solver_arguments.get("electrostatic_feedback_model", MAGNETIC_ONLY_FBIS))
    if feedback_model == ZERO_POTENTIAL_MAGNETIC_REFERENCE:
        raise ValueError("zero_potential_magnetic_reference is not a coupled fast D and fast T operating point")
    working_arguments = dict(solver_arguments)
    runtime_profile = working_arguments.get("_runtime_profile_accumulator")
    if bool(working_arguments.get("_runtime_eq70_confined_interpolation_audit_enabled", False)) and isinstance(runtime_profile, dict):
        audit_private = runtime_profile.setdefault("_eq70_confined_interpolation_audit_private", {})
        if not isinstance(audit_private, dict):
            raise TypeError("Eq 70 confined interpolation audit private state must be a mapping")
        fast_system_index = int(audit_private.get("fast_ion_system_solve_count", 0)) + 1
        audit_private["fast_ion_system_solve_count"] = fast_system_index
        working_arguments["_runtime_eq70_fast_system_solve_index"] = fast_system_index
    initial_potential = working_arguments.get("_eq70_initial_potential_energy_J")
    run_assessments = bool(working_arguments.get("_run_numerical_convergence_assessments", True))
    reference_lambda1 = None
    throat_drop = None
    feedback_history: list[float] = []
    feedback_converged = feedback_model == MAGNETIC_ONLY_FBIS
    warm_starts = dict(eq59_warm_start_states_by_species or {})
    if feedback_model == EGEDAL_2022_PUBLISHED_ELECTROSTATIC_FBIS_APPROXIMATION:
        basis = working_arguments.get("basis")
        if basis is None:
            basis = build_modal_fbis_basis(mirror_ratio=float(working_arguments["mirror_ratio"]), B_tilde_function=working_arguments["B_tilde_function"], zeta_faces=working_arguments["zeta_faces"], B_tilde_midpoints=working_arguments["B_tilde_midpoints"], numerics=working_arguments["numerics"])
            working_arguments["basis"] = basis
        reference_profile = rescaled_magnetic_profile_for_reference_ratio(working_arguments["B_tilde_function"], float(working_arguments["mirror_ratio"]))
        active_midpoints = np.asarray(working_arguments["B_tilde_midpoints"], dtype=float)
        reference_midpoints = 1.0 + 0.5 * (active_midpoints - 1.0) / (float(working_arguments["mirror_ratio"]) - 1.0)
        reference_basis = build_modal_fbis_basis(mirror_ratio=1.5, B_tilde_function=reference_profile, zeta_faces=working_arguments["zeta_faces"], B_tilde_midpoints=reference_midpoints, numerics=working_arguments["numerics"])
        reference_lambda1 = float(reference_basis.physical_basis.eigenvalues[0])
    states, current_balance, electrostatic = _solve_mixed_system_once(requests=request_map, collision_states=collision_map, prompt_only_loss_rates_s=prompt_rates, prompt_only_midplane_power_W=prompt_powers, ion_loss_rates_s_by_component_override=current_override, solver_arguments=working_arguments, throat_drop_override_J=None, published_reference_lambda1=reference_lambda1, initial_warm_starts=warm_starts, warm_start_compatibility_policies=warm_start_policies, initial_potential_energy_J=initial_potential, run_numerical_convergence_assessments=False if feedback_model == EGEDAL_2022_PUBLISHED_ELECTROSTATIC_FBIS_APPROXIMATION else run_assessments, runtime_mixed_pass_kind='initial')
    if current_balance is None:
        return FastIonSystemState(species_states=states, shared_basis=working_arguments.get("basis"), status="inactive_zero_source")
    if electrostatic is not None:
        initial_potential = np.asarray(electrostatic.profile.potential_energy_J, dtype=float)
        throat_drop = float(electrostatic.profile.throat_potential_energy_J)
    if feedback_model == EGEDAL_2022_PUBLISHED_ELECTROSTATIC_FBIS_APPROXIMATION:
        limit = max(int(working_arguments["numerics"].electrostatic_feedback_iterations), 1)
        tolerance = float(working_arguments["numerics"].electrostatic_feedback_relative_tolerance)
        for feedback_iteration_index in range(1, limit + 1):
            warm_starts = {key: state.modal_result.eq59_warm_start_state for key, state in states.items() if state.modal_result is not None}
            previous = throat_drop
            states, current_balance, electrostatic = _solve_mixed_system_once(requests=request_map, collision_states=collision_map, prompt_only_loss_rates_s=prompt_rates, prompt_only_midplane_power_W=prompt_powers, ion_loss_rates_s_by_component_override=current_override, solver_arguments=working_arguments, throat_drop_override_J=throat_drop, published_reference_lambda1=reference_lambda1, initial_warm_starts=warm_starts, warm_start_compatibility_policies={key: 'closure_continuation' for key in warm_starts}, initial_potential_energy_J=initial_potential, run_numerical_convergence_assessments=False, runtime_mixed_pass_kind=f'electrostatic_feedback_iteration_{feedback_iteration_index}')
            if electrostatic is None:
                break
            initial_potential = np.asarray(electrostatic.profile.potential_energy_J, dtype=float)
            throat_drop = float(electrostatic.profile.throat_potential_energy_J)
            feedback_history.append(throat_drop)
            change = abs(throat_drop - previous) / max(abs(throat_drop), abs(previous), np.finfo(float).tiny)
            if change <= tolerance:
                feedback_converged = True
                break
        if run_assessments:
            warm_starts = {key: state.modal_result.eq59_warm_start_state for key, state in states.items() if state.modal_result is not None}
            states, current_balance, electrostatic = _solve_mixed_system_once(requests=request_map, collision_states=collision_map, prompt_only_loss_rates_s=prompt_rates, prompt_only_midplane_power_W=prompt_powers, ion_loss_rates_s_by_component_override=current_override, solver_arguments=working_arguments, throat_drop_override_J=throat_drop, published_reference_lambda1=reference_lambda1, initial_warm_starts=warm_starts, warm_start_compatibility_policies={key: 'closure_continuation' for key in warm_starts}, initial_potential_energy_J=initial_potential, run_numerical_convergence_assessments=True, runtime_mixed_pass_kind='electrostatic_feedback_final_assessment')
    global_current_claimed = False
    if current_override is not None:
        current_scope = "explicit_component_ion_current_override_candidate"
    else:
        current_scope = 'fast_D_fast_T_electron_FBIS_subsystem'
    working_arguments["_global_current_balance_claimed"] = global_current_claimed
    working_arguments["_current_balance_scope"] = current_scope
    runtime_profile = working_arguments.get("_runtime_profile_accumulator")
    finalized = {}
    for key, state in states.items():
        finalize_runtime_started = perf_counter()
        finalized[key] = _finalize_species_result(state=state, current_balance=current_balance, electrostatic=electrostatic, solver_arguments=working_arguments)
        record_runtime(runtime_profile, f"system_finalize_{key}", perf_counter() - finalize_runtime_started, nested=True)
    shared_basis = next(state.modal_result.basis for state in finalized.values() if state.modal_result is not None)
    electron_balance = current_balance.electron_balance
    electron_power = float(electron_balance.electron_wall_power_loss_W)
    fast_cross_pair_ids = tuple(sorted({pair_id for state in finalized.values() for pair_id in _active_fast_cross_pair_ids(state)}))
    fast_cross_available = len(active_ids) < 2 or all(bool(_active_fast_cross_pair_ids(finalized[species_id])) for species_id in active_ids)
    metadata = {'fast_ion_operating_point_model': 'coupled_fast_D_fast_T_shared_electrostatic_state', 'active_fast_species': active_ids, 'source_fast_species': tuple(sorted((key for key, request in request_map.items() if float(request.attenuated_source.total_birth_rate_s) > 0.0))), 'shared_eq42_basis': True, 'shared_eq68_wall_barrier': True, 'shared_eq70_potential': electrostatic is not None, 'fast_D_fast_T_cross_collisions_available': bool(fast_cross_available), 'fast_D_fast_T_cross_collision_pair_ids': fast_cross_pair_ids, 'fast_D_fast_T_cross_collision_model': None if not fast_cross_pair_ids else 'reduced_isotropic_first_mode_Rosenbluth_cross_species', 'global_current_balance_claimed': global_current_claimed, 'total_current_balance_candidate_active': False, 'physical_electron_energy_equation_active': False, 'fast_T_fusion_channels_available': True, 'current_balance_scope': current_scope, 'ion_loss_rates_s_by_component_override': current_override, 'electrostatic_feedback_model': feedback_model, 'electrostatic_feedback_iteration_history_throat_energy_J': feedback_history, 'electrostatic_feedback_iteration_converged': feedback_converged, 'system_total_ion_loss_current_A': float(current_balance.ion_current_loss.total_ion_current_A), 'system_electron_loss_current_A': float(electron_balance.electron_current_A), 'system_current_residual_A': float(current_balance.current_residual_A), 'system_electron_wall_power_loss_W': electron_power, 'electron_tail_integral_production_model': 'exact_integrals_preceding_Egedal_Eq68_Eq69', 'electron_tail_integral_asymptotic_role': 'diagnostic_only', 'eq68_exact_tail_kernel': None if electron_balance.exact_tail_kernel is None else float(electron_balance.exact_tail_kernel), 'eq68_asymptotic_tail_kernel': None if electron_balance.eq68_asymptotic_tail_kernel is None else float(electron_balance.eq68_asymptotic_tail_kernel), 'eq68_asymptotic_relative_correction': None if electron_balance.eq68_asymptotic_relative_correction is None else float(electron_balance.eq68_asymptotic_relative_correction), 'eq69_asymptotic_relative_correction': None if electron_balance.eq69_asymptotic_relative_correction is None else float(electron_balance.eq69_asymptotic_relative_correction), 'electron_mean_wall_kinetic_energy_J': None if electron_balance.mean_wall_kinetic_energy_J is None else float(electron_balance.mean_wall_kinetic_energy_J)}
   
    return FastIonSystemState(species_states=finalized, shared_basis=shared_basis, shared_current_balance=current_balance, shared_electrostatic_profile=None if electrostatic is None else electrostatic.profile, electron_wall_power_loss_W=electron_power, status="active_coupled_species", metadata=metadata)
