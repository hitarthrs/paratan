"""
Electron temperature closure around the fixed temperature operating point and expander state

Fixed temperature mode evaluates the represented energy residual as a diagnostic
Self consistent mode searches for a temperature whose closed beam, kinetic, terminal current, expander, and electron energy state satisfies the configured residual tolerance
"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import replace
from math import isfinite
from time import perf_counter
from typing import Any, Callable
import numpy as np
from source_model_revamp.coupling.electron_temperature_power_balance import ElectronTemperatureFinalQualificationComparison, ElectronTemperatureSolveConfig, ElectronTemperatureTrialEvaluation, ElectronTemperatureTrialRejectedError, solve_electron_temperature_power_balance
from source_model_revamp.coupling.electron_temperature_preflight import ElectronTemperaturePreflightDisposition, evaluate_self_consistent_temperature_preflight
from source_model_revamp.fbis.collision_parameters import energy_J_from_keV
from source_model_revamp.fbis.modal.types import ModalPhysicalDistributionError
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.evaluation import OperatingPointEvaluationMode, OperatingPointProgressCallback, OperatingPointProgressEvent, emit_progress
from source_model_revamp.integration.modal_stage import build_reusable_modal_basis
from source_model_revamp.integration.assessment import kinetic_numerical_convergence
from source_model_revamp.integration.pipeline_stages.electron_energy_balance_stage import build_electron_energy_balance_stage
from source_model_revamp.integration.pipeline_stages.expander_stage import build_expander_stage
from source_model_revamp.integration.pipeline_stages.operating_point_stage import OperatingPointWarmStartIncompatibleError, build_fixed_temperature_operating_point
from source_model_revamp.integration.pipeline_stages.ion_current_helpers import mark_pre_expander_current_candidate
from source_model_revamp.integration.pipeline_types import GeometryStageResult, OperatingPointStageResult, OperatingPointWarmStartState
from source_model_revamp.pipeline_errors import InvalidPipelineStateError
from source_model_revamp.runtime_profile import finalize_runtime_profile, new_runtime_profile, record_runtime

def _electron_temperature_solve_config(config: SourceModelRunConfig) -> ElectronTemperatureSolveConfig:
    """Build scalar electron temperature root search controls from the power balance configuration"""
    pb = config.power_balance
   
    return ElectronTemperatureSolveConfig(bracket_min_keV=pb.bracket_min_keV, bracket_max_keV=pb.bracket_max_keV, bracket_scan_points=pb.bracket_scan_points, max_iterations=pb.max_iterations, residual_tolerance_W=pb.residual_tolerance_W, relative_residual_tolerance=pb.relative_residual_tolerance, root_temperature_relative_tolerance=pb.root_temperature_relative_tolerance)

def _temperature_scope_metadata() -> dict[str, object]:
    """Return metadata defining the represented electron energy balance scope and its non global interpretation"""
    return {'electron_temperature_power_balance_scope': 'classical_represented_confined_plasma_electron_energy_balance', 'global_plasma_power_balance_available': False}

def _with_operating_point_metadata(operating_point: OperatingPointStageResult, metadata: Mapping[str, Any]) -> OperatingPointStageResult:
    """Attach temperature closure metadata to the operating point, kinetic state, and expander stage"""
    combined = {**operating_point.metadata, **metadata}
    operating_point.kinetic.metadata.update(metadata)
    expander = operating_point.expander
    if expander is not None:
        expander_system = expander.system_state
        if expander_system is not None:
            expander_system = replace(expander_system, metadata={**expander_system.metadata, **metadata})
        expander = replace(expander, system_state=expander_system, metadata={**expander.metadata, **metadata})
    assignment = combined.get("electron_temperature_assignment_valid") is True
    balance_value = combined.get("electron_temperature_power_balance_converged")
    balance = bool(balance_value) if isinstance(balance_value, bool) else None
   
    return replace(operating_point, metadata=combined, expander=expander, electron_temperature_assignment_valid=assignment, electron_temperature_power_balance_converged=balance)

def _valid_warm_start(state: OperatingPointWarmStartState | None, geometry: GeometryStageResult) -> OperatingPointWarmStartState | None:
    """Return a warm state only when its confined grid matches the current geometry"""
    if state is None:
        return None
    target = np.asarray(state.target_density_profile_m3, dtype=float)
    mask = np.asarray(state.target_density_defined_mask, dtype=bool)
    expected = geometry.zeta_centers.shape
    if target.shape != expected or mask.shape != expected:
        return None
    if np.any(~np.isfinite(target)) or np.any(target < 0.0) or not np.any(mask):
        return None
  
    return state

def _warm_start_from_current_closed_state(operating_point: OperatingPointStageResult) -> OperatingPointWarmStartState | None:
    """Build a continuation warm state from the current closed operating point when one is available"""
    kinetic = operating_point.kinetic
    density = kinetic.operating_point_density_state
    if density is None:
        return None
  
    return OperatingPointWarmStartState(
        target_density_profile_m3=np.asarray(density.electron_cell_density_m3, dtype=float).copy(),
        target_density_defined_mask=np.asarray(density.profile_support_mask, dtype=bool).copy(),
        kinetic_target_state=kinetic,
        eq59_warm_start_state=kinetic.eq59_warm_start_state,
        density_state=density,
        modal_basis_warm_start=kinetic.modal_basis_warm_start,
        eq70_potential_energy_warm_start_J=kinetic.eq70_potential_energy_warm_start_J,
        eq42_density_profile_warm_start=kinetic.eq42_density_profile_warm_start,
        eq59_warm_start_state_by_species=kinetic.eq59_warm_start_state_by_species,
        scalar_collision_operator_warm_start_by_species=kinetic.scalar_collision_operator_warm_start_by_species,
        scalar_external_fast_field_warm_start_by_test_species=kinetic.scalar_external_fast_field_warm_start_by_test_species,
        scalar_wall_barrier_energy_warm_start_J=kinetic.scalar_wall_barrier_energy_warm_start_J,
    )

def _compact_trial_record(temperature_keV: float, *, valid: bool, residual_W: float | None = None, relative_residual: float | None = None, failure_reason: str | None = None, warm_start_used: bool = False) -> dict[str, object]:
    """Return the compact temperature trial record used in scalar search metadata"""
    return {
        "electron_temperature_keV": float(temperature_keV),
        "valid": bool(valid),
        "residual_W": residual_W,
        "relative_residual": relative_residual,
        "failure_reason": failure_reason,
        "warm_start_used": bool(warm_start_used),
    }

def _final_kinetic_evidence_converged(metadata: Mapping[str, Any]) -> bool:
    """Return whether the final kinetic convergence evidence required by temperature qualification is present and passed"""
    return kinetic_numerical_convergence(metadata)[0]

def _expander_closure_failure_reason(config: SourceModelRunConfig, expander: Any) -> str:
    """Return a specific expander closure failure reason from terminal state and side potential convergence diagnostics"""
    system = expander.system_state
   
    if system is None:
        return "expander_state_unavailable"
    if expander.metadata.get("expander_potential_converged") is not True:
        reason = "expander_center_nodal_Eq70_potential_not_converged"
    else:
        reason = "expander_terminal_state_unavailable"
  
    details = []
    iteration_limit = int(config.kinetic_electrostatic.expander_potential_iterations)
    tolerance = float(config.kinetic_electrostatic.expander_potential_relative_tolerance)
  
    for side in ("left", "right"):
        state = system.potential_by_side.get(side)
        if state is None:
            continue
        details.append(f"{side}[quasineutrality_relative_error={state.quasineutrality_relative_error:.17g}, tolerance={tolerance:.17g}, iterations={state.iterations}/{iteration_limit}, potential_solver_failure_reason={state.potential_solver_failure_reason}]")
   
    return reason if not details else reason + "; " + "; ".join(details)

def _coupled_temperature_operating_point(config: SourceModelRunConfig, geometry: GeometryStageResult, *, electron_temperature_J: float, modal_basis: Any, warm_start_state: OperatingPointWarmStartState | None, collision_state_recompute_policy: str, collision_state_is_self_consistent_electron_temperature: bool, evaluation_mode: OperatingPointEvaluationMode, operating_point_builder: Callable[..., OperatingPointStageResult], progress_callback: OperatingPointProgressCallback | None, temperature_trial_index: int | None, clock: Callable[[], float]) -> OperatingPointStageResult:
    """
    Close the beam density fixed point, terminal current expander iteration, and represented electron energy balance at one trial temperature
    
    The result is admissible to the scalar temperature search only when the coupled operating point, current closure, expander state, and energy terms are all evaluable
    """
    coupled_runtime_started = clock()
    runtime_profile = new_runtime_profile()
    runtime_iterations: list[dict[str, object]] = []
    active_warm_start = _valid_warm_start(warm_start_state, geometry)
    history: list[dict[str, object]] = []
    last: OperatingPointStageResult | None = None
    limit = max(2, int(config.kinetic_electrostatic.total_current_balance_iterations))
   
    for iteration in range(1, limit + 1):
        iteration_runtime_started = clock()
        base_runtime_started = clock()
        base = operating_point_builder(
            config,
            geometry,
            electron_temperature_J=float(electron_temperature_J),
            modal_basis=modal_basis,
            warm_start_state=active_warm_start,
            collision_state_recompute_policy=collision_state_recompute_policy,
            collision_state_is_self_consistent_electron_temperature=collision_state_is_self_consistent_electron_temperature,
            evaluation_mode=evaluation_mode,
            progress_callback=progress_callback,
            temperature_trial_index=temperature_trial_index,
            clock=clock,
        )
        base_runtime_s = max(0.0, float(clock() - base_runtime_started))
        record_runtime(runtime_profile, "fixed_temperature_operating_point", base_runtime_s)

        kinetic = mark_pre_expander_current_candidate(config, base.kinetic)

        expander_runtime_started = clock()
        kinetic, expander = build_expander_stage(config, geometry, base.beam, kinetic, electron_temperature_J=float(electron_temperature_J), evaluation_mode=evaluation_mode)
        expander_runtime_s = max(0.0, float(clock() - expander_runtime_started))
        record_runtime(runtime_profile, "expander", expander_runtime_s)

        energy_runtime_started = clock()
        kinetic, energy = build_electron_energy_balance_stage(config, kinetic, electron_temperature_J=float(electron_temperature_J), expander=expander)
        energy_runtime_s = max(0.0, float(clock() - energy_runtime_started))
        record_runtime(runtime_profile, "electron_energy_balance", energy_runtime_s)

        beam_current_consistent = kinetic.metadata.get("total_current_balance_postclosure_beam_density_consistency_passed") is True
        current_numerical = kinetic.metadata.get("total_current_numerical_fixed_point_converged", kinetic.metadata.get("terminal_current_balance_numerical_fixed_point_converged")) is True
        expander_available = expander.system_state is not None
        expander_potential_converged = expander.metadata.get("expander_potential_converged") is True
        expander_center_nodal_eq70_active = expander.metadata.get("expander_center_nodal_potential_active") is True
        expander_terminal_current_authoritative = expander.metadata.get("terminal_current_balance_authoritative") is True
        expander_terminal_state_unavailable = bool(not expander_available or not expander_potential_converged)
        energy_available = energy is not None and energy.represented_terms_available
        energy_evaluable = energy is not None and energy.represented_terms_evaluable
        coupled_numerical = bool(base.beam_density_fixed_point_converged and current_numerical and beam_current_consistent and expander_available and expander_center_nodal_eq70_active and not expander_terminal_state_unavailable and expander_terminal_current_authoritative and energy_evaluable)
        history.append({
            "iteration": iteration,
            "beam_density_fixed_point_converged": base.beam_density_fixed_point_converged,
            "total_current_numerical_fixed_point_converged": current_numerical,
            "postclosure_beam_density_consistency_passed": beam_current_consistent,
            "expander_state_available": expander_available,
            "expander_center_nodal_Eq70_active": expander_center_nodal_eq70_active,
            "expander_potential_converged": expander_potential_converged,
            "expander_terminal_current_authoritative": expander_terminal_current_authoritative,
            "expander_terminal_state_unavailable": expander_terminal_state_unavailable,
            "represented_electron_energy_terms_available": energy_available,
            "represented_electron_energy_terms_evaluable": energy_evaluable,
            "coupled_iteration_converged": coupled_numerical,
        })
        iteration_runtime_s = max(0.0, float(clock() - iteration_runtime_started))
        runtime_iterations.append({'iteration': int(iteration), 'total_s': iteration_runtime_s, 'fixed_temperature_operating_point_s': base_runtime_s, 'expander_s': expander_runtime_s, 'electron_energy_balance_s': energy_runtime_s, 'base_runtime_profile': base.metadata.get('runtime_fixed_temperature_operating_point_profile'), 'kinetic_evaluations': base.metadata.get('runtime_kinetic_evaluations', ()), 'beam_density_fixed_point_converged': base.metadata.get('beam_density_fixed_point_converged'), 'beam_density_coupling_iteration_count': base.metadata.get('beam_density_coupling_iteration_count'), 'beam_density_coupling_iteration_limit': base.metadata.get('beam_density_coupling_iteration_limit'), 'beam_density_coupling_reached_iteration_limit': base.metadata.get('beam_density_coupling_reached_iteration_limit'), 'beam_density_coupling_exhausted_iteration_limit': base.metadata.get('beam_density_coupling_exhausted_iteration_limit'), 'beam_density_species_change_metrics_by_species': base.metadata.get('species_density_change_metrics_by_species'), 'beam_density_coupling_history': base.metadata.get('beam_density_coupling_history', ()) if config.output.write_iteration_histories else (), 'beam_density_final_consistency_attempt_count': base.metadata.get('beam_density_final_consistency_attempt_count'), 'beam_density_final_consistency_beam_metrics_passed': base.metadata.get('beam_density_final_consistency_beam_metrics_passed'), 'beam_density_final_consistency_stationary_particle_balance_required': base.metadata.get('beam_density_final_consistency_stationary_particle_balance_required'), 'beam_density_final_consistency_stationary_particle_balance_passed': base.metadata.get('beam_density_final_consistency_stationary_particle_balance_passed'), 'beam_density_final_consistency_overall_passed': base.metadata.get('beam_density_final_consistency_overall_passed'), 'beam_density_final_consistency_metric_gates': base.metadata.get('beam_density_final_consistency_metric_gates'), 'beam_density_final_consistency_stationary_particle_balance': base.metadata.get('beam_density_final_consistency_stationary_particle_balance'), 'beam_density_final_consistency_attempt_history': base.metadata.get('beam_density_final_consistency_attempt_history', ()) if config.output.write_iteration_histories else ()})
      
        if coupled_numerical:
            failure_reason = None
        elif expander_terminal_state_unavailable or not expander_center_nodal_eq70_active:
            failure_reason = _expander_closure_failure_reason(config, expander)
        elif not current_numerical or not expander_terminal_current_authoritative:
            failure_reason = "terminal_current_fixed_point_not_converged"
        elif not base.beam_density_fixed_point_converged or not beam_current_consistent:
            failure_reason = "beam_current_outer_iteration_not_converged"
        elif not energy_evaluable:
            if energy is not None and energy.failure_reason:
                failure_reason = f"electron_energy_balance_not_evaluable: {energy.failure_reason}"
            else:
                failure_reason = "electron_energy_balance_not_evaluable"
        else:
            failure_reason = "coupled_operating_point_not_converged"
       
        runtime_coupled_profile = finalize_runtime_profile(runtime_profile, total_s=clock() - coupled_runtime_started)
        metadata = {
            **base.metadata,
            **kinetic.metadata,
            
            **expander.metadata,
            "beam_current_outer_iteration_count": iteration,
            "beam_current_outer_iteration_limit": limit,
            "beam_current_outer_iteration_converged": coupled_numerical,
            "beam_current_outer_iteration_stopped_due_to_unavailable_expander_closure": expander_terminal_state_unavailable,
            "beam_current_outer_iteration_history": tuple(history),
            "runtime_coupled_temperature_operating_point_profile": runtime_coupled_profile,
            "runtime_operating_point_iterations": tuple(runtime_iterations),
        }
        last = replace(base, kinetic=kinetic, expander=expander, electron_energy_balance=energy, density_state=kinetic.operating_point_density_state, converged=coupled_numerical, failure_reason=failure_reason, warm_start_state=_warm_start_from_current_closed_state(replace(base, kinetic=kinetic)) or base.warm_start_state, metadata=metadata, state_finite=bool(base.state_finite and kinetic.metadata.get('state_finite', True)), state_physically_evaluable=bool(base.state_physically_evaluable and energy_evaluable), beam_density_fixed_point_converged=bool(base.beam_density_fixed_point_converged and beam_current_consistent), kinetic_numerical_convergence_passed=bool(kinetic.metadata.get('kinetic_numerical_convergence_passed')), end_loss_power_terms_available=bool(kinetic.metadata.get('modal_end_loss_power_terms_available')))
      
        if coupled_numerical:
            return last
        if expander_terminal_state_unavailable:
            break
        active_warm_start = _warm_start_from_current_closed_state(last)
        if active_warm_start is None:
            break
    if last is None:
        raise RuntimeError("coupled beam current temperature evaluation produced no state")
   
    return last

def _fixed_temperature_operating_point(config: SourceModelRunConfig, geometry: GeometryStageResult, *, modal_basis: Any, operating_point_builder: Callable[..., OperatingPointStageResult], progress_callback: OperatingPointProgressCallback | None, clock: Callable[[], float]) -> OperatingPointStageResult:
    """Build the configured fixed electron temperature operating point and retain the represented energy residual as a diagnostic"""
    temperature_keV = float(config.plasma_closure.electron_temperature_initial_guess_keV)
    temperature_J = float(energy_J_from_keV(temperature_keV))
    operating_point = _coupled_temperature_operating_point(
        config,
        geometry,
        electron_temperature_J=temperature_J,
        modal_basis=modal_basis,
        warm_start_state=None,
        collision_state_recompute_policy="computed_for_fixed_electron_temperature_operating_point",
        collision_state_is_self_consistent_electron_temperature=False,
        evaluation_mode=OperatingPointEvaluationMode.FIXED_TEMPERATURE_FINAL,
        operating_point_builder=operating_point_builder,
        progress_callback=progress_callback,
        temperature_trial_index=None,
        clock=clock,
    )
    energy = operating_point.electron_energy_balance
    residual_within_tolerance = bool(energy is not None and (abs(energy.residual_W) <= config.power_balance.residual_tolerance_W or abs(energy.relative_residual) <= config.power_balance.relative_residual_tolerance))
    metadata = {
        "electron_temperature_mode": "fixed_closure",
        "electron_temperature_assignment_valid": bool(isfinite(temperature_keV) and temperature_keV > 0.0),
        "electron_temperature_fixed_input_keV": temperature_keV,
        "electron_temperature_solved_keV": temperature_keV,
        "electron_temperature_diagnostic_root_keV": None,
        "electron_temperature_diagnostic_root_converged": False,
        "electron_temperature_solve_converged": None,
        "electron_temperature_power_balance_converged": None,
        "electron_temperature_power_residual_W": None if energy is None else energy.residual_W,
        "electron_temperature_relative_power_residual": None if energy is None else energy.relative_residual,
        "electron_temperature_power_balance_terms_available": energy is not None,
        "electron_energy_balance_residual_within_tolerance": residual_within_tolerance,
        "electron_temperature_final_qualification_passed": False,
        "electron_temperature_final_residual_consistent_with_trial": False,
        **_temperature_scope_metadata(),
    }
 
    return _with_operating_point_metadata(operating_point, metadata)

def _self_consistent_temperature_operating_point(config: SourceModelRunConfig, geometry: GeometryStageResult, *, modal_basis: Any, operating_point_builder: Callable[..., OperatingPointStageResult], progress_callback: OperatingPointProgressCallback | None, clock: Callable[[], float]) -> OperatingPointStageResult:
    """
    Solve for the self consistent electron temperature
    
    Each scalar trial evaluates the full coupled operating point and retries without warm state when continuation is incompatible
    After root convergence, a clean final qualification is rebuilt without warm continuation and compared with the accepted trial
    """
    preflight = evaluate_self_consistent_temperature_preflight(config)
    emit_progress(progress_callback, OperatingPointProgressEvent(phase="temperature_preflight_completed", evaluation_mode=OperatingPointEvaluationMode.TEMPERATURE_TRIAL, message=preflight.disposition.value))
    preflight_metadata = preflight.as_metadata()
    initial_temperature_keV = float(config.plasma_closure.electron_temperature_initial_guess_keV)
  
    if preflight.disposition is ElectronTemperaturePreflightDisposition.UNAVAILABLE:
        operating_point = _coupled_temperature_operating_point(
            config,
            geometry,
            electron_temperature_J=float(energy_J_from_keV(initial_temperature_keV)),
            modal_basis=modal_basis,
            warm_start_state=None,
            collision_state_recompute_policy="recomputed_for_diagnostic_initial_electron_temperature",
            collision_state_is_self_consistent_electron_temperature=True,
            evaluation_mode=OperatingPointEvaluationMode.TEMPERATURE_TRIAL,
            operating_point_builder=operating_point_builder,
            progress_callback=progress_callback,
            temperature_trial_index=1,
            clock=clock,
        )
        metadata = {
            "electron_temperature_mode": "self_consistent_electron_energy",
            "electron_temperature_assignment_valid": False,
            "electron_temperature_solved_keV": None,
            "electron_temperature_diagnostic_root_keV": None,
            "electron_temperature_diagnostic_root_converged": False,
            "electron_temperature_solve_converged": False,
            "electron_temperature_power_balance_converged": False,
            "electron_temperature_solve_failure_reason": "static_temperature_model_applicability_preflight_failed",
            "electron_temperature_search_skipped": True,
            "electron_temperature_search_skip_reason": "static_required_physics_unavailable",
            "electron_temperature_diagnostic_initial_guess_state_built": True,
            "electron_temperature_search_bracket_found": False,
            "electron_temperature_search_converged": False,
            "electron_temperature_intermediate_trial_base_physics_valid": False,
            "electron_temperature_intermediate_trial_residual_admitted": False,
            "electron_temperature_final_qualification_passed": False,
            "electron_temperature_final_residual_consistent_with_trial": False,
            "electron_temperature_power_balance_terms_available": False,
            **preflight_metadata,
            **_temperature_scope_metadata(),
        }
     
        return _with_operating_point_metadata(operating_point, metadata)

    solve_config = _electron_temperature_solve_config(config)
    trial_records: list[dict[str, object]] = []
    previous_warm_start: OperatingPointWarmStartState | None = None

    def evaluate(temperature_keV: float) -> ElectronTemperatureTrialEvaluation:
        """Evaluate one temperature trial, using compatible warm continuation when possible and returning only an admissible energy residual to the scalar solver"""
        nonlocal previous_warm_start
        warm = _valid_warm_start(previous_warm_start, geometry)
        used_warm = warm is not None
        try:
            operating_point = _coupled_temperature_operating_point(
                config,
                geometry,
                electron_temperature_J=float(energy_J_from_keV(temperature_keV)),
                modal_basis=modal_basis,
                warm_start_state=warm,
                collision_state_recompute_policy="recomputed_for_each_physical_electron_energy_temperature_trial",
                collision_state_is_self_consistent_electron_temperature=True,
                evaluation_mode=OperatingPointEvaluationMode.TEMPERATURE_TRIAL,
                operating_point_builder=operating_point_builder,
                progress_callback=progress_callback,
                temperature_trial_index=len(trial_records) + 1,
                clock=clock,
            )
        except (OperatingPointWarmStartIncompatibleError, ModalPhysicalDistributionError, ElectronTemperatureTrialRejectedError):
            if warm is None:
                raise
            operating_point = _coupled_temperature_operating_point(
                config,
                geometry,
                electron_temperature_J=float(energy_J_from_keV(temperature_keV)),
                modal_basis=modal_basis,
                warm_start_state=None,
                collision_state_recompute_policy="recomputed_for_each_physical_electron_energy_temperature_trial",
                collision_state_is_self_consistent_electron_temperature=True,
                evaluation_mode=OperatingPointEvaluationMode.TEMPERATURE_TRIAL,
                operating_point_builder=operating_point_builder,
                progress_callback=progress_callback,
                temperature_trial_index=len(trial_records) + 1,
                clock=clock,
            )
            used_warm = False
        energy = operating_point.electron_energy_balance
        current_numerical = operating_point.kinetic.metadata.get("total_current_numerical_fixed_point_converged", operating_point.kinetic.metadata.get("total_current_balance_numerical_fixed_point_converged")) is True
        residual_admitted = bool(operating_point.converged and current_numerical and energy is not None and energy.represented_terms_evaluable)
      
        if not residual_admitted:
            reasons = [
                operating_point.failure_reason,
                None if energy is not None else "represented_electron_energy_state_unavailable",
                None if energy is None or energy.represented_terms_evaluable else energy.failure_reason or "represented_electron_energy_terms_not_evaluable",
                None if current_numerical else "total_current_numerical_fixed_point_not_converged",
            ]
            reason = next((str(item) for item in reasons if item), "temperature_trial_residual_not_admitted")
            trial_records.append(_compact_trial_record(temperature_keV, valid=False, failure_reason=reason, warm_start_used=used_warm))
            raise ElectronTemperatureTrialRejectedError(f"temperature trial represented electron energy terms were not admitted; reason={reason}", rejection_kind=reason)
      
        trial_records.append(_compact_trial_record(temperature_keV, valid=True, residual_W=float(energy.residual_W), relative_residual=float(energy.relative_residual), warm_start_used=used_warm))
        previous_warm_start = operating_point.warm_start_state
        return ElectronTemperatureTrialEvaluation(electron_temperature_keV=float(temperature_keV), terms=energy, state=operating_point, evaluation_mode=OperatingPointEvaluationMode.TEMPERATURE_TRIAL)

    solve_result = solve_electron_temperature_power_balance(evaluate, solve_config, progress_callback=progress_callback, clock=clock)
   
    if not solve_result.converged:
        if solve_result.final_evaluation is None:
            attempted = ", ".join(f"{temperature:.6g}" for temperature, _ in solve_result.failed_trials)
            counts: dict[str, int] = {}
            for _, error in solve_result.failed_trials:
                kind = error.split(":", 1)[0]
                counts[kind] = counts.get(kind, 0) + 1
            failure_counts = ", ".join(f"{kind}={count}" for kind, count in counts.items())
            representative = solve_result.failed_trials[-1][1] if solve_result.failed_trials else "unavailable"
            raise InvalidPipelineStateError(
                "no physically valid electron temperature trial was found in configured bracket "
                f"[{solve_config.bracket_min_keV:.6g}, {solve_config.bracket_max_keV:.6g}] keV; "
                f"attempted_temperatures_keV=[{attempted}]; failure_counts=[{failure_counts}]; "
                f"representative_failure={representative}",
                stage="electron_temperature",
            )
        operating_point = solve_result.final_evaluation.state
        metadata = {
            "electron_temperature_mode": "self_consistent_electron_energy",
            "electron_temperature_assignment_valid": False,
            "electron_temperature_solved_keV": None,
            "electron_temperature_diagnostic_root_keV": solve_result.electron_temperature_keV,
            "electron_temperature_diagnostic_root_converged": False,
            "electron_temperature_best_diagnostic_trial_keV": solve_result.electron_temperature_keV,
            "electron_temperature_solve_converged": False,
            "electron_temperature_power_balance_converged": False,
            "electron_temperature_solve_failure_reason": solve_result.failure_reason,
            "electron_temperature_power_residual_W": solve_result.residual_W,
            "electron_temperature_relative_power_residual": solve_result.relative_residual,
            "electron_temperature_solve_numerics_settings": config.power_balance.solve_numerics_metadata(),
            "electron_temperature_search_skipped": False,
            "electron_temperature_search_skip_reason": None,
            "electron_temperature_diagnostic_initial_guess_state_built": False,
            "electron_temperature_search_bracket_found": solve_result.valid_bracket_keV is not None,
            "electron_temperature_search_converged": False,
            "electron_temperature_intermediate_trial_base_physics_valid": bool(solve_result.evaluations),
            "electron_temperature_intermediate_trial_residual_admitted": bool(solve_result.evaluations),
            "electron_temperature_final_qualification_passed": False,
            "electron_temperature_final_residual_consistent_with_trial": False,
            "electron_temperature_power_balance_terms_available": solve_result.final_evaluation is not None,
            **preflight_metadata,
            **_temperature_scope_metadata(),
        }
        if config.output.write_iteration_histories:
            metadata["electron_temperature_iteration_history"] = trial_records
       
        return _with_operating_point_metadata(operating_point, metadata)

    if solve_result.final_evaluation is None:
        raise RuntimeError("converged electron temperature solve has no final operating point")
  
    candidate_temperature_keV = float(solve_result.final_evaluation.electron_temperature_keV)
    emit_progress(progress_callback, OperatingPointProgressEvent(phase="final_qualification_started", evaluation_mode=OperatingPointEvaluationMode.FINAL_QUALIFICATION, temperature_keV=candidate_temperature_keV))
    final_started = clock()
    operating_point = _coupled_temperature_operating_point(
        config,
        geometry,
        electron_temperature_J=float(energy_J_from_keV(candidate_temperature_keV)),
        modal_basis=modal_basis,
        warm_start_state=None,
        collision_state_recompute_policy="recomputed_from_clean_state_for_final_physical_electron_energy_qualification",
        collision_state_is_self_consistent_electron_temperature=True,
        evaluation_mode=OperatingPointEvaluationMode.FINAL_QUALIFICATION,
        operating_point_builder=operating_point_builder,
        progress_callback=progress_callback,
        temperature_trial_index=None,
        clock=clock,
    )
    final_runtime_s = max(0.0, float(clock() - final_started))
    final_energy = operating_point.electron_energy_balance
    trial_energy = solve_result.final_evaluation.terms
    trial_residual_W = float(solve_result.final_evaluation.residual_W)
    final_residual_W = None if final_energy is None else float(final_energy.residual_W)
    difference_W = None if final_residual_W is None else abs(final_residual_W - trial_residual_W)
    difference_normalization_power_W = None if final_energy is None else max(
        abs(float(trial_energy.total_electron_heating_W)),
        abs(float(trial_energy.total_electron_loss_W)),
        abs(float(final_energy.total_electron_heating_W)),
        abs(float(final_energy.total_electron_loss_W)),
        1.0,
    )
    difference_relative = None if difference_W is None or difference_normalization_power_W is None else difference_W / difference_normalization_power_W
    residual_passed = bool(final_energy is not None and (abs(final_energy.residual_W) <= solve_config.residual_tolerance_W or abs(final_energy.relative_residual) <= solve_config.relative_residual_tolerance))
    residual_consistent = bool(difference_W is not None and difference_relative is not None and (difference_W <= 2.0 * solve_config.residual_tolerance_W or difference_relative <= 2.0 * solve_config.relative_residual_tolerance))
    kinetic_evidence_converged = _final_kinetic_evidence_converged(operating_point.kinetic.metadata)
    current_qualified = operating_point.kinetic.metadata.get("total_current_balance_converged") is True
    expander_qualified = operating_point.expander is not None and operating_point.expander.metadata.get("expander_qualified") is True
    energy_qualified = final_energy is not None and final_energy.represented_terms_qualified
    beam_current_qualified = operating_point.metadata.get("beam_current_outer_iteration_converged") is True
    final_comparison = ElectronTemperatureFinalQualificationComparison(electron_temperature_keV=candidate_temperature_keV, trial_residual_W=trial_residual_W, final_residual_W=final_residual_W, absolute_difference_W=difference_W, relative_difference=difference_relative, relative_difference_normalization_power_W=difference_normalization_power_W, final_residual_within_tolerance=residual_passed, residual_consistent_with_trial=residual_consistent, required_final_evidence_passed=bool(kinetic_evidence_converged and current_qualified and expander_qualified and energy_qualified and beam_current_qualified), operating_point_numerically_converged=operating_point.converged)
    final_passed = final_comparison.passed
    emit_progress(
        progress_callback,
        OperatingPointProgressEvent(
            phase="final_qualification_completed",
            evaluation_mode=OperatingPointEvaluationMode.FINAL_QUALIFICATION,
            temperature_keV=candidate_temperature_keV,
            elapsed_stage_s=final_runtime_s,
            current_residual_W=final_residual_W,
            current_relative_residual=None if final_energy is None else float(final_energy.relative_residual),
            residual_tolerance_W=float(solve_config.residual_tolerance_W),
            relative_residual_tolerance=float(solve_config.relative_residual_tolerance),
            power_balance_converged=bool(residual_passed),
            message="passed" if final_passed else "failed",
        ),
    )
    metadata = {
        "electron_temperature_mode": "self_consistent_electron_energy",
        "electron_temperature_assignment_valid": final_passed,
        "electron_temperature_solved_keV": candidate_temperature_keV if final_passed else None,
        "electron_temperature_diagnostic_root_keV": candidate_temperature_keV,
        "electron_temperature_diagnostic_root_converged": True,
        "electron_temperature_solve_converged": final_passed,
        "electron_temperature_power_balance_converged": final_passed,
        "electron_temperature_solve_failure_reason": None if final_passed else "final_qualification_failed",
        "electron_temperature_power_residual_W": final_residual_W,
        "electron_temperature_relative_power_residual": None if final_energy is None else final_energy.relative_residual,
        "electron_temperature_solve_iterations": solve_result.iterations,
        "electron_temperature_solve_numerics_settings": config.power_balance.solve_numerics_metadata(),
        "electron_temperature_valid_bracket_keV": None if solve_result.valid_bracket_keV is None else list(solve_result.valid_bracket_keV),
        "electron_temperature_search_skipped": False,
        "electron_temperature_search_skip_reason": None,
        "electron_temperature_diagnostic_initial_guess_state_built": False,
        "electron_temperature_search_bracket_found": solve_result.valid_bracket_keV is not None,
        "electron_temperature_search_converged": solve_result.converged,
        "electron_temperature_intermediate_trial_base_physics_valid": bool(solve_result.evaluations),
        "electron_temperature_intermediate_trial_residual_admitted": bool(solve_result.evaluations),
        "electron_temperature_final_qualification_passed": final_passed,
        "electron_temperature_final_residual_consistent_with_trial": residual_consistent,
        "electron_temperature_power_balance_terms_available": final_energy is not None,
        **preflight_metadata,
        **final_comparison.as_metadata(),
        **_temperature_scope_metadata(),
    }
  
    if final_energy is not None:
        metadata["electron_energy_balance_converged"] = final_passed
    if config.output.write_iteration_histories:
        metadata["electron_temperature_iteration_history"] = trial_records
  
    return _with_operating_point_metadata(operating_point, metadata)

def build_temperature_resolved_operating_point(config: SourceModelRunConfig, geometry: GeometryStageResult, *, modal_basis: Any | None = None, operating_point_builder: Callable[..., OperatingPointStageResult] = build_fixed_temperature_operating_point, progress_callback: OperatingPointProgressCallback | None = None, clock: Callable[[], float] = perf_counter) -> OperatingPointStageResult:
    """
    Build the operating point for the configured fixed or self consistent electron temperature closure
    
    A reusable modal basis is built once and shared across the selected temperature closure path
    """
    temperature_runtime_started = clock()
    runtime_profile = new_runtime_profile()
    mode = config.power_balance.electron_temperature_mode
    basis_started = clock()
    active_modal_basis = build_reusable_modal_basis(config, geometry) if modal_basis is None else modal_basis
    basis_runtime_s = max(0.0, float(clock() - basis_started))
   
    record_runtime(runtime_profile, "reusable_modal_basis", basis_runtime_s)
    emit_progress(progress_callback, OperatingPointProgressEvent(phase="basis_construction_completed", evaluation_mode=OperatingPointEvaluationMode.FIXED_TEMPERATURE_FINAL if mode == "fixed_closure" else OperatingPointEvaluationMode.TEMPERATURE_TRIAL, elapsed_stage_s=basis_runtime_s))
  
    operating_runtime_started = clock()
  
    if mode == "fixed_closure":
        result = _fixed_temperature_operating_point(config, geometry, modal_basis=active_modal_basis, operating_point_builder=operating_point_builder, progress_callback=progress_callback, clock=clock)
    elif mode == "self_consistent_electron_energy":
        result = _self_consistent_temperature_operating_point(config, geometry, modal_basis=active_modal_basis, operating_point_builder=operating_point_builder, progress_callback=progress_callback, clock=clock)
    else:
        raise ValueError(f"Unsupported power_balance.electron_temperature_mode {mode!r}")
  
    record_runtime(runtime_profile, "temperature_operating_point", clock() - operating_runtime_started)
    runtime_temperature_closure_profile = finalize_runtime_profile(runtime_profile, total_s=clock() - temperature_runtime_started)
    
    return _with_operating_point_metadata(result, {"runtime_temperature_closure_profile": runtime_temperature_closure_profile})

__all__ = ["build_temperature_resolved_operating_point"]
