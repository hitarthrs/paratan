"""Drive the nonlinear fixed point solution of the Egedal Eq 59 modal equations"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState, rebuild_fbis_collision_parameter_state_with_fast_density
from source_model_revamp.fbis.modal.types import ModalRosenbluthCoefficients
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid
from source_model_revamp.fbis.modal.eq59.types import Eq59ConvergenceDiagnostics, Eq59WarmStartState, _Eq59IterationMetrics
from source_model_revamp.fbis.modal.eq59.maxwellian_pairs import maxwellian_rosenbluth_coefficients
from source_model_revamp.fbis.modal.eq59.collision_operator_state import Eq59CollisionOperatorState, build_eq59_collision_operator_state
from source_model_revamp.fbis.modal.eq59.cross_species import ExternalFastFieldCollisionState
from source_model_revamp.fbis.modal.eq59.operator import _density_normalized_rosenbluth_coefficients, _solve_eq59_mode_finite_difference
from source_model_revamp.fbis.modal.eq59.state import _fixed_point_pattern_flags, _iteration_metrics, _modal_inventory_and_effective_energy, _prepare_eq59_warm_start, _validated_optional_physical_measures

def _positive_finite_scalar(value: float, name: str) -> float:
    """Validate one strictly positive finite scalar"""
    scalar = float(value)
    if not np.isfinite(scalar) or scalar <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
   
    return scalar

def _operator_from_distribution(*, collision_state: FBISCollisionParameterState, speed_grid: SpeedGrid, modal_distribution: np.ndarray, mode_eta_integrals: np.ndarray, particle_mass_kg: float, external_fast_field_states: tuple[ExternalFastFieldCollisionState, ...]) -> tuple[ModalRosenbluthCoefficients, Eq59CollisionOperatorState, float, float]:
    """Rebuild the nonlinear fast self collision state from one modal iterate
    
    The first mode supplies the density normalized Rosenbluth shape while the full modal state supplies physical fast density and mean energy
    """
    coefficients = _density_normalized_rosenbluth_coefficients(speed_grid, modal_distribution[0])
    density, energy = _modal_inventory_and_effective_energy(speed_grid=speed_grid, modal_distribution=modal_distribution, mode_eta_integrals=mode_eta_integrals, particle_mass_kg=particle_mass_kg)
    if density is None or not np.isfinite(density) or density <= 0.0:
        raise ValueError("Eq 59 modal distribution must have positive physical fast density")
    updated_collision_state = rebuild_fbis_collision_parameter_state_with_fast_density(collision_state, float(density))
    operator = build_eq59_collision_operator_state(pairwise_collision_state=updated_collision_state.pairwise_collision_state, speed_grid=speed_grid, spitzer_slowing_down_time_s=updated_collision_state.spitzer_slowing_down_time_s, fast_self_rosenbluth_coefficients=coefficients, fast_self_physical_density_m3=float(density), external_fast_field_states=external_fast_field_states)
  
    return coefficients, operator, float(density), float(energy if energy is not None else 0.0)

def _initial_pairwise_operator(*, collision_state: FBISCollisionParameterState, speed_grid: SpeedGrid, particle_mass_kg: float, source_mean_birth_energy_J: float, external_fast_field_states: tuple[ExternalFastFieldCollisionState, ...]) -> tuple[ModalRosenbluthCoefficients, Eq59CollisionOperatorState]:
    """Build the initial Eq 59 collision operator from a Maxwellian fast self seed
    
    The seed uses T_seed = 2 E_birth / 3 with the configured fast ion density
    """
    fast_density = _positive_finite_scalar(collision_state.fast_ion_density_m3, "collision_state.fast_ion_density_m3")
    mean_birth_energy = _positive_finite_scalar(source_mean_birth_energy_J, "source_mean_birth_energy_J")
    seed_temperature = 2.0 * mean_birth_energy / 3.0
    seed_coefficients = maxwellian_rosenbluth_coefficients(speed_m_s=speed_grid.centers_m_s, field_temperature_J=seed_temperature, field_mass_kg=particle_mass_kg)
    seeded_collision_state = rebuild_fbis_collision_parameter_state_with_fast_density(collision_state, fast_density)
    operator = build_eq59_collision_operator_state(pairwise_collision_state=seeded_collision_state.pairwise_collision_state, speed_grid=speed_grid, spitzer_slowing_down_time_s=seeded_collision_state.spitzer_slowing_down_time_s, fast_self_rosenbluth_coefficients=seed_coefficients, fast_self_physical_density_m3=fast_density, external_fast_field_states=external_fast_field_states)
  
    return seed_coefficients, operator

def _linear_eq59_solve(*, speed_grid: SpeedGrid, source_coefficients: np.ndarray, source_speeds: np.ndarray, eigenvalues_by_speed: np.ndarray, collision_operator_state: Eq59CollisionOperatorState, charge_exchange_loss_frequency_s: np.ndarray | None = None) -> np.ndarray:
    """Solve every retained modal speed equation for one fixed collision operator"""
    solution = np.zeros((eigenvalues_by_speed.shape[0], speed_grid.centers_m_s.size), dtype=float)
    for mode_index, eigenvalue in enumerate(eigenvalues_by_speed):
        solution[mode_index] = _solve_eq59_mode_finite_difference(speed_grid=speed_grid, source_coefficients_by_component=source_coefficients[:, mode_index], source_speeds_m_s=source_speeds, eigenvalue=eigenvalue, collision_operator_state=collision_operator_state, charge_exchange_loss_frequency_s=charge_exchange_loss_frequency_s)
   
    return solution

def hot_ion_rosenbluth_modal_distribution_eq59(*, speed_grid: SpeedGrid, source_coefficients_by_component_j: ArrayLike, source_speeds_m_s: ArrayLike, source_mean_birth_energy_J: float, eigenvalues: ArrayLike, collision_state: FBISCollisionParameterState, eigenvalues_by_speed: ArrayLike | None = None, iterations: int = 4, max_iterations: int | None = None, relative_tolerance: float = 1.0e-6, absolute_tolerance: float = 0.0, relaxation: float = 1.0, min_iterations: int | None = None, mode_eta_integrals: ArrayLike | None = None, physical_eigenfunctions_eta: ArrayLike | None = None, eta_grid: ArrayLike | None = None, particle_mass_kg: float | None = None, warm_start_state: Eq59WarmStartState | None = None, warm_start_compatibility_policy: str = "exact_physics", species_label: str = "D", charge_exchange_loss_frequency_s: ArrayLike | None = None, external_fast_field_states: tuple[ExternalFastFieldCollisionState, ...] = (), return_diagnostics: bool = False) -> tuple:
    """Iteratively solve Egedal Eq 59 with pair resolved ion collisions and analytic electron drag
    
    The modal array has shape (n_mode, n_speed) and source coefficients have shape (n_component, n_mode)
    Each accepted iterate rebuilds the first mode Rosenbluth coefficients physical fast density and collision operator
    Convergence requires the fixed point update nonlinear residual density mean energy and collision coefficient changes to satisfy tolerance
    """
    source_coefficients = np.asarray(source_coefficients_by_component_j, dtype=float)
    speeds = np.asarray(source_speeds_m_s, dtype=float)
    eigenvalue_array = np.asarray(eigenvalues, dtype=float)
    if source_coefficients.ndim != 2 or source_coefficients.shape[1] != eigenvalue_array.size:
        raise ValueError("source_coefficients_by_component_j must have shape (n_component, n_mode)")
    if speeds.shape != (source_coefficients.shape[0],):
        raise ValueError("source_speeds_m_s must have one value per component")
    if not np.all(np.isfinite(source_coefficients)):
        raise ValueError("source coefficients must be finite")
    if not np.all(np.isfinite(speeds)) or np.any(speeds <= 0.0):
        raise ValueError("source speeds must be positive and finite")
    if eigenvalue_array.ndim != 1 or eigenvalue_array.size == 0:
        raise ValueError("eigenvalues must be a nonempty one-dimensional array")
    if not np.all(np.isfinite(eigenvalue_array)) or np.any(eigenvalue_array <= 0.0):
        raise ValueError("eigenvalues must be positive and finite")
    if collision_state.fast_ion_species.symbol.strip().lower() != str(species_label).strip().lower():
        raise ValueError("collision state species must match species_label")
    if charge_exchange_loss_frequency_s is None:
        charge_exchange_frequency = None
    else:
        charge_exchange_frequency = np.asarray(charge_exchange_loss_frequency_s, dtype=float)
        if charge_exchange_frequency.shape != speed_grid.centers_m_s.shape or np.any(~np.isfinite(charge_exchange_frequency)) or np.any(charge_exchange_frequency < 0.0):
            raise ValueError("charge_exchange_loss_frequency_s must be finite and nonnegative on the Eq 59 speed grid")
    if eigenvalues_by_speed is None:
        eigenvalue_speed_array = np.broadcast_to(eigenvalue_array[:, None], (eigenvalue_array.size, speed_grid.centers_m_s.size)).copy()
    else:
        eigenvalue_speed_array = np.asarray(eigenvalues_by_speed, dtype=float)
        expected_eigenvalue_shape = (eigenvalue_array.size, speed_grid.centers_m_s.size)
        if eigenvalue_speed_array.shape != expected_eigenvalue_shape:
            raise ValueError("eigenvalues_by_speed must have shape (n_mode, n_speed_cell)")
        if np.any(~np.isfinite(eigenvalue_speed_array)) or np.any(eigenvalue_speed_array <= 0.0):
            raise ValueError("eigenvalues_by_speed must be positive and finite")
    if max_iterations is None:
        iteration_limit = max(int(iterations), 1)
        minimum_iterations = iteration_limit if min_iterations is None else int(min_iterations)
    else:
        iteration_limit = int(max_iterations)
        minimum_iterations = min(2, iteration_limit) if min_iterations is None else int(min_iterations)
    if iteration_limit < 1:
        raise ValueError("max_iterations must be at least one")
    if minimum_iterations < 1 or minimum_iterations > iteration_limit:
        raise ValueError("min_iterations must satisfy 1 <= min_iterations <= max_iterations")
    relative_tolerance_value = float(relative_tolerance)
    absolute_tolerance_value = float(absolute_tolerance)
    relaxation_value = float(relaxation)
    if not np.isfinite(relative_tolerance_value) or relative_tolerance_value <= 0.0:
        raise ValueError("relative_tolerance must be positive and finite")
    if not np.isfinite(absolute_tolerance_value) or absolute_tolerance_value < 0.0:
        raise ValueError("absolute_tolerance must be finite and nonnegative")
    if not np.isfinite(relaxation_value) or not 0.0 < relaxation_value <= 1.0:
        raise ValueError("relaxation must satisfy 0 < relaxation <= 1")
    eta_integrals, eigenfunctions, mass = _validated_optional_physical_measures(mode_count=eigenvalue_array.size, mode_eta_integrals=mode_eta_integrals, physical_eigenfunctions_eta=physical_eigenfunctions_eta, eta_grid=eta_grid, particle_mass_kg=particle_mass_kg)
    if eta_integrals is None or mass is None:
        raise ValueError("Eq 59 pairwise fast self scaling requires mode eta integrals and particle mass")
    external_fields = tuple(external_fast_field_states)
    seed_coefficients, seed_operator = _initial_pairwise_operator(collision_state=collision_state, speed_grid=speed_grid, particle_mass_kg=mass, source_mean_birth_energy_J=source_mean_birth_energy_J, external_fast_field_states=external_fields)
    warm_start = _prepare_eq59_warm_start(warm_start_state=warm_start_state, speed_grid=speed_grid, source_coefficients=source_coefficients, source_speeds=speeds, eigenvalues=eigenvalue_array, active_eigenvalues_by_speed=eigenvalue_speed_array, collision_operator_state=seed_operator, particle_mass_kg=mass, species_label=species_label, mode_eta_integrals=eta_integrals, physical_eigenfunctions_eta=eigenfunctions, eta_grid=(None if eta_grid is None else np.asarray(eta_grid, dtype=float)), charge_exchange_loss_frequency_s=charge_exchange_frequency, compatibility_policy=warm_start_compatibility_policy)
    if warm_start.used:
        current = np.asarray(warm_start.distribution, dtype=float)
    else:
        current = _linear_eq59_solve(speed_grid=speed_grid, source_coefficients=source_coefficients, source_speeds=speeds, eigenvalues_by_speed=eigenvalue_speed_array, collision_operator_state=seed_operator, charge_exchange_loss_frequency_s=charge_exchange_frequency)
    coefficients, operator, previous_inventory, previous_energy = _operator_from_distribution(collision_state=collision_state, speed_grid=speed_grid, modal_distribution=current, mode_eta_integrals=eta_integrals, particle_mass_kg=mass, external_fast_field_states=external_fields)
    history: list[dict[str, float]] = []
    metrics_history: list[_Eq59IterationMetrics] = []
    update_cosines: list[float] = []
    previous_update: np.ndarray | None = None
    converged = False
    diverged = False
    stagnated = False
    oscillatory = False
    for iteration_index in range(iteration_limit):
        fixed_point_solution = _linear_eq59_solve(speed_grid=speed_grid, source_coefficients=source_coefficients, source_speeds=speeds, eigenvalues_by_speed=eigenvalue_speed_array, collision_operator_state=operator, charge_exchange_loss_frequency_s=charge_exchange_frequency)
        accepted = current + relaxation_value * (fixed_point_solution - current)
        if not np.all(np.isfinite(accepted)):
            diverged = True
            break
        # Rebuild the fast self Rosenbluth state from each accepted nonlinear iterate
        accepted_coefficients, accepted_operator, accepted_inventory, accepted_energy = _operator_from_distribution(collision_state=collision_state, speed_grid=speed_grid, modal_distribution=accepted, mode_eta_integrals=eta_integrals, particle_mass_kg=mass, external_fast_field_states=external_fields)
        metrics = _iteration_metrics(speed_grid=speed_grid, previous=current, current=accepted, source_coefficients_by_component_j=source_coefficients, source_speeds_m_s=speeds, eigenvalues_by_speed=eigenvalue_speed_array, collision_operator_state=accepted_operator, previous_collision_operator_state=operator, mode_eta_integrals=eta_integrals, particle_mass_kg=mass, previous_inventory=previous_inventory, previous_effective_energy_J=previous_energy, charge_exchange_loss_frequency_s=charge_exchange_frequency)
        metrics_history.append(metrics)
        update = accepted - current
        if previous_update is not None:
            norm_product = float(np.linalg.norm(update) * np.linalg.norm(previous_update))
            cosine = float(np.vdot(update, previous_update).real / norm_product) if norm_product > 0.0 else 1.0
            update_cosines.append(cosine)
        previous_update = update
        history.append({
            "iteration": float(iteration_index + 1),
            "relative_change": metrics.relative_change,
            "absolute_change": metrics.absolute_change,
            "relative_residual": metrics.relative_residual,
            "absolute_residual": metrics.absolute_residual,
            "fast_self_physical_density_m3": metrics.fast_self_density_m3,
            "fast_self_g_scale_m3_s3": metrics.fast_self_g_scale_m3_s3,
            "fast_self_net_drag_scale_m3_s3": metrics.fast_self_net_drag_scale_m3_s3,
            "ion_drag_operator_relative_change": float(metrics.ion_drag_operator_relative_change or 0.0),
            "ion_energy_diffusion_operator_relative_change": float(metrics.ion_energy_diffusion_operator_relative_change or 0.0),
            "ion_pitch_scattering_operator_relative_change": float(metrics.ion_pitch_scattering_operator_relative_change or 0.0),
            "f1_density_normalization_m3": float(accepted_coefficients.density_normalization_m3),
            "h_tilde_max": float(np.max(accepted_coefficients.h_tilde)),
            "g_tilde_1_max": float(np.max(accepted_coefficients.g_tilde_1)),
            "g_tilde_2_max": float(np.max(accepted_coefficients.g_tilde_2)),
        })
        current = accepted
        coefficients = accepted_coefficients
        operator = accepted_operator
        previous_inventory = accepted_inventory
        previous_energy = accepted_energy
        relative_changes = [item.relative_change for item in metrics_history]
        residuals = [item.relative_residual for item in metrics_history]
        stagnated, oscillatory, detected_divergence = _fixed_point_pattern_flags(relative_changes, residuals, update_cosines, relative_tolerance_value)
        diverged = diverged or detected_divergence
        absolute_and_relative_change_converged = metrics.absolute_change <= absolute_tolerance_value + relative_tolerance_value * metrics.solution_scale
        residual_converged = metrics.relative_residual <= relative_tolerance_value
        inventory_converged = metrics.inventory_relative_change is None or metrics.inventory_relative_change <= relative_tolerance_value
        energy_converged = metrics.effective_energy_relative_change is None or metrics.effective_energy_relative_change <= relative_tolerance_value
        operator_changes = (metrics.ion_drag_operator_relative_change, metrics.ion_energy_diffusion_operator_relative_change, metrics.ion_pitch_scattering_operator_relative_change)
        operator_converged = all(value is None or value <= relative_tolerance_value for value in operator_changes)
        converged = bool(iteration_index + 1 >= minimum_iterations and absolute_and_relative_change_converged and residual_converged and inventory_converged and energy_converged and operator_converged)
        if converged or diverged:
            break
    final_metrics = metrics_history[-1] if metrics_history else None
    fixed_point_converged = bool(final_metrics is not None and final_metrics.absolute_change <= absolute_tolerance_value + relative_tolerance_value * final_metrics.solution_scale)
    residual_converged = bool(final_metrics is not None and final_metrics.relative_residual <= relative_tolerance_value)
    inventory_converged = bool(final_metrics is not None and (final_metrics.inventory_relative_change is None or final_metrics.inventory_relative_change <= relative_tolerance_value))
    effective_energy_converged = bool(final_metrics is not None and (final_metrics.effective_energy_relative_change is None or final_metrics.effective_energy_relative_change <= relative_tolerance_value))
    operator_converged = bool(final_metrics is not None and all(value is None or value <= relative_tolerance_value for value in (final_metrics.ion_drag_operator_relative_change, final_metrics.ion_energy_diffusion_operator_relative_change, final_metrics.ion_pitch_scattering_operator_relative_change)))
    if converged:
        status = "converged"
        failure_reason = ""
    elif diverged:
        status = "diverged"
        failure_reason = "Eq 59 nonlinear iteration diverged"
    elif oscillatory:
        status = "oscillatory"
        failure_reason = "Eq 59 nonlinear iteration is oscillatory"
    elif stagnated:
        status = "stagnated"
        failure_reason = "Eq 59 nonlinear iteration stagnated above tolerance"
    else:
        status = "maximum_iterations_reached"
        failure_reason = "Eq 59 nonlinear iteration did not converge within max_iterations"
    diagnostics = Eq59ConvergenceDiagnostics(
        assessed=True,
        converged=converged,
        iterations=len(metrics_history),
        status=status,
        failure_reason=failure_reason,
        relative_change_history=tuple(item.relative_change for item in metrics_history),
        absolute_change_history=tuple(item.absolute_change for item in metrics_history),
        residual_history=tuple(item.relative_residual for item in metrics_history),
        absolute_residual_history=tuple(item.absolute_residual for item in metrics_history),
        inventory_history=tuple(item.inventory for item in metrics_history if item.inventory is not None),
        inventory_relative_change_history=tuple(item.inventory_relative_change for item in metrics_history if item.inventory_relative_change is not None),
        effective_energy_history_J=tuple(item.effective_energy_J for item in metrics_history if item.effective_energy_J is not None),
        effective_energy_relative_change_history=tuple(item.effective_energy_relative_change for item in metrics_history if item.effective_energy_relative_change is not None),
        fast_self_density_history_m3=tuple(item.fast_self_density_m3 for item in metrics_history),
        fast_self_density_relative_change_history=tuple(item.fast_self_density_relative_change for item in metrics_history if item.fast_self_density_relative_change is not None),
        fast_self_g_scale_history_m3_s3=tuple(item.fast_self_g_scale_m3_s3 for item in metrics_history),
        fast_self_net_drag_scale_history_m3_s3=tuple(item.fast_self_net_drag_scale_m3_s3 for item in metrics_history),
        ion_drag_operator_relative_change_history=tuple(item.ion_drag_operator_relative_change for item in metrics_history if item.ion_drag_operator_relative_change is not None),
        ion_energy_diffusion_operator_relative_change_history=tuple(item.ion_energy_diffusion_operator_relative_change for item in metrics_history if item.ion_energy_diffusion_operator_relative_change is not None),
        ion_pitch_scattering_operator_relative_change_history=tuple(item.ion_pitch_scattering_operator_relative_change for item in metrics_history if item.ion_pitch_scattering_operator_relative_change is not None),
        stagnated=stagnated,
        oscillatory=oscillatory,
        diverged=diverged,
        fixed_point_iteration_converged=fixed_point_converged,
        nonlinear_residual_converged=residual_converged,
        inventory_converged=inventory_converged,
        effective_energy_converged=effective_energy_converged,
        fast_self_density_converged=inventory_converged,
        collision_operator_converged=operator_converged,
        speed_resolution_assessed=False,
        speed_resolution_converged=False,
        speed_domain_assessed=False,
        speed_domain_converged=False,
        warm_start_used=warm_start.used,
        warm_start_remapped=warm_start.remapped,
        warm_start_status=warm_start.status,
        warm_start_rejection_reason=warm_start.rejection_reason,
        solve_method="scipy.linalg.solve_banded_tridiagonal_pairwise_operator",
        collision_operator_model=operator.operator_model,
    )
    if warm_start.used and not diagnostics.converged:
        cold_distribution, cold_coefficients, cold_operator, cold_history, cold_diagnostics = hot_ion_rosenbluth_modal_distribution_eq59(speed_grid=speed_grid, source_coefficients_by_component_j=source_coefficients_by_component_j, source_speeds_m_s=source_speeds_m_s, source_mean_birth_energy_J=source_mean_birth_energy_J, eigenvalues=eigenvalues, collision_state=collision_state, eigenvalues_by_speed=eigenvalues_by_speed, iterations=iterations, max_iterations=max_iterations, relative_tolerance=relative_tolerance, absolute_tolerance=absolute_tolerance, relaxation=relaxation, min_iterations=min_iterations, mode_eta_integrals=mode_eta_integrals, physical_eigenfunctions_eta=physical_eigenfunctions_eta, eta_grid=eta_grid, particle_mass_kg=particle_mass_kg, warm_start_state=None, warm_start_compatibility_policy="exact_physics", species_label=species_label, charge_exchange_loss_frequency_s=charge_exchange_frequency, external_fast_field_states=external_fields, return_diagnostics=True)
        fallback_diagnostics = replace(cold_diagnostics, warm_start_status="warm_start_fallback_to_pairwise_seed_after_nonconvergence", warm_start_fallback_to_pairwise_seed=True, discarded_warm_start_iterations=diagnostics.iterations, discarded_warm_start_status=diagnostics.status, discarded_warm_start_remapped=warm_start.remapped)
      
        if return_diagnostics:
            return cold_distribution, cold_coefficients, cold_operator, cold_history, fallback_diagnostics
     
        return cold_distribution, cold_coefficients, cold_operator, cold_history
  
    if return_diagnostics:
        return current, coefficients, operator, history, diagnostics

    return current, coefficients, operator, history
