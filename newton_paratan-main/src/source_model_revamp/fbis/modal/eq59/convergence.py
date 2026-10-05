"""Assess nonlinear speed resolution and speed domain convergence for Egedal Eq 59"""
from __future__ import annotations
from time import perf_counter
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState
from source_model_revamp.fbis.modal.types import ModalRosenbluthCoefficients
from source_model_revamp.fbis.modal.eq59.collision_operator_state import Eq59CollisionOperatorState
from source_model_revamp.fbis.modal.eq59.energy_moments import build_eq59_energy_moment_state
from source_model_revamp.fbis.modal.utils import _EPS
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid, source_aligned_stretched_speed_grid
from source_model_revamp.fbis.modal.eq59.types import Eq59ConvergenceDiagnostics, Eq59SpeedConvergenceDiagnostics, Eq59WarmStartState
from source_model_revamp.fbis.modal.eq59.state import _conservative_shell_remap, _modal_inventory_and_effective_energy, _relative_scalar_change, build_eq59_warm_start_state
from source_model_revamp.fbis.modal.eq59.solve import hot_ion_rosenbluth_modal_distribution_eq59
from source_model_revamp.fbis.modal.eq59.cross_species import ExternalFastFieldCollisionState

def _upper_tail_fraction_from_cell_weights(speed_grid: SpeedGrid, weights: np.ndarray, start_fraction: float = 0.95, radial_measure_power: int = 3) -> float:
    """Integrate a cellwise moment above a fractional upper speed boundary"""
    values = np.asarray(weights, dtype=float)
    if values.shape != speed_grid.centers_m_s.shape or np.any(~np.isfinite(values)):
        raise ValueError("tail weights must contain one finite value per speed cell")
    total = float(np.sum(values))
    if total <= 0.0:
        return 0.0
    power = int(radial_measure_power)
    if power not in {3, 5}:
        raise ValueError("radial_measure_power must be 3 for population or 5 for kinetic energy")
    boundary = float(start_fraction) * float(speed_grid.faces_m_s[-1])
    lower = np.asarray(speed_grid.faces_m_s[:-1], dtype=float)
    upper = np.asarray(speed_grid.faces_m_s[1:], dtype=float)
    clipped_lower = np.maximum(lower, boundary)
    fraction = np.divide(np.maximum(upper**power - clipped_lower**power, 0.0), upper**power - lower**power, out=np.zeros_like(lower), where=upper > lower)

    return float(np.sum(values * fraction) / total)

def _eq59_modal_loss_moments(*, speed_grid: SpeedGrid, modal_distribution: np.ndarray, mode_eta_integrals: np.ndarray, eigenvalues_by_speed: np.ndarray, collision_operator_state: Eq59CollisionOperatorState, particle_mass_kg: float) -> tuple[float, float]:
    """Return modal pitch scattering loss rate and kinetic power density
    
    The speed dependent loss frequency is λ_j P(v) / (τ_s v³)
    """
    v = np.asarray(speed_grid.centers_m_s, dtype=float)
    shell = np.asarray(speed_grid.shell_volumes_m3_s3, dtype=float)
    pitch = np.asarray(collision_operator_state.ion_pitch_scattering_velocity_cubed_m3_s3, dtype=float)
    if pitch.shape != v.shape or np.any(~np.isfinite(pitch)) or np.any(pitch < 0.0):
        raise ValueError("Eq 59 pitch scattering array must match the speed grid")
    loss_frequency = pitch / (float(collision_operator_state.spitzer_slowing_down_time_s) * np.maximum(v, _EPS) ** 3)
    rate_by_mode_and_speed = (modal_distribution * eigenvalues_by_speed * mode_eta_integrals[:, None] * loss_frequency[None, :] * shell[None, :])
    rate = float(np.sum(rate_by_mode_and_speed))
    power = float(np.sum(rate_by_mode_and_speed * (0.5 * float(particle_mass_kg) * v[None, :] ** 2)))

    return rate, power

def build_eq59_speed_convergence_grid_sequence(*, production_speed_grid: SpeedGrid, source_speeds_m_s: ArrayLike, core_cell_fraction: float = 0.7, tail_stretch_power: float = 2.0) -> tuple[tuple[SpeedGrid, ...], int]:
    """Build four distinct speed grids around the configured solve grid
    
    The sequence separates resolution changes from upper speed domain changes while retaining source aligned cells
    """
    production = production_speed_grid
    source_speeds = np.asarray(source_speeds_m_s, dtype=float).reshape(-1)
    if (source_speeds.size == 0 or np.any(~np.isfinite(source_speeds)) or np.any(source_speeds <= 0.0)):
        raise ValueError("source_speeds_m_s must contain positive finite values")
    production_cells = int(production.centers_m_s.size)
    production_upper = float(production.faces_m_s[-1])
    if not np.isfinite(production_upper) or production_upper <= 0.0:
        raise ValueError("production_speed_grid must have a positive finite upper face")
    largest_source = float(np.max(source_speeds))
    roundoff = 256.0 * np.finfo(float).eps * max(production_upper, largest_source, 1.0)
    smaller_upper = max(0.75 * production_upper, largest_source + roundoff)
    if smaller_upper >= production_upper - roundoff:
        smaller_upper = max(0.9 * production_upper, largest_source + roundoff)
    if smaller_upper >= production_upper - roundoff:
        raise ValueError("production speed domain is too close to the largest source speed for an independent domain convergence assessment")
    core_fraction = float(core_cell_fraction)
    tail_power = float(tail_stretch_power)
    if not np.isfinite(core_fraction) or not 0.0 < core_fraction < 1.0:
        raise ValueError("core_cell_fraction must lie inside (0, 1)")
    if not np.isfinite(tail_power) or tail_power < 1.0:
        raise ValueError("tail_stretch_power must be finite and at least one")
    minimum_core_cells = int(np.unique(source_speeds).size) + 2
    coarse_cells = max(production_cells // 2, minimum_core_cells, 4)
    coarse_cells = min(coarse_cells, production_cells - 1)
    smaller_coarse = source_aligned_stretched_speed_grid(max_speed_m_s=smaller_upper, num_cells=coarse_cells, source_speeds_m_s=source_speeds, core_cell_fraction=core_fraction, tail_stretch_power=tail_power)
    smaller_requested = source_aligned_stretched_speed_grid(max_speed_m_s=smaller_upper, num_cells=production_cells, source_speeds_m_s=source_speeds, core_cell_fraction=core_fraction, tail_stretch_power=tail_power)
    production_refined = source_aligned_stretched_speed_grid(max_speed_m_s=production_upper, num_cells=2 * production_cells, source_speeds_m_s=source_speeds, core_cell_fraction=core_fraction, tail_stretch_power=tail_power)
    grids = (smaller_coarse, smaller_requested, production, production_refined)
    signatures = tuple((int(grid.centers_m_s.size), float(grid.faces_m_s[-1])) for grid in grids)
    if len(signatures) != len(set(signatures)):
        raise RuntimeError("Eq 59 speed grid sequence contains duplicate levels")
    
    return grids, 2

def _mapped_charge_exchange_loss_frequency(*, source_speed_grid: SpeedGrid | None, source_loss_frequency_s: ArrayLike | None, target_speed_grid: SpeedGrid) -> np.ndarray | None:
    """Map a speed resolved charge exchange loss frequency onto one comparison grid"""
    if (source_speed_grid is None) != (source_loss_frequency_s is None):
        raise ValueError("charge exchange source grid and loss frequency must be supplied together")
    if source_speed_grid is None or source_loss_frequency_s is None:
        return None
    source_speed = np.asarray(source_speed_grid.centers_m_s, dtype=float)
    source_frequency = np.asarray(source_loss_frequency_s, dtype=float)
    target_speed = np.asarray(target_speed_grid.centers_m_s, dtype=float)
    if source_frequency.shape != source_speed.shape or np.any(~np.isfinite(source_frequency)) or np.any(source_frequency < 0.0):
        raise ValueError("charge exchange loss frequency must be finite and nonnegative on its source grid")
    if source_speed.shape == target_speed.shape and np.allclose(source_speed, target_speed, rtol=1.0e-13, atol=0.0):
        return source_frequency.copy()

    return np.interp(target_speed, source_speed, source_frequency, left=0.0, right=0.0)

def assess_eq59_speed_convergence(*, speed_grid_sequence: tuple[SpeedGrid, ...], source_coefficients_by_component_j: ArrayLike, source_speeds_m_s: ArrayLike, source_mean_birth_energy_J: float, eigenvalues: ArrayLike, collision_state: FBISCollisionParameterState, mode_eta_integrals: ArrayLike, physical_eigenfunctions_eta: ArrayLike, eta_grid: ArrayLike, particle_mass_kg: float, eigenvalues_by_speed_sequence: tuple[ArrayLike, ...] | None = None, max_iterations: int = 32, min_iterations: int = 2, relative_tolerance: float = 1.0e-4, absolute_tolerance: float = 0.0, relaxation: float = 1.0, speed_convergence_relative_tolerance: float | None = None, high_speed_tail_population_tolerance: float = 1.0e-4, high_speed_tail_energy_tolerance: float | None = None, production_grid_index: int = 0, precomputed_modal_distribution: np.ndarray | None = None, precomputed_rosenbluth_coefficients: ModalRosenbluthCoefficients | None = None, precomputed_collision_operator_state: Eq59CollisionOperatorState | None = None, precomputed_nonlinear_diagnostics: Eq59ConvergenceDiagnostics | None = None, external_fast_field_states: tuple[ExternalFastFieldCollisionState, ...] = (), charge_exchange_loss_frequency_source_grid: SpeedGrid | None = None, charge_exchange_loss_frequency_s: ArrayLike | None = None) -> Eq59SpeedConvergenceDiagnostics:
    """Compare Eq 59 solutions across nested speed resolution and domain levels
    
    Physical f(v,η) is compared on common domains after conservative shell remapping
    The assessment also tracks density energy modal losses electron heating source to loss balance and upper speed tails
    """
    grids = tuple(speed_grid_sequence)
    selected_index = int(production_grid_index)
    if selected_index < 0 or selected_index >= len(grids):
        raise ValueError("production_grid_index lies outside speed_grid_sequence")
    precomputed_items = (precomputed_modal_distribution, precomputed_rosenbluth_coefficients, precomputed_collision_operator_state, precomputed_nonlinear_diagnostics)
    if any(item is not None for item in precomputed_items) and not all(item is not None for item in precomputed_items):
        raise ValueError("precomputed Eq 59 state inputs must be supplied together")
    signatures = tuple((grid.centers_m_s.size, float(grid.faces_m_s[-1])) for grid in grids)
    if len(set(signatures)) != len(signatures):
        raise ValueError("speed_grid_sequence levels must be distinct")
    for previous, current in zip(grids[:-1], grids[1:], strict=True):
        if current.centers_m_s.size < previous.centers_m_s.size or current.faces_m_s[-1] < previous.faces_m_s[-1]:
            raise ValueError("speed grids must be nested in nondecreasing resolution and extent")
    if eigenvalues_by_speed_sequence is None:
        eigenvalue_levels: tuple[ArrayLike | None, ...] = (None,) * len(grids)
    else:
        eigenvalue_levels = tuple(eigenvalues_by_speed_sequence)
        if len(eigenvalue_levels) != len(grids):
            raise ValueError("eigenvalues_by_speed_sequence must match speed_grid_sequence")
    source_coefficients = np.asarray(source_coefficients_by_component_j, dtype=float)
    source_speeds = np.asarray(source_speeds_m_s, dtype=float)
    magnetic_eigenvalues = np.asarray(eigenvalues, dtype=float)
    mode_integrals = np.asarray(mode_eta_integrals, dtype=float)
    eigenfunctions = np.asarray(physical_eigenfunctions_eta, dtype=float)
    eta = np.asarray(eta_grid, dtype=float)
    comparison_tolerance = (float(relative_tolerance) if speed_convergence_relative_tolerance is None else float(speed_convergence_relative_tolerance))
    energy_tail_tolerance = (float(high_speed_tail_population_tolerance) if high_speed_tail_energy_tolerance is None else float(high_speed_tail_energy_tolerance))
    if not np.isfinite(comparison_tolerance) or comparison_tolerance <= 0.0:
        raise ValueError("speed_convergence_relative_tolerance must be positive and finite")
    if (not np.isfinite(high_speed_tail_population_tolerance) or high_speed_tail_population_tolerance < 0.0):
        raise ValueError("high_speed_tail_population_tolerance must be finite and nonnegative")
    if not np.isfinite(energy_tail_tolerance) or energy_tail_tolerance < 0.0:
        raise ValueError("high_speed_tail_energy_tolerance must be finite and nonnegative")

    solutions: list[np.ndarray] = []
    coefficients_history: list[ModalRosenbluthCoefficients] = []
    operator_history: list[Eq59CollisionOperatorState] = []
    diagnostics_history: list[Eq59ConvergenceDiagnostics] = []
    runtimes: list[float] = []
    inventories: list[float] = []
    energies: list[float] = []
    losses: list[float] = []
    powers: list[float] = []
    electron_heating_powers: list[float] = []
    tail_population: list[float] = []
    tail_energy: list[float] = []
    source_to_loss: list[float] = []
    warm_state: Eq59WarmStartState | None = None
    source_rate = float(4.0 * np.pi * np.sum(source_coefficients * mode_integrals[None, :]))
    for level_index, (grid, eigenvalue_level) in enumerate(zip(grids, eigenvalue_levels, strict=True)):
        charge_exchange_frequency = _mapped_charge_exchange_loss_frequency(
            source_speed_grid=charge_exchange_loss_frequency_source_grid,
            source_loss_frequency_s=charge_exchange_loss_frequency_s,
            target_speed_grid=grid,
        )
        if level_index == selected_index and precomputed_modal_distribution is not None:
            solution = np.asarray(precomputed_modal_distribution, dtype=float)
            coefficients = precomputed_rosenbluth_coefficients
            operator_state = precomputed_collision_operator_state
            diagnostics = precomputed_nonlinear_diagnostics
            expected_shape = (magnetic_eigenvalues.size, grid.centers_m_s.size)
            if solution.shape != expected_shape or np.any(~np.isfinite(solution)):
                raise ValueError("precomputed_modal_distribution does not match the selected production grid")
            if coefficients is None or operator_state is None or diagnostics is None:
                raise RuntimeError("validated precomputed Eq 59 state is incomplete")
            runtimes.append(0.0)
        else:
            start = perf_counter()
            solution, coefficients, operator_state, _, diagnostics = hot_ion_rosenbluth_modal_distribution_eq59(
                speed_grid=grid,
                source_coefficients_by_component_j=source_coefficients,
                source_speeds_m_s=source_speeds,
                source_mean_birth_energy_J=source_mean_birth_energy_J,
                eigenvalues=magnetic_eigenvalues,
                eigenvalues_by_speed=eigenvalue_level,
                collision_state=collision_state,
                max_iterations=max_iterations,
                min_iterations=min_iterations,
                relative_tolerance=relative_tolerance,
                absolute_tolerance=absolute_tolerance,
                relaxation=relaxation,
                mode_eta_integrals=mode_integrals,
                physical_eigenfunctions_eta=eigenfunctions,
                eta_grid=eta,
                particle_mass_kg=particle_mass_kg,
                warm_start_state=warm_state,
                warm_start_compatibility_policy=("exact_physics" if warm_state is None else "closure_continuation"),
                species_label=collision_state.fast_ion_species.symbol,
                external_fast_field_states=tuple(external_fast_field_states),
                charge_exchange_loss_frequency_s=charge_exchange_frequency,
                return_diagnostics=True,
            )
            runtimes.append(perf_counter() - start)
        solutions.append(solution)
        coefficients_history.append(coefficients)
        operator_history.append(operator_state)
        diagnostics_history.append(diagnostics)
        inventory, energy = _modal_inventory_and_effective_energy(speed_grid=grid, modal_distribution=solution, mode_eta_integrals=mode_integrals, particle_mass_kg=particle_mass_kg,)
        if inventory is None or energy is None:
            raise RuntimeError("nested Eq 59 speed assessment requires inventory and energy moments")
        inventories.append(float(inventory))
        energies.append(float(energy))
        eigenvalue_array = (np.broadcast_to(magnetic_eigenvalues[:, None], solution.shape) if eigenvalue_level is None else np.asarray(eigenvalue_level, dtype=float))
        loss, power = _eq59_modal_loss_moments(speed_grid=grid, modal_distribution=solution, mode_eta_integrals=mode_integrals, eigenvalues_by_speed=eigenvalue_array, collision_operator_state=operator_state, particle_mass_kg=particle_mass_kg,)
        losses.append(loss)
        powers.append(power)
        energy_moment_state = build_eq59_energy_moment_state(
            speed_grid=grid,
            modal_distribution_f_j_v=solution,
            source_coefficients_by_component_j=source_coefficients,
            source_speeds_m_s=source_speeds,
            active_eigenvalues_by_speed=eigenvalue_array,
            collision_operator_state=operator_state,
            mode_eta_integrals=mode_integrals,
            particle_mass_kg=particle_mass_kg,
            volume_m3=1.0,
            exact_confined_source_power_W=source_rate * float(source_mean_birth_energy_J),
            charge_exchange_loss_frequency_s=charge_exchange_frequency,
            energy_residual_relative_tolerance=max(float(relative_tolerance), np.finfo(float).eps),
            source_power_relative_tolerance=max(comparison_tolerance, np.finfo(float).eps),
        )
        electron_heating_powers.append(float(energy_moment_state.electron_heating_power_W))
        eta_integrated = mode_integrals @ solution
        population_weights = np.asarray(grid.shell_volumes_m3_s3, dtype=float) * eta_integrated
        lower_faces = np.asarray(grid.faces_m_s[:-1], dtype=float)
        upper_faces = np.asarray(grid.faces_m_s[1:], dtype=float)
        energy_weights = (eta_integrated * (2.0 * np.pi * float(particle_mass_kg) / 5.0) * (upper_faces**5 - lower_faces**5))
        tail_population.append(_upper_tail_fraction_from_cell_weights(grid, population_weights))
        tail_energy.append(_upper_tail_fraction_from_cell_weights(grid, energy_weights, radial_measure_power=5))
        source_to_loss.append(source_rate / loss if loss > 0.0 else np.inf)
        warm_state = build_eq59_warm_start_state(
            speed_grid=grid,
            modal_distribution_f_j_v=solution,
            magnetic_eigenvalues=magnetic_eigenvalues,
            active_eigenvalues_by_speed=eigenvalue_array,
            source_coefficients_by_component_j=source_coefficients,
            source_speeds_m_s=source_speeds,
            collision_operator_state=operator_state,
            particle_mass_kg=particle_mass_kg,
            charge_exchange_loss_frequency_s=charge_exchange_frequency,
            species_label=collision_state.fast_ion_species.symbol,
            nonlinear_converged=bool(diagnostics.converged),
            mode_eta_integrals=mode_integrals,
            physical_eigenfunctions_eta=eigenfunctions,
            eta_grid=eta,
        )

    distribution_changes: list[float] = []
    inventory_changes: list[float] = []
    energy_changes: list[float] = []
    loss_changes: list[float] = []
    power_changes: list[float] = []
    electron_heating_changes: list[float] = []
    source_to_loss_changes: list[float] = []
    for index in range(1, len(grids)):
        comparison_grid = grids[index - 1]
        # Compare consecutive solutions after remapping onto the smaller common speed domain
        current_on_common = _conservative_shell_remap(grids[index].faces_m_s, solutions[index], comparison_grid)
        previous_physical = np.einsum("jv,je->ve", solutions[index - 1], eigenfunctions, optimize=True)
        current_physical = np.einsum("jv,je->ve", current_on_common, eigenfunctions, optimize=True)
        difference = current_physical - previous_physical
        shell = np.asarray(comparison_grid.shell_volumes_m3_s3, dtype=float)
        difference_norm = np.sqrt(float(np.sum(shell * np.trapezoid(difference**2, eta, axis=1))))
        current_norm = np.sqrt(float(np.sum(shell * np.trapezoid(current_physical**2, eta, axis=1))))
        distribution_changes.append(difference_norm / max(current_norm, np.finfo(float).tiny))
        inventory_changes.append(_relative_scalar_change(inventories[index], inventories[index - 1]))
        energy_changes.append(_relative_scalar_change(energies[index], energies[index - 1]))
        loss_changes.append(_relative_scalar_change(losses[index], losses[index - 1]))
        power_changes.append(_relative_scalar_change(powers[index], powers[index - 1]))
        electron_heating_changes.append(_relative_scalar_change(electron_heating_powers[index], electron_heating_powers[index - 1]))
        source_to_loss_changes.append(_relative_scalar_change(source_to_loss[index], source_to_loss[index - 1]))
    resolution_indices = [index for index in range(1, len(grids)) if np.isclose(grids[index].faces_m_s[-1], grids[index - 1].faces_m_s[-1], rtol=1.0e-13, atol=0.0) and grids[index].centers_m_s.size > grids[index - 1].centers_m_s.size]
    domain_indices = [index for index in range(1, len(grids)) if grids[index].faces_m_s[-1] > grids[index - 1].faces_m_s[-1]]
    selected_resolution_indices = [index for index in resolution_indices if index == selected_index or index - 1 == selected_index]
    selected_domain_indices = [index for index in domain_indices if index == selected_index or index - 1 == selected_index]
    comparison_arrays = (distribution_changes, inventory_changes, energy_changes, loss_changes, power_changes, electron_heating_changes, source_to_loss_changes)
    resolution_assessed = bool(selected_resolution_indices)
    domain_assessed = bool(selected_domain_indices)
    electron_heating_resolution_assessed = bool(selected_resolution_indices)
    electron_heating_domain_assessed = bool(selected_domain_indices)
    electron_heating_resolution_converged = bool(electron_heating_resolution_assessed and all(electron_heating_changes[index - 1] <= comparison_tolerance for index in selected_resolution_indices))
    electron_heating_domain_converged = bool(electron_heating_domain_assessed and all(electron_heating_changes[index - 1] <= comparison_tolerance for index in selected_domain_indices))
    tail_moments_valid = bool(np.all(np.isfinite(tail_population)) and np.all(np.isfinite(tail_energy)) and np.all(np.asarray(tail_population, dtype=float) >= 0.0) and np.all(np.asarray(tail_population, dtype=float) <= 1.0) and np.all(np.asarray(tail_energy, dtype=float) >= 0.0) and np.all(np.asarray(tail_energy, dtype=float) <= 1.0))
    resolution_converged = bool(resolution_assessed and all(all(values[index - 1] <= comparison_tolerance for values in comparison_arrays) for index in selected_resolution_indices))
    domain_converged = bool(
        domain_assessed
        and all(all(values[index - 1] <= comparison_tolerance for values in comparison_arrays) for index in selected_domain_indices)
        and tail_moments_valid
        and tail_population[selected_index]
        <= high_speed_tail_population_tolerance
        and tail_energy[selected_index] <= energy_tail_tolerance
        and tail_population[-1] <= high_speed_tail_population_tolerance
        and tail_energy[-1] <= energy_tail_tolerance
    )
    nonlinear_converged = np.asarray([item.converged for item in diagnostics_history], dtype=bool)
    required_nonlinear_levels = {selected_index}
    for transition_index in (selected_resolution_indices + selected_domain_indices):
        required_nonlinear_levels.add(transition_index - 1)
        required_nonlinear_levels.add(transition_index)
    required_nonlinear_converged = bool(all(nonlinear_converged[index] for index in required_nonlinear_levels))
    converged = bool(required_nonlinear_converged and resolution_converged and domain_converged)
    failures: list[str] = []
    if not required_nonlinear_converged:
        failures.append("selected_or_comparison_nonlinear_level_not_converged")
    if not resolution_assessed:
        failures.append("speed_resolution_not_assessed")
    elif not resolution_converged:
        failures.append("speed_resolution_not_converged")
    if not domain_assessed:
        failures.append("speed_domain_not_assessed")
    elif not domain_converged:
        failures.append("speed_domain_not_converged")
    if not tail_moments_valid:
        failures.append("one_or_more_tail_moments_outside_unit_interval")
    elif (tail_population[selected_index] > high_speed_tail_population_tolerance or tail_population[-1] > high_speed_tail_population_tolerance):
        failures.append("high_speed_particle_tail_above_tolerance")
    elif (tail_energy[selected_index] > energy_tail_tolerance or tail_energy[-1] > energy_tail_tolerance):
        failures.append("high_speed_energy_tail_above_tolerance")

    return Eq59SpeedConvergenceDiagnostics(
        assessed=True,
        converged=converged,
        speed_cell_count_history=np.asarray([grid.centers_m_s.size for grid in grids], dtype=int),
        speed_domain_upper_m_s_history=np.asarray([grid.faces_m_s[-1] for grid in grids], dtype=float),
        runtime_s_history=np.asarray(runtimes, dtype=float),
        nonlinear_converged_history=nonlinear_converged,
        common_domain_distribution_relative_change_history=np.asarray(distribution_changes, dtype=float),
        inventory_history=np.asarray(inventories, dtype=float),
        inventory_relative_change_history=np.asarray(inventory_changes, dtype=float),
        effective_energy_history_J=np.asarray(energies, dtype=float),
        effective_energy_relative_change_history=np.asarray(energy_changes, dtype=float),
        modal_particle_loss_rate_history=np.asarray(losses, dtype=float),
        modal_particle_loss_relative_change_history=np.asarray(loss_changes, dtype=float),
        modal_power_loss_history_W_per_m3=np.asarray(powers, dtype=float),
        modal_power_loss_relative_change_history=np.asarray(power_changes, dtype=float),
        electron_heating_power_history_W_per_m3=np.asarray(electron_heating_powers, dtype=float),
        electron_heating_power_relative_change_history=np.asarray(electron_heating_changes, dtype=float),
        high_speed_tail_population_fraction_history=np.asarray(tail_population, dtype=float),
        high_speed_tail_energy_fraction_history=np.asarray(tail_energy, dtype=float),
        source_to_modal_loss_ratio_history=np.asarray(source_to_loss, dtype=float),
        source_to_modal_loss_ratio_relative_change_history=np.asarray(source_to_loss_changes, dtype=float),
        speed_resolution_assessed=resolution_assessed,
        speed_resolution_converged=resolution_converged,
        speed_domain_assessed=domain_assessed,
        speed_domain_converged=domain_converged,
        electron_heating_speed_resolution_assessed=electron_heating_resolution_assessed,
        electron_heating_speed_resolution_converged=electron_heating_resolution_converged,
        electron_heating_speed_domain_assessed=electron_heating_domain_assessed,
        electron_heating_speed_domain_converged=electron_heating_domain_converged,
        production_grid_index=selected_index,
        production_grid_matched=bool(precomputed_modal_distribution is not None),
        selected_speed_cell_count=int(grids[selected_index].centers_m_s.size),
        selected_speed_domain_upper_m_s=float(grids[selected_index].faces_m_s[-1]),
        failure_reason=";".join(failures),
    )
