"""Warm start compatibility remapping and nonlinear iteration metrics for Eq 59"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid
from source_model_revamp.fbis.modal.eq59.collision_operator_state import Eq59CollisionOperatorState
from source_model_revamp.fbis.modal.eq59.types import Eq59WarmStartState, _Eq59WarmStartDecision, _Eq59IterationMetrics
from source_model_revamp.fbis.modal.eq59.operator import _eq59_nonlinear_residual

def _maximum_retained_mode_relative_change(previous: np.ndarray, current: np.ndarray) -> tuple[float, float, float]:
    """Return the largest retained mode relative update absolute update and solution scale"""
    difference = np.abs(current - previous)
    mode_scale = np.maximum(np.max(np.abs(previous), axis=1), np.max(np.abs(current), axis=1))
    mode_change = np.max(difference, axis=1)
    zero_scale = mode_scale == 0.0
    relative_by_mode = np.zeros_like(mode_change)
    relative_by_mode[~zero_scale] = mode_change[~zero_scale] / mode_scale[~zero_scale]
    relative_by_mode[zero_scale] = np.where(mode_change[zero_scale] == 0.0, 0.0, np.inf)

    return (float(np.max(relative_by_mode)), float(np.max(difference)), float(np.max(mode_scale)))

def _relative_scalar_change(current: float, previous: float) -> float:
    """Return a symmetric scale relative change for two finite scalars"""
    if not np.isfinite(current) or not np.isfinite(previous):
        return float("inf")
    scale = max(abs(current), abs(previous), np.finfo(float).tiny)
  
    return abs(current - previous) / scale

def _relative_array_change(current: ArrayLike, previous: ArrayLike) -> float:
    """Return the maximum scale relative change between two arrays"""
    current_array = np.asarray(current, dtype=float)
    previous_array = np.asarray(previous, dtype=float)
    if current_array.shape != previous_array.shape or np.any(~np.isfinite(current_array)) or np.any(~np.isfinite(previous_array)):
        return float("inf")
    difference = float(np.max(np.abs(current_array - previous_array)))
    scale = max(float(np.max(np.abs(current_array))), float(np.max(np.abs(previous_array))), np.finfo(float).tiny)
   
    return difference / scale

def _modal_inventory_and_effective_energy(*, speed_grid: SpeedGrid, modal_distribution: np.ndarray, mode_eta_integrals: np.ndarray | None, particle_mass_kg: float | None) -> tuple[float | None, float | None]:
    """Return volume averaged fast density and mean kinetic energy from the modal state
    
    The density measure is Σ_j ∫ I_j dη Σ_v f_j Δ(4πv³/3)
    """
    if mode_eta_integrals is None:
        return None, None
    eta_integrated_distribution = mode_eta_integrals @ modal_distribution
    shell = np.asarray(speed_grid.shell_volumes_m3_s3, dtype=float)
    inventory = float(np.sum(shell * eta_integrated_distribution))
    if particle_mass_kg is None or not np.isfinite(inventory) or inventory <= 0.0:
        return inventory, None
    kinetic_energy = (0.5 * particle_mass_kg * np.asarray(speed_grid.centers_m_s, dtype=float) ** 2)
    energy = float(np.sum(shell * eta_integrated_distribution * kinetic_energy) / inventory)

    return inventory, energy

def _validated_optional_physical_measures(*, mode_count: int, mode_eta_integrals: ArrayLike | None, physical_eigenfunctions_eta: ArrayLike | None, eta_grid: ArrayLike | None, particle_mass_kg: float | None) -> tuple[np.ndarray | None, np.ndarray | None, float | None]:
    """Validate or derive modal η integrals eigenfunctions η grid and particle mass"""
    eigenfunctions: np.ndarray | None = None
    eta: np.ndarray | None = None
    if physical_eigenfunctions_eta is not None:
        eigenfunctions = np.asarray(physical_eigenfunctions_eta, dtype=float)
        if eigenfunctions.ndim != 2 or eigenfunctions.shape[0] != mode_count:
            raise ValueError("physical_eigenfunctions_eta must have shape (n_mode, n_eta)")
        if not np.all(np.isfinite(eigenfunctions)):
            raise ValueError("physical_eigenfunctions_eta must contain only finite values")
    if eta_grid is not None:
        eta = np.asarray(eta_grid, dtype=float)
        if eta.ndim != 1 or eta.size < 2 or not np.all(np.isfinite(eta)):
            raise ValueError("eta_grid must be a finite 1D with at least two points")
        if np.any(np.diff(eta) <= 0.0):
            raise ValueError("eta_grid must be strictly increasing")
        if eigenfunctions is not None and eigenfunctions.shape[1] != eta.size:
            raise ValueError("eta_grid must match the eigenfunction eta dimension")
    elif eigenfunctions is not None and mode_eta_integrals is None:
        raise ValueError("eta_grid is required to derive mode eta integrals")
    integrals: np.ndarray | None
    if mode_eta_integrals is None:
        if eigenfunctions is not None and eta is not None:
            integrals = np.trapezoid(eigenfunctions, eta, axis=1)
        else:
            integrals = None
    else:
        integrals = np.asarray(mode_eta_integrals, dtype=float)
        if integrals.shape != (mode_count,) or not np.all(np.isfinite(integrals)):
            raise ValueError("mode_eta_integrals must contain one finite value per mode")
    mass: float | None
    if particle_mass_kg is None:
        mass = None
    else:
        mass = float(particle_mass_kg)
        if not np.isfinite(mass) or mass <= 0.0:
            raise ValueError("particle_mass_kg must be positive and finite")
        
    return integrals, eigenfunctions, mass

def build_eq59_warm_start_state(*, speed_grid: SpeedGrid, modal_distribution_f_j_v: ArrayLike, magnetic_eigenvalues: ArrayLike, active_eigenvalues_by_speed: ArrayLike | None, source_coefficients_by_component_j: ArrayLike, source_speeds_m_s: ArrayLike, collision_operator_state: Eq59CollisionOperatorState, particle_mass_kg: float | None, charge_exchange_loss_frequency_s: ArrayLike | None = None, species_label: str = "D", nonlinear_converged: bool = True, mode_eta_integrals: ArrayLike | None = None, physical_eigenfunctions_eta: ArrayLike | None = None, eta_grid: ArrayLike | None = None) -> Eq59WarmStartState:
    """Build a compatibility record from a completed Eq 59 nonlinear solve"""
    modal = np.asarray(modal_distribution_f_j_v, dtype=float).copy()
    eigenvalues = np.asarray(magnetic_eigenvalues, dtype=float).copy()
    active_eigenvalues = (np.broadcast_to(eigenvalues[:, None], modal.shape).copy() if active_eigenvalues_by_speed is None else np.asarray(active_eigenvalues_by_speed, dtype=float).copy())
    coefficients = np.asarray(source_coefficients_by_component_j, dtype=float).copy()
    speeds = np.asarray(source_speeds_m_s, dtype=float).copy()
    if modal.shape != (eigenvalues.size, speed_grid.centers_m_s.size):
        raise ValueError("modal_distribution_f_j_v has incompatible mode or speed dimensions")
    if coefficients.shape != (speeds.size, eigenvalues.size):
        raise ValueError("source coefficient shape is incompatible with components and modes")
    if active_eigenvalues.shape != modal.shape:
        raise ValueError("active_eigenvalues_by_speed must match the modal distribution")
    if any(np.any(~np.isfinite(values)) for values in (modal, eigenvalues, active_eigenvalues, coefficients, speeds)):
        raise ValueError("Eq 59 warm-start arrays must contain only finite values")
    if np.any(eigenvalues <= 0.0) or np.any(active_eigenvalues <= 0.0) or np.any(speeds <= 0.0):
        raise ValueError("Eq 59 warm-start eigenvalues and source speeds must be positive")
    tau_s = float(collision_operator_state.spitzer_slowing_down_time_s)
    if not np.isfinite(tau_s) or tau_s <= 0.0:
        raise ValueError("Eq 59 warm-start slowing time must be finite and positive")
    if collision_operator_state.test_species.symbol.strip().lower() != str(species_label).strip().lower():
        raise ValueError("Eq 59 warm-start operator species must match species_label")
    mass = None if particle_mass_kg is None else float(particle_mass_kg)
    if mass is not None and (not np.isfinite(mass) or mass <= 0.0):
        raise ValueError("particle_mass_kg must be finite and positive")
    species = str(species_label).strip()
    if not species:
        raise ValueError("species_label must not be empty")
    mode_integrals, eigenfunctions, _ = _validated_optional_physical_measures(mode_count=eigenvalues.size, mode_eta_integrals=mode_eta_integrals, physical_eigenfunctions_eta=physical_eigenfunctions_eta, eta_grid=eta_grid, particle_mass_kg=mass)
    eta = None if eta_grid is None else np.asarray(eta_grid, dtype=float).copy()
    charge_exchange = None if charge_exchange_loss_frequency_s is None else np.asarray(charge_exchange_loss_frequency_s, dtype=float).copy()
    if charge_exchange is not None and (charge_exchange.shape != speed_grid.centers_m_s.shape or np.any(~np.isfinite(charge_exchange)) or np.any(charge_exchange < 0.0)):
        raise ValueError("charge_exchange_loss_frequency_s must be finite and nonnegative on the Eq 59 speed grid")
    return Eq59WarmStartState(
        speed_faces_m_s=np.asarray(speed_grid.faces_m_s, dtype=float).copy(),
        modal_distribution_f_j_v=modal,
        magnetic_eigenvalues=eigenvalues,
        active_eigenvalues_by_speed=active_eigenvalues,
        source_coefficients_by_component_j=coefficients,
        source_speeds_m_s=speeds,
        spitzer_slowing_down_time_s=tau_s,
        collision_operator_model=collision_operator_state.operator_model,
        collision_operator_structure_fingerprint=collision_operator_state.structure_fingerprint(),
        collision_operator_closure_fingerprint=collision_operator_state.closure_fingerprint(),
        particle_mass_kg=mass,
        charge_exchange_loss_frequency_s=charge_exchange,
        species_label=species,
        nonlinear_converged=bool(nonlinear_converged),
        mode_eta_integrals=(None if mode_integrals is None else np.asarray(mode_integrals, dtype=float).copy()),
        physical_eigenfunctions_eta=(None if eigenfunctions is None else np.asarray(eigenfunctions, dtype=float).copy()),
        eta_grid=eta,
    )

def _conservative_shell_remap(source_faces_m_s: np.ndarray, source_values: np.ndarray, target_grid: SpeedGrid) -> np.ndarray:
    """Conservatively remap signed cell averages in exact spherical shell measure
    
    Each overlap uses ΔV_v = 4π(v_R³ − v_L³) / 3
    """
    source_faces = np.asarray(source_faces_m_s, dtype=float)
    values = np.asarray(source_values, dtype=float)
    target_faces = np.asarray(target_grid.faces_m_s, dtype=float)
    if source_faces.ndim != 1 or values.ndim != 2 or values.shape[1] != source_faces.size - 1:
        raise ValueError("warm start source grid and modal distribution shapes are inconsistent")
    shell_factor = 4.0 * np.pi / 3.0
    remapped_integrals = np.zeros((values.shape[0], target_faces.size - 1), dtype=float)
    for source_index, (source_left, source_right) in enumerate(zip(source_faces[:-1], source_faces[1:], strict=True)):
        left = np.maximum(target_faces[:-1], source_left)
        right = np.minimum(target_faces[1:], source_right)
        overlap = shell_factor * np.maximum(right**3 - left**3, 0.0)
        remapped_integrals += values[:, source_index, None] * overlap[None, :]
  
    return remapped_integrals / np.asarray(target_grid.shell_volumes_m3_s3, dtype=float)[None, :]

def _prepare_eq59_warm_start(*, warm_start_state: Eq59WarmStartState | None, speed_grid: SpeedGrid, source_coefficients: np.ndarray, source_speeds: np.ndarray, eigenvalues: np.ndarray, active_eigenvalues_by_speed: np.ndarray, collision_operator_state: Eq59CollisionOperatorState, particle_mass_kg: float | None, species_label: str, mode_eta_integrals: np.ndarray | None, physical_eigenfunctions_eta: np.ndarray | None, eta_grid: np.ndarray | None, charge_exchange_loss_frequency_s: np.ndarray | None, compatibility_policy: str) -> _Eq59WarmStartDecision:
    """Validate a prior Eq 59 state for reuse on the requested solve
    
    Exact physics requires matching closure inputs while closure continuation permits selected closure changes
    A larger target speed domain may receive a conservative shell remap
    """
    if warm_start_state is None:
        return _Eq59WarmStartDecision(None, False, False, "pairwise_seed_no_prior_state", "")
    state = warm_start_state
    reasons: list[str] = []
    if state.model != "egedal_eq59_hot_ion_rosenbluth":
        reasons.append("model_mismatch")
    if state.species_label.strip().lower() != str(species_label).strip().lower():
        reasons.append("species_mismatch")
    if state.modal_distribution_f_j_v.shape[0] != eigenvalues.size:
        reasons.append("retained_mode_count_mismatch")
    if state.source_coefficients_by_component_j.shape != source_coefficients.shape:
        reasons.append("source_component_shape_mismatch")
    elif not np.allclose(state.source_coefficients_by_component_j, source_coefficients, rtol=1.0e-12, atol=0.0):
        reasons.append("source_component_coefficients_mismatch")
    if state.source_speeds_m_s.shape != source_speeds.shape or not np.allclose(state.source_speeds_m_s, source_speeds, rtol=1.0e-12, atol=0.0):
        reasons.append("source_component_speeds_mismatch")
    if state.magnetic_eigenvalues.shape != eigenvalues.shape or not np.allclose(state.magnetic_eigenvalues, eigenvalues, rtol=1.0e-12, atol=0.0):
        reasons.append("magnetic_eigenbasis_mismatch")
    active_profiles_match = bool(state.active_eigenvalues_by_speed.shape == active_eigenvalues_by_speed.shape and np.allclose(state.active_eigenvalues_by_speed, active_eigenvalues_by_speed, rtol=1.0e-12, atol=0.0))
    both_profiles_are_magnetic = bool(np.allclose(state.active_eigenvalues_by_speed, state.magnetic_eigenvalues[:, None], rtol=1.0e-12, atol=0.0) and np.allclose(active_eigenvalues_by_speed, eigenvalues[:, None], rtol=1.0e-12, atol=0.0))
    if not active_profiles_match and not both_profiles_are_magnetic:
        reasons.append("active_eigenvalue_profile_mismatch")
    if state.collision_operator_model != collision_operator_state.operator_model:
        reasons.append("collision_operator_model_mismatch")
    if state.collision_operator_structure_fingerprint != collision_operator_state.structure_fingerprint():
        reasons.append("collision_operator_structure_mismatch")
    if state.collision_operator_closure_fingerprint != collision_operator_state.closure_fingerprint():
        reasons.append("collision_operator_closure_mismatch")
    previous_charge_exchange = None if state.charge_exchange_loss_frequency_s is None else np.asarray(state.charge_exchange_loss_frequency_s, dtype=float)
    current_charge_exchange = None if charge_exchange_loss_frequency_s is None else np.asarray(charge_exchange_loss_frequency_s, dtype=float)
    if (previous_charge_exchange is None) != (current_charge_exchange is None):
        reasons.append("charge_exchange_loss_frequency_mismatch")
    elif previous_charge_exchange is not None and current_charge_exchange is not None and (previous_charge_exchange.shape != current_charge_exchange.shape or not np.allclose(previous_charge_exchange, current_charge_exchange, rtol=1.0e-12, atol=0.0)):
        reasons.append("charge_exchange_loss_frequency_mismatch")
    target_mass = None if particle_mass_kg is None else float(particle_mass_kg)
    if (state.particle_mass_kg is None) != (target_mass is None):
        reasons.append("particle_mass_availability_mismatch")
    elif target_mass is not None and not np.isclose(float(state.particle_mass_kg), target_mass, rtol=1.0e-12, atol=0.0):
        reasons.append("particle_mass_mismatch")
    faces = np.asarray(state.speed_faces_m_s, dtype=float)
    if faces.ndim != 1 or faces.size != state.modal_distribution_f_j_v.shape[1] + 1 or np.any(~np.isfinite(faces)) or np.any(np.diff(faces) <= 0.0):
        reasons.append("invalid_speed_grid")
    if np.any(~np.isfinite(state.modal_distribution_f_j_v)):
        reasons.append("nonfinite_distribution")
    for name, previous, current in (("mode_eta_integrals", state.mode_eta_integrals, mode_eta_integrals), ("physical_eigenfunctions_eta", state.physical_eigenfunctions_eta, physical_eigenfunctions_eta), ("eta_grid", state.eta_grid, eta_grid)):
        if (previous is None) != (current is None):
            reasons.append(f"{name}_availability_mismatch")
        elif previous is not None and current is not None:
            previous_array = np.asarray(previous, dtype=float)
            current_array = np.asarray(current, dtype=float)
            if (previous_array.shape != current_array.shape or not np.allclose(previous_array, current_array, rtol=1.0e-12, atol=0.0)):
                reasons.append(f"{name}_mismatch")
    if not bool(state.nonlinear_converged):
        reasons.append("prior_nonlinear_state_not_converged")
    policy = str(compatibility_policy).strip().lower()
    if policy not in {"exact_physics", "closure_continuation"}:
        raise ValueError("warm_start_compatibility_policy must be exact_physics or closure_continuation")
    continuation_mismatches = {"active_eigenvalue_profile_mismatch", "collision_operator_closure_mismatch", "charge_exchange_loss_frequency_mismatch", "source_component_coefficients_mismatch"}
    continuation_used = bool(reasons)
    hard_reasons = list(reasons)
    if policy == "closure_continuation":
        hard_reasons = [reason for reason in reasons if reason not in continuation_mismatches]
    if hard_reasons:
        return _Eq59WarmStartDecision(None, False, False, "warm_start_rejected", ";".join(hard_reasons))
    target_faces = np.asarray(speed_grid.faces_m_s, dtype=float)
    if faces.shape == target_faces.shape and np.allclose(faces, target_faces, rtol=1.0e-13, atol=0.0):
        candidate = np.asarray(state.modal_distribution_f_j_v, dtype=float).copy()
        remapped = False
        status = "warm_start_used_exact_grid"
    elif np.isclose(faces[0], target_faces[0], rtol=0.0, atol=1.0e-14) and target_faces[-1] >= faces[-1]:
        candidate = _conservative_shell_remap(faces, state.modal_distribution_f_j_v, speed_grid)
        original_integrals = state.modal_distribution_f_j_v @ ((4.0 * np.pi / 3.0) * (faces[1:] ** 3 - faces[:-1] ** 3))
        remapped_integrals = candidate @ np.asarray(speed_grid.shell_volumes_m3_s3, dtype=float)
        if not np.allclose(remapped_integrals, original_integrals, rtol=5.0e-13, atol=1.0e-300):
            return _Eq59WarmStartDecision(None, False, False, "warm_start_rejected", "conservative_remap_failed")
        remapped = True
        status = "warm_start_used_conservative_speed_remap"
    else:
        return _Eq59WarmStartDecision(None, False, False, "warm_start_rejected", "speed_grid_not_nested_or_target_domain_smaller")
    if continuation_used:
        status += "_closure_continuation"

    return _Eq59WarmStartDecision(candidate, True, remapped, status, "")

def _iteration_metrics(*, speed_grid: SpeedGrid, previous: np.ndarray, current: np.ndarray, source_coefficients_by_component_j: np.ndarray, source_speeds_m_s: np.ndarray, eigenvalues_by_speed: np.ndarray, collision_operator_state: Eq59CollisionOperatorState, previous_collision_operator_state: Eq59CollisionOperatorState | None, mode_eta_integrals: np.ndarray | None, particle_mass_kg: float | None, previous_inventory: float | None, previous_effective_energy_J: float | None, charge_exchange_loss_frequency_s: np.ndarray | None = None) -> _Eq59IterationMetrics:
    """Evaluate solution residual density energy and collision operator changes for one accepted iterate"""
    relative_change, absolute_change, solution_scale = _maximum_retained_mode_relative_change(previous, current)
    absolute_residual, relative_residual = _eq59_nonlinear_residual(speed_grid=speed_grid, modal_distribution=current, source_coefficients_by_component_j=source_coefficients_by_component_j, source_speeds_m_s=source_speeds_m_s, eigenvalues=eigenvalues_by_speed, collision_operator_state=collision_operator_state, charge_exchange_loss_frequency_s=charge_exchange_loss_frequency_s)
    inventory, effective_energy_J = _modal_inventory_and_effective_energy(speed_grid=speed_grid, modal_distribution=current, mode_eta_integrals=mode_eta_integrals, particle_mass_kg=particle_mass_kg)
    inventory_change = (None if inventory is None or previous_inventory is None else _relative_scalar_change(inventory, previous_inventory))
    energy_change = (None if effective_energy_J is None or previous_effective_energy_J is None else _relative_scalar_change(effective_energy_J, previous_effective_energy_J))
    fast_self = next((contribution for contribution in collision_operator_state.pair_contributions_by_id.values() if contribution.field_population_id == "fast_self"), None)
    fast_self_g_scale = 0.0 if fast_self is None else float(fast_self.g_scale_velocity_cubed_m3_s3)
    fast_self_drag_scale = 0.0 if fast_self is None else float(fast_self.net_drag_scale_velocity_cubed_m3_s3)
    if previous_collision_operator_state is None:
        drag_change = None
        diffusion_change = None
        pitch_change = None
    else:
        drag_change = _relative_array_change(collision_operator_state.ion_drag_velocity_cubed_m3_s3, previous_collision_operator_state.ion_drag_velocity_cubed_m3_s3)
        diffusion_change = _relative_array_change(collision_operator_state.ion_energy_diffusion_velocity_fourth_m4_s4, previous_collision_operator_state.ion_energy_diffusion_velocity_fourth_m4_s4)
        pitch_change = _relative_array_change(collision_operator_state.ion_pitch_scattering_velocity_cubed_m3_s3, previous_collision_operator_state.ion_pitch_scattering_velocity_cubed_m3_s3)
  
    return _Eq59IterationMetrics(
        relative_change=relative_change,
        absolute_change=absolute_change,
        solution_scale=solution_scale,
        relative_residual=relative_residual,
        absolute_residual=absolute_residual,
        inventory=inventory,
        inventory_relative_change=inventory_change,
        effective_energy_J=effective_energy_J,
        effective_energy_relative_change=energy_change,
        fast_self_density_m3=float(collision_operator_state.fast_self_physical_density_m3),
        fast_self_density_relative_change=inventory_change,
        fast_self_g_scale_m3_s3=fast_self_g_scale,
        fast_self_net_drag_scale_m3_s3=fast_self_drag_scale,
        ion_drag_operator_relative_change=drag_change,
        ion_energy_diffusion_operator_relative_change=diffusion_change,
        ion_pitch_scattering_operator_relative_change=pitch_change,
    )

def _fixed_point_pattern_flags(relative_changes: list[float], residuals: list[float], update_cosines: list[float], relative_tolerance: float) -> tuple[bool, bool, bool]:
    """Detect sustained stagnation oscillation or residual growth in recent nonlinear iterations"""
    stagnated = False
    oscillatory = False
    diverged = False
    if len(residuals) >= 4:
        recent_residuals = np.asarray(residuals[-4:], dtype=float)
        recent_changes = np.asarray(relative_changes[-4:], dtype=float)
        if np.all(np.isfinite(recent_residuals)) and np.all(np.isfinite(recent_changes)):
            residual_floor = max(relative_tolerance, np.finfo(float).eps)
            residual_plateau = (float(np.min(recent_residuals)) > 5.0 * residual_floor and float(np.max(recent_residuals)) <= 1.02 * float(np.min(recent_residuals)))
            change_plateau = (float(np.min(recent_changes)) > 5.0 * relative_tolerance and float(np.max(recent_changes)) <= 1.02 * float(np.min(recent_changes)))
            stagnated = bool(residual_plateau and change_plateau)
            diverged = bool(recent_residuals[-1] > max(1.0, 100.0 * relative_tolerance) and np.all(recent_residuals[1:] > 1.5 * recent_residuals[:-1]))
    if len(update_cosines) >= 3 and len(relative_changes) >= 4:
        recent_cosines = np.asarray(update_cosines[-3:], dtype=float)
        recent_changes = np.asarray(relative_changes[-4:], dtype=float)
        oscillatory = bool(np.all(recent_cosines < -0.9) and recent_changes[-1] > 5.0 * relative_tolerance and recent_changes[-1] >= 0.8 * recent_changes[-2])
        
    return stagnated, oscillatory, diverged
