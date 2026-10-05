"""Particle loss current and power balance assembly for one modal fast ion species"""
from __future__ import annotations
import numpy as np
from scipy.constants import elementary_charge
from source_model_revamp.electrostatic.current_balance import AmbipolarCurrentBalance, build_ion_current_loss, egedal_tail_refilling_current_balance_for_actual_midplane_density
from source_model_revamp.fbis.modal.local.mapping import _density_from_v_eta, _effective_temperature_from_v_eta, _modal_distribution_to_lambda_grid
from source_model_revamp.fbis.modal.utils import _EPS
from collections.abc import Mapping
from typing import Any
from source_model_revamp.fbis.modal.solver.solver_loss_rates import _conservative_eq61_mode_slopes_dI_dlambda, _eq61_boundary_reconstruction_diagnostics, _ion_boundary_flux_from_modal, _modal_eigenvalue_sink_audit, _select_global_loss_rate
from source_model_revamp.fbis.modal.solver.solver_convergence import _projected_modal_source_rate_from_coefficients, _safe_ratio

def solve_modal_loss_and_current_balance(context: Mapping[str, Any]) -> dict[str, Any]:
    """Assemble Eq 61 ion losses, Eq 68 electron current balance, confinement, and consistency diagnostics"""
    active_eigenvalues_by_speed = context['active_eigenvalues_by_speed']
    active_global_loss_rate_model = context['active_global_loss_rate_model']
    auxiliary_input_power_W = context['auxiliary_input_power_W']
    basis = context['basis']
    collision_state = context['collision_state']
    component_results = context['component_results']
    electron_midplane_density_m3 = context['electron_midplane_density_m3']
    electron_collision_density_m3 = context['electron_collision_density_m3']
    electron_temperature_J = context['electron_temperature_J']
    eq59_collision_operator_state = context['eq59_collision_operator_state']
    f_eta = context['f_eta']
    fast_ion_species = context['fast_ion_species']
    half_length_m = context['half_length_m']
    lambda_grid = context['lambda_grid']
    midplane_area_m2 = context['midplane_area_m2']
    mirror_ratio = context['mirror_ratio']
    modal_prompt_loss_birth_power_W = context['modal_prompt_loss_birth_power_W']
    modal_prompt_loss_birth_rate_s = context['modal_prompt_loss_birth_rate_s']
    modal_source_birth_rate_s = context['modal_source_birth_rate_s']
    modal_total = context['modal_total']
    modal_total_birth_current_A = context['modal_total_birth_current_A']
    modal_total_birth_rate_s = context['modal_total_birth_rate_s']
    modal_total_deposited_birth_power_W = context['modal_total_deposited_birth_power_W']
    numerics = context['numerics']
    published_heuristic = context['published_heuristic']
    speed_grid = context['speed_grid']
    volume_m3 = context['volume_m3']
    f_lambda = _modal_distribution_to_lambda_grid(speed_grid=speed_grid, lambda_grid=lambda_grid, basis=basis, distribution_v_eta=f_eta)
    density = _density_from_v_eta(speed_grid, basis.physical_basis.eta_grid, f_eta)
    effective_ion_temperature_J = _effective_temperature_from_v_eta(speed_grid, basis.physical_basis.eta_grid, f_eta, fast_ion_species.mass_kg)
    inventory = density * float(volume_m3)
    dfdlambda, boundary_flux_single_density_v, boundary_flux_single_loss_rate_density, boundary_flux_single_midplane_power_density, mean_loss_energy = _ion_boundary_flux_from_modal(speed_grid=speed_grid, basis=basis, modal_f_j_v=modal_total, collision_state=collision_state, fast_ion_mass_kg=fast_ion_species.mass_kg, collision_operator_state=eq59_collision_operator_state, species_id=fast_ion_species.species_id)
    eq61_boundary_reconstruction_diagnostics = _eq61_boundary_reconstruction_diagnostics(
        speed_grid=speed_grid,
        basis=basis,
        modal_f_j_v=modal_total,
        collision_state=collision_state,
        fast_ion_mass_kg=fast_ion_species.mass_kg,
        collision_operator_state=eq59_collision_operator_state,
        species_id=fast_ion_species.species_id,
    )
    # Eq 61 reconstructs one directed end before the symmetric device total is formed
    boundary_flux_single_particle_loss_rate_s = boundary_flux_single_loss_rate_density * float(volume_m3)
    boundary_flux_single_midplane_power_loss_W = boundary_flux_single_midplane_power_density * float(volume_m3)
    boundary_flux_two_boundary_particle_loss_rate_s = 2.0 * boundary_flux_single_particle_loss_rate_s
    boundary_flux_two_boundary_midplane_power_loss_W = 2.0 * boundary_flux_single_midplane_power_loss_W
    projected_source_particle_rate_s, projected_source_component_rates_s = _projected_modal_source_rate_from_coefficients(basis=basis, component_results=component_results, volume_m3=volume_m3)
    eigen_audit = _modal_eigenvalue_sink_audit(speed_grid=speed_grid, basis=basis, modal_f_j_v=modal_total, collision_state=collision_state, fast_ion_mass_kg=fast_ion_species.mass_kg, collision_operator_state=eq59_collision_operator_state, eigenvalues_by_speed=active_eigenvalues_by_speed)
    eigenvalue_sink_particle_rate_s = float(eigen_audit["total_particle_rate_density_m3_s"]) * float(volume_m3)
    eigenvalue_sink_current_A = float(elementary_charge * fast_ion_species.charge_number * eigenvalue_sink_particle_rate_s)
    eigenvalue_sink_midplane_power_W = float(eigen_audit["total_midplane_power_density_W_m3"]) * float(volume_m3)
    eigenvalue_sink_mean_energy_J = float(eigen_audit["mean_energy_J"])
    eigenvalue_first_mode_particle_rate_s = float(eigen_audit["first_mode_particle_rate_density_m3_s"]) * float(volume_m3)
    eigenvalue_first_mode_current_A = float(elementary_charge * fast_ion_species.charge_number * eigenvalue_first_mode_particle_rate_s)
    eigenvalue_first_mode_midplane_power_W = float(eigen_audit["first_mode_midplane_power_density_W_m3"]) * float(volume_m3)
    eigenvalue_first_mode_mean_energy_J = float(eigen_audit["first_mode_mean_energy_J"])
    eigenvalue_mode_particle_rates_s = [float(x * float(volume_m3)) for x in np.asarray(eigen_audit["mode_particle_rate_density_m3_s"], dtype=float)]
    eigenvalue_mode_midplane_power_W = [float(x * float(volume_m3)) for x in np.asarray(eigen_audit["mode_midplane_power_density_W_m3"], dtype=float)]
    active_loss = _select_global_loss_rate(model=active_global_loss_rate_model, source_particle_rate_s=modal_source_birth_rate_s, eigenvalue_sink_particle_rate_s=eigenvalue_sink_particle_rate_s, eigenvalue_sink_midplane_power_W=eigenvalue_sink_midplane_power_W, boundary_single_particle_rate_s=boundary_flux_single_particle_loss_rate_s, boundary_single_midplane_power_W=boundary_flux_single_midplane_power_loss_W, boundary_mean_loss_energy_J=mean_loss_energy, volume_m3=volume_m3, boundary_flux_density_v=boundary_flux_single_density_v)
    confined_total_device_loss_rate_s = float(active_loss["particle_loss_rate_s"])
    confined_total_device_midplane_power_W = float(active_loss["midplane_power_loss_W"])
    # Prompt magnetic losses bypass the confined modal inventory and are added after the collisional loss closure
    ion_particle_loss_rate_s = confined_total_device_loss_rate_s + modal_prompt_loss_birth_rate_s
    ion_midplane_power_loss_W = confined_total_device_midplane_power_W + modal_prompt_loss_birth_power_W
    loss_rate_density = ion_particle_loss_rate_s / max(float(volume_m3), _EPS)
    flux = np.asarray(active_loss["boundary_flux_density_v"], dtype=float)
    conservative_mode_slopes = _conservative_eq61_mode_slopes_dI_dlambda(basis)
    direct_cell_average_mode_slopes = np.asarray(basis.mode_slopes_dI_dlambda_at_boundary, dtype=float)
    first_slope = float(conservative_mode_slopes[0])
    if first_slope <= 0.0 or not np.isfinite(first_slope):
        raise ValueError("first modal loss cone slope must be positive and finite")
    ion_current_loss = build_ion_current_loss(ion_particle_loss_rates_s=np.asarray([ion_particle_loss_rate_s], dtype=float), ion_charge_numbers=np.asarray([fast_ion_species.charge_number], dtype=float))
    electron_balance = egedal_tail_refilling_current_balance_for_actual_midplane_density(target_current_A=float(ion_current_loss.total_ion_current_A), electron_midplane_density_m3=electron_midplane_density_m3, electron_collision_density_m3=electron_collision_density_m3, electron_temperature_J=electron_temperature_J, midplane_area_m2=midplane_area_m2, half_length_m=half_length_m, mirror_ratio=mirror_ratio, loss_geometry_factor=basis.geometry_factor_G, loss_cone_slope_dI_dlambda=first_slope, electron_collision_frequency_s=collision_state.electron_electron_collision_frequency_s, number_of_ends=2.0)
    current_balance = AmbipolarCurrentBalance(ion_current_loss=ion_current_loss, electron_balance=electron_balance, current_residual_A=electron_balance.current_residual_A, electron_wall_potential_relative_to_midplane_V=electron_balance.wall_potential_relative_to_midplane_V)
    barrier = float(np.asarray(current_balance.electron_balance.barrier_energy_J, dtype=float))
    if ion_particle_loss_rate_s <= 0.0:
        ion_barrier_power_W = 0.0
    elif np.isfinite(barrier):
        ion_barrier_power_W = ion_particle_loss_rate_s * barrier
    else:
        ion_barrier_power_W = float("inf")
    ion_wall_power_W = ion_midplane_power_loss_W + ion_barrier_power_W
    electron_power_W = float(np.asarray(current_balance.electron_balance.electron_wall_power_loss_W, dtype=float))
    confinement = inventory / ion_particle_loss_rate_s if ion_particle_loss_rate_s > 0.0 else None
    confined_collisional_confinement = (inventory / confined_total_device_loss_rate_s if confined_total_device_loss_rate_s > 0.0 else None)
    ion_loss_current_A = float(elementary_charge * fast_ion_species.charge_number * ion_particle_loss_rate_s)
    electron_loss_current_A = float(np.asarray(current_balance.electron_balance.electron_current_A, dtype=float))
    electron_current_residual_A = float(np.asarray(current_balance.electron_balance.current_residual_A, dtype=float))
    electron_current_relative_error = electron_current_residual_A / max(abs(ion_loss_current_A), 1.0e-300)
    electron_loss_prefactor_s = float(np.asarray(current_balance.electron_balance.electron_loss_prefactor_s, dtype=float))
    left_ion_loss_rate_s = 0.5 * ion_particle_loss_rate_s
    right_ion_loss_rate_s = 0.5 * ion_particle_loss_rate_s
    left_electron_loss_rate_s = 0.5 * float(np.asarray(current_balance.electron_balance.electron_particle_loss_rate_s, dtype=float))
    right_electron_loss_rate_s = left_electron_loss_rate_s
    target_electron_particle_loss_rate_s = ion_particle_loss_rate_s
    prefactor_over_target = electron_loss_prefactor_s / max(target_electron_particle_loss_rate_s, 1.0e-300)
    normalized_barrier = float(barrier / electron_temperature_J) if electron_temperature_J > 0.0 else float("nan")
    tail_factor = float(np.exp(-normalized_barrier)) if np.isfinite(normalized_barrier) else 0.0
    y_exp_y = float(normalized_barrier * np.exp(normalized_barrier)) if np.isfinite(normalized_barrier) and normalized_barrier < 700.0 else float("inf")
    source_particle_rate = max(modal_total_birth_rate_s, 0.0)
    loss_particle_rate = max(ion_particle_loss_rate_s, 0.0)
    beam_source_power = max(modal_total_deposited_birth_power_W, 0.0)
    auxiliary_power = float(auxiliary_input_power_W)
    if not np.isfinite(auxiliary_power) or auxiliary_power < 0.0:
        raise ValueError("auxiliary_input_power_W must be finite and nonnegative")
    total_input_power = beam_source_power + auxiliary_power
    source_mean_energy = beam_source_power / source_particle_rate if source_particle_rate > 0.0 else 0.0
    loss_equivalent_birth_power_W = loss_particle_rate * source_mean_energy
    source_to_loss_particle_ratio = source_particle_rate / loss_particle_rate if loss_particle_rate > 0.0 else None
    loss_to_source_particle_ratio = loss_particle_rate / source_particle_rate if source_particle_rate > 0.0 else None
    source_to_loss_current_ratio = modal_total_birth_current_A / ion_loss_current_A if ion_loss_current_A > 0.0 else None
    loss_to_source_current_ratio = ion_loss_current_A / modal_total_birth_current_A if modal_total_birth_current_A > 0.0 else None
    total_wall_power_W = ion_wall_power_W + electron_power_W
    source_power_to_loss_power_ratio = total_input_power / total_wall_power_W if total_wall_power_W > 0.0 else None
    loss_power_to_source_power_ratio = total_wall_power_W / total_input_power if total_input_power > 0.0 else None
    particle_balance_relative_error = (loss_particle_rate - source_particle_rate) / source_particle_rate if source_particle_rate > 0.0 else None
    projected_source_to_birth_ratio = _safe_ratio(projected_source_particle_rate_s, modal_source_birth_rate_s)
    eigenvalue_sink_to_birth_ratio = _safe_ratio(eigenvalue_sink_particle_rate_s, modal_source_birth_rate_s)
    two_end_eq61_to_birth_ratio = _safe_ratio(boundary_flux_two_boundary_particle_loss_rate_s, modal_source_birth_rate_s)
    two_end_eq61_to_eigenvalue_ratio = _safe_ratio(boundary_flux_two_boundary_particle_loss_rate_s, eigenvalue_sink_particle_rate_s)
    current_balance_relative_error = electron_current_relative_error
    power_balance_relative_error = (total_wall_power_W - total_input_power) / total_input_power if total_input_power > 0.0 else None
    loss_validation_tolerance = float(numerics.loss_convention_relative_tolerance)
    loss_validation_failures: list[str] = []
    ratio_checks = [("projected_source_to_birth", projected_source_to_birth_ratio), ("eigenvalue_sink_to_birth", eigenvalue_sink_to_birth_ratio)]
    if published_heuristic is None:
        ratio_checks.extend([("two_end_eq61_to_birth", two_end_eq61_to_birth_ratio), ("two_end_eq61_to_eigenvalue", two_end_eq61_to_eigenvalue_ratio),])
    for name, ratio in ratio_checks:
        if ratio is None or abs(float(ratio) - 1.0) > loss_validation_tolerance:
            loss_validation_failures.append(f"{name}_outside_tolerance")
    for name, residual in (("particle_balance", particle_balance_relative_error), ("current_balance", current_balance_relative_error)):
        if residual is None or abs(float(residual)) > loss_validation_tolerance:
            loss_validation_failures.append(f"{name}_outside_tolerance")
    loss_convention_validation_passed = not loss_validation_failures
    if loss_particle_rate > 0.0 and source_mean_energy > 0.0:
        modal_per_lost_ion_midplane_energy_over_birth_energy = ion_midplane_power_loss_W / (loss_particle_rate * source_mean_energy)
        modal_per_lost_ion_barrier_energy_over_birth_energy = ion_barrier_power_W / (loss_particle_rate * source_mean_energy)
        modal_per_lost_ion_wall_energy_over_birth_energy = ion_wall_power_W / (loss_particle_rate * source_mean_energy)
        modal_per_lost_ion_electron_energy_over_birth_energy = electron_power_W / (loss_particle_rate * source_mean_energy)
        modal_per_lost_ion_total_wall_energy_over_birth_energy = total_wall_power_W / (loss_particle_rate * source_mean_energy)
    else:
        modal_per_lost_ion_midplane_energy_over_birth_energy = None
        modal_per_lost_ion_barrier_energy_over_birth_energy = None
        modal_per_lost_ion_wall_energy_over_birth_energy = None
        modal_per_lost_ion_electron_energy_over_birth_energy = None
        modal_per_lost_ion_total_wall_energy_over_birth_energy = None
  
    return {
        'active_loss': active_loss,
        'auxiliary_power': auxiliary_power,
        'barrier': barrier,
        'boundary_flux_single_density_v': boundary_flux_single_density_v,
        'boundary_flux_single_midplane_power_loss_W': boundary_flux_single_midplane_power_loss_W,
        'boundary_flux_single_particle_loss_rate_s': boundary_flux_single_particle_loss_rate_s,
        'boundary_flux_two_boundary_midplane_power_loss_W': boundary_flux_two_boundary_midplane_power_loss_W,
        'boundary_flux_two_boundary_particle_loss_rate_s': boundary_flux_two_boundary_particle_loss_rate_s,
        'confined_collisional_confinement': confined_collisional_confinement,
        'confined_total_device_loss_rate_s': confined_total_device_loss_rate_s,
        'confined_total_device_midplane_power_W': confined_total_device_midplane_power_W,
        'confinement': confinement,
        'conservative_mode_slopes': conservative_mode_slopes,
        'current_balance': current_balance,
        'current_balance_relative_error': current_balance_relative_error,
        'density': density,
        'dfdlambda': dfdlambda,
        'direct_cell_average_mode_slopes': direct_cell_average_mode_slopes,
        'eq61_boundary_reconstruction_diagnostics': eq61_boundary_reconstruction_diagnostics,
        'effective_ion_temperature_J': effective_ion_temperature_J,
        'eigen_audit': eigen_audit,
        'eigenvalue_first_mode_current_A': eigenvalue_first_mode_current_A,
        'eigenvalue_first_mode_mean_energy_J': eigenvalue_first_mode_mean_energy_J,
        'eigenvalue_first_mode_midplane_power_W': eigenvalue_first_mode_midplane_power_W,
        'eigenvalue_first_mode_particle_rate_s': eigenvalue_first_mode_particle_rate_s,
        'eigenvalue_mode_midplane_power_W': eigenvalue_mode_midplane_power_W,
        'eigenvalue_mode_particle_rates_s': eigenvalue_mode_particle_rates_s,
        'eigenvalue_sink_current_A': eigenvalue_sink_current_A,
        'eigenvalue_sink_mean_energy_J': eigenvalue_sink_mean_energy_J,
        'eigenvalue_sink_midplane_power_W': eigenvalue_sink_midplane_power_W,
        'eigenvalue_sink_particle_rate_s': eigenvalue_sink_particle_rate_s,
        'eigenvalue_sink_to_birth_ratio': eigenvalue_sink_to_birth_ratio,
        'electron_current_relative_error': electron_current_relative_error,
        'electron_current_residual_A': electron_current_residual_A,
        'electron_loss_current_A': electron_loss_current_A,
        'electron_loss_prefactor_s': electron_loss_prefactor_s,
        'electron_power_W': electron_power_W,
        'f_lambda': f_lambda,
        'first_slope': first_slope,
        'flux': flux,
        'inventory': inventory,
        'ion_barrier_power_W': ion_barrier_power_W,
        'ion_loss_current_A': ion_loss_current_A,
        'ion_midplane_power_loss_W': ion_midplane_power_loss_W,
        'ion_particle_loss_rate_s': ion_particle_loss_rate_s,
        'ion_wall_power_W': ion_wall_power_W,
        'left_electron_loss_rate_s': left_electron_loss_rate_s,
        'left_ion_loss_rate_s': left_ion_loss_rate_s,
        'loss_convention_validation_passed': loss_convention_validation_passed,
        'loss_equivalent_birth_power_W': loss_equivalent_birth_power_W,
        'loss_power_to_source_power_ratio': loss_power_to_source_power_ratio,
        'loss_rate_density': loss_rate_density,
        'loss_to_source_current_ratio': loss_to_source_current_ratio,
        'loss_to_source_particle_ratio': loss_to_source_particle_ratio,
        'loss_validation_failures': loss_validation_failures,
        'loss_validation_tolerance': loss_validation_tolerance,
        'mean_loss_energy': mean_loss_energy,
        'modal_per_lost_ion_barrier_energy_over_birth_energy': modal_per_lost_ion_barrier_energy_over_birth_energy,
        'modal_per_lost_ion_electron_energy_over_birth_energy': modal_per_lost_ion_electron_energy_over_birth_energy,
        'modal_per_lost_ion_midplane_energy_over_birth_energy': modal_per_lost_ion_midplane_energy_over_birth_energy,
        'modal_per_lost_ion_total_wall_energy_over_birth_energy': modal_per_lost_ion_total_wall_energy_over_birth_energy,
        'modal_per_lost_ion_wall_energy_over_birth_energy': modal_per_lost_ion_wall_energy_over_birth_energy,
        'normalized_barrier': normalized_barrier,
        'particle_balance_relative_error': particle_balance_relative_error,
        'power_balance_relative_error': power_balance_relative_error,
        'prefactor_over_target': prefactor_over_target,
        'projected_source_component_rates_s': projected_source_component_rates_s,
        'projected_source_particle_rate_s': projected_source_particle_rate_s,
        'projected_source_to_birth_ratio': projected_source_to_birth_ratio,
        'right_electron_loss_rate_s': right_electron_loss_rate_s,
        'right_ion_loss_rate_s': right_ion_loss_rate_s,
        'source_mean_energy': source_mean_energy,
        'source_power_to_loss_power_ratio': source_power_to_loss_power_ratio,
        'source_to_loss_current_ratio': source_to_loss_current_ratio,
        'source_to_loss_particle_ratio': source_to_loss_particle_ratio,
        'tail_factor': tail_factor,
        'target_electron_particle_loss_rate_s': target_electron_particle_loss_rate_s,
        'total_input_power': total_input_power,
        'total_wall_power_W': total_wall_power_W,
        'two_end_eq61_to_birth_ratio': two_end_eq61_to_birth_ratio,
        'two_end_eq61_to_eigenvalue_ratio': two_end_eq61_to_eigenvalue_ratio,
        'y_exp_y': y_exp_y,
    }
