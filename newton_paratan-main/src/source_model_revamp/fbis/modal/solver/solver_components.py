"""Beam source projection and velocity solve assembly for one modal fast ion species"""
from __future__ import annotations
import numpy as np
from scipy.constants import elementary_charge
from source_model_revamp.fbis.modal.basis import assess_retained_mode_reconstruction_convergence
from source_model_revamp.fbis.modal.electrostatic_feedback import PublishedLambda1Heuristic, egedal_2022_published_low_energy_lambda1_heuristic
from source_model_revamp.fbis.modal.local.mapping import _modal_distribution_to_eta
from source_model_revamp.fbis.modal.models import COLD_ION_EQ14, HOT_ION_ROSENBLUTH_EQ59, velocity_solution_model
from source_model_revamp.fbis.modal.source_projection import _component_modal_rhs_from_spatial_birth, _solve_modal_source_coefficients
from source_model_revamp.fbis.modal.types import ModalBeamComponentResult, ModalPhysicalDistributionError, ModalRosenbluthCoefficients
from source_model_revamp.fbis.modal.velocity_solve_cold import cold_ion_modal_distribution_eq14
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid, _distribution_nonnegativity_metrics
from source_model_revamp.fbis.modal.eq59.types import Eq59ConvergenceDiagnostics, Eq59SpeedConvergenceDiagnostics
from source_model_revamp.fbis.modal.eq59.collision_operator_state import Eq59CollisionOperatorState
from source_model_revamp.fbis.modal.eq59.convergence import assess_eq59_speed_convergence, build_eq59_speed_convergence_grid_sequence
from source_model_revamp.fbis.modal.eq59.solve import hot_ion_rosenbluth_modal_distribution_eq59
from collections.abc import Mapping
from typing import Any
from source_model_revamp.fbis.modal.solver.solver_loss_rates import _modal_mode_eta_integrals
from source_model_revamp.fbis.modal.solver.solver_convergence import _upper_speed_tail_energy_fraction, _upper_speed_tail_population_fraction

def _correct_reconstructed_distribution(*, distribution: np.ndarray, speed_grid: SpeedGrid, eta_grid: np.ndarray, particle_mass_kg: float, particle_tolerance: float, energy_tolerance: float) -> tuple[np.ndarray, dict[str, float | bool]]:
    """Remove roundoff scale negative reconstruction cells while preserving the signed particle inventory"""
    f = np.asarray(distribution, dtype=float)
    eta = np.asarray(eta_grid, dtype=float)
    scale, _, minimum, negative_cell_fraction, _ = _distribution_nonnegativity_metrics(f, name="modal physical reconstruction")
    if f.ndim != 2 or eta.ndim != 1 or f.shape != (speed_grid.centers_m_s.size, eta.size) or eta.size < 2 or np.any(np.diff(eta) <= 0.0):
        raise ValueError("distribution and eta grid shapes are inconsistent")
    if not np.isfinite(particle_mass_kg) or particle_mass_kg <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    if not np.isfinite(particle_tolerance) or particle_tolerance < 0.0 or not np.isfinite(energy_tolerance) or energy_tolerance < 0.0:
        raise ValueError("reconstruction materiality tolerances must be finite and nonnegative")
    positive = np.maximum(f, 0.0)
    negative = np.maximum(-f, 0.0)
    positive_eta = np.trapezoid(positive, eta, axis=1)
    negative_eta = np.trapezoid(negative, eta, axis=1)
    particle_measure = np.asarray(speed_grid.shell_volumes_m3_s3, dtype=float)
    energy_measure = particle_measure * (0.5 * float(particle_mass_kg) * np.asarray(speed_grid.centers_m_s, dtype=float) ** 2)
    positive_particles = float(np.sum(particle_measure * positive_eta))
    negative_particles = float(np.sum(particle_measure * negative_eta))
    positive_energy = float(np.sum(energy_measure * positive_eta))
    negative_energy = float(np.sum(energy_measure * negative_eta))
    if positive_particles <= 0.0:
        if negative_particles > 0.0:
            raise ModalPhysicalDistributionError("modal reconstruction has no positive particle inventory")
        return f.copy(), {"applied": False, "material": False, "relative_minimum": 0.0, "negative_cell_fraction": 0.0, "negative_particle_fraction": 0.0, "energy_correction_relative": 0.0, "renormalization_factor": 1.0}
    signed_particles = positive_particles - negative_particles
    signed_energy = positive_energy - negative_energy
    if signed_particles <= 0.0 or signed_energy <= 0.0:
        raise ModalPhysicalDistributionError("modal reconstruction has nonpositive signed particle or energy inventory")
    renormalization = signed_particles / positive_particles
    corrected = positive * renormalization
    negative_particle_fraction = negative_particles / positive_particles
    energy_correction = abs(renormalization * positive_energy - signed_energy) / signed_energy
    return corrected, {"applied": bool(negative_particles > 0.0), "material": bool(negative_particle_fraction > float(particle_tolerance) or energy_correction > float(energy_tolerance)), "relative_minimum": minimum / scale, "negative_cell_fraction": negative_cell_fraction, "negative_particle_fraction": negative_particle_fraction, "energy_correction_relative": energy_correction, "renormalization_factor": renormalization}

def solve_modal_components_and_velocity(context: Mapping[str, Any]) -> dict[str, Any]:
    """Project beam births, solve Eq 14 or Eq 59 in speed, reconstruct f(v, η), and assess convergence"""
    _eq59_warm_start_compatibility_policy = context['_eq59_warm_start_compatibility_policy']
    _eq59_warm_start_state = context['_eq59_warm_start_state']
    _published_reference_lambda1 = context['_published_reference_lambda1']
    _published_throat_drop_override_J = context['_published_throat_drop_override_J']
    _run_numerical_convergence_assessments = context['_run_numerical_convergence_assessments']
    attenuated_source = context['attenuated_source']
    basis = context['basis']
    collision_state = context['collision_state']
    electron_temperature_J = context['electron_temperature_J']
    feedback_model = context['feedback_model']
    fast_ion_species = context['fast_ion_species']
    fixed_boundary_production = context['fixed_boundary_production']
    mirror_ratio = context['mirror_ratio']
    numerics = context['numerics']
    published_feedback_names = context['published_feedback_names']
    speed_grid = context['speed_grid']
    volume_m3 = context['volume_m3']
    charge_exchange_sink_state = context.get('_charge_exchange_sink_state')
    external_fast_field_states = tuple(context.get('_external_fast_field_collision_states') or ())
    charge_exchange_loss_frequency_s = None if charge_exchange_sink_state is None else np.asarray(charge_exchange_sink_state.loss_frequency_s, dtype=float)
    charge_exchange_sink_active = bool(charge_exchange_sink_state is not None and float(charge_exchange_sink_state.reference_event_rate_s) > 0.0)
    charge_exchange_sink_rate_s = 0.0 if charge_exchange_sink_state is None else float(charge_exchange_sink_state.reference_event_rate_s)
    charge_exchange_sink_energy_W = 0.0 if charge_exchange_sink_state is None else float(charge_exchange_sink_state.reference_energy_removal_W)
    charge_exchange_sink_model = "inactive" if charge_exchange_sink_state is None else str(charge_exchange_sink_state.model)
    charge_exchange_sink_reference_identity_error = 0.0 if charge_exchange_sink_state is None else float(charge_exchange_sink_state.reference_rate_identity_relative_error)
    n_modes = basis.physical_basis.eigenvalues.size
    component_results: list[ModalBeamComponentResult] = []
    component_coeffs: list[np.ndarray] = []
    component_speeds: list[float] = []
    cold_components: list[np.ndarray] = []
    for comp in attenuated_source.component_sources:
        rhs, active_cells, lambda_mean, eta_mean, first_mode_mean, confined_birth_rate_s = _component_modal_rhs_from_spatial_birth(component=comp, basis=basis, volume_m3=volume_m3)
        coeff = _solve_modal_source_coefficients(rhs, basis.physical_basis)
        cold_component = cold_ion_modal_distribution_eq14(speed_grid=speed_grid, source_coefficients_j=coeff, source_speed_m_s=float(comp.spatial_source.speed_m_s), eigenvalues=basis.physical_basis.eigenvalues, spitzer_slowing_down_time_s=collision_state.spitzer_slowing_down_time_s, critical_velocity_m_s=collision_state.critical_velocity_m_s, beta_m=collision_state.beta_m)
        component_coeffs.append(coeff)
        component_speeds.append(float(comp.spatial_source.speed_m_s))
        cold_components.append(cold_component)
        # Births outside the confined η interval remain prompt losses and do not enter the modal right hand side
        total_component_birth_rate_s = float(comp.birth_rate_s)
        prompt_difference = total_component_birth_rate_s - confined_birth_rate_s
        prompt_roundoff_tolerance = 1.0e-12 * max(total_component_birth_rate_s, 1.0)
        if abs(prompt_difference) <= prompt_roundoff_tolerance:
            confined_birth_rate_s = total_component_birth_rate_s
            prompt_loss_birth_rate_s = 0.0
        elif prompt_difference < 0.0:
            raise ValueError("projected confined component birth exceeds the physical component birth rate")
        else:
            prompt_loss_birth_rate_s = prompt_difference
        confined_birth_power_W = float(confined_birth_rate_s * comp.spatial_source.energy_J)
        prompt_loss_birth_power_W = float(prompt_loss_birth_rate_s * comp.spatial_source.energy_J)
        component_results.append(ModalBeamComponentResult(
            energy_J=float(comp.spatial_source.energy_J),
            speed_m_s=float(comp.spatial_source.speed_m_s),
            birth_rate_s=float(comp.birth_rate_s),
            deposited_birth_power_W=float(comp.deposited_birth_power_W),
            confined_birth_rate_s=float(confined_birth_rate_s),
            confined_birth_power_W=confined_birth_power_W,
            prompt_loss_birth_rate_s=prompt_loss_birth_rate_s,
            prompt_loss_birth_power_W=prompt_loss_birth_power_W,
            modal_rhs=rhs,
            modal_source_coefficients=coeff,
            modal_distribution_f_j_v=cold_component,
            contributing_axial_cells=active_cells,
            birth_lambda_weighted_mean=lambda_mean,
            birth_eta_weighted_mean=eta_mean,
            first_mode_weighted_mean=first_mode_mean,
        ))
    modal_total_birth_rate_s = float(sum(c.birth_rate_s for c in component_results))
    modal_total_deposited_birth_power_W = float(sum(c.deposited_birth_power_W for c in component_results))
    modal_source_birth_rate_s = float(sum(c.confined_birth_rate_s for c in component_results))
    modal_source_deposited_birth_power_W = float(sum(c.confined_birth_power_W for c in component_results))
    modal_prompt_loss_birth_rate_s = float(sum(c.prompt_loss_birth_rate_s for c in component_results))
    modal_prompt_loss_birth_power_W = float(sum(c.prompt_loss_birth_power_W for c in component_results))
    modal_source_birth_current_A = float(elementary_charge * fast_ion_species.charge_number * modal_source_birth_rate_s)
    modal_total_birth_current_A = float(elementary_charge * fast_ion_species.charge_number * modal_total_birth_rate_s)
    if modal_source_birth_rate_s > 0.0:
        modal_source_mean_birth_energy_J = modal_source_deposited_birth_power_W / modal_source_birth_rate_s
    else:
        modal_source_mean_birth_energy_J = 0.0
    velocity_model = velocity_solution_model(numerics.velocity_solution_model)
    rosen_coeffs: ModalRosenbluthCoefficients | None = None
    rosen_history: list[dict[str, float]] = []
    eq59_diagnostics: Eq59ConvergenceDiagnostics | None = None
    eq59_collision_operator_state: Eq59CollisionOperatorState | None = None
    eq59_speed_assessment: Eq59SpeedConvergenceDiagnostics | None = None
    eq59_production_grid_index: int | None = None
    published_heuristic: PublishedLambda1Heuristic | None = None
    active_eigenvalues_by_speed: np.ndarray | None = None
    if feedback_model in published_feedback_names and not fixed_boundary_production:
        if _published_throat_drop_override_J is None or _published_reference_lambda1 is None:
            raise RuntimeError("feedback internal state is incomplete")
        published_heuristic = egedal_2022_published_low_energy_lambda1_heuristic(speed_centers_m_s=speed_grid.centers_m_s, particle_mass_kg=fast_ion_species.mass_kg, mirror_ratio=mirror_ratio, throat_potential_drop_magnitude_J=_published_throat_drop_override_J, magnetic_eigenvalues=basis.physical_basis.eigenvalues, reference_lambda1_at_mirror_ratio_1p5=_published_reference_lambda1, geometry_extrapolation=True)
        active_eigenvalues_by_speed = published_heuristic.eigenvalues_by_speed
    if velocity_model == HOT_ION_ROSENBLUTH_EQ59:
        modal_total, rosen_coeffs, eq59_collision_operator_state, rosen_history, eq59_diagnostics = hot_ion_rosenbluth_modal_distribution_eq59(
            speed_grid=speed_grid,
            source_coefficients_by_component_j=np.asarray(component_coeffs, dtype=float),
            source_speeds_m_s=np.asarray(component_speeds, dtype=float),
            source_mean_birth_energy_J=modal_source_mean_birth_energy_J,
            eigenvalues=basis.physical_basis.eigenvalues,
            eigenvalues_by_speed=active_eigenvalues_by_speed,
            collision_state=collision_state,
            max_iterations=numerics.rosenbluth_max_iterations,
            relative_tolerance=numerics.rosenbluth_relative_tolerance,
            absolute_tolerance=numerics.rosenbluth_absolute_tolerance,
            relaxation=numerics.rosenbluth_relaxation,
            min_iterations=numerics.rosenbluth_min_iterations,
            mode_eta_integrals=_modal_mode_eta_integrals(basis),
            physical_eigenfunctions_eta=basis.physical_basis.eigenfunctions,
            eta_grid=basis.physical_basis.eta_grid,
            particle_mass_kg=fast_ion_species.mass_kg,
            warm_start_state=_eq59_warm_start_state,
            warm_start_compatibility_policy=_eq59_warm_start_compatibility_policy,
            species_label=fast_ion_species.symbol,
            charge_exchange_loss_frequency_s=charge_exchange_loss_frequency_s,
            external_fast_field_states=external_fast_field_states,
            return_diagnostics=True,
        )
        velocity_solution_equation = "egedal_eq59_hot_ion_rosenbluth_energy_diffusion"
        rosen_status = "enabled"
    elif velocity_model == COLD_ION_EQ14:
        if charge_exchange_sink_active:
            raise ValueError("charge exchange target redistribution requires the hot ion Eq 59 solver")
        if published_heuristic is not None:
            raise ValueError("low energy lambda1 feedback requires the hot ion Eq 59 solver")
        modal_total = np.sum(np.asarray(cold_components, dtype=float), axis=0) if cold_components else np.zeros((n_modes, speed_grid.centers_m_s.size), dtype=float)
        velocity_solution_equation = "egedal_eq14_cold_ion_slowing_down"
        rosen_status = "disabled_cold_ion_reference_mode"
    else:
        raise ValueError(f"Unsupported modal_velocity_solution_model {numerics.velocity_solution_model!r}")
    # The signed finite modal sum is checked before any roundoff scale negativity correction is accepted
    raw_f_eta = _modal_distribution_to_eta(modal_total, basis.physical_basis)
    f_eta, reconstruction_correction = _correct_reconstructed_distribution(distribution=raw_f_eta, speed_grid=speed_grid, eta_grid=basis.physical_basis.eta_grid, particle_mass_kg=fast_ion_species.mass_kg, particle_tolerance=numerics.reconstruction_negative_particle_fraction_tolerance, energy_tolerance=numerics.reconstruction_energy_correction_relative_tolerance)
    if reconstruction_correction["material"]:
        eq59_status = "not_applicable" if eq59_diagnostics is None else eq59_diagnostics.status
        eq59_residual = None if eq59_diagnostics is None or not eq59_diagnostics.residual_history else float(eq59_diagnostics.residual_history[-1])
        speed_index, eta_index = np.unravel_index(int(np.argmin(raw_f_eta)), raw_f_eta.shape)
        raise ModalPhysicalDistributionError(
            "modal reconstruction correction exceeds materiality tolerances; "
            f"species={fast_ion_species.species_id}; minimum={float(raw_f_eta[speed_index, eta_index]):.17g}; relative_minimum={reconstruction_correction['relative_minimum']:.17g}; negative_cell_fraction={reconstruction_correction['negative_cell_fraction']:.17g}; "
            f"negative_particle_fraction={reconstruction_correction['negative_particle_fraction']:.17g}; particle_tolerance={numerics.reconstruction_negative_particle_fraction_tolerance:.17g}; "
            f"energy_correction_relative={reconstruction_correction['energy_correction_relative']:.17g}; energy_tolerance={numerics.reconstruction_energy_correction_relative_tolerance:.17g}; "
            f"speed_index={speed_index}; speed_m_s={float(speed_grid.centers_m_s[speed_index]):.17g}; eta_index={eta_index}; eta={float(basis.physical_basis.eta_grid[eta_index]):.17g}; "
            f"n_physical_modes={n_modes}; electron_temperature_J={float(electron_temperature_J):.17g}; Eq59_status={eq59_status}; Eq59_final_relative_residual={eq59_residual}"
        )
    eta_integral_by_speed = np.trapezoid(np.maximum(f_eta, 0.0), basis.physical_basis.eta_grid, axis=1)
    invariant_energy_cell_population_weights = (np.asarray(speed_grid.shell_volumes_m3_s3, dtype=float) * eta_integral_by_speed)
    if published_heuristic is not None:
        published_heuristic = egedal_2022_published_low_energy_lambda1_heuristic(speed_centers_m_s=speed_grid.centers_m_s, particle_mass_kg=fast_ion_species.mass_kg, mirror_ratio=mirror_ratio, throat_potential_drop_magnitude_J=float(_published_throat_drop_override_J), magnetic_eigenvalues=basis.physical_basis.eigenvalues, reference_lambda1_at_mirror_ratio_1p5=float(_published_reference_lambda1), invariant_energy_population_weights=invariant_energy_cell_population_weights, geometry_extrapolation=True)
        active_eigenvalues_by_speed = published_heuristic.eigenvalues_by_speed
    eq59_speed_domain_applicable = eq59_diagnostics is not None
    eq59_high_speed_tail_population_fraction = (_upper_speed_tail_population_fraction(speed_grid, invariant_energy_cell_population_weights) if eq59_speed_domain_applicable else 0.0)
    eq59_high_speed_tail_energy_fraction = (_upper_speed_tail_energy_fraction(speed_grid, eta_integral_by_speed, particle_mass_kg=fast_ion_species.mass_kg) if eq59_speed_domain_applicable else 0.0)
    eq59_single_grid_tail_indicator_passed = bool(not eq59_speed_domain_applicable or (eq59_high_speed_tail_population_fraction <= float(numerics.rosenbluth_high_speed_tail_population_tolerance) and eq59_high_speed_tail_energy_fraction <= float(numerics.rosenbluth_high_speed_tail_energy_tolerance)))
    retained_mode_assessment = (None if n_modes < 2 or not _run_numerical_convergence_assessments else assess_retained_mode_reconstruction_convergence(basis=basis, modal_distribution_f_j_v=modal_total, speed_cell_weights=speed_grid.shell_volumes_m3_s3, relative_tolerance=float(numerics.retained_mode_relative_tolerance)))
    if (_run_numerical_convergence_assessments and eq59_diagnostics is not None and eq59_diagnostics.converged):
        speed_sequence, eq59_production_grid_index = (build_eq59_speed_convergence_grid_sequence(production_speed_grid=speed_grid, source_speeds_m_s=np.asarray(component_speeds, dtype=float), core_cell_fraction=numerics.speed_core_cell_fraction, tail_stretch_power=numerics.speed_tail_stretch_power))
        if active_eigenvalues_by_speed is None:
            eigenvalue_sequence = None
        else:
            if (_published_throat_drop_override_J is None or _published_reference_lambda1 is None):
                raise RuntimeError("speed dependent Eq 59 eigenvalues require the published throat drop and reference lambda1")
            eigenvalue_sequence = tuple(egedal_2022_published_low_energy_lambda1_heuristic(speed_centers_m_s=grid.centers_m_s, particle_mass_kg=fast_ion_species.mass_kg, mirror_ratio=mirror_ratio, throat_potential_drop_magnitude_J=float(_published_throat_drop_override_J), magnetic_eigenvalues=basis.physical_basis.eigenvalues, reference_lambda1_at_mirror_ratio_1p5=float(_published_reference_lambda1), geometry_extrapolation=True).eigenvalues_by_speed for grid in speed_sequence)
        eq59_speed_assessment = assess_eq59_speed_convergence(
            speed_grid_sequence=speed_sequence,
            source_coefficients_by_component_j=np.asarray(component_coeffs, dtype=float),
            source_speeds_m_s=np.asarray(component_speeds, dtype=float),
            eigenvalues=basis.physical_basis.eigenvalues,
            eigenvalues_by_speed_sequence=eigenvalue_sequence,
            source_mean_birth_energy_J=modal_source_mean_birth_energy_J,
            collision_state=collision_state,
            mode_eta_integrals=_modal_mode_eta_integrals(basis),
            physical_eigenfunctions_eta=(basis.physical_basis.eigenfunctions),
            eta_grid=basis.physical_basis.eta_grid,
            particle_mass_kg=fast_ion_species.mass_kg,
            max_iterations=numerics.rosenbluth_max_iterations,
            min_iterations=numerics.rosenbluth_min_iterations,
            relative_tolerance=numerics.rosenbluth_relative_tolerance,
            absolute_tolerance=numerics.rosenbluth_absolute_tolerance,
            relaxation=numerics.rosenbluth_relaxation,
            speed_convergence_relative_tolerance=(numerics.rosenbluth_speed_convergence_relative_tolerance),
            high_speed_tail_population_tolerance=(numerics.rosenbluth_high_speed_tail_population_tolerance),
            high_speed_tail_energy_tolerance=(numerics.rosenbluth_high_speed_tail_energy_tolerance),
            production_grid_index=eq59_production_grid_index,
            precomputed_modal_distribution=modal_total,
            precomputed_rosenbluth_coefficients=rosen_coeffs,
            precomputed_collision_operator_state=eq59_collision_operator_state,
            precomputed_nonlinear_diagnostics=eq59_diagnostics,
            external_fast_field_states=external_fast_field_states,
            charge_exchange_loss_frequency_source_grid=(speed_grid if charge_exchange_sink_active else None),
            charge_exchange_loss_frequency_s=(charge_exchange_loss_frequency_s if charge_exchange_sink_active else None),
        )

    return {
        'active_eigenvalues_by_speed': active_eigenvalues_by_speed,
        'charge_exchange_loss_frequency_s': charge_exchange_loss_frequency_s,
        'charge_exchange_sink_active': charge_exchange_sink_active,
        'charge_exchange_sink_rate_s': charge_exchange_sink_rate_s,
        'charge_exchange_sink_energy_W': charge_exchange_sink_energy_W,
        'charge_exchange_sink_model': charge_exchange_sink_model,
        'charge_exchange_sink_reference_identity_error': charge_exchange_sink_reference_identity_error,
        'external_fast_field_collision_states': external_fast_field_states,
        'external_fast_cross_collision_active': bool(external_fast_field_states),
        'external_fast_cross_collision_pair_ids': tuple(state.pair_id for state in external_fast_field_states),
        'component_coeffs': component_coeffs,
        'component_results': component_results,
        'component_speeds': component_speeds,
        'eq59_diagnostics': eq59_diagnostics,
        'eq59_collision_operator_state': eq59_collision_operator_state,
        'eq59_high_speed_tail_energy_fraction': eq59_high_speed_tail_energy_fraction,
        'eq59_high_speed_tail_population_fraction': eq59_high_speed_tail_population_fraction,
        'eq59_production_grid_index': eq59_production_grid_index,
        'eq59_single_grid_tail_indicator_passed': eq59_single_grid_tail_indicator_passed,
        'eq59_speed_assessment': eq59_speed_assessment,
        'f_eta': f_eta,
        'invariant_energy_cell_population_weights': invariant_energy_cell_population_weights,
        'modal_prompt_loss_birth_power_W': modal_prompt_loss_birth_power_W,
        'modal_prompt_loss_birth_rate_s': modal_prompt_loss_birth_rate_s,
        'modal_source_birth_current_A': modal_source_birth_current_A,
        'modal_source_birth_rate_s': modal_source_birth_rate_s,
        'modal_source_deposited_birth_power_W': modal_source_deposited_birth_power_W,
        'modal_source_mean_birth_energy_J': modal_source_mean_birth_energy_J,
        'modal_total': modal_total,
        'modal_total_birth_current_A': modal_total_birth_current_A,
        'modal_total_birth_rate_s': modal_total_birth_rate_s,
        'modal_total_deposited_birth_power_W': modal_total_deposited_birth_power_W,
        'n_modes': n_modes,
        'published_heuristic': published_heuristic,
        'reconstruction_correction': reconstruction_correction,
        'retained_mode_assessment': retained_mode_assessment,
        'rosen_coeffs': rosen_coeffs,
        'rosen_history': rosen_history,
        'rosen_status': rosen_status,
        'velocity_solution_equation': velocity_solution_equation,
    }
