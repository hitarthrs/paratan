"""Metadata assembly for modal FBIS solver results and numerical diagnostics"""
from __future__ import annotations
import numpy as np
from source_model_revamp.fusion.reactions import KEV_TO_J
from source_model_revamp.fbis.modal.lost_ion_distribution import BALDWIN_1972_THROAT_DENSITY_CLOSURE
from source_model_revamp.fbis.modal.utils import _EPS
from collections.abc import Mapping
from typing import Any
from source_model_revamp.fbis.modal.solver.solver_convergence import _safe_ratio
from source_model_revamp.fbis.modal.eq59.collision_operator_state import eq59_collision_operator_metadata
from source_model_revamp.fbis.modal.models import COLD_ION_EQ14, HOT_ION_ROSENBLUTH_EQ59

def build_modal_solver_metadata(context: Mapping[str, Any]) -> dict[str, Any]:
    """Collect solver choices, Eq 59 and Eq 70 convergence, loss closure, and reconstruction diagnostics"""
    _run_numerical_convergence_assessments = context['_run_numerical_convergence_assessments']
    active_loss = context['active_loss']
    barrier = context['barrier']
    basis = context['basis']
    basis_convergence_assessment = context['basis_convergence_assessment']
    boundary_flux_single_particle_loss_rate_s = context['boundary_flux_single_particle_loss_rate_s']
    boundary_flux_two_boundary_particle_loss_rate_s = context['boundary_flux_two_boundary_particle_loss_rate_s']
    component_results = context['component_results']
    confined_collisional_confinement = context['confined_collisional_confinement']
    confinement = context['confinement']
    current_balance = context['current_balance']
    current_balance_relative_error = context['current_balance_relative_error']
    density = context['density']
    directed_eq63_throat_rate_v_lambda_s = context['directed_eq63_throat_rate_v_lambda_s']
    effective_ion_temperature_J = context['effective_ion_temperature_J']
    eigen_audit = context['eigen_audit']
    electron_current_relative_error = context['electron_current_relative_error']
    electron_current_residual_A = context['electron_current_residual_A']
    electron_collision_density_m3 = context['electron_collision_density_m3']
    electron_midplane_density_m3 = context['electron_midplane_density_m3']
    eq70_prescribed_electron_midplane_density_m3 = context['eq70_prescribed_electron_midplane_density_m3']
    electron_power_W = context['electron_power_W']
    electron_temperature_J = context['electron_temperature_J']
    electrostatic_profile = context['electrostatic_profile']
    eq61_boundary_reconstruction_diagnostics = context['eq61_boundary_reconstruction_diagnostics']
    eq59_diagnostics = context['eq59_diagnostics']
    eq59_collision_operator_state = context['eq59_collision_operator_state']
    eq59_high_speed_tail_energy_fraction = context['eq59_high_speed_tail_energy_fraction']
    eq59_high_speed_tail_population_fraction = context['eq59_high_speed_tail_population_fraction']
    eq59_single_grid_tail_indicator_passed = context['eq59_single_grid_tail_indicator_passed']
    eq59_speed_assessment = context['eq59_speed_assessment']
    f_lambda = context['f_lambda']
    feedback_model = context['feedback_model']
    fixed_boundary_production = context['fixed_boundary_production']
    fast_ion_species = context['fast_ion_species']
    fixed_boundary_reconstruction_failure_reason = context['fixed_boundary_reconstruction_failure_reason']
    fixed_boundary_reconstruction_policy = context['fixed_boundary_reconstruction_policy']
    fixed_boundary_reconstruction_status = context['fixed_boundary_reconstruction_status']
    inventory = context['inventory']
    ion_barrier_power_W = context['ion_barrier_power_W']
    ion_midplane_power_loss_W = context['ion_midplane_power_loss_W']
    ion_particle_loss_rate_s = context['ion_particle_loss_rate_s']
    ion_wall_power_W = context['ion_wall_power_W']
    lambda_grid = context['lambda_grid']
    left_electron_loss_rate_s = context['left_electron_loss_rate_s']
    left_eq63_branch = context['left_eq63_branch']
    left_ion_loss_rate_s = context['left_ion_loss_rate_s']
    left_throat_temperature = context['left_throat_temperature']
    local_reconstruction = context['local_reconstruction']
    loss_convention_validation_passed = context['loss_convention_validation_passed']
    loss_model = context['loss_model']
    loss_validation_failures = context['loss_validation_failures']
    loss_validation_tolerance = context['loss_validation_tolerance']
    modal_per_lost_ion_barrier_energy_over_birth_energy = context['modal_per_lost_ion_barrier_energy_over_birth_energy']
    modal_per_lost_ion_midplane_energy_over_birth_energy = context['modal_per_lost_ion_midplane_energy_over_birth_energy']
    modal_per_lost_ion_total_wall_energy_over_birth_energy = context['modal_per_lost_ion_total_wall_energy_over_birth_energy']
    modal_prompt_loss_birth_power_W = context['modal_prompt_loss_birth_power_W']
    modal_prompt_loss_birth_rate_s = context['modal_prompt_loss_birth_rate_s']
    modal_source_birth_current_A = context['modal_source_birth_current_A']
    modal_source_birth_rate_s = context['modal_source_birth_rate_s']
    modal_source_deposited_birth_power_W = context['modal_source_deposited_birth_power_W']
    modal_source_mean_birth_energy_J = context['modal_source_mean_birth_energy_J']
    modal_total_birth_current_A = context['modal_total_birth_current_A']
    modal_total_birth_rate_s = context['modal_total_birth_rate_s']
    modal_total_deposited_birth_power_W = context['modal_total_deposited_birth_power_W']
    n_modes = context['n_modes']
    numerics = context['numerics']
    particle_balance_relative_error = context['particle_balance_relative_error']
    pitch_grid = context['pitch_grid']
    power_balance_relative_error = context['power_balance_relative_error']
    projected_source_component_rates_s = context['projected_source_component_rates_s']
    projected_source_particle_rate_s = context['projected_source_particle_rate_s']
    projected_source_to_birth_ratio = context['projected_source_to_birth_ratio']
    published_feedback_names = context['published_feedback_names']
    published_heuristic = context['published_heuristic']
    reconstruction_correction = context['reconstruction_correction']
    retained_mode_assessment = context['retained_mode_assessment']
    right_electron_loss_rate_s = context['right_electron_loss_rate_s']
    right_eq63_branch = context['right_eq63_branch']
    right_ion_loss_rate_s = context['right_ion_loss_rate_s']
    right_throat_temperature = context['right_throat_temperature']
    rosen_history = context['rosen_history']
    source_mean_energy = context['source_mean_energy']
    speed_grid = context['speed_grid']
    total_input_power = context['total_input_power']
    two_end_eq61_to_eigenvalue_ratio = context['two_end_eq61_to_eigenvalue_ratio']
    total_wall_power_W = context['total_wall_power_W']
    velocity_solution_equation = context['velocity_solution_equation']
    # Metadata keeps applicability separate from convergence so deferred checks are not reported as successful checks
    physical_basis = basis.physical_basis
    basis_assessed = basis_convergence_assessment is not None
    basis_converged = bool(basis_convergence_assessment.converged if basis_convergence_assessment is not None else False)
    eq59_applicable = eq59_diagnostics is not None
    eq59_speed_resolution_assessed = bool(eq59_speed_assessment is not None and eq59_speed_assessment.speed_resolution_assessed)
    eq59_speed_resolution_converged = bool(eq59_speed_assessment is not None and eq59_speed_assessment.speed_resolution_converged)
    eq59_speed_domain_assessed = bool(eq59_speed_assessment is not None and eq59_speed_assessment.speed_domain_assessed)
    eq59_speed_domain_converged = bool(eq59_speed_assessment is not None and eq59_speed_assessment.speed_domain_converged)
    eq59_converged = bool(eq59_diagnostics.converged and eq59_single_grid_tail_indicator_passed and eq59_speed_resolution_assessed and eq59_speed_resolution_converged and eq59_speed_domain_assessed and eq59_speed_domain_converged) if eq59_diagnostics is not None else True
    if eq59_diagnostics is None:
        eq59_overall_status = "not_applicable_cold_ion_eq14"
        eq59_failure_reason = None
    elif not eq59_diagnostics.converged:
        eq59_overall_status = eq59_diagnostics.status
        eq59_failure_reason = eq59_diagnostics.failure_reason
    elif not eq59_single_grid_tail_indicator_passed:
        eq59_overall_status = "selected_speed_grid_tail_not_converged"
        eq59_failure_reason = ("selected production speed grid exceeds the particle or energy tail tolerance")
    elif not _run_numerical_convergence_assessments:
        eq59_overall_status = "nested_speed_convergence_deferred"
        eq59_failure_reason = "intermediate_closure_state"
    elif eq59_speed_assessment is None:
        eq59_overall_status = "nested_speed_convergence_not_assessed"
        eq59_failure_reason = "speed_assessment_unavailable"
    elif not eq59_speed_assessment.converged:
        eq59_overall_status = "nested_speed_convergence_failed"
        eq59_failure_reason = eq59_speed_assessment.failure_reason
    else:
        eq59_overall_status = "converged"
        eq59_failure_reason = None
    eq59_operator_metadata = {} if eq59_collision_operator_state is None else eq59_collision_operator_metadata(eq59_collision_operator_state)
    electron_balance = current_balance.electron_balance
    left_T_L_J = (None if left_throat_temperature is None else float(left_throat_temperature.eq63_parallel_temperature_J))
    right_T_L_J = (None if right_throat_temperature is None else float(right_throat_temperature.eq63_parallel_temperature_J))
    left_eq63_rate_s = (None if left_eq63_branch is None else float(left_eq63_branch.reconstructed_total_rate_s))
    right_eq63_rate_s = (None if right_eq63_branch is None else float(right_eq63_branch.reconstructed_total_rate_s))
    metadata = {
        "kinetic_backend": "egedal_modal",
        "kinetic_model": "egedal_modal_fbis",
        "modal_fbis_model": "egedal_orbit_averaged_eigenbasis",
        "ion_loss_closure_model": loss_model,
        "modal_ion_loss_closure_model": loss_model,
        "fixed_boundary_lost_reconstruction_policy": (fixed_boundary_reconstruction_policy if fixed_boundary_production else None),
        "fixed_boundary_lost_reconstruction_status": (fixed_boundary_reconstruction_status),
        "fixed_boundary_lost_reconstruction_failure_reason": (fixed_boundary_reconstruction_failure_reason),
        "fixed_boundary_lost_reconstruction_complete": bool(not fixed_boundary_production or fixed_boundary_reconstruction_status == "represented"),
        "active_ion_loss_boundary_model": ("fixed_magnetic_boundary_Lambda_M_equals_1_over_RM" if fixed_boundary_production else str(active_loss["model"])),
        "active_ion_loss_uses_Eq71_moving_boundary": False if fixed_boundary_production else bool(published_heuristic is not None),
        "Eq71_72_loss_boundary_role": ("diagnostic_local_accessibility_and_quasineutral_reconstruction_only" if fixed_boundary_production else "published_approximate_electrostatic_loss_feedback"),
        "lost_ion_parallel_temperature_model": (BALDWIN_1972_THROAT_DENSITY_CLOSURE if fixed_boundary_production else None),
        "lost_ion_parallel_temperature_left_J": left_T_L_J,
        "lost_ion_parallel_temperature_right_J": right_T_L_J,
        "lost_ion_parallel_temperature_left_keV": None if left_T_L_J is None else float(left_T_L_J / KEV_TO_J),
        "lost_ion_parallel_temperature_right_keV": None if right_T_L_J is None else float(right_T_L_J / KEV_TO_J),
        "lost_ion_parallel_temperature_left_diagnostics": None if left_throat_temperature is None else dict(left_throat_temperature.diagnostics or {}),
        "lost_ion_parallel_temperature_right_diagnostics": None if right_throat_temperature is None else dict(right_throat_temperature.diagnostics or {}),
        "modal_eq61_left_one_end_particle_loss_rate_s": float(boundary_flux_single_particle_loss_rate_s),
        "modal_eq61_right_one_end_particle_loss_rate_s": float(boundary_flux_single_particle_loss_rate_s),
        "modal_eq61_boundary_reconstruction_diagnostics": dict(eq61_boundary_reconstruction_diagnostics),
        "modal_eq63_left_throat_crossing_rate_s": left_eq63_rate_s,
        "modal_eq63_right_throat_crossing_rate_s": right_eq63_rate_s,
        "modal_eq63_left_rate_relative_error": None if left_eq63_branch is None else float(left_eq63_branch.relative_rate_error),
        "modal_eq63_right_rate_relative_error": None if right_eq63_branch is None else float(right_eq63_branch.relative_rate_error),
        "modal_eq63_left_H_U": None if left_eq63_branch is None else left_eq63_branch.H_U,
        "modal_eq63_right_H_U": None if right_eq63_branch is None else right_eq63_branch.H_U,
        "modal_eq63_left_throat_distribution_v_lambda": None if left_eq63_branch is None else left_eq63_branch.throat_distribution_v_lambda,
        "modal_eq63_right_throat_distribution_v_lambda": None if right_eq63_branch is None else right_eq63_branch.throat_distribution_v_lambda,
        "modal_eq63_left_throat_rate_v_lambda_s": None if left_eq63_branch is None else left_eq63_branch.throat_rate_v_lambda_s,
        "modal_eq63_right_throat_rate_v_lambda_s": None if right_eq63_branch is None else right_eq63_branch.throat_rate_v_lambda_s,
        "modal_eq64_matching_model": None if left_eq63_branch is None else left_eq63_branch.matching_model,
        "modal_eq63_pitch_integration_model": "analytic_exact_magnetic_loss_cone_boundary" if fixed_boundary_production else None,
        "electrostatic_feedback_model": (feedback_model if fixed_boundary_production else ("magnetic_only_fbis" if published_heuristic is None else published_heuristic.model)),
        "electrostatic_feedback_reference": ("Egedal_et_al_Nuclear_Fusion_62_126053_2022"),
        "electrostatic_feedback_is_exact": False,
        "published_lambda1_heuristic_population_fraction": (0.0 if published_heuristic is None else float(published_heuristic.active_population_fraction)),
        "published_lambda1_heuristic_higher_modes_modified": False,
        "electron_midplane_density_m3": (float(electrostatic_profile.electron_midplane_density_m3) if electrostatic_profile is not None else float(electron_midplane_density_m3)),
        "electron_midplane_density_used_by_Eq68_m3": float(electron_midplane_density_m3),
        "electron_midplane_density_prescribed_to_Eq70_m3": (None if eq70_prescribed_electron_midplane_density_m3 is None else float(eq70_prescribed_electron_midplane_density_m3)),
        "electron_parent_maxwellian_n0_m3": float(current_balance.electron_balance.electron_parent_maxwellian_n0_m3),
        "electron_collision_density_m3": float(electron_collision_density_m3),
        "electron_confined_volume_average_density_m3": (None if electrostatic_profile is None else float(electrostatic_profile.electron_volume_average_density_m3)),
        "electron_density_profile_m3": (None if electrostatic_profile is None else [float(x) for x in electrostatic_profile.electron_density_m3]),
        "electron_profile_scope": (None if electrostatic_profile is None else "confined_throat_to_throat"),
        "expander_electron_profile_available": False,
        "modal_velocity_solution_model": str(numerics.velocity_solution_model),
        "modal_velocity_solution_equation": velocity_solution_equation,
        "modal_component_distribution_role": ("cold_Eq14_reference_diagnostic" if numerics.velocity_solution_model == HOT_ION_ROSENBLUTH_EQ59 else "active_cold_Eq14_component_solution"),
        "modal_component_distributions_sum_to_active_solution": bool(numerics.velocity_solution_model == COLD_ION_EQ14),
        "modal_legacy_critical_velocity_authority": ("diagnostic_only" if numerics.velocity_solution_model == HOT_ION_ROSENBLUTH_EQ59 else "active_cold_Eq14_reference_coefficient"),
        "modal_legacy_beta_m_authority": ("diagnostic_only" if numerics.velocity_solution_model == HOT_ION_ROSENBLUTH_EQ59 else "active_cold_Eq14_reference_coefficient"),
        **eq59_operator_metadata,
        "modal_rosenbluth_iterations": len(rosen_history),
        "modal_rosenbluth_iteration_history": rosen_history,
        "modal_rosenbluth_convergence_applicable": eq59_applicable,
        "modal_rosenbluth_convergence_assessed": eq59_applicable,
        "modal_rosenbluth_converged": eq59_converged,
        "modal_rosenbluth_convergence_status": eq59_overall_status,
        "modal_rosenbluth_convergence_failure_reason": eq59_failure_reason,
        "modal_rosenbluth_fixed_point_iteration_converged": False if eq59_diagnostics is None else bool(eq59_diagnostics.fixed_point_iteration_converged),
        "modal_rosenbluth_nonlinear_residual_converged": False if eq59_diagnostics is None else bool(eq59_diagnostics.nonlinear_residual_converged),
        "modal_rosenbluth_inventory_converged": False if eq59_diagnostics is None else bool(eq59_diagnostics.inventory_converged),
        "modal_rosenbluth_effective_energy_converged": False if eq59_diagnostics is None else bool(eq59_diagnostics.effective_energy_converged),
        "modal_rosenbluth_fast_self_density_converged": False if eq59_diagnostics is None else bool(eq59_diagnostics.fast_self_density_converged),
        "modal_rosenbluth_collision_operator_converged": False if eq59_diagnostics is None else bool(eq59_diagnostics.collision_operator_converged),
        "modal_rosenbluth_fast_self_density_history_m3": [] if eq59_diagnostics is None else [float(x) for x in eq59_diagnostics.fast_self_density_history_m3],
        "modal_rosenbluth_fast_self_density_relative_change_history": [] if eq59_diagnostics is None else [float(x) for x in eq59_diagnostics.fast_self_density_relative_change_history],
        "modal_rosenbluth_fast_self_g_scale_history_m3_s3": [] if eq59_diagnostics is None else [float(x) for x in eq59_diagnostics.fast_self_g_scale_history_m3_s3],
        "modal_rosenbluth_fast_self_net_drag_scale_history_m3_s3": [] if eq59_diagnostics is None else [float(x) for x in eq59_diagnostics.fast_self_net_drag_scale_history_m3_s3],
        "modal_rosenbluth_ion_drag_operator_relative_change_history": [] if eq59_diagnostics is None else [float(x) for x in eq59_diagnostics.ion_drag_operator_relative_change_history],
        "modal_rosenbluth_ion_energy_diffusion_operator_relative_change_history": [] if eq59_diagnostics is None else [float(x) for x in eq59_diagnostics.ion_energy_diffusion_operator_relative_change_history],
        "modal_rosenbluth_ion_pitch_scattering_operator_relative_change_history": [] if eq59_diagnostics is None else [float(x) for x in eq59_diagnostics.ion_pitch_scattering_operator_relative_change_history],
        "modal_rosenbluth_warm_start_used": False if eq59_diagnostics is None else bool(eq59_diagnostics.warm_start_used),
        "modal_rosenbluth_warm_start_remapped": False if eq59_diagnostics is None else bool(eq59_diagnostics.warm_start_remapped),
        "modal_rosenbluth_warm_start_status": None if eq59_diagnostics is None else eq59_diagnostics.warm_start_status,
        "modal_rosenbluth_warm_start_rejection_reason": None if eq59_diagnostics is None else eq59_diagnostics.warm_start_rejection_reason,
        "modal_rosenbluth_warm_start_fallback_to_pairwise_seed": False if eq59_diagnostics is None else bool(eq59_diagnostics.warm_start_fallback_to_pairwise_seed),
        "modal_reconstruction_correction_model": "clip_and_particle_inventory_renormalization",
        "modal_reconstruction_correction_applied": bool(reconstruction_correction["applied"]),
        "modal_reconstruction_pointwise_relative_minimum": float(reconstruction_correction["relative_minimum"]),
        "modal_reconstruction_negative_cell_fraction": float(reconstruction_correction["negative_cell_fraction"]),
        "modal_reconstruction_negative_particle_fraction": float(reconstruction_correction["negative_particle_fraction"]),
        "modal_reconstruction_negative_particle_fraction_tolerance": float(numerics.reconstruction_negative_particle_fraction_tolerance),
        "modal_reconstruction_energy_correction_relative": float(reconstruction_correction["energy_correction_relative"]),
        "modal_reconstruction_energy_correction_relative_tolerance": float(numerics.reconstruction_energy_correction_relative_tolerance),
        "modal_reconstruction_particle_renormalization_factor": float(reconstruction_correction["renormalization_factor"]),
        "modal_rosenbluth_relative_tolerance": float(numerics.rosenbluth_relative_tolerance),
        "modal_rosenbluth_absolute_tolerance": float(numerics.rosenbluth_absolute_tolerance),
        "modal_rosenbluth_relaxation": float(numerics.rosenbluth_relaxation),
        "modal_rosenbluth_min_iterations": int(numerics.rosenbluth_min_iterations),
        "modal_rosenbluth_high_speed_tail_population_fraction": float(eq59_high_speed_tail_population_fraction),
        "modal_rosenbluth_high_speed_tail_population_tolerance": float(numerics.rosenbluth_high_speed_tail_population_tolerance),
        "modal_rosenbluth_high_speed_tail_energy_fraction": float(eq59_high_speed_tail_energy_fraction),
        "modal_rosenbluth_high_speed_tail_energy_tolerance": float(numerics.rosenbluth_high_speed_tail_energy_tolerance),
        "modal_rosenbluth_single_grid_tail_indicator_passed": bool(eq59_single_grid_tail_indicator_passed),
        "modal_rosenbluth_speed_resolution_assessed": bool(eq59_speed_resolution_assessed),
        "modal_rosenbluth_speed_resolution_converged": bool(eq59_speed_resolution_converged),
        "modal_rosenbluth_speed_domain_assessed": bool(eq59_speed_domain_assessed),
        "modal_rosenbluth_speed_domain_converged": bool(eq59_speed_domain_converged),
        "modal_rosenbluth_speed_convergence_relative_tolerance": float(numerics.rosenbluth_speed_convergence_relative_tolerance),
        "modal_rosenbluth_speed_selected_grid_matched": False if eq59_speed_assessment is None else bool(eq59_speed_assessment.production_grid_matched),
        "modal_rosenbluth_selected_speed_cell_count": int(speed_grid.centers_m_s.size) if eq59_speed_assessment is None else int(eq59_speed_assessment.selected_speed_cell_count),
        "modal_rosenbluth_selected_speed_domain_upper_m_s": float(speed_grid.faces_m_s[-1]) if eq59_speed_assessment is None else float(eq59_speed_assessment.selected_speed_domain_upper_m_s),
        "modal_rosenbluth_speed_common_domain_distribution_relative_change_history": [] if eq59_speed_assessment is None else [float(x) for x in eq59_speed_assessment.common_domain_distribution_relative_change_history],
        "modal_rosenbluth_speed_inventory_relative_change_history": [] if eq59_speed_assessment is None else [float(x) for x in eq59_speed_assessment.inventory_relative_change_history],
        "modal_rosenbluth_speed_effective_energy_relative_change_history": [] if eq59_speed_assessment is None else [float(x) for x in eq59_speed_assessment.effective_energy_relative_change_history],
        "modal_rosenbluth_speed_particle_loss_relative_change_history": [] if eq59_speed_assessment is None else [float(x) for x in eq59_speed_assessment.modal_particle_loss_relative_change_history],
        "modal_rosenbluth_speed_power_loss_relative_change_history": [] if eq59_speed_assessment is None else [float(x) for x in eq59_speed_assessment.modal_power_loss_relative_change_history],
        "modal_rosenbluth_speed_electron_heating_power_history_W_per_m3": [] if eq59_speed_assessment is None else [float(x) for x in eq59_speed_assessment.electron_heating_power_history_W_per_m3],
        "modal_rosenbluth_speed_electron_heating_relative_change_history": [] if eq59_speed_assessment is None else [float(x) for x in eq59_speed_assessment.electron_heating_power_relative_change_history],
        "modal_rosenbluth_speed_electron_heating_resolution_assessed": False if eq59_speed_assessment is None else bool(eq59_speed_assessment.electron_heating_speed_resolution_assessed),
        "modal_rosenbluth_speed_electron_heating_resolution_converged": False if eq59_speed_assessment is None else bool(eq59_speed_assessment.electron_heating_speed_resolution_converged),
        "modal_rosenbluth_speed_electron_heating_domain_assessed": False if eq59_speed_assessment is None else bool(eq59_speed_assessment.electron_heating_speed_domain_assessed),
        "modal_rosenbluth_speed_electron_heating_domain_converged": False if eq59_speed_assessment is None else bool(eq59_speed_assessment.electron_heating_speed_domain_converged),
        "modal_rosenbluth_speed_source_to_loss_ratio_relative_change_history": [] if eq59_speed_assessment is None else [float(x) for x in eq59_speed_assessment.source_to_modal_loss_ratio_relative_change_history],
        "modal_rosenbluth_final_relative_residual": None if eq59_diagnostics is None or not eq59_diagnostics.residual_history else float(eq59_diagnostics.residual_history[-1]),
        "modal_electrostatic_phi_z_status": "egedal_eq70_72_quasineutral_iteration_enabled" if electrostatic_profile is not None else "not_evaluated_no_pitch_grid",
        "modal_electron_loss_equation": "egedal_eq68_eq69_tail_refilling",
        "modal_local_reconstruction_model": "egedal_eq71_72_conserved_total_energy_and_magnetic_moment",
        "modal_lost_ion_distribution_model": ("egedal_eq63_directed_left_and_right_fixed_magnetic_boundary" if fixed_boundary_production else "unavailable_active_hot_eq59_eq63_eq64_population"),
        "modal_eq63_eq64_published_scope": "hot_eq59_fixed_magnetic_boundary_with_U_equals_E_plus_ePhi_energy_shift",
        "modal_eq64_g_tilde_1_printed_form_status": ("active_pitch_array_is_stored_once_and_used_by_Eq59_Eq61_Eq63_and_Eq64_matching" if fixed_boundary_production else "published_g_tilde_1_prefactor_not_activated_for_the_unavailable_moving_boundary_population"),
        "modal_eq59_maxwellian_g_tilde_subscript_interpretation": "v_Psi_over_x_is_g_tilde_2_following_the_defining_derivatives_and_dimensions",
        "modal_eq63_eq64_active_implementation_status": ("implemented_for_selected_hot_fixed_magnetic_boundary_model" if fixed_boundary_production else "not_implemented_or_validated_for_active_hot_electrostatic_population"),
        "modal_n_lambda_grid": int(basis.eta_lambda_map.lambda_grid.size),
        "modal_n_eta_grid": int(physical_basis.eta_grid.size),
        "modal_n_square_basis_modes": int(numerics.n_square_basis_modes),
        "modal_n_physical_modes": int(numerics.n_physical_modes),
        "modal_basis_convergence_assessed": basis_assessed,
        "modal_basis_converged": basis_converged,
        "modal_basis_convergence_failure_reason": None if basis_convergence_assessment is None or basis_convergence_assessment.converged else basis_convergence_assessment.failure_reason,
        "modal_basis_convergence_active_basis_role": "configured_production_grid",
        "modal_basis_convergence_reference_basis_role": None if basis_convergence_assessment is None else "independent_finer_grid_assessment_only",
        "modal_basis_convergence_grid_sequence": [] if basis_convergence_assessment is None else [[int(a), int(b)] for a, b in basis_convergence_assessment.grid_sequence],
        "modal_basis_convergence_active_grid": [int(basis.eta_lambda_map.lambda_grid.size), int(physical_basis.eta_grid.size)],
        "modal_basis_convergence_finest_reference_grid": None if basis_convergence_assessment is None else [int(basis_convergence_assessment.grid_sequence[-1][0]), int(basis_convergence_assessment.grid_sequence[-1][1])],
        "modal_basis_convergence_eigenvalue_relative_change_history": [] if basis_convergence_assessment is None else [[float(x) for x in values] for values in basis_convergence_assessment.eigenvalue_relative_change_history],
        "modal_basis_convergence_maximum_eigenvalue_relative_change_history": [] if basis_convergence_assessment is None else [float(np.max(values)) for values in basis_convergence_assessment.eigenvalue_relative_change_history],
        "modal_basis_convergence_eigenfunction_overlap_history": [] if basis_convergence_assessment is None else [[float(x) for x in values] for values in basis_convergence_assessment.eigenfunction_overlap_history],
        "modal_basis_convergence_minimum_eigenfunction_overlap_history": [] if basis_convergence_assessment is None else [float(np.min(values)) for values in basis_convergence_assessment.eigenfunction_overlap_history],
        "modal_basis_convergence_max_eigenpair_residual": None if basis_convergence_assessment is None else float(basis_convergence_assessment.max_eigenpair_residual),
        "modal_basis_convergence_gram_condition_number": None if basis_convergence_assessment is None else float(basis_convergence_assessment.gram_condition_number),
        "modal_retained_mode_count_assessed": retained_mode_assessment is not None,
        "modal_retained_mode_count_converged": None if retained_mode_assessment is None else retained_mode_assessment.converged,
        "modal_retained_mode_relative_tolerance": float(numerics.retained_mode_relative_tolerance),
        "modal_retained_mode_selected_count": int(n_modes),
        "modal_retained_mode_final_distribution_relative_change": None if retained_mode_assessment is None else float(retained_mode_assessment.weighted_relative_change_history[-1]),
        "modal_retained_mode_final_inventory_relative_change": None if retained_mode_assessment is None else float(retained_mode_assessment.inventory_relative_change_history[-1]),
        "modal_basis_gram_condition_number": float(physical_basis.gram_condition_number) if basis_convergence_assessment is None else float(basis_convergence_assessment.gram_condition_number),
        "modal_eq42_density_weighting_model": basis.eq42_density_weighting_model,
        "modal_eq42_nonuniform_density_weighting_active": basis.eq42_nonuniform_density_weighting_active,
        "modal_eq42_nonuniform_density_weighting_required": feedback_model in published_feedback_names,
        "modal_eigenvalues": [float(x) for x in basis.physical_basis.eigenvalues],
        "modal_eta_trapped_passing_boundary": float(basis.eta_boundary),
        "modal_lambda_boundary": float(basis.lambda_boundary),
        "left_ion_loss_rate_s": float(left_ion_loss_rate_s),
        "right_ion_loss_rate_s": float(right_ion_loss_rate_s),
        "left_electron_loss_rate_s": float(left_electron_loss_rate_s),
        "right_electron_loss_rate_s": float(right_electron_loss_rate_s),
        "modal_single_end_eq61_loss_rate_s": float(boundary_flux_single_particle_loss_rate_s),
        "modal_two_end_eq61_loss_rate_s": float(boundary_flux_two_boundary_particle_loss_rate_s),
        "modal_two_end_eq61_to_eigenvalue_ratio": None if two_end_eq61_to_eigenvalue_ratio is None else float(two_end_eq61_to_eigenvalue_ratio),
        "modal_eq59_pitch_scattering_model": str(eigen_audit["pitch_scattering_model"]),
        "modal_projected_source_to_birth_ratio": projected_source_to_birth_ratio,
        "modal_current_balance_relative_error": current_balance_relative_error,
        "modal_power_balance_relative_error": power_balance_relative_error,
        "modal_loss_convention_relative_tolerance": loss_validation_tolerance,
        "modal_source_projection_relative_tolerance": float(numerics.source_projection_relative_tolerance),
        "modal_loss_convention_validation_passed": bool(loss_convention_validation_passed),
        "modal_loss_convention_failure_reason": None if loss_convention_validation_passed else ";".join(loss_validation_failures),
        "modal_component_confined_birth_power_W": [float(c.confined_birth_power_W) for c in component_results],
        "modal_component_prompt_loss_birth_power_W": [float(c.prompt_loss_birth_power_W) for c in component_results],
        "modal_total_beam_birth_rate_s": float(modal_total_birth_rate_s),
        "modal_total_beam_deposited_birth_power_W": float(modal_total_deposited_birth_power_W),
        "modal_total_input_power_W": float(total_input_power),
        "modal_source_birth_rate_s": float(modal_source_birth_rate_s),
        "modal_source_birth_current_A": float(modal_source_birth_current_A),
        "modal_total_birth_current_A": float(modal_total_birth_current_A),
        "modal_source_deposited_birth_power_W": float(modal_source_deposited_birth_power_W),
        "modal_source_mean_birth_energy_J": float(modal_source_mean_birth_energy_J),
        "modal_source_mean_birth_energy_keV": float(modal_source_mean_birth_energy_J / KEV_TO_J) if modal_source_mean_birth_energy_J > 0.0 else 0.0,
        "modal_prompt_loss_birth_rate_s": float(modal_prompt_loss_birth_rate_s),
        "modal_prompt_loss_birth_power_W": float(modal_prompt_loss_birth_power_W),
        "modal_prompt_loss_birth_fraction": _safe_ratio(modal_prompt_loss_birth_rate_s, modal_total_birth_rate_s),
        "modal_confined_birth_fraction": _safe_ratio(modal_source_birth_rate_s, modal_total_birth_rate_s),
        "modal_projected_source_particle_rate_from_coefficients_s": float(projected_source_particle_rate_s),
        "modal_projected_source_component_rates_s": [float(x) for x in projected_source_component_rates_s],
        "modal_component_source_coefficients": [[float(x) for x in c.modal_source_coefficients] for c in component_results],
        "modal_effective_ion_temperature_J": float(effective_ion_temperature_J),
        "modal_effective_ion_temperature_keV": float(effective_ion_temperature_J / KEV_TO_J),
        "modal_effective_ion_temperature_over_beam_energy": float(effective_ion_temperature_J / max(sum(c.confined_birth_power_W for c in component_results) / max(sum(c.confined_birth_rate_s for c in component_results), _EPS), _EPS)) if component_results else None,
        "modal_ion_particle_loss_rate_s": float(ion_particle_loss_rate_s),
        "modal_ion_midplane_kinetic_power_loss_W": float(ion_midplane_power_loss_W),
        "modal_ion_wall_power_loss_W": float(ion_wall_power_W),
        "modal_ion_barrier_power_loss_W": float(ion_barrier_power_W),
        "modal_particle_balance_relative_error_loss_minus_source_over_source": None if particle_balance_relative_error is None else float(particle_balance_relative_error),
        "modal_total_wall_power_loss_W": float(total_wall_power_W),
        "modal_per_lost_ion_midplane_energy_over_birth_energy": None if modal_per_lost_ion_midplane_energy_over_birth_energy is None else float(modal_per_lost_ion_midplane_energy_over_birth_energy),
        "modal_per_lost_ion_barrier_energy_over_birth_energy": None if modal_per_lost_ion_barrier_energy_over_birth_energy is None else float(modal_per_lost_ion_barrier_energy_over_birth_energy),
        "modal_per_lost_ion_total_wall_energy_over_birth_energy": None if modal_per_lost_ion_total_wall_energy_over_birth_energy is None else float(modal_per_lost_ion_total_wall_energy_over_birth_energy),
        "modal_electron_wall_power_loss_W": float(electron_power_W),
        "modal_electron_tail_integral_production_model": "exact_integrals_preceding_Egedal_Eq68_Eq69",
        "modal_electron_tail_integral_asymptotic_role": "diagnostic_only",
        "modal_electron_tail_integral_exact": True,
        "modal_eq68_exact_tail_kernel": None if electron_balance.exact_tail_kernel is None else float(electron_balance.exact_tail_kernel),
        "modal_eq68_asymptotic_tail_kernel": None if electron_balance.eq68_asymptotic_tail_kernel is None else float(electron_balance.eq68_asymptotic_tail_kernel),
        "modal_eq68_asymptotic_relative_correction": None if electron_balance.eq68_asymptotic_relative_correction is None else float(electron_balance.eq68_asymptotic_relative_correction),
        "modal_eq69_asymptotic_relative_correction": None if electron_balance.eq69_asymptotic_relative_correction is None else float(electron_balance.eq69_asymptotic_relative_correction),
        "modal_electron_mean_wall_kinetic_energy_J": None if electron_balance.mean_wall_kinetic_energy_J is None else float(electron_balance.mean_wall_kinetic_energy_J),
        "modal_electron_number_of_loss_ends": float(electron_balance.number_of_ends),
        "modal_electron_current_balance_residual_A": float(electron_current_residual_A),
        "modal_electron_current_balance_relative_error": float(electron_current_relative_error),
        "modal_wall_barrier_energy_J": float(barrier),
        "modal_wall_barrier_over_Te": float(barrier / electron_temperature_J) if electron_temperature_J > 0.0 else None,
        "modal_electron_wall_potential_relative_to_midplane_V": float(np.asarray(current_balance.electron_wall_potential_relative_to_midplane_V, dtype=float)),
        "modal_ion_confinement_time_s": float(confinement) if confinement is not None else None,
        "modal_confined_ion_collisional_confinement_time_s": float(confined_collisional_confinement) if confined_collisional_confinement is not None else None,
        "fast_ion_species_id": fast_ion_species.species_id,
        "fast_ion_species_name": fast_ion_species.name,
        "fast_ion_species_symbol": fast_ion_species.symbol,
        "fast_ion_species_mass_kg": float(fast_ion_species.mass_kg),
        "fast_ion_species_charge_number": float(fast_ion_species.charge_number),
        "electron_loss_model": "egedal_tail_refilling",
        "electrostatic_model": "egedal_modal_wall_current_balance_and_phi_z_quasineutrality",
    }
    if fast_ion_species.species_id == "deuterium":
        metadata.update({
            "solved_fast_D_inventory_particles": float(inventory),
            "solved_fast_D_ion_particle_loss_rate_s": float(ion_particle_loss_rate_s),
            "solved_fast_D_ion_confinement_time_s": float(confinement) if confinement is not None else None,
            "solved_fast_D_volume_average_density_m3_from_inventory_sum": float(density),
            "final_fast_D_distribution_v_lambda": f_lambda,
        })
    if eq70_prescribed_electron_midplane_density_m3 is not None:
        metadata.update({'background_ion_temperature_keV': None, 'background_profile_model': None, 'background_profile_midplane_normalized': None})
    metadata.update({
        "speed_centers_m_s": speed_grid.centers_m_s,
        "speed_faces_m_s": speed_grid.faces_m_s,
        "lambda_centers": lambda_grid.centers,
        "lambda_faces": lambda_grid.faces,
        "pitch_centers": None if pitch_grid is None else pitch_grid.centers,
        "pitch_faces": None if pitch_grid is None else pitch_grid.faces,
        "modal_eta_lambda_grid": basis.eta_lambda_map.lambda_grid,
        "modal_eta_of_lambda_grid": basis.eta_lambda_map.eta_grid,
        "modal_physical_eta_grid": basis.physical_basis.eta_grid,
        "modal_physical_eigenfunctions_eta": basis.physical_basis.eigenfunctions,
        "modal_source_mean_birth_speed_m_s": float(np.sqrt(2.0 * source_mean_energy / fast_ion_species.mass_kg)) if source_mean_energy > 0.0 else None,
    })
    if electrostatic_profile is not None:
        metadata.update({
            "modal_phi_z_iterations": int(electrostatic_profile.iterations),
            "modal_phi_z_converged": bool(electrostatic_profile.converged),
            "modal_phi_z_max_relative_quasineutrality_error": float(electrostatic_profile.max_relative_quasineutrality_error),
            "modal_phi_z_min_V": float(np.min(electrostatic_profile.potential_relative_to_midplane_V)),
            "modal_phi_z_max_V": float(np.max(electrostatic_profile.potential_relative_to_midplane_V)),
            "modal_phi_z_profile_V": [float(x) for x in electrostatic_profile.potential_relative_to_midplane_V],
            "modal_phi_z_electron_density_profile_m3": [float(x) for x in electrostatic_profile.electron_density_m3],
            "modal_phi_z_ion_density_profile_m3": [float(x) for x in electrostatic_profile.ion_density_m3],
            "electron_parent_maxwellian_n0_m3": float(electrostatic_profile.electron_parent_maxwellian_n0_m3),
            "electron_collision_density_m3": float(electrostatic_profile.electron_collision_density_m3),
            "electron_volume_average_density_m3": float(electrostatic_profile.electron_volume_average_density_m3),
            "electron_confined_volume_average_density_m3": float(electrostatic_profile.electron_volume_average_density_m3),
            "electron_density_profile_m3": [float(x) for x in electrostatic_profile.electron_density_m3],
            "electron_profile_scope": "confined_throat_to_throat",
            "expander_electron_profile_available": False,
            "electron_density_closure_model": electrostatic_profile.electron_density_closure_model,
            "electron_density_closure_converged": bool(electrostatic_profile.converged),
            "electron_density_closure_failure_reason": electrostatic_profile.failure_reason,
            "modal_phi_z_relative_residual_profile": [float(x) for x in electrostatic_profile.relative_residual_profile],
            "modal_phi_z_absolute_residual_profile_m3": [float(x) for x in electrostatic_profile.absolute_residual_profile_m3],
            "modal_phi_z_residual_history": [float(x) for x in electrostatic_profile.residual_history],
            "modal_phi_z_failure_reason": electrostatic_profile.failure_reason,
            "modal_phi_z_throat_potential_energy_J": float(electrostatic_profile.throat_potential_energy_J),
            "modal_phi_z_throat_potential_left_energy_J": float(electrostatic_profile.throat_potential_left_energy_J),
            "modal_phi_z_throat_potential_right_energy_J": float(electrostatic_profile.throat_potential_right_energy_J),
            "modal_phi_z_low_energy_approximation_max_weight_fraction": float(electrostatic_profile.low_energy_approximation_max_weight_fraction),
            "modal_phi_z_low_energy_approximation_valid": bool(electrostatic_profile.low_energy_approximation_valid),
            "modal_phi_z_eq71_closed_interval_intersects_distribution_support": bool(electrostatic_profile.eq71_closed_interval_intersects_distribution_support),
            "modal_phi_z_maximum_absolute_density_residual_normalized_to_reference": float(electrostatic_profile.maximum_absolute_density_residual_normalized_to_reference),
            "modal_phi_z_volume_integrated_absolute_particle_mismatch_fraction": float(electrostatic_profile.volume_integrated_absolute_particle_mismatch_fraction),
            "modal_phi_z_electrostatic_node_potential_V": None if electrostatic_profile.electrostatic_node_potential_relative_to_midplane_V is None else [float(x) for x in electrostatic_profile.electrostatic_node_potential_relative_to_midplane_V],
            "modal_phi_z_exact_midplane_quasineutrality_residual": float(electrostatic_profile.exact_midplane_quasineutrality_residual),
            "modal_phi_z_exact_left_throat_quasineutrality_residual": float(electrostatic_profile.exact_left_throat_quasineutrality_residual),
            "modal_phi_z_exact_right_throat_quasineutrality_residual": float(electrostatic_profile.exact_right_throat_quasineutrality_residual),
            "modal_phi_z_direct_throat_roots_valid": bool(electrostatic_profile.direct_throat_roots_valid),
            "modal_phi_z_direct_central_roots_valid": bool(electrostatic_profile.direct_central_roots_valid),
            "modal_phi_z_direct_central_no_root_cell_count": int(electrostatic_profile.direct_central_no_root_cell_count),
            "modal_phi_z_direct_central_interior_scan_cell_count": int(electrostatic_profile.direct_central_interior_scan_cell_count),
            "modal_phi_z_direct_central_multiple_root_cell_count": int(electrostatic_profile.direct_central_multiple_root_cell_count),
            "modal_phi_z_direct_central_maximum_root_count": int(electrostatic_profile.direct_central_maximum_root_count),
            "modal_phi_z_direct_central_root_scan_points": int(electrostatic_profile.direct_central_root_scan_points),
            "modal_phi_z_direct_central_root_selection_model": "bounded_Egedal_interval_previous_profile_continuation",
            "modal_phi_z_interior_converged": bool(electrostatic_profile.converged),
            "modal_phi_z_exact_throat_population_qualification_status": "unavailable_accepted_references_do_not_define_stationary_exact_throat_ion_density",
            "modal_phi_z_exact_throat_population_qualification_available": False,
            "modal_phi_z_exact_throat_population_classes_included": ["confined_fast_plus_optional_prescribed_reference_charge"],
            "modal_phi_z_exact_throat_population_classes_excluded": ["directed_outbound_Eq63_population", "reflected_population", "local_well_population", "population_requiring_residence_time_closure"],
            "modal_phi_z_midplane_reference_V": float(electrostatic_profile.midplane_reference_V),
            "modal_phi_z_throat_extrapolation_valid": bool(electrostatic_profile.throat_extrapolation_valid),
            "modal_phi_z_effective_potential_throat_boundary_valid": bool(electrostatic_profile.effective_potential_throat_boundary_valid),
        })
    if local_reconstruction is not None:
        active_fixed_eq63_available = bool(fixed_boundary_production and left_eq63_branch is not None and right_eq63_branch is not None)
        metadata.update({
            "modal_lost_ion_parallel_temperature_J": local_reconstruction.lost_ion_parallel_temperature_J,
            "modal_lost_ion_parallel_temperature_keV": (None if local_reconstruction.lost_ion_parallel_temperature_J is None else float(local_reconstruction.lost_ion_parallel_temperature_J / KEV_TO_J)),
            "modal_lost_ion_parallel_temperature_model": (BALDWIN_1972_THROAT_DENSITY_CLOSURE if active_fixed_eq63_available else "unavailable_for_selected_loss_model"),
            "modal_lost_ion_distribution_normalization_reference": None,
            "modal_lost_ion_distribution_normalization_factor": 1.0,
            "modal_lost_ion_distribution_forced_normalization_applied": False,
            "modal_local_distribution_z_v_lambda": local_reconstruction.local_distribution_z_v_lambda,
            "modal_local_distribution_z_v_pitch": local_reconstruction.local_distribution_z_v_pitch,
            "modal_lost_ion_distribution_z_v_lambda": None,
            "modal_lost_ion_geometry_profile_G_z": None,
            "modal_local_density_profile_m3": local_reconstruction.local_density_m3,
            "modal_local_speed_grid_faces_m_s": np.asarray(local_reconstruction.local_speed_grid.faces_m_s, dtype=float),
            "modal_local_speed_grid_diagnostics": local_reconstruction.local_speed_grid_diagnostics,
            "modal_local_speed_refinement_assessed": bool((local_reconstruction.local_speed_grid_diagnostics or {}).get("local_grid_refinement_assessed", False)),
            "modal_local_speed_refinement_converged": bool((local_reconstruction.local_speed_grid_diagnostics or {}).get("local_grid_refinement_converged", False)),
            "modal_local_speed_refinement_failure_reason": (local_reconstruction.local_speed_grid_diagnostics or {}).get("local_grid_refinement_failure_reason"),
            "modal_local_speed_refinement_relative_tolerance": float((local_reconstruction.local_speed_grid_diagnostics or {}).get("local_grid_refinement_relative_tolerance", numerics.local_speed_refinement_relative_tolerance)),
            "modal_local_speed_final_density_relative_change": (local_reconstruction.local_speed_grid_diagnostics or {}).get("local_grid_final_density_relative_change"),
            "modal_local_speed_final_energy_density_relative_change": (local_reconstruction.local_speed_grid_diagnostics or {}).get("local_grid_final_energy_density_relative_change"),
            "modal_local_speed_final_particle_inventory_relative_change": (local_reconstruction.local_speed_grid_diagnostics or {}).get("local_grid_final_particle_inventory_relative_change"),
            "modal_local_speed_final_energy_inventory_relative_change": (local_reconstruction.local_speed_grid_diagnostics or {}).get("local_grid_final_energy_inventory_relative_change"),
            "modal_local_speed_domain_sufficient": bool((local_reconstruction.local_speed_grid_diagnostics or {}).get("local_speed_domain_sufficient", False)),
            "modal_local_speed_overflow_particle_fraction": (local_reconstruction.local_speed_grid_diagnostics or {}).get("local_speed_overflow_particle_fraction"),
            "modal_local_speed_overflow_energy_fraction": (local_reconstruction.local_speed_grid_diagnostics or {}).get("local_speed_overflow_energy_fraction"),
            "modal_local_speed_tail_clipped": bool((local_reconstruction.local_speed_grid_diagnostics or {}).get("local_speed_tail_clipped", False)),
            "modal_directed_eq63_throat_rate_v_lambda_s": directed_eq63_throat_rate_v_lambda_s,
            "modal_directed_eq63_throat_rate_model": ("separate_left_right_per_speed_eq64_matched_fixed_magnetic_boundary" if active_fixed_eq63_available else "unavailable_for_selected_loss_model"),
            "modal_directed_eq63_throat_rate_per_end_s": right_eq63_rate_s,
            "modal_eq63_reference_compatible": active_fixed_eq63_available,
            "modal_active_hot_electrostatic_eq63_population_available": active_fixed_eq63_available,
            "modal_eq63_reference_may_substitute_for_active_population": False,
            "modal_active_material_deposition_rate_available": False,
            "modal_loss_and_global_population_conservation_status": ("fixed_boundary_throat_chain_evaluable_material_deposition_pending_pipeline" if active_fixed_eq63_available else "unavailable_physics_not_a_failed_numerical_residual"),
        })
        
    return {'metadata': metadata}
