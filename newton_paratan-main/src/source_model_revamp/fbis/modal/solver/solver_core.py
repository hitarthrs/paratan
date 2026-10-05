"""One species modal FBIS state assembly across basis, kinetics, losses, electrostatics, and diagnostics"""
from __future__ import annotations
from collections.abc import Callable
from dataclasses import replace
from time import perf_counter
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.beam.source_from_attenuation import AttenuatedMultiEnergyBeamSource
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState, rebuild_fbis_collision_parameter_state_with_fast_density
from source_model_revamp.fbis.species import IonSpecies
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid
from source_model_revamp.fusion.reactions import KEV_TO_J
from source_model_revamp.fbis.modal.basis import ModalBasisConvergenceAssessment, assess_modal_fbis_basis_convergence, build_modal_fbis_basis
from source_model_revamp.fbis.modal.electrostatic_feedback import rescaled_magnetic_profile_for_reference_ratio
from source_model_revamp.fbis.modal.local.pitch import reconstruct_zero_potential_magnetic_reference
from source_model_revamp.fbis.modal.models import COLD_ION_EQ14, EGEDAL_2022_PUBLISHED_ELECTROSTATIC_FBIS_APPROXIMATION, MAGNETIC_ONLY_FBIS, ZERO_POTENTIAL_MAGNETIC_REFERENCE, electrostatic_feedback_model as canonical_electrostatic_feedback_model
from source_model_revamp.fbis.modal.lost_ion_distribution import BALDWIN_1972_THROAT_DENSITY_CLOSURE, EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY, _lost_ion_distribution_eq63, _lost_ion_distribution_eq63_exact_right_throat, ion_loss_model
from source_model_revamp.fbis.modal.types import ModalFBISNumerics, ModalFBISResult, ModalLocalReconstruction
from source_model_revamp.fbis.modal.eq59.types import Eq59WarmStartState
from source_model_revamp.fbis.modal.eq59.state import build_eq59_warm_start_state
from source_model_revamp.fbis.modal.eq59.energy_moments import build_fast_ion_electron_heating_state, fast_ion_electron_heating_metadata
from source_model_revamp.fbis.modal.types import ModalFBISBasis
from source_model_revamp.fbis.modal.solver.solver_loss_rates import _raw_directed_eq63_throat_rate, _modal_mode_eta_integrals
from source_model_revamp.fbis.modal.solver.solver_components import solve_modal_components_and_velocity
from source_model_revamp.fbis.modal.solver.solver_balance import solve_modal_loss_and_current_balance
from source_model_revamp.fbis.modal.solver.solver_fixed_boundary import solve_modal_fixed_boundary_reconstruction
from source_model_revamp.fbis.modal.solver.solver_metadata import build_modal_solver_metadata
from source_model_revamp.integration.evaluation import OperatingPointEvaluationMode
from source_model_revamp.runtime_profile import record_runtime

_FIXED_BOUNDARY_RECONSTRUCTION_REQUIRED = "required"
_FIXED_BOUNDARY_RECONSTRUCTION_DEFERRED = "deferred_eq42_iteration"
_FIXED_BOUNDARY_RECONSTRUCTION_BEST_EFFORT = "best_effort_diagnostic"
_FIXED_BOUNDARY_RECONSTRUCTION_POLICIES = {_FIXED_BOUNDARY_RECONSTRUCTION_REQUIRED, _FIXED_BOUNDARY_RECONSTRUCTION_DEFERRED, _FIXED_BOUNDARY_RECONSTRUCTION_BEST_EFFORT}

def _solve_modal_fbis_state(*, fast_ion_species: IonSpecies, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, pitch_grid: PitchGrid | None, attenuated_source: AttenuatedMultiEnergyBeamSource, volume_m3: float, midplane_area_m2: float, half_length_m: float, mirror_ratio: float, B_tilde_function: Callable[[float], float], zeta_faces: ArrayLike, B_tilde_midpoints: ArrayLike, electron_temperature_J: float, electron_midplane_density_m3: float, electron_collision_density_m3: float, eq70_prescribed_electron_midplane_density_m3: float | None, ion_temperature_J: float, cell_volumes_m3: ArrayLike | None, collision_state: FBISCollisionParameterState, numerics: ModalFBISNumerics, auxiliary_input_power_W: float = 0.0, global_loss_rate_model: str = EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY, electrostatic_feedback_model: str = "magnetic_only_fbis", ion_loss_closure_model: str = EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY, lost_ion_parallel_temperature_model: str = BALDWIN_1972_THROAT_DENSITY_CLOSURE, background_positive_charge_density_m3: ArrayLike | None = None, background_midplane_positive_charge_density_m3: float = 0.0, background_left_throat_positive_charge_density_m3: float = 0.0, background_right_throat_positive_charge_density_m3: float = 0.0, basis: ModalFBISBasis | None = None, basis_convergence_assessment: ModalBasisConvergenceAssessment | None = None, _eq59_warm_start_state: Eq59WarmStartState | None = None, _eq59_warm_start_compatibility_policy: str = "exact_physics", _eq70_initial_potential_energy_J: np.ndarray | None = None, _run_numerical_convergence_assessments: bool = True, _published_throat_drop_override_J: float | None = None, _published_reference_lambda1: float | None = None, _fixed_boundary_reconstruction_policy: str = (_FIXED_BOUNDARY_RECONSTRUCTION_REQUIRED), _evaluation_mode: OperatingPointEvaluationMode = OperatingPointEvaluationMode.FIXED_TEMPERATURE_FINAL, _charge_exchange_sink_state: object | None = None, _external_fast_field_collision_states: tuple[object, ...] = (), _runtime_profile_accumulator: dict[str, object] | None = None) -> ModalFBISResult:
    """Solve one modal fast ion state through basis, Eq 14 or Eq 59, Eq 61, Eq 68, Eq 70 to Eq 72, and loss reconstruction"""
    _evaluation_mode = OperatingPointEvaluationMode.parse(_evaluation_mode)
    runtime_species_id = fast_ion_species.species_id
    electron_midplane_density_m3 = float(electron_midplane_density_m3)
    electron_collision_density_m3 = float(electron_collision_density_m3)
    if not np.isfinite(electron_midplane_density_m3) or electron_midplane_density_m3 <= 0.0:
        raise ValueError("electron_midplane_density_m3 must be finite and positive")
    if not np.isfinite(electron_collision_density_m3) or electron_collision_density_m3 <= 0.0:
        raise ValueError("electron_collision_density_m3 must be finite and positive")
    represented_collision_density = float(collision_state.coulomb_logs.electron_density_m3)
    density_scale = max(abs(electron_collision_density_m3), abs(represented_collision_density), 1.0)
    if abs(electron_collision_density_m3 - represented_collision_density) > 1.0e-12 * density_scale:
        raise ValueError("electron collision density must match the supplied collision state")
    if eq70_prescribed_electron_midplane_density_m3 is not None:
        eq70_prescribed_electron_midplane_density_m3 = float(eq70_prescribed_electron_midplane_density_m3)
        if not np.isfinite(eq70_prescribed_electron_midplane_density_m3) or eq70_prescribed_electron_midplane_density_m3 <= 0.0:
            raise ValueError("Eq 70 prescribed electron midplane density must be finite and positive")
    basis_reused = basis is not None
    basis_runtime_started = perf_counter()
    if basis is None and _run_numerical_convergence_assessments and int(numerics.basis_convergence_levels) >= 3:
        basis_convergence_assessment = assess_modal_fbis_basis_convergence(mirror_ratio=mirror_ratio, B_tilde_function=B_tilde_function, zeta_faces=zeta_faces, B_tilde_midpoints=B_tilde_midpoints, numerics=numerics, eigenvalue_relative_tolerance=numerics.basis_eigenvalue_relative_tolerance, eigenfunction_minimum_overlap=numerics.basis_eigenfunction_overlap_tolerance)
        basis = basis_convergence_assessment.basis_history[0]
    elif basis is None:
        basis = build_modal_fbis_basis(mirror_ratio=mirror_ratio, B_tilde_function=B_tilde_function, zeta_faces=zeta_faces, B_tilde_midpoints=B_tilde_midpoints, numerics=numerics)
    record_runtime(_runtime_profile_accumulator, f"modal_basis_resolution_{runtime_species_id}", perf_counter() - basis_runtime_started, nested=True)
    feedback_model = canonical_electrostatic_feedback_model(electrostatic_feedback_model)
    magnetic_feedback_names = {MAGNETIC_ONLY_FBIS}
    zero_potential_feedback_names = {ZERO_POTENTIAL_MAGNETIC_REFERENCE}
    published_feedback_names = {EGEDAL_2022_PUBLISHED_ELECTROSTATIC_FBIS_APPROXIMATION}
    loss_model = ion_loss_model(ion_loss_closure_model)
    if str(lost_ion_parallel_temperature_model).strip().lower() != BALDWIN_1972_THROAT_DENSITY_CLOSURE:
        raise ValueError("lost_ion_parallel_temperature_model must be baldwin_1972_throat_density_closure")
    fixed_boundary_production = loss_model == EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY
    fixed_boundary_reconstruction_policy = str(_fixed_boundary_reconstruction_policy).strip().lower()
    if (fixed_boundary_reconstruction_policy not in _FIXED_BOUNDARY_RECONSTRUCTION_POLICIES):
        raise ValueError("_fixed_boundary_reconstruction_policy must be required, deferred_eq42_iteration, or best_effort_diagnostic")
    active_global_loss_rate_model = (EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY if fixed_boundary_production else str(global_loss_rate_model))
    # The literal zero potential reference reuses the same magnetic and velocity solve with local electrostatics disabled
    if feedback_model in zero_potential_feedback_names:
        if pitch_grid is None or cell_volumes_m3 is None:
            raise ValueError("zero_potential_magnetic_reference requires pitch and axial-volume grids")
        magnetic_reference = _solve_modal_fbis_state(
            fast_ion_species=fast_ion_species,
            speed_grid=speed_grid,
            lambda_grid=lambda_grid,
            pitch_grid=None,
            attenuated_source=attenuated_source,
            volume_m3=volume_m3,
            midplane_area_m2=midplane_area_m2,
            half_length_m=half_length_m,
            mirror_ratio=mirror_ratio,
            B_tilde_function=B_tilde_function,
            zeta_faces=zeta_faces,
            B_tilde_midpoints=B_tilde_midpoints,
            electron_temperature_J=electron_temperature_J,
            electron_midplane_density_m3=electron_midplane_density_m3,
            electron_collision_density_m3=electron_collision_density_m3,
            eq70_prescribed_electron_midplane_density_m3=eq70_prescribed_electron_midplane_density_m3,
            ion_temperature_J=ion_temperature_J,
            cell_volumes_m3=None,
            collision_state=collision_state,
            numerics=numerics,
            auxiliary_input_power_W=auxiliary_input_power_W,
            global_loss_rate_model=active_global_loss_rate_model,
            electrostatic_feedback_model=MAGNETIC_ONLY_FBIS,
            ion_loss_closure_model=loss_model,
            lost_ion_parallel_temperature_model=lost_ion_parallel_temperature_model,
            background_positive_charge_density_m3=background_positive_charge_density_m3,
            background_midplane_positive_charge_density_m3=background_midplane_positive_charge_density_m3,
            background_left_throat_positive_charge_density_m3=background_left_throat_positive_charge_density_m3,
            background_right_throat_positive_charge_density_m3=background_right_throat_positive_charge_density_m3,
            basis=basis,
            basis_convergence_assessment=basis_convergence_assessment,
            _eq59_warm_start_state=_eq59_warm_start_state,
            _eq59_warm_start_compatibility_policy=_eq59_warm_start_compatibility_policy,
            _eq70_initial_potential_energy_J=None,
            _run_numerical_convergence_assessments=(_run_numerical_convergence_assessments),
            _fixed_boundary_reconstruction_policy=(fixed_boundary_reconstruction_policy),
            _evaluation_mode=_evaluation_mode,
            _charge_exchange_sink_state=_charge_exchange_sink_state,
            _external_fast_field_collision_states=_external_fast_field_collision_states,
            _runtime_profile_accumulator=_runtime_profile_accumulator,
        )
        eta_to_local_normalization = float(basis.volume_geometry_integral / basis.eta_lambda_map.lambda_normalization)
        zero_local_lambda, zero_local_pitch, zero_density, zero_diagnostics = (reconstruct_zero_potential_magnetic_reference(
                speed_grid=speed_grid,
                lambda_grid=lambda_grid,
                pitch_grid=pitch_grid,
                base_distribution_v_lambda=magnetic_reference.distribution_v_lambda,
                B_tilde=np.asarray(B_tilde_midpoints, dtype=float),
                cell_volumes_m3=np.asarray(cell_volumes_m3, dtype=float),
                mirror_ratio=mirror_ratio,
                eta_to_local_phase_space_normalization=eta_to_local_normalization,
                expected_global_inventory_particles=magnetic_reference.inventory_particles,
                quadrature_order=numerics.local_velocity_quadrature_order,
                particle_mass_kg=fast_ion_species.mass_kg,
            ))
        cold_eq63_compatible = str(numerics.velocity_solution_model).strip().lower() == COLD_ION_EQ14
        lost_distribution: np.ndarray | None = None
        lost_geometry: np.ndarray | None = None
        lost_temperature: float | None = None
        raw_eq63_rate: np.ndarray | None = None
        raw_eq63_per_end_rate: float | None = None
        eq63_to_eq61_relative_difference: float | None = None
        if cold_eq63_compatible:
            lost_temperature = float(ion_temperature_J)
            lost_distribution, lost_geometry = _lost_ion_distribution_eq63(speed_grid=speed_grid, lambda_grid=lambda_grid, dfdlambda_boundary_v=magnetic_reference.dF_dlambda_boundary_v, zeta_faces=np.asarray(zeta_faces, dtype=float), B_tilde_midpoints=np.asarray(B_tilde_midpoints, dtype=float), mirror_ratio=mirror_ratio, collision_state=collision_state, half_length_m=half_length_m, lost_ion_parallel_temperature_J=lost_temperature, particle_mass_kg=fast_ion_species.mass_kg, collision_operator_state=magnetic_reference.eq59_collision_operator_state)
            exact_throat_distribution = _lost_ion_distribution_eq63_exact_right_throat(speed_grid=speed_grid, lambda_grid=lambda_grid, dfdlambda_boundary_v=magnetic_reference.dF_dlambda_boundary_v, zeta_faces=np.asarray(zeta_faces, dtype=float), B_tilde_midpoints=np.asarray(B_tilde_midpoints, dtype=float), mirror_ratio=mirror_ratio, collision_state=collision_state, half_length_m=half_length_m, lost_ion_parallel_temperature_J=lost_temperature, particle_mass_kg=fast_ion_species.mass_kg, collision_operator_state=magnetic_reference.eq59_collision_operator_state)
            raw_eq63_rate = _raw_directed_eq63_throat_rate(speed_grid=speed_grid, lambda_grid=lambda_grid, lost_ion_distribution_z_v_lambda=exact_throat_distribution, mirror_ratio=mirror_ratio, midplane_area_m2=midplane_area_m2)
            raw_eq63_per_end_rate = float(np.sum(raw_eq63_rate))
            eq61_per_end = float(magnetic_reference.metadata["modal_single_end_eq61_loss_rate_s"])
            eq63_to_eq61_relative_difference = ((raw_eq63_per_end_rate - eq61_per_end) / eq61_per_end if eq61_per_end > 0.0 else None)
        local_reconstruction = ModalLocalReconstruction(
            local_distribution_z_v_lambda=zero_local_lambda,
            local_distribution_z_v_pitch=zero_local_pitch,
            lost_ion_distribution_z_v_lambda=lost_distribution,
            lost_ion_parallel_temperature_J=lost_temperature,
            lost_ion_geometry_profile_G_z=lost_geometry,
            local_density_m3=zero_density,
            local_speed_grid=speed_grid,
            local_speed_grid_diagnostics={
                "invariant_speed_max_m_s": float(speed_grid.faces_m_s[-1]),
                "required_local_physical_speed_max_m_s": float(speed_grid.faces_m_s[-1]),
                "local_physical_speed_max_m_s": float(speed_grid.faces_m_s[-1]),
                "local_speed_appended_cell_count": 0,
                "local_speed_domain_sufficient": True,
                "local_speed_grid_model": "literal_zero_potential_invariant_grid_identity",
                "local_speed_overflow_particle_fraction": 0.0,
                "local_speed_overflow_energy_fraction": 0.0,
                "local_speed_tail_clipped": False,
            },
        )
        zero_invalid_prefixes = ("modal_phi_z_", "modal_electron_", "modal_current_", "modal_power_balance_", "modal_per_lost_ion_", "modal_egedal_reduced_power_balance_", "modal_loss_power_to_", "modal_source_deposited_power_to_loss_", "modal_total_input_power_to_loss_")
        zero_invalid_keys = {"electron_parent_maxwellian_n0_m3", "electron_volume_average_density_m3", "electron_density_profile_m3", "electron_density_closure_model", "electron_density_closure_converged", "electron_density_closure_failure_reason", "modal_total_wall_power_loss_W", "modal_total_input_power_W", "modal_power_balance_relative_error"}
        metadata = {key: value for key, value in magnetic_reference.metadata.items() if key not in zero_invalid_keys and not key.startswith(zero_invalid_prefixes)}
        metadata.update({
                "electrostatic_feedback_model": "zero_potential_magnetic_reference",
                "electrostatic_feedback_reference": "literal_Phi_zero_fixed_magnetic_boundary_reference",
                "electrostatic_feedback_is_exact": True,
                "electrostatic_feedback_iteration_converged": None,
                "electrostatic_feedback_iteration_history_throat_energy_J": [],
                "kinetic_eq70_profile_applicable": False,
                "kinetic_current_balance_applicable": False,
                "modal_electrostatic_phi_z_status": "not_applicable_literal_zero_potential_reference",
                "modal_phi_z_iterations": 0,
                "modal_phi_z_converged": None,
                "modal_phi_z_failure_reason": None,
                "modal_phi_z_profile_V": np.zeros(np.asarray(B_tilde_midpoints, dtype=float).shape),
                "modal_phi_z_min_V": 0.0,
                "modal_phi_z_max_V": 0.0,
                "modal_phi_z_residual_history": [],
                "modal_phi_z_relative_residual_profile": None,
                "modal_phi_z_absolute_residual_profile_m3": None,
                "modal_phi_z_electron_density_profile_m3": None,
                "modal_phi_z_ion_density_profile_m3": None,
                "modal_phi_z_electrostatic_node_potential_V": None,
                "modal_phi_z_throat_potential_energy_J": 0.0,
                "modal_phi_z_throat_potential_left_energy_J": 0.0,
                "modal_phi_z_throat_potential_right_energy_J": 0.0,
                "modal_phi_z_eq70_evaluated": False,
                "modal_phi_z_eq71_evaluated": False,
                "modal_phi_z_eq72_evaluated": False,
                "modal_phi_z_potential_identically_zero": True,
                "modal_local_reconstruction_model": "literal_zero_potential_fixed_magnetic_moment_mapping",
                "modal_local_distribution_z_v_lambda": zero_local_lambda,
                "modal_local_distribution_z_v_pitch": zero_local_pitch,
                "modal_local_density_profile_m3": zero_density,
                "modal_local_speed_grid_faces_m_s": np.asarray(speed_grid.faces_m_s, dtype=float),
                "modal_local_speed_grid_diagnostics": local_reconstruction.local_speed_grid_diagnostics,
                "modal_zero_potential_local_inventory_particles": zero_diagnostics["local_inventory_particles"],
                "modal_zero_potential_global_inventory_particles": zero_diagnostics["expected_global_inventory_particles"],
                "modal_zero_potential_local_to_global_inventory_relative_error": zero_diagnostics["local_to_global_inventory_relative_error"],
                "modal_wall_barrier_energy_J": 0.0,
                "modal_wall_barrier_over_Te": 0.0,
                "modal_electron_wall_potential_relative_to_midplane_V": None,
                "modal_electron_current_balance_relative_error": None,
                "modal_electron_current_balance_residual_A": None,
                "modal_current_balance_relative_error": None,
                "modal_electron_wall_power_loss_W": None,
                "modal_ion_wall_power_loss_W": None,
                "modal_ion_barrier_power_loss_W": 0.0,
                "modal_total_wall_power_loss_W": None,
                "modal_power_balance_relative_error": None,
                "modal_end_loss_power_terms_available": False,
                "modal_end_loss_power_unavailability_reason": "zero_potential_reference_has_no_electron_current_or_wall_power_closure",
                "modal_electron_loss_equation": "not_applicable_literal_zero_potential_reference",
                "modal_lost_ion_distribution_model": ("egedal_eq63_eq64_cold_lorentz_fixed_magnetic_reference" if cold_eq63_compatible else "unavailable_conservative_scope_hot_fixed_boundary_eq63_eq64_not_implemented"),
                "modal_lost_ion_distribution_z_v_lambda": lost_distribution,
                "modal_lost_ion_geometry_profile_G_z": lost_geometry,
                "modal_lost_ion_parallel_temperature_J": lost_temperature,
                "modal_lost_ion_parallel_temperature_keV": (None if lost_temperature is None else float(lost_temperature / KEV_TO_J)),
                "modal_lost_ion_parallel_temperature_model": ("prescribed_background_ion_temperature_as_nonpredictive_eq63_reference_parameter" if cold_eq63_compatible else "unavailable_hot_eq59_eq63_eq64_active_mapping_not_implemented"),
                "modal_lost_ion_distribution_normalization_factor": 1.0,
                "modal_lost_ion_distribution_forced_normalization_applied": False,
                "modal_lost_ion_distribution_normalization_reference": None,
                "modal_directed_eq63_throat_rate_v_lambda_s": raw_eq63_rate,
                "modal_directed_eq63_throat_rate_per_end_s": raw_eq63_per_end_rate,
                "modal_directed_eq63_throat_geometry_sample": ("exact_right_throat_face_zeta_plus_one" if cold_eq63_compatible else None),
                "modal_directed_eq63_throat_rate_model": ("raw_exact_throat_face_outbound_one_sign_flux_times_one_throat_area_no_normalization" if cold_eq63_compatible else "not_evaluable_hot_eq59"),
                "modal_eq63_reference_compatible": cold_eq63_compatible,
                "modal_eq63_eq64_published_scope": "hot_eq59_fixed_magnetic_boundary_with_U_equals_E_plus_ePhi_energy_shift",
                "modal_eq64_g_tilde_1_printed_form_status": "published_g_tilde_1_prefactor_not_activated_for_the_unavailable_moving_boundary_population",
                "modal_eq63_eq64_active_implementation_status": "cold_zero_potential_diagnostic_subset_only",
                "modal_eq63_to_eq61_single_end_relative_difference": eq63_to_eq61_relative_difference,
                "modal_active_hot_electrostatic_eq63_population_available": False,
                "modal_eq63_reference_may_substitute_for_active_population": False,
                "modal_active_material_deposition_rate_available": False,
                "modal_loss_convention_validation_passed": False,
                "modal_loss_convention_failure_reason": "electron_current_and_material_deposition_closures_not_applicable_to_literal_zero_reference",
                "modal_loss_and_global_population_conservation_status": ("cold_magnetic_eq61_eq63_reference_comparison_evaluable" if cold_eq63_compatible else "unavailable_physics_hot_eq59_eq63_eq64_active_mapping_not_implemented"),
            })
        return replace(magnetic_reference, pitch_grid=pitch_grid, ion_wall_power_loss_W=None, electron_wall_power_loss_W=None, wall_barrier_energy_J=0.0, normalized_wall_barrier=0.0, electron_wall_potential_relative_to_midplane_V=None, electrostatic_profile=None, local_reconstruction=local_reconstruction, metadata=metadata)
    # The low energy first eigenvalue feedback path iterates on the throat drop using the previous state as a continuation seed
    if feedback_model in published_feedback_names and _published_throat_drop_override_J is None:
        magnetic_reference = _solve_modal_fbis_state(
            fast_ion_species=fast_ion_species,
            speed_grid=speed_grid,
            lambda_grid=lambda_grid,
            pitch_grid=pitch_grid,
            attenuated_source=attenuated_source,
            volume_m3=volume_m3,
            midplane_area_m2=midplane_area_m2,
            half_length_m=half_length_m,
            mirror_ratio=mirror_ratio,
            B_tilde_function=B_tilde_function,
            zeta_faces=zeta_faces,
            B_tilde_midpoints=B_tilde_midpoints,
            electron_temperature_J=electron_temperature_J,
            electron_midplane_density_m3=electron_midplane_density_m3,
            electron_collision_density_m3=electron_collision_density_m3,
            eq70_prescribed_electron_midplane_density_m3=eq70_prescribed_electron_midplane_density_m3,
            ion_temperature_J=ion_temperature_J,
            cell_volumes_m3=cell_volumes_m3,
            collision_state=collision_state,
            numerics=numerics,
            auxiliary_input_power_W=auxiliary_input_power_W,
            global_loss_rate_model=active_global_loss_rate_model,
            electrostatic_feedback_model=MAGNETIC_ONLY_FBIS,
            ion_loss_closure_model=loss_model,
            lost_ion_parallel_temperature_model=lost_ion_parallel_temperature_model,
            background_positive_charge_density_m3=background_positive_charge_density_m3,
            background_midplane_positive_charge_density_m3=background_midplane_positive_charge_density_m3,
            background_left_throat_positive_charge_density_m3=background_left_throat_positive_charge_density_m3,
            background_right_throat_positive_charge_density_m3=background_right_throat_positive_charge_density_m3,
            basis=basis,
            basis_convergence_assessment=basis_convergence_assessment,
            _eq59_warm_start_state=_eq59_warm_start_state,
            _eq59_warm_start_compatibility_policy=_eq59_warm_start_compatibility_policy,
            _eq70_initial_potential_energy_J=_eq70_initial_potential_energy_J,
            _run_numerical_convergence_assessments=False,
            _fixed_boundary_reconstruction_policy=(fixed_boundary_reconstruction_policy),
            _evaluation_mode=_evaluation_mode,
            _charge_exchange_sink_state=_charge_exchange_sink_state,
            _external_fast_field_collision_states=_external_fast_field_collision_states,
            _runtime_profile_accumulator=_runtime_profile_accumulator,
        )
        if magnetic_reference.electrostatic_profile is None:
            raise ValueError("published electrostatic feedback requires an Eq 70-72 profile")
        if _published_reference_lambda1 is None:
            reference_profile = rescaled_magnetic_profile_for_reference_ratio(B_tilde_function, mirror_ratio)
            active_midpoints = np.asarray(B_tilde_midpoints, dtype=float)
            reference_midpoints = 1.0 + 0.5 * (active_midpoints - 1.0) / (float(mirror_ratio) - 1.0)
            reference_basis = build_modal_fbis_basis(mirror_ratio=1.5, B_tilde_function=reference_profile, zeta_faces=zeta_faces, B_tilde_midpoints=reference_midpoints, numerics=numerics)
            reference_lambda1 = float(reference_basis.physical_basis.eigenvalues[0])
        else:
            reference_lambda1 = float(_published_reference_lambda1)
        seed = float(magnetic_reference.electrostatic_profile.throat_potential_energy_J)
        feedback_history: list[float] = []
        latest: ModalFBISResult | None = None
        feedback_converged = False
        iteration_limit = max(int(numerics.electrostatic_feedback_iterations), 1)
        feedback_tolerance = float(numerics.electrostatic_feedback_relative_tolerance)
        feedback_warm_state = magnetic_reference.eq59_warm_start_state
        feedback_potential_energy_J = np.asarray(magnetic_reference.electrostatic_profile.potential_energy_J, dtype=float)
        for _ in range(iteration_limit):
            latest = _solve_modal_fbis_state(
                fast_ion_species=fast_ion_species,
                speed_grid=speed_grid,
                lambda_grid=lambda_grid,
                pitch_grid=pitch_grid,
                attenuated_source=attenuated_source,
                volume_m3=volume_m3,
                midplane_area_m2=midplane_area_m2,
                half_length_m=half_length_m,
                mirror_ratio=mirror_ratio,
                B_tilde_function=B_tilde_function,
                zeta_faces=zeta_faces,
                B_tilde_midpoints=B_tilde_midpoints,
                electron_temperature_J=electron_temperature_J,
                electron_midplane_density_m3=electron_midplane_density_m3,
                electron_collision_density_m3=electron_collision_density_m3,
                eq70_prescribed_electron_midplane_density_m3=eq70_prescribed_electron_midplane_density_m3,
                ion_temperature_J=ion_temperature_J,
                cell_volumes_m3=cell_volumes_m3,
                collision_state=collision_state,
                numerics=numerics,
                auxiliary_input_power_W=auxiliary_input_power_W,
                global_loss_rate_model=(active_global_loss_rate_model if fixed_boundary_production else "egedal_2022_published_electrostatic_modal_sink"),
                electrostatic_feedback_model=feedback_model,
                ion_loss_closure_model=loss_model,
                lost_ion_parallel_temperature_model=lost_ion_parallel_temperature_model,
                background_positive_charge_density_m3=background_positive_charge_density_m3,
                background_midplane_positive_charge_density_m3=background_midplane_positive_charge_density_m3,
                background_left_throat_positive_charge_density_m3=background_left_throat_positive_charge_density_m3,
                background_right_throat_positive_charge_density_m3=background_right_throat_positive_charge_density_m3,
                basis=basis,
                basis_convergence_assessment=basis_convergence_assessment,
                _eq59_warm_start_state=feedback_warm_state,
                _eq59_warm_start_compatibility_policy="closure_continuation",
                _eq70_initial_potential_energy_J=feedback_potential_energy_J,
                _run_numerical_convergence_assessments=False,
                _published_throat_drop_override_J=seed,
                _published_reference_lambda1=reference_lambda1,
                _fixed_boundary_reconstruction_policy=(fixed_boundary_reconstruction_policy),
                _evaluation_mode=_evaluation_mode,
                _charge_exchange_sink_state=_charge_exchange_sink_state,
                _external_fast_field_collision_states=_external_fast_field_collision_states,
                _runtime_profile_accumulator=_runtime_profile_accumulator,
            )
            feedback_warm_state = latest.eq59_warm_start_state
            if latest.electrostatic_profile is not None:
                feedback_potential_energy_J = np.asarray(latest.electrostatic_profile.potential_energy_J, dtype=float)
            if latest.electrostatic_profile is None:
                break
            updated = float(latest.electrostatic_profile.throat_potential_energy_J)
            feedback_history.append(updated)
            change = abs(updated - seed) / max(abs(updated), abs(seed), np.finfo(float).tiny)
            seed = updated
            if change <= feedback_tolerance:
                feedback_converged = True
                break
        if latest is None:
            raise RuntimeError("published electrostatic feedback iteration did not produce a state")
        if _run_numerical_convergence_assessments:
            latest = _solve_modal_fbis_state(
                fast_ion_species=fast_ion_species,
                speed_grid=speed_grid,
                lambda_grid=lambda_grid,
                pitch_grid=pitch_grid,
                attenuated_source=attenuated_source,
                volume_m3=volume_m3,
                midplane_area_m2=midplane_area_m2,
                half_length_m=half_length_m,
                mirror_ratio=mirror_ratio,
                B_tilde_function=B_tilde_function,
                zeta_faces=zeta_faces,
                B_tilde_midpoints=B_tilde_midpoints,
                electron_temperature_J=electron_temperature_J,
                electron_midplane_density_m3=electron_midplane_density_m3,
                electron_collision_density_m3=electron_collision_density_m3,
                eq70_prescribed_electron_midplane_density_m3=eq70_prescribed_electron_midplane_density_m3,
                ion_temperature_J=ion_temperature_J,
                cell_volumes_m3=cell_volumes_m3,
                collision_state=collision_state,
                numerics=numerics,
                auxiliary_input_power_W=auxiliary_input_power_W,
                global_loss_rate_model=(active_global_loss_rate_model if fixed_boundary_production else "egedal_2022_published_electrostatic_modal_sink"),
                electrostatic_feedback_model=feedback_model,
                ion_loss_closure_model=loss_model,
                lost_ion_parallel_temperature_model=(lost_ion_parallel_temperature_model),
                background_positive_charge_density_m3=(background_positive_charge_density_m3),
                background_midplane_positive_charge_density_m3=(background_midplane_positive_charge_density_m3),
                background_left_throat_positive_charge_density_m3=(background_left_throat_positive_charge_density_m3),
                background_right_throat_positive_charge_density_m3=(background_right_throat_positive_charge_density_m3),
                basis=basis,
                basis_convergence_assessment=basis_convergence_assessment,
                _eq59_warm_start_state=latest.eq59_warm_start_state,
                _eq59_warm_start_compatibility_policy="closure_continuation",
                _eq70_initial_potential_energy_J=(None if latest.electrostatic_profile is None else np.asarray(latest.electrostatic_profile.potential_energy_J, dtype=float)),
                _run_numerical_convergence_assessments=True,
                _published_throat_drop_override_J=seed,
                _published_reference_lambda1=reference_lambda1,
                _fixed_boundary_reconstruction_policy=(fixed_boundary_reconstruction_policy),
                _evaluation_mode=_evaluation_mode,
                _charge_exchange_sink_state=_charge_exchange_sink_state,
                _external_fast_field_collision_states=_external_fast_field_collision_states,
                _runtime_profile_accumulator=_runtime_profile_accumulator,
            )
        metadata = dict(latest.metadata)
        metadata.update({
                "electrostatic_feedback_model": "egedal_2022_published_electrostatic_fbis_approximation",
                "electrostatic_feedback_reference": "Egedal_et_al_Nuclear_Fusion_62_126053_2022_pp15_16",
                "electrostatic_feedback_is_exact": False,
                "electrostatic_feedback_iteration_history_throat_energy_J": feedback_history,
                "electrostatic_feedback_iteration_converged": feedback_converged,
                "published_lambda1_reference_basis_reused": bool(_published_reference_lambda1 is not None),
                "magnetic_only_comparison_available": True,
                "magnetic_only_particle_inventory": float(magnetic_reference.inventory_particles),
                "magnetic_only_loss_rate": float(magnetic_reference.ion_particle_loss_rate_s),
                "electrostatic_approximation_particle_inventory": float(latest.inventory_particles),
                "electrostatic_approximation_loss_rate": float(latest.ion_particle_loss_rate_s),
            })
        return replace(latest, metadata=metadata)
    reconstructed_flux_tube_volume_m3 = (float(midplane_area_m2) * float(half_length_m) * float(basis.volume_geometry_integral))
    geometry_volume_relative_error = (reconstructed_flux_tube_volume_m3 - float(volume_m3)) / max(abs(float(volume_m3)), np.finfo(float).tiny)
    if abs(geometry_volume_relative_error) > 1.0e-10:
        raise ValueError("modal geometry violates V = A0 * half_length * integral(dzeta/B_tilde)")
    # The common solve path assembles the source and velocity solution before loss and electrostatic reconstruction
    components_runtime_started = perf_counter()
    _solve_modal_components_and_velocity_result = solve_modal_components_and_velocity(locals())
    record_runtime(_runtime_profile_accumulator, f"eq59_components_and_velocity_{runtime_species_id}", perf_counter() - components_runtime_started, nested=True)
    active_eigenvalues_by_speed = _solve_modal_components_and_velocity_result['active_eigenvalues_by_speed']
    charge_exchange_loss_frequency_s = _solve_modal_components_and_velocity_result['charge_exchange_loss_frequency_s']
    charge_exchange_sink_active = _solve_modal_components_and_velocity_result['charge_exchange_sink_active']
    charge_exchange_sink_rate_s = _solve_modal_components_and_velocity_result['charge_exchange_sink_rate_s']
    charge_exchange_sink_energy_W = _solve_modal_components_and_velocity_result['charge_exchange_sink_energy_W']
    charge_exchange_sink_model = _solve_modal_components_and_velocity_result['charge_exchange_sink_model']
    charge_exchange_sink_reference_identity_error = _solve_modal_components_and_velocity_result['charge_exchange_sink_reference_identity_error']
    component_coeffs = _solve_modal_components_and_velocity_result['component_coeffs']
    component_results = _solve_modal_components_and_velocity_result['component_results']
    component_speeds = _solve_modal_components_and_velocity_result['component_speeds']
    eq59_diagnostics = _solve_modal_components_and_velocity_result['eq59_diagnostics']
    eq59_collision_operator_state = _solve_modal_components_and_velocity_result['eq59_collision_operator_state']
    if eq59_collision_operator_state is not None:
        collision_state = rebuild_fbis_collision_parameter_state_with_fast_density(collision_state, eq59_collision_operator_state.fast_self_physical_density_m3)
    eq59_high_speed_tail_energy_fraction = _solve_modal_components_and_velocity_result['eq59_high_speed_tail_energy_fraction']
    eq59_high_speed_tail_population_fraction = _solve_modal_components_and_velocity_result['eq59_high_speed_tail_population_fraction']
    eq59_production_grid_index = _solve_modal_components_and_velocity_result['eq59_production_grid_index']
    eq59_single_grid_tail_indicator_passed = _solve_modal_components_and_velocity_result['eq59_single_grid_tail_indicator_passed']
    eq59_speed_assessment = _solve_modal_components_and_velocity_result['eq59_speed_assessment']
    f_eta = _solve_modal_components_and_velocity_result['f_eta']
    invariant_energy_cell_population_weights = _solve_modal_components_and_velocity_result['invariant_energy_cell_population_weights']
    modal_prompt_loss_birth_power_W = _solve_modal_components_and_velocity_result['modal_prompt_loss_birth_power_W']
    modal_prompt_loss_birth_rate_s = _solve_modal_components_and_velocity_result['modal_prompt_loss_birth_rate_s']
    modal_source_birth_current_A = _solve_modal_components_and_velocity_result['modal_source_birth_current_A']
    modal_source_birth_rate_s = _solve_modal_components_and_velocity_result['modal_source_birth_rate_s']
    modal_source_deposited_birth_power_W = _solve_modal_components_and_velocity_result['modal_source_deposited_birth_power_W']
    modal_source_mean_birth_energy_J = _solve_modal_components_and_velocity_result['modal_source_mean_birth_energy_J']
    modal_total = _solve_modal_components_and_velocity_result['modal_total']
    modal_total_birth_current_A = _solve_modal_components_and_velocity_result['modal_total_birth_current_A']
    modal_total_birth_rate_s = _solve_modal_components_and_velocity_result['modal_total_birth_rate_s']
    modal_total_deposited_birth_power_W = _solve_modal_components_and_velocity_result['modal_total_deposited_birth_power_W']
    n_modes = _solve_modal_components_and_velocity_result['n_modes']
    published_heuristic = _solve_modal_components_and_velocity_result['published_heuristic']
    reconstruction_correction = _solve_modal_components_and_velocity_result['reconstruction_correction']
    retained_mode_assessment = _solve_modal_components_and_velocity_result['retained_mode_assessment']
    rosen_coeffs = _solve_modal_components_and_velocity_result['rosen_coeffs']
    rosen_history = _solve_modal_components_and_velocity_result['rosen_history']
    rosen_status = _solve_modal_components_and_velocity_result['rosen_status']
    velocity_solution_equation = _solve_modal_components_and_velocity_result['velocity_solution_equation']
    balance_runtime_started = perf_counter()
    _solve_modal_loss_and_current_balance_result = solve_modal_loss_and_current_balance(locals())
    record_runtime(_runtime_profile_accumulator, f"modal_loss_and_current_balance_{runtime_species_id}", perf_counter() - balance_runtime_started, nested=True)
    active_loss = _solve_modal_loss_and_current_balance_result['active_loss']
    auxiliary_power = _solve_modal_loss_and_current_balance_result['auxiliary_power']
    barrier = _solve_modal_loss_and_current_balance_result['barrier']
    boundary_flux_single_density_v = _solve_modal_loss_and_current_balance_result['boundary_flux_single_density_v']
    boundary_flux_single_midplane_power_loss_W = _solve_modal_loss_and_current_balance_result['boundary_flux_single_midplane_power_loss_W']
    boundary_flux_single_particle_loss_rate_s = _solve_modal_loss_and_current_balance_result['boundary_flux_single_particle_loss_rate_s']
    boundary_flux_two_boundary_midplane_power_loss_W = _solve_modal_loss_and_current_balance_result['boundary_flux_two_boundary_midplane_power_loss_W']
    boundary_flux_two_boundary_particle_loss_rate_s = _solve_modal_loss_and_current_balance_result['boundary_flux_two_boundary_particle_loss_rate_s']
    confined_collisional_confinement = _solve_modal_loss_and_current_balance_result['confined_collisional_confinement']
    confined_total_device_loss_rate_s = _solve_modal_loss_and_current_balance_result['confined_total_device_loss_rate_s']
    confined_total_device_midplane_power_W = _solve_modal_loss_and_current_balance_result['confined_total_device_midplane_power_W']
    confinement = _solve_modal_loss_and_current_balance_result['confinement']
    conservative_mode_slopes = _solve_modal_loss_and_current_balance_result['conservative_mode_slopes']
    current_balance = _solve_modal_loss_and_current_balance_result['current_balance']
    current_balance_relative_error = _solve_modal_loss_and_current_balance_result['current_balance_relative_error']
    density = _solve_modal_loss_and_current_balance_result['density']
    dfdlambda = _solve_modal_loss_and_current_balance_result['dfdlambda']
    direct_cell_average_mode_slopes = _solve_modal_loss_and_current_balance_result['direct_cell_average_mode_slopes']
    eq61_boundary_reconstruction_diagnostics = _solve_modal_loss_and_current_balance_result['eq61_boundary_reconstruction_diagnostics']
    effective_ion_temperature_J = _solve_modal_loss_and_current_balance_result['effective_ion_temperature_J']
    eigen_audit = _solve_modal_loss_and_current_balance_result['eigen_audit']
    eigenvalue_first_mode_current_A = _solve_modal_loss_and_current_balance_result['eigenvalue_first_mode_current_A']
    eigenvalue_first_mode_mean_energy_J = _solve_modal_loss_and_current_balance_result['eigenvalue_first_mode_mean_energy_J']
    eigenvalue_first_mode_midplane_power_W = _solve_modal_loss_and_current_balance_result['eigenvalue_first_mode_midplane_power_W']
    eigenvalue_first_mode_particle_rate_s = _solve_modal_loss_and_current_balance_result['eigenvalue_first_mode_particle_rate_s']
    eigenvalue_mode_midplane_power_W = _solve_modal_loss_and_current_balance_result['eigenvalue_mode_midplane_power_W']
    eigenvalue_mode_particle_rates_s = _solve_modal_loss_and_current_balance_result['eigenvalue_mode_particle_rates_s']
    eigenvalue_sink_current_A = _solve_modal_loss_and_current_balance_result['eigenvalue_sink_current_A']
    eigenvalue_sink_mean_energy_J = _solve_modal_loss_and_current_balance_result['eigenvalue_sink_mean_energy_J']
    eigenvalue_sink_midplane_power_W = _solve_modal_loss_and_current_balance_result['eigenvalue_sink_midplane_power_W']
    eigenvalue_sink_particle_rate_s = _solve_modal_loss_and_current_balance_result['eigenvalue_sink_particle_rate_s']
    eigenvalue_sink_to_birth_ratio = _solve_modal_loss_and_current_balance_result['eigenvalue_sink_to_birth_ratio']
    electron_current_relative_error = _solve_modal_loss_and_current_balance_result['electron_current_relative_error']
    electron_current_residual_A = _solve_modal_loss_and_current_balance_result['electron_current_residual_A']
    electron_loss_current_A = _solve_modal_loss_and_current_balance_result['electron_loss_current_A']
    electron_loss_prefactor_s = _solve_modal_loss_and_current_balance_result['electron_loss_prefactor_s']
    electron_power_W = _solve_modal_loss_and_current_balance_result['electron_power_W']
    f_lambda = _solve_modal_loss_and_current_balance_result['f_lambda']
    first_slope = _solve_modal_loss_and_current_balance_result['first_slope']
    flux = _solve_modal_loss_and_current_balance_result['flux']
    inventory = _solve_modal_loss_and_current_balance_result['inventory']
    ion_barrier_power_W = _solve_modal_loss_and_current_balance_result['ion_barrier_power_W']
    ion_loss_current_A = _solve_modal_loss_and_current_balance_result['ion_loss_current_A']
    ion_midplane_power_loss_W = _solve_modal_loss_and_current_balance_result['ion_midplane_power_loss_W']
    ion_particle_loss_rate_s = _solve_modal_loss_and_current_balance_result['ion_particle_loss_rate_s']
    ion_wall_power_W = _solve_modal_loss_and_current_balance_result['ion_wall_power_W']
    left_electron_loss_rate_s = _solve_modal_loss_and_current_balance_result['left_electron_loss_rate_s']
    left_ion_loss_rate_s = _solve_modal_loss_and_current_balance_result['left_ion_loss_rate_s']
    loss_convention_validation_passed = _solve_modal_loss_and_current_balance_result['loss_convention_validation_passed']
    loss_equivalent_birth_power_W = _solve_modal_loss_and_current_balance_result['loss_equivalent_birth_power_W']
    loss_power_to_source_power_ratio = _solve_modal_loss_and_current_balance_result['loss_power_to_source_power_ratio']
    loss_rate_density = _solve_modal_loss_and_current_balance_result['loss_rate_density']
    loss_to_source_current_ratio = _solve_modal_loss_and_current_balance_result['loss_to_source_current_ratio']
    loss_to_source_particle_ratio = _solve_modal_loss_and_current_balance_result['loss_to_source_particle_ratio']
    loss_validation_failures = _solve_modal_loss_and_current_balance_result['loss_validation_failures']
    loss_validation_tolerance = _solve_modal_loss_and_current_balance_result['loss_validation_tolerance']
    mean_loss_energy = _solve_modal_loss_and_current_balance_result['mean_loss_energy']
    modal_per_lost_ion_barrier_energy_over_birth_energy = _solve_modal_loss_and_current_balance_result['modal_per_lost_ion_barrier_energy_over_birth_energy']
    modal_per_lost_ion_electron_energy_over_birth_energy = _solve_modal_loss_and_current_balance_result['modal_per_lost_ion_electron_energy_over_birth_energy']
    modal_per_lost_ion_midplane_energy_over_birth_energy = _solve_modal_loss_and_current_balance_result['modal_per_lost_ion_midplane_energy_over_birth_energy']
    modal_per_lost_ion_total_wall_energy_over_birth_energy = _solve_modal_loss_and_current_balance_result['modal_per_lost_ion_total_wall_energy_over_birth_energy']
    modal_per_lost_ion_wall_energy_over_birth_energy = _solve_modal_loss_and_current_balance_result['modal_per_lost_ion_wall_energy_over_birth_energy']
    normalized_barrier = _solve_modal_loss_and_current_balance_result['normalized_barrier']
    particle_balance_relative_error = _solve_modal_loss_and_current_balance_result['particle_balance_relative_error']
    power_balance_relative_error = _solve_modal_loss_and_current_balance_result['power_balance_relative_error']
    prefactor_over_target = _solve_modal_loss_and_current_balance_result['prefactor_over_target']
    projected_source_component_rates_s = _solve_modal_loss_and_current_balance_result['projected_source_component_rates_s']
    projected_source_particle_rate_s = _solve_modal_loss_and_current_balance_result['projected_source_particle_rate_s']
    projected_source_to_birth_ratio = _solve_modal_loss_and_current_balance_result['projected_source_to_birth_ratio']
    right_electron_loss_rate_s = _solve_modal_loss_and_current_balance_result['right_electron_loss_rate_s']
    right_ion_loss_rate_s = _solve_modal_loss_and_current_balance_result['right_ion_loss_rate_s']
    source_mean_energy = _solve_modal_loss_and_current_balance_result['source_mean_energy']
    source_power_to_loss_power_ratio = _solve_modal_loss_and_current_balance_result['source_power_to_loss_power_ratio']
    source_to_loss_current_ratio = _solve_modal_loss_and_current_balance_result['source_to_loss_current_ratio']
    source_to_loss_particle_ratio = _solve_modal_loss_and_current_balance_result['source_to_loss_particle_ratio']
    tail_factor = _solve_modal_loss_and_current_balance_result['tail_factor']
    target_electron_particle_loss_rate_s = _solve_modal_loss_and_current_balance_result['target_electron_particle_loss_rate_s']
    total_input_power = _solve_modal_loss_and_current_balance_result['total_input_power']
    total_wall_power_W = _solve_modal_loss_and_current_balance_result['total_wall_power_W']
    two_end_eq61_to_birth_ratio = _solve_modal_loss_and_current_balance_result['two_end_eq61_to_birth_ratio']
    two_end_eq61_to_eigenvalue_ratio = _solve_modal_loss_and_current_balance_result['two_end_eq61_to_eigenvalue_ratio']
    y_exp_y = _solve_modal_loss_and_current_balance_result['y_exp_y']
    boundary_runtime_started = perf_counter()
    _solve_modal_fixed_boundary_reconstruction_result = solve_modal_fixed_boundary_reconstruction(locals())
    record_runtime(_runtime_profile_accumulator, f"fixed_boundary_reconstruction_{runtime_species_id}", perf_counter() - boundary_runtime_started, nested=True)
    directed_eq63_throat_rate_v_lambda_s = _solve_modal_fixed_boundary_reconstruction_result['directed_eq63_throat_rate_v_lambda_s']
    electrostatic_profile = _solve_modal_fixed_boundary_reconstruction_result['electrostatic_profile']
    fixed_boundary_reconstruction_failure_reason = _solve_modal_fixed_boundary_reconstruction_result['fixed_boundary_reconstruction_failure_reason']
    fixed_boundary_reconstruction_status = _solve_modal_fixed_boundary_reconstruction_result['fixed_boundary_reconstruction_status']
    left_eq63_branch = _solve_modal_fixed_boundary_reconstruction_result['left_eq63_branch']
    left_throat_temperature = _solve_modal_fixed_boundary_reconstruction_result['left_throat_temperature']
    local_reconstruction = _solve_modal_fixed_boundary_reconstruction_result['local_reconstruction']
    right_eq63_branch = _solve_modal_fixed_boundary_reconstruction_result['right_eq63_branch']
    right_throat_temperature = _solve_modal_fixed_boundary_reconstruction_result['right_throat_temperature']
    metadata_runtime_started = perf_counter()
    _build_modal_solver_metadata_result = build_modal_solver_metadata(locals())
    record_runtime(_runtime_profile_accumulator, f"modal_metadata_{runtime_species_id}", perf_counter() - metadata_runtime_started, nested=True)
    metadata = _build_modal_solver_metadata_result['metadata']
    metadata["operating_point_evaluation_mode"] = _evaluation_mode.value
    metadata.update({
        "modal_charge_exchange_target_sink_active": charge_exchange_sink_active,
        "modal_charge_exchange_target_sink_rate_s": charge_exchange_sink_rate_s,
        "modal_charge_exchange_target_sink_energy_W": charge_exchange_sink_energy_W,
        "modal_charge_exchange_target_sink_model": charge_exchange_sink_model,
        "modal_charge_exchange_target_sink_reference_identity_relative_error": charge_exchange_sink_reference_identity_error,
        "modal_charge_exchange_target_sink_pitch_resolution": "pitch_averaged_speed_resolved",
    })
    warm_start_runtime_started = perf_counter()
    eq59_final_warm_start_state = (None if eq59_diagnostics is None or eq59_collision_operator_state is None else build_eq59_warm_start_state(speed_grid=speed_grid, modal_distribution_f_j_v=modal_total, magnetic_eigenvalues=basis.physical_basis.eigenvalues, active_eigenvalues_by_speed=active_eigenvalues_by_speed, source_coefficients_by_component_j=np.asarray(component_coeffs, dtype=float), source_speeds_m_s=np.asarray(component_speeds, dtype=float), collision_operator_state=eq59_collision_operator_state, particle_mass_kg=fast_ion_species.mass_kg, charge_exchange_loss_frequency_s=charge_exchange_loss_frequency_s, species_label=fast_ion_species.symbol, nonlinear_converged=bool(eq59_diagnostics.converged), mode_eta_integrals=_modal_mode_eta_integrals(basis), physical_eigenfunctions_eta=basis.physical_basis.eigenfunctions, eta_grid=basis.physical_basis.eta_grid))
    record_runtime(_runtime_profile_accumulator, f"eq59_warm_start_packaging_{runtime_species_id}", perf_counter() - warm_start_runtime_started, nested=True)
    
    result = ModalFBISResult(
        species=fast_ion_species,
        speed_grid=speed_grid,
        lambda_grid=lambda_grid,
        pitch_grid=pitch_grid,
        basis=basis,
        component_results=tuple(component_results),
        modal_distribution_f_j_v=modal_total,
        distribution_v_eta=f_eta,
        distribution_v_lambda=f_lambda,
        density_m3=float(density),
        inventory_particles=float(inventory),
        ion_particle_loss_rate_s=float(ion_particle_loss_rate_s),
        ion_midplane_kinetic_power_loss_W=float(ion_midplane_power_loss_W),
        ion_wall_power_loss_W=float(ion_wall_power_W),
        electron_wall_power_loss_W=float(electron_power_W),
        wall_barrier_energy_J=float(barrier),
        normalized_wall_barrier=float(barrier / electron_temperature_J) if electron_temperature_J > 0.0 else np.nan,
        electron_wall_potential_relative_to_midplane_V=float(np.asarray(current_balance.electron_wall_potential_relative_to_midplane_V, dtype=float)),
        dF_dlambda_boundary_v=dfdlambda,
        boundary_flux_density_v=flux,
        confinement_time_s=confinement,
        rosenbluth_coefficients=rosen_coeffs,
        eq59_collision_operator_state=eq59_collision_operator_state,
        electrostatic_profile=electrostatic_profile,
        local_reconstruction=local_reconstruction,
        metadata=metadata,
        eq59_warm_start_state=eq59_final_warm_start_state,
        current_balance=current_balance,
    )
 
    if eq59_collision_operator_state is None:
        return result
 
    energy_moment_runtime_started = perf_counter()
    fast_ion_electron_heating_state = build_fast_ion_electron_heating_state(modal_result=result, volume_m3=volume_m3, energy_residual_relative_tolerance=numerics.rosenbluth_relative_tolerance, source_power_relative_tolerance=numerics.source_projection_relative_tolerance)
    record_runtime(_runtime_profile_accumulator, f"eq59_energy_moment_{runtime_species_id}", perf_counter() - energy_moment_runtime_started, nested=True)
    result_metadata = dict(result.metadata)
    result_metadata.update(fast_ion_electron_heating_metadata(fast_ion_electron_heating_state, prefix=f"fast_{fast_ion_species.symbol}"))
    if charge_exchange_sink_active:
        discrete_charge_exchange_sink_power_W = float(fast_ion_electron_heating_state.charge_exchange_sink_power_W)
        charge_exchange_energy_scale_W = max(abs(charge_exchange_sink_energy_W), abs(discrete_charge_exchange_sink_power_W), 1.0)
        result_metadata[f"fast_{fast_ion_species.symbol}_energy_identity_external_charge_exchange_sink_W"] = charge_exchange_sink_energy_W
        result_metadata[f"fast_{fast_ion_species.symbol}_energy_identity_discrete_charge_exchange_sink_W"] = discrete_charge_exchange_sink_power_W
        result_metadata[f"fast_{fast_ion_species.symbol}_energy_identity_charge_exchange_sink_relative_error"] = abs(discrete_charge_exchange_sink_power_W - charge_exchange_sink_energy_W) / charge_exchange_energy_scale_W
        result_metadata[f"fast_{fast_ion_species.symbol}_energy_identity_includes_external_charge_exchange_sink"] = bool(fast_ion_electron_heating_state.energy_identity_includes_charge_exchange_sink)
  
    return replace(result, metadata=result_metadata, fast_ion_electron_heating_state=fast_ion_electron_heating_state)
