"""
Single species modal density closure and Eq 42 density weighted basis iteration

This path supports the deuterium benchmark closure and couples the solved quasineutral density shape back into the orbit averaged Eq 42 scattering basis
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState
from source_model_revamp.fbis.species import DEUTERON, TRITON
from source_model_revamp.fbis.modal.basis import match_eigenmodes_on_common_eta
from source_model_revamp.fbis.modal.density_weighting import EQ42_PRESCRIBED_AXIAL_DENSITY_PROFILE, EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, Eq42BasisChange, Eq42DensityProfile, build_eq42_density_profile, eq42_density_weighting_model, compare_eq42_density_profiles, relax_eq42_density_profile
from source_model_revamp.fbis.modal import ModalFBISBasis, ModalFBISNumerics, ModalFBISResult
from source_model_revamp.fbis.modal.species_state import FastIonSpeciesRequest
from source_model_revamp.fbis.modal.solver import solve_fast_ion_species
from source_model_revamp.integration.pipeline_types import BeamEnsembleResult, GeometryStageResult
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.modal_stage.types import OperatingPointDensityState, _ModalDensityClosureResult, _ReusableModalBasisBundle
from source_model_revamp.integration.modal_stage.basis_density import _background_positive_charge_on_confined_grid, _build_modal_basis_bundle_for_eq42_profile, _collision_state_for_density, _initial_eq42_density_profile
from source_model_revamp.integration.evaluation import OperatingPointEvaluationMode

_EPS_DENSITY_M3 = 1.0e6

def _relative_scalar_change(current: float, target: float, floor: float = 1.0e-300) -> float:
    """Return the symmetric relative change between two positive scalar values"""
  
    return abs(float(target) - float(current)) / max(abs(float(target)), abs(float(current)), float(floor))

def _relative_array_change(current: np.ndarray, target: np.ndarray) -> float:
    """Return the maximum normalized absolute change between two finite arrays"""
    current_values = np.asarray(current, dtype=float)
    target_values = np.asarray(target, dtype=float)
    if current_values.shape != target_values.shape or np.any(~np.isfinite(current_values)) or np.any(~np.isfinite(target_values)):
        return float("inf")
    difference = float(np.max(np.abs(target_values - current_values)))
    scale = max(float(np.max(np.abs(target_values))), float(np.max(np.abs(current_values))), np.finfo(float).tiny)
  
    return difference / scale

def _seed_reference_species_on_confined_grid(geometry: GeometryStageResult) -> tuple[np.ndarray, np.ndarray]:
    """Return configured D and T reference density profiles on confined cell centers in m⁻³"""
    background = geometry.background_profiles
    if background is None:
        raise ValueError("seed profile evaluation requires geometry owned ion profiles")
    z_m = np.asarray(geometry.z_centers_m, dtype=float)
    deuterium = np.asarray(background.deuterium_density_m3_at(z_m), dtype=float)
    tritium = np.asarray(background.tritium_density_m3_at(z_m), dtype=float)
  
    return deuterium, tritium


def _solve_modal_with_collision_state(*, config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, electron_temperature_J: float, ion_temperature_J: float, electron_midplane_density_m3: float, electron_collision_density_m3: float, eq70_prescribed_electron_midplane_density_m3: float | None, collision_state: FBISCollisionParameterState, numerics: ModalFBISNumerics, modal_basis: _ReusableModalBasisBundle | ModalFBISBasis | None, eq59_warm_start_state: object | None = None, eq59_warm_start_compatibility_policy: str = "exact_physics", eq70_initial_potential_energy_J: np.ndarray | None = None, run_numerical_convergence_assessments: bool = False, fixed_boundary_reconstruction_policy: str = "required", evaluation_mode: OperatingPointEvaluationMode = OperatingPointEvaluationMode.FIXED_TEMPERATURE_FINAL) -> ModalFBISResult:
    """
    Run one deuterium modal FBIS solve with a fixed collision state
    
    The helper forwards the active basis, Eq 59 warm state, Eq 70 potential seed, loss model, and evaluation mode into solve_fast_ion_species
    """
    if isinstance(modal_basis, _ReusableModalBasisBundle):
        basis = modal_basis.basis
        assessment = modal_basis.convergence_assessment
        published_reference_lambda1 = modal_basis.published_reference_lambda1
    else:
        basis = modal_basis
        assessment = None
        published_reference_lambda1 = None
    background_positive = None
    background_midplane = 0.0
    background_left = 0.0
    background_right = 0.0
    feedback_model = str(getattr(config.kinetic_electrostatic, "electrostatic_feedback_model", "magnetic_only_fbis"))

    species_state = solve_fast_ion_species(
        request=FastIonSpeciesRequest(species=DEUTERON, speed_grid=beam.speed_grid_by_species[DEUTERON.species_id], attenuated_source=beam.attenuated_source_by_species[DEUTERON.species_id]),
        lambda_grid=beam.lambda_grid,
        volume_m3=geometry.volume_m3,
        midplane_area_m2=geometry.midplane_area_m2,
        half_length_m=geometry.half_length_m,
        mirror_ratio=geometry.mirror_ratio,
        B_tilde_function=geometry.B_tilde_function,
        zeta_faces=geometry.zeta_edges,
        B_tilde_midpoints=geometry.B_tilde_midpoints,
        electron_temperature_J=float(electron_temperature_J),
        electron_midplane_density_m3=float(electron_midplane_density_m3),
        electron_collision_density_m3=float(electron_collision_density_m3),
        eq70_prescribed_electron_midplane_density_m3=eq70_prescribed_electron_midplane_density_m3,
        ion_temperature_J=float(ion_temperature_J),
        cell_volumes_m3=geometry.cell_volumes_m3,
        pitch_grid=beam.pitch_grid,
        collision_state=collision_state,
        numerics=numerics,
        auxiliary_input_power_W=0.0,
        global_loss_rate_model=config.kinetic_electrostatic.ion_loss_closure_model,
        electrostatic_feedback_model=feedback_model,
        ion_loss_closure_model=config.kinetic_electrostatic.ion_loss_closure_model,
        lost_ion_parallel_temperature_model=(config.kinetic_electrostatic.lost_ion_parallel_temperature_model),
        background_positive_charge_density_m3=background_positive,
        background_midplane_positive_charge_density_m3=background_midplane,
        background_left_throat_positive_charge_density_m3=background_left,
        background_right_throat_positive_charge_density_m3=background_right,
        basis=basis,
        basis_convergence_assessment=assessment,
        _eq59_warm_start_state=eq59_warm_start_state,
        _eq59_warm_start_compatibility_policy=eq59_warm_start_compatibility_policy,
        _eq70_initial_potential_energy_J=eq70_initial_potential_energy_J,
        _run_numerical_convergence_assessments=run_numerical_convergence_assessments,
        _published_reference_lambda1=published_reference_lambda1,
        _fixed_boundary_reconstruction_policy=(fixed_boundary_reconstruction_policy),
        _evaluation_mode=evaluation_mode,
    )
    if species_state.modal_result is None:
        raise RuntimeError("the active deuterium source did not produce a modal result")

    return species_state.modal_result


def _egedal_beam_plasma_density_closure(*, config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, electron_temperature_J: float, ion_temperature_J: float, numerics: ModalFBISNumerics, modal_basis: _ReusableModalBasisBundle | ModalFBISBasis | None, initial_eq59_warm_start_state: object | None = None, initial_eq70_potential_energy_J: np.ndarray | None = None, fixed_boundary_reconstruction_policy: str = "required", evaluation_mode: OperatingPointEvaluationMode = OperatingPointEvaluationMode.FIXED_TEMPERATURE_FINAL) -> _ModalDensityClosureResult:
    """
    Iterate the deuterium benchmark density closure at the prescribed electron density scale
    
    Each iteration rebuilds the collision state, solves the modal system, compares the solved density and Eq 59 collision operator to the previous iterate, and relaxes the scalar density when needed
    """
    k = config.kinetic_electrostatic
    benchmark_density = config.plasma_closure.benchmark_prescribed_electron_density_m3
    if benchmark_density is None:
        raise ValueError("benchmark prescribed electron density is unavailable")
    density = max(float(benchmark_density), _EPS_DENSITY_M3)
    iterations = int(k.modal_density_closure_iterations)
    tolerance = float(k.modal_density_closure_relative_tolerance)
    relaxation = float(k.modal_density_closure_relaxation)
    if iterations < 1:
        raise ValueError("modal_density_closure_iterations must be at least one")
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("modal_density_closure_relative_tolerance must be nonnegative and finite")
    if not np.isfinite(relaxation) or relaxation <= 0.0 or relaxation > 1.0:
        raise ValueError("modal_density_closure_relaxation must be in the interval (0, 1]")
    ion_charges = np.asarray([1.0], dtype=float)
    ion_masses = np.asarray([DEUTERON.mass_kg], dtype=float)
    history: list[dict[str, float]] = []
    modal: ModalFBISResult | None = None
    collision_state: FBISCollisionParameterState | None = None
    relative_error = float("nan")
    converged = False
    eq59_warm_start_state: object | None = initial_eq59_warm_start_state
    eq70_potential_energy_J = initial_eq70_potential_energy_J
    previous_operator_drag: np.ndarray | None = None
    previous_operator_diffusion: np.ndarray | None = None
    previous_operator_pitch: np.ndarray | None = None
    for iteration in range(1, iterations + 1):
        ion_densities = np.asarray([density], dtype=float)
        collision_state = _collision_state_for_density(fast_ion_species=DEUTERON, electron_density_m3=density, ion_densities_m3=ion_densities, ion_charge_numbers=ion_charges, ion_masses_kg=ion_masses, electron_temperature_J=electron_temperature_J, ion_temperature_J=ion_temperature_J, fast_ion_density_m3=density, collision_scope="pure_fast_D_benchmark")
        modal = _solve_modal_with_collision_state(config=config, geometry=geometry, beam=beam, electron_temperature_J=electron_temperature_J, ion_temperature_J=ion_temperature_J, electron_midplane_density_m3=density, electron_collision_density_m3=density, eq70_prescribed_electron_midplane_density_m3=density, collision_state=collision_state, numerics=numerics, modal_basis=modal_basis, eq59_warm_start_state=eq59_warm_start_state, eq59_warm_start_compatibility_policy=("exact_physics" if eq59_warm_start_state is None else "closure_continuation"), eq70_initial_potential_energy_J=eq70_potential_energy_J, fixed_boundary_reconstruction_policy=(fixed_boundary_reconstruction_policy), evaluation_mode=evaluation_mode)
        eq59_warm_start_state = modal.eq59_warm_start_state
        eq70_potential_energy_J = (None if modal.electrostatic_profile is None else np.asarray(modal.electrostatic_profile.potential_energy_J, dtype=float))
        solved_density = float(modal.density_m3)
        scale = max(abs(density), abs(solved_density), _EPS_DENSITY_M3)
        density_relative_error = (solved_density - density) / scale
        operator_state = modal.eq59_collision_operator_state
        if operator_state is None:
            operator_drag_change = 0.0
            operator_diffusion_change = 0.0
            operator_pitch_change = 0.0
            operator_model = "legacy_scalar_cold_Eq14_reference"
        else:
            operator_drag = np.asarray(operator_state.ion_drag_velocity_cubed_m3_s3, dtype=float)
            operator_diffusion = np.asarray(operator_state.ion_energy_diffusion_velocity_fourth_m4_s4, dtype=float)
            operator_pitch = np.asarray(operator_state.ion_pitch_scattering_velocity_cubed_m3_s3, dtype=float)
            operator_drag_change = float("inf") if previous_operator_drag is None else _relative_array_change(previous_operator_drag, operator_drag)
            operator_diffusion_change = float("inf") if previous_operator_diffusion is None else _relative_array_change(previous_operator_diffusion, operator_diffusion)
            operator_pitch_change = float("inf") if previous_operator_pitch is None else _relative_array_change(previous_operator_pitch, operator_pitch)
            operator_model = operator_state.operator_model
        relative_error = max(abs(float(density_relative_error)), operator_drag_change, operator_diffusion_change, operator_pitch_change)
        history.append({"iteration": float(iteration), "collision_electron_density_m3": float(density), "modal_density_m3": float(solved_density), "density_relative_error": float(density_relative_error), "relative_error": float(relative_error), "beta_m_diagnostic": float(collision_state.beta_m), "spitzer_slowing_down_time_s": float(collision_state.spitzer_slowing_down_time_s), "critical_velocity_m_s_diagnostic": float(collision_state.critical_velocity_m_s), "eq59_collision_operator_model": operator_model, "eq59_ion_drag_operator_relative_change": operator_drag_change, "eq59_ion_energy_diffusion_operator_relative_change": operator_diffusion_change, "eq59_ion_pitch_scattering_operator_relative_change": operator_pitch_change, "modal_ion_particle_loss_rate_s": float(modal.ion_particle_loss_rate_s), "modal_wall_barrier_over_Te": float(modal.metadata.get("modal_wall_barrier_over_Te", float("nan")))})
        if operator_state is not None:
            previous_operator_drag = np.asarray(operator_state.ion_drag_velocity_cubed_m3_s3, dtype=float).copy()
            previous_operator_diffusion = np.asarray(operator_state.ion_energy_diffusion_velocity_fourth_m4_s4, dtype=float).copy()
            previous_operator_pitch = np.asarray(operator_state.ion_pitch_scattering_velocity_cubed_m3_s3, dtype=float).copy()
        if relative_error <= tolerance:
            converged = True
            break
        if not np.isfinite(solved_density) or solved_density <= 0.0:
            raise ValueError("modal density closure requires a positive solved modal density")
        target_density = float(np.sqrt(max(density, _EPS_DENSITY_M3) * solved_density))
        log_next = (1.0 - relaxation) * np.log(density) + relaxation * np.log(target_density)
        density = max(float(np.exp(log_next)), _EPS_DENSITY_M3)
    if modal is None or collision_state is None:
        raise RuntimeError("modal density closure did not evaluate any states")
    final_ion_densities = np.asarray([float(collision_state.coulomb_logs.electron_density_m3)], dtype=float)

    return _ModalDensityClosureResult(
        modal=modal,
        collision_state=collision_state,
        electron_midplane_density_m3=float(collision_state.coulomb_logs.electron_density_m3),
        electron_collision_density_m3=float(collision_state.coulomb_logs.electron_density_m3),
        ion_densities_m3=final_ion_densities,
        ion_charge_numbers=ion_charges,
        ion_masses_kg=ion_masses,
        converged=bool(converged),
        iterations=len(history),
        relative_error=float(relative_error),
        history=history,
    )

def _solve_scalar_modal_density_closure(*, config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, electron_temperature_J: float, ion_temperature_J: float, numerics: ModalFBISNumerics, modal_basis: _ReusableModalBasisBundle | ModalFBISBasis | None, initial_eq59_warm_start_state: object | None = None, initial_representative_fast_density_m3: float | None = None, initial_fast_D_midpoint_density_m3: float | None = None, initial_fast_D_confined_average_density_m3: float | None = None, initial_electron_midpoint_density_m3: float | None = None, initial_electron_confined_average_density_m3: float | None = None, initial_electron_parent_n0_m3: float | None = None, initial_electron_collision_density_m3: float | None = None, initial_eq70_potential_energy_J: np.ndarray | None = None, fixed_boundary_reconstruction_policy: str = "required", evaluation_mode: OperatingPointEvaluationMode = OperatingPointEvaluationMode.FIXED_TEMPERATURE_FINAL) -> _ModalDensityClosureResult:
    """Dispatch the configured single species scalar density closure"""
    model = str(config.kinetic_electrostatic.modal_density_closure_model).strip().lower()
    pass
    if model == "egedal_beam_plasma_quasineutral":
        return _egedal_beam_plasma_density_closure(
            config=config,
            geometry=geometry,
            beam=beam,
            electron_temperature_J=electron_temperature_J,
            ion_temperature_J=ion_temperature_J,
            numerics=numerics,
            modal_basis=modal_basis,
            initial_eq59_warm_start_state=initial_eq59_warm_start_state,
            initial_eq70_potential_energy_J=initial_eq70_potential_energy_J,
            fixed_boundary_reconstruction_policy=(fixed_boundary_reconstruction_policy),
            evaluation_mode=evaluation_mode,
        )
    raise ValueError(f"Unsupported modal_density_closure_model {config.kinetic_electrostatic.modal_density_closure_model!r}")

def _eq42_profile_from_modal_state(*, config: SourceModelRunConfig, geometry: GeometryStageResult, modal: ModalFBISResult, numerics: ModalFBISNumerics, density_state: OperatingPointDensityState | None = None) -> tuple[Eq42DensityProfile, bool, str]:
    """
    Build the next Eq 42 profile from the represented quasineutral state
    
    When a typed operating point density state is available its total positive charge profile is authoritative
    Otherwise the converged Eq 70 ion profile supplies the density shape
    """
    zeta_cells = np.asarray(geometry.zeta_centers, dtype=float)
    cell_volumes = np.asarray(geometry.cell_volumes_m3, dtype=float)
    electrostatic = modal.electrostatic_profile
    if density_state is not None:
        density_cells = np.asarray(density_state.total_positive_charge_cell_density_m3, dtype=float)
        density_nodes = np.asarray(density_state.total_positive_charge_node_density_m3, dtype=float)
        zeta_nodes = np.asarray(density_state.zeta_nodes, dtype=float)
        if density_cells.shape != zeta_cells.shape or density_nodes.shape != zeta_nodes.shape:
            raise ValueError("operating point density state does not match the Eq 42 grids")
        profile = build_eq42_density_profile(model=EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, source="operating_point_density_state_total_positive_charge", zeta_cells=zeta_cells, density_cells_m3=density_cells, cell_volumes_m3=cell_volumes, zeta_nodes=zeta_nodes, density_nodes_m3=density_nodes, symmetry_tolerance=numerics.eq42_density_symmetry_tolerance)
        return profile, bool(electrostatic is not None and electrostatic.converged), "typed_operating_point_total_positive_charge_profile"
    if electrostatic is not None:
        source_zeta = np.asarray(electrostatic.zeta, dtype=float)
        source_density = np.asarray(electrostatic.ion_density_m3, dtype=float)
        density_cells = np.interp(zeta_cells, source_zeta, source_density)
        node_zeta = electrostatic.electrostatic_node_zeta
        node_density = electrostatic.electrostatic_node_ion_density_m3
        if node_zeta is not None and node_density is not None:
            zeta_nodes = np.asarray(node_zeta, dtype=float)
            density_nodes = np.asarray(node_density, dtype=float)
        else:
            zeta_nodes = np.asarray(geometry.zeta_edges, dtype=float)
            density_nodes = np.interp(zeta_nodes, source_zeta, source_density, left=float(source_density[0]), right=float(source_density[-1]))
        profile = build_eq42_density_profile(model=EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, source=("Eq70_ion_side_quasineutral_target_density_profile"), zeta_cells=zeta_cells, density_cells_m3=density_cells, cell_volumes_m3=cell_volumes, zeta_nodes=zeta_nodes, density_nodes_m3=density_nodes, symmetry_tolerance=numerics.eq42_density_symmetry_tolerance)
        return profile, bool(electrostatic.converged), "Eq70_ion_side_quasineutral_target_profile"
    if modal.local_reconstruction is None:
        raise ValueError("self consistent Eq 42 weighting requires an Eq 70 profile or local ion density")
    local_fast = np.asarray(modal.local_reconstruction.local_density_m3, dtype=float)
    if local_fast.shape != zeta_cells.shape:
        raise ValueError("local fast ion density must match the confined Eq 42 grid")
    background_positive, _, _, background_source = _background_positive_charge_on_confined_grid(geometry)
    if background_positive is None or config.plasma_closure.model == "nbi_supported_stationary":
        density_cells = local_fast
    else:
        density_cells = local_fast + np.asarray(background_positive, dtype=float)
    zeta_nodes = np.asarray(geometry.zeta_edges, dtype=float)
    density_nodes = np.interp(zeta_nodes, zeta_cells, density_cells, left=float(density_cells[0]), right=float(density_cells[-1]))
    feedback_model = str(config.kinetic_electrostatic.electrostatic_feedback_model).strip().lower()
    reference_profile_valid = feedback_model in {"magnetic_only_fbis", "zero_potential_magnetic_reference",}
    profile = build_eq42_density_profile(model=EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, source=f"local_positive_charge_density_without_Eq70:{background_source}", zeta_cells=zeta_cells, density_cells_m3=density_cells, cell_volumes_m3=cell_volumes, zeta_nodes=zeta_nodes, density_nodes_m3=density_nodes, symmetry_tolerance=numerics.eq42_density_symmetry_tolerance)
    
    return profile, reference_profile_valid, "local_positive_charge_profile"

def _compare_eq42_bases(current: ModalFBISBasis, target: ModalFBISBasis,) -> Eq42BasisChange:
    """Compare two density dependent physical eigenbases by eigenvalue change and matched eigenfunction overlap on their common η interval"""
    match = match_eigenmodes_on_common_eta(current, target)
    current_eigenvalues = np.asarray(current.physical_basis.eigenvalues, dtype=float)
    target_eigenvalues = np.asarray(target.physical_basis.eigenvalues, dtype=float)[match.fine_mode_indices]
    count = min(current_eigenvalues.size, target_eigenvalues.size)
    relative = np.abs(target_eigenvalues[:count] - current_eigenvalues[:count]) / np.maximum(np.maximum(np.abs(target_eigenvalues[:count]), np.abs(current_eigenvalues[:count])), np.finfo(float).tiny,)
    overlaps = np.asarray(match.weighted_overlaps[:count], dtype=float)
    
    return Eq42BasisChange(
        maximum_eigenvalue_relative_change=float(np.max(relative)) if relative.size else 0.0,
        minimum_eigenfunction_overlap=float(np.min(overlaps)) if overlaps.size else 1.0,
        eigenvalue_relative_changes=np.asarray(relative, dtype=float),
        eigenfunction_overlaps=overlaps,
    )

def _eq42_state_relative_change(previous_distribution: np.ndarray | None, current_distribution: np.ndarray) -> float | None:
    """Return the relative change in the physical invariant distribution `F(v, Λ)` between two Eq 42 iterations"""
    if previous_distribution is None:
        return None
    previous = np.asarray(previous_distribution, dtype=float)
    current = np.asarray(current_distribution, dtype=float)
    if previous.shape != current.shape:
        return None
    numerator = float(np.linalg.norm(current - previous))
    denominator = max(float(np.linalg.norm(current)), np.finfo(float).tiny)

    return numerator / denominator

def _eq42_metadata(*, model: str, profile: Eq42DensityProfile, numerics: ModalFBISNumerics, converged: bool, iterations: int, history: list[dict[str, object]], quasineutral_profile_converged: bool, final_profile_change: float | None, final_basis_change: Eq42BasisChange | None, outer_iteration_status: str, outer_iteration_failure_reason: str | None, last_valid_iteration: int, candidate_rejection_count: int) -> dict[str, object]:
    """Return Eq 42 density profile identity, symmetry, basis change, iteration history, and convergence metadata"""
    outer_applicable = model == EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY
    return {
        "modal_eq42_density_weighting_model": model,
        "modal_eq42_density_weighting_coupled": True,
        "modal_eq42_nonuniform_density_weighting_active": profile.nonuniform_weighting_active,
        "modal_eq42_nonuniform_density_weighting_required": True,
        "modal_eq42_density_profile_source": profile.source,
        "modal_eq42_density_profile_volume_average_density_m3": profile.volume_average_density_m3,
        "modal_eq42_density_profile_raw_volume_average_density_m3": profile.raw_volume_average_density_m3,
        "modal_eq42_density_ratio_zeta": [float(value) for value in profile.zeta_positive_nodes],
        "modal_eq42_density_ratio_values": [float(value) for value in profile.density_ratio_positive_nodes],
        "modal_eq42_density_profile_identity": profile.profile_identity,
        "modal_eq42_density_profile_symmetry_relative_error": profile.symmetry_relative_error,
        "modal_eq42_density_profile_symmetry_tolerance": profile.symmetry_tolerance,
        "modal_eq42_density_profile_symmetric": profile.symmetric_within_tolerance,
        "modal_eq42_density_normalization_measure": "physical_flux_tube_cell_volume_average_n_equals_N_over_V",
        "modal_eq42_eta_mapping_density_weighted": False,
        "modal_eq42_eta_mapping_model": "Egedal_2022_equations_47_48_bounce_time_phase_space_measure",
        "modal_eq42_basis_cache_identity_includes_normalized_density_shape": True,
        "modal_eq42_density_amplitude_changes_basis": False,
        "modal_eq42_outer_iteration_applicable": outer_applicable,
        "modal_eq42_outer_iteration_converged": bool(converged),
        "modal_eq42_outer_iteration_iterations": int(iterations),
        "modal_eq42_outer_iteration_history": history,
        "modal_eq42_outer_iteration_status": str(outer_iteration_status),
        "modal_eq42_outer_iteration_failure_reason": outer_iteration_failure_reason,
        "modal_eq42_outer_iteration_last_valid_iteration": int(last_valid_iteration),
        "modal_eq42_outer_iteration_candidate_rejection_count": int(candidate_rejection_count),
        "modal_eq42_quasineutral_profile_source_converged": bool(quasineutral_profile_converged),
        "modal_eq42_density_shape_final_relative_change": final_profile_change,
        "modal_eq42_density_shape_relative_tolerance": (numerics.eq42_density_shape_relative_tolerance),
        "modal_eq42_density_shape_relaxation": (numerics.eq42_density_shape_relaxation),
        "modal_eq42_basis_eigenvalue_final_relative_change": (None if final_basis_change is None else final_basis_change.maximum_eigenvalue_relative_change),
        "modal_eq42_basis_eigenfunction_final_minimum_overlap": (None if final_basis_change is None else final_basis_change.minimum_eigenfunction_overlap),
        "modal_eq42_basis_eigenvalue_relative_tolerance": (numerics.eq42_basis_eigenvalue_relative_tolerance),
        "modal_eq42_basis_eigenfunction_overlap_tolerance": (numerics.eq42_basis_eigenfunction_overlap_tolerance),
        "modal_eq42_reference": "Egedal_2022_equation_42",
    }

def _solve_modal_density_closure(*, config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, electron_temperature_J: float, ion_temperature_J: float, numerics: ModalFBISNumerics, modal_basis: _ReusableModalBasisBundle | ModalFBISBasis | None, initial_eq59_warm_start_state: object | None = None, initial_representative_fast_density_m3: float | None = None, initial_fast_D_midpoint_density_m3: float | None = None, initial_fast_D_confined_average_density_m3: float | None = None, initial_electron_midpoint_density_m3: float | None = None, initial_electron_confined_average_density_m3: float | None = None, initial_electron_parent_n0_m3: float | None = None, initial_electron_collision_density_m3: float | None = None, initial_eq42_density_profile: Eq42DensityProfile | None = None, initial_eq70_potential_energy_J: np.ndarray | None = None, evaluation_mode: OperatingPointEvaluationMode = OperatingPointEvaluationMode.FIXED_TEMPERATURE_FINAL) -> _ModalDensityClosureResult:
    """
    Couple the scalar density closure to the Eq 42 density weighted modal basis
    
    For self consistent Eq 42 weighting the loop alternates scalar kinetic solves, quasineutral density profile reconstruction, density relaxation, and physical basis rebuilding before the final basis assessment
    """
    model = eq42_density_weighting_model(config.kinetic_electrostatic.modal_eq42_density_weighting_model)
    if isinstance(modal_basis, _ReusableModalBasisBundle):
        initial_bundle = modal_basis
        if initial_eq42_density_profile is not None:
            bundle_profile = initial_bundle.eq42_density_profile
            if bundle_profile is None or (bundle_profile.profile_identity != initial_eq42_density_profile.profile_identity):
                raise ValueError("Eq42 warm profile is incompatible with the warm basis")
    else:
        profile = (_initial_eq42_density_profile(config, geometry, numerics) if initial_eq42_density_profile is None else initial_eq42_density_profile)
        if eq42_density_weighting_model(profile.model) != model:
            raise ValueError("Eq42 warm profile model is incompatible")
        if modal_basis is None:
            initial_bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=profile, assess_grid_convergence=(model != EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY),)
        else:
            initial_bundle = _ReusableModalBasisBundle(basis=modal_basis, convergence_assessment=None, published_reference_lambda1=None, eq42_density_profile=profile)
    if initial_bundle.eq42_density_profile is None:
        raise RuntimeError("Eq 42 basis bundle is missing its density profile")
    if model != EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY:
        bundle = initial_bundle
        if (evaluation_mode.runs_final_qualification and bundle.convergence_assessment is None and numerics.basis_convergence_levels >= 3):
            bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=bundle.eq42_density_profile, assess_grid_convergence=True, published_reference_lambda1=(bundle.published_reference_lambda1), basis_cache=bundle.basis_cache, bypass_cache=True,)
        scalar = _solve_scalar_modal_density_closure(config=config, geometry=geometry, beam=beam, electron_temperature_J=electron_temperature_J, ion_temperature_J=ion_temperature_J, numerics=numerics, modal_basis=bundle, initial_eq59_warm_start_state=initial_eq59_warm_start_state, initial_representative_fast_density_m3=initial_representative_fast_density_m3, initial_fast_D_midpoint_density_m3=initial_fast_D_midpoint_density_m3, initial_fast_D_confined_average_density_m3=initial_fast_D_confined_average_density_m3, initial_electron_midpoint_density_m3=initial_electron_midpoint_density_m3, initial_electron_confined_average_density_m3=initial_electron_confined_average_density_m3, initial_electron_parent_n0_m3=initial_electron_parent_n0_m3, initial_electron_collision_density_m3=initial_electron_collision_density_m3, initial_eq70_potential_energy_J=initial_eq70_potential_energy_J, fixed_boundary_reconstruction_policy="required", evaluation_mode=evaluation_mode)
        profile = bundle.eq42_density_profile
        if profile is None:
            raise RuntimeError("Eq 42 profile disappeared from the basis bundle")
        eq42_converged = bool(profile.symmetric_within_tolerance)
        history = [{"iteration": 1, "profile_model": model, "profile_source": profile.source, "profile_symmetry_relative_error": profile.symmetry_relative_error, "profile_symmetric": profile.symmetric_within_tolerance, "nonuniform_weighting_active": profile.nonuniform_weighting_active, "status": ("prescribed_profile" if model == EQ42_PRESCRIBED_AXIAL_DENSITY_PROFILE else "uniform_reference")}]
        modal = replace(scalar.modal,metadata={**scalar.modal.metadata, **_eq42_metadata(
                    model=model,
                    profile=profile,
                    numerics=numerics,
                    converged=eq42_converged,
                    iterations=1,
                    history=history,
                    quasineutral_profile_converged=True,
                    final_profile_change=0.0,
                    final_basis_change=None,
                    outer_iteration_status=("not_applicable_prescribed_profile" if model == EQ42_PRESCRIBED_AXIAL_DENSITY_PROFILE else "not_applicable_uniform_reference"),
                    outer_iteration_failure_reason=None,
                    last_valid_iteration=1,
                    candidate_rejection_count=0,
                ),},)
        
        return replace(scalar, modal=modal, modal_basis=bundle, eq42_density_profile=profile, eq42_density_basis_converged=eq42_converged, eq42_density_basis_iterations=1, eq42_density_basis_history=history,)

    current_profile = initial_bundle.eq42_density_profile
    current_bundle = initial_bundle
    if current_bundle.convergence_assessment is not None:
        current_bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=current_profile, assess_grid_convergence=False, published_reference_lambda1=(current_bundle.published_reference_lambda1), basis_cache=current_bundle.basis_cache,)
    history: list[dict[str, object]] = []
    previous_distribution: np.ndarray | None = None
    final_scalar: _ModalDensityClosureResult | None = None
    final_target: Eq42DensityProfile | None = None
    final_basis_change: Eq42BasisChange | None = None
    final_profile_change: float | None = None
    quasineutral_profile_converged = False
    shape_basis_converged = False
    eq59_warm_start_state = initial_eq59_warm_start_state
    eq70_potential_energy_J: np.ndarray | None = initial_eq70_potential_energy_J
    representative_fast_density_seed_m3: float | None = initial_representative_fast_density_m3
    fast_midpoint_seed_m3 = initial_fast_D_midpoint_density_m3
    fast_average_seed_m3 = initial_fast_D_confined_average_density_m3
    electron_midpoint_seed_m3 = initial_electron_midpoint_density_m3
    electron_average_seed_m3 = initial_electron_confined_average_density_m3
    electron_parent_n0_seed_m3 = initial_electron_parent_n0_m3
    electron_collision_seed_m3 = initial_electron_collision_density_m3

    # Alternate the scalar kinetic solve with an Eq 42 density shape and basis update
    for iteration in range(1, numerics.eq42_density_basis_iterations + 1):
        scalar = _solve_scalar_modal_density_closure(
            config=config,
            geometry=geometry,
            beam=beam,
            electron_temperature_J=electron_temperature_J,
            ion_temperature_J=ion_temperature_J,
            numerics=numerics,
            modal_basis=current_bundle,
            initial_eq59_warm_start_state=eq59_warm_start_state,
            initial_representative_fast_density_m3=(representative_fast_density_seed_m3),
            initial_fast_D_midpoint_density_m3=fast_midpoint_seed_m3,
            initial_fast_D_confined_average_density_m3=fast_average_seed_m3,
            initial_electron_midpoint_density_m3=electron_midpoint_seed_m3,
            initial_electron_confined_average_density_m3=electron_average_seed_m3,
            initial_electron_parent_n0_m3=electron_parent_n0_seed_m3,
            initial_electron_collision_density_m3=electron_collision_seed_m3,
            initial_eq70_potential_energy_J=eq70_potential_energy_J,
            fixed_boundary_reconstruction_policy=("deferred_eq42_iteration"),
            evaluation_mode=evaluation_mode,
        )
        eq59_warm_start_state = scalar.modal.eq59_warm_start_state
        if scalar.modal.electrostatic_profile is not None:
            eq70_potential_energy_J = np.asarray(scalar.modal.electrostatic_profile.potential_energy_J, dtype=float,)
            profile_zeta = np.asarray(scalar.modal.electrostatic_profile.zeta, dtype=float,)
            profile_fast_density = np.asarray(scalar.modal.electrostatic_profile.fast_ion_density_m3,dtype=float,)
            representative_fast_density_seed_m3 = float(profile_fast_density[int(np.argmin(np.abs(profile_zeta)))])
        else:
            representative_fast_density_seed_m3 = float(scalar.modal.density_m3)
        if scalar.density_state is not None:
            fast_midpoint_seed_m3 = scalar.density_state.fast_deuterium_midplane_density_m3
            fast_average_seed_m3 = scalar.density_state.fast_deuterium_confined_volume_average_density_m3
            electron_midpoint_seed_m3 = scalar.density_state.electron_midplane_density_m3
            electron_average_seed_m3 = scalar.density_state.electron_confined_volume_average_density_m3
            electron_parent_n0_seed_m3 = scalar.density_state.electron_parent_maxwellian_n0_m3
            electron_collision_seed_m3 = scalar.density_state.electron_collision_density_m3
        target_profile, profile_source_converged, profile_source_status = (_eq42_profile_from_modal_state(config=config, geometry=geometry, modal=scalar.modal, numerics=numerics, density_state=scalar.density_state,))
        target_bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=target_profile, assess_grid_convergence=False, published_reference_lambda1=(current_bundle.published_reference_lambda1), basis_cache=current_bundle.basis_cache,)
        profile_change = compare_eq42_density_profiles(current_profile, target_profile)
        basis_change = _compare_eq42_bases(current_bundle.basis, target_bundle.basis,)
        state_change = _eq42_state_relative_change(previous_distribution, scalar.modal.distribution_v_lambda)
        shape_basis_converged = bool(profile_change.volume_weighted_l2_relative_change <= numerics.eq42_density_shape_relative_tolerance and basis_change.maximum_eigenvalue_relative_change <= numerics.eq42_basis_eigenvalue_relative_tolerance and basis_change.minimum_eigenfunction_overlap >= numerics.eq42_basis_eigenfunction_overlap_tolerance and target_profile.symmetric_within_tolerance)
        iteration_converged = bool(shape_basis_converged and profile_source_converged and scalar.converged)
        history.append({
                "iteration": int(iteration),
                "status": "evaluated",
                "fixed_boundary_reconstruction_status": (scalar.modal.metadata.get("fixed_boundary_lost_reconstruction_status")),
                "profile_source_status": profile_source_status,
                "quasineutral_profile_source_converged": bool(profile_source_converged),
                "density_shape_volume_weighted_l2_relative_change": (profile_change.volume_weighted_l2_relative_change),
                "density_shape_maximum_absolute_ratio_change": (profile_change.maximum_absolute_ratio_change),
                "density_shape_maximum_supported_relative_change": (profile_change.maximum_supported_relative_change),
                "basis_maximum_eigenvalue_relative_change": (basis_change.maximum_eigenvalue_relative_change),
                "basis_minimum_eigenfunction_overlap": (basis_change.minimum_eigenfunction_overlap),
                "physical_distribution_relative_change": state_change,
                "profile_symmetry_relative_error": (target_profile.symmetry_relative_error),
                "profile_symmetric": target_profile.symmetric_within_tolerance,
                "scalar_density_closure_converged": bool(scalar.converged),
                "shape_and_basis_converged": shape_basis_converged,
                "outer_iteration_converged": iteration_converged,
                "current_first_eigenvalue": float(current_bundle.basis.physical_basis.eigenvalues[0]),
                "target_first_eigenvalue": float(target_bundle.basis.physical_basis.eigenvalues[0]),
            })
        final_scalar = scalar
        final_target = target_profile
        final_basis_change = basis_change
        final_profile_change = (profile_change.volume_weighted_l2_relative_change)
        quasineutral_profile_converged = bool(profile_source_converged)
        if iteration_converged:
            break
        if iteration == numerics.eq42_density_basis_iterations:
            break
        previous_distribution = np.asarray(scalar.modal.distribution_v_lambda,dtype=float,)
        current_profile = relax_eq42_density_profile(current=current_profile, target=target_profile, relaxation=numerics.eq42_density_shape_relaxation, model=EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, source=f"relaxed_Eq42_quasineutral_profile_iteration_{iteration}",)
        current_bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=current_profile, assess_grid_convergence=False, published_reference_lambda1=(current_bundle.published_reference_lambda1), basis_cache=current_bundle.basis_cache)
    if final_scalar is None or final_target is None:
        raise RuntimeError("self consistent Eq 42 iteration did not evaluate a modal state")
    assessment_profile = final_target
    assessed_bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=assessment_profile, assess_grid_convergence=evaluation_mode.runs_final_qualification, published_reference_lambda1=(current_bundle.published_reference_lambda1), basis_cache=current_bundle.basis_cache, bypass_cache=evaluation_mode.runs_final_qualification,)
    final_scalar = _solve_scalar_modal_density_closure(config=config, geometry=geometry, beam=beam, electron_temperature_J=electron_temperature_J, ion_temperature_J=ion_temperature_J, numerics=numerics, modal_basis=assessed_bundle, initial_eq59_warm_start_state=eq59_warm_start_state, initial_representative_fast_density_m3=(representative_fast_density_seed_m3), initial_fast_D_midpoint_density_m3=fast_midpoint_seed_m3, initial_fast_D_confined_average_density_m3=fast_average_seed_m3, initial_electron_midpoint_density_m3=electron_midpoint_seed_m3, initial_electron_confined_average_density_m3=electron_average_seed_m3, initial_electron_parent_n0_m3=electron_parent_n0_seed_m3, initial_electron_collision_density_m3=electron_collision_seed_m3, initial_eq70_potential_energy_J=eq70_potential_energy_J, fixed_boundary_reconstruction_policy=("best_effort_diagnostic" if evaluation_mode.runs_final_qualification else "required"), evaluation_mode=evaluation_mode,)
    final_target, final_source_converged, final_source_status = (_eq42_profile_from_modal_state(config=config, geometry=geometry, modal=final_scalar.modal, numerics=numerics, density_state=final_scalar.density_state))
    final_profile_assessment = compare_eq42_density_profiles(assessment_profile, final_target,)
    final_target_bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=final_target, assess_grid_convergence=False, published_reference_lambda1=(assessed_bundle.published_reference_lambda1), basis_cache=assessed_bundle.basis_cache,)
    final_basis_change = _compare_eq42_bases(assessed_bundle.basis, final_target_bundle.basis,)
    final_profile_change = (final_profile_assessment.volume_weighted_l2_relative_change)
    shape_basis_converged = bool(final_profile_change <= numerics.eq42_density_shape_relative_tolerance and final_basis_change.maximum_eigenvalue_relative_change <= numerics.eq42_basis_eigenvalue_relative_tolerance and final_basis_change.minimum_eigenfunction_overlap >= numerics.eq42_basis_eigenfunction_overlap_tolerance and final_target.symmetric_within_tolerance)
    quasineutral_profile_converged = bool(final_source_converged)
    eq42_converged = bool(shape_basis_converged and quasineutral_profile_converged and final_scalar.converged)
    if eq42_converged:
        outer_iteration_status = "converged"
        outer_iteration_failure_reason = None
    elif not quasineutral_profile_converged:
        outer_iteration_status = "not_converged_quasineutral_profile"
        outer_iteration_failure_reason = "Eq70_quasineutral_profile_not_converged"
    elif not shape_basis_converged:
        outer_iteration_status = "not_converged_density_shape_or_basis"
        outer_iteration_failure_reason = ("Eq42_density_shape_or_basis_change_above_tolerance")
    else:
        outer_iteration_status = "not_converged_scalar_density_closure"
        outer_iteration_failure_reason = "scalar_density_closure_not_converged"
    history.append({
            "iteration": "final_grid_assessment",
            "status": outer_iteration_status,
            "fixed_boundary_reconstruction_status": (final_scalar.modal.metadata.get("fixed_boundary_lost_reconstruction_status")),
            "fixed_boundary_reconstruction_failure_reason": (final_scalar.modal.metadata.get("fixed_boundary_lost_reconstruction_failure_reason")),
            "profile_source_status": final_source_status,
            "quasineutral_profile_source_converged": (quasineutral_profile_converged),
            "density_shape_volume_weighted_l2_relative_change": (final_profile_change),
            "density_shape_maximum_absolute_ratio_change": (final_profile_assessment.maximum_absolute_ratio_change),
            "density_shape_maximum_supported_relative_change": (final_profile_assessment.maximum_supported_relative_change),
            "basis_maximum_eigenvalue_relative_change": (final_basis_change.maximum_eigenvalue_relative_change),
            "basis_minimum_eigenfunction_overlap": (final_basis_change.minimum_eigenfunction_overlap),
            "profile_symmetry_relative_error": (final_target.symmetry_relative_error),
            "profile_symmetric": final_target.symmetric_within_tolerance,
            "scalar_density_closure_converged": bool(final_scalar.converged),
            "shape_and_basis_converged": shape_basis_converged,
            "overall_eq42_outer_iteration_converged": eq42_converged,
            "failure_reason": outer_iteration_failure_reason,
        })
    modal = replace(final_scalar.modal, metadata={**final_scalar.modal.metadata, **_eq42_metadata(
                model=model,
                profile=assessment_profile,
                numerics=numerics,
                converged=eq42_converged,
                iterations=len(history),
                history=history,
                quasineutral_profile_converged=(quasineutral_profile_converged),
                final_profile_change=final_profile_change,
                final_basis_change=final_basis_change,
                outer_iteration_status=outer_iteration_status,
                outer_iteration_failure_reason=(outer_iteration_failure_reason),
                last_valid_iteration=int(min(numerics.eq42_density_basis_iterations, len([ item for item in history if isinstance(item.get("iteration"), int)]),)),
                candidate_rejection_count=0,
            ),},)

    return replace(final_scalar, modal=modal, modal_basis=assessed_bundle, eq42_density_profile=assessment_profile, eq42_density_basis_converged=eq42_converged, eq42_density_basis_iterations=len(history), eq42_density_basis_history=history)
