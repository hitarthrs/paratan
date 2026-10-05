"""
Coupled fast D and fast T modal kinetic integration stage

This file iterates species densities, Eq 59 collision operators, reduced fast cross species Rosenbluth fields, shared Eq 68 current balance, shared Eq 70 quasineutrality, Eq 42 density weighting, and final convergence metadata
"""
from __future__ import annotations
from dataclasses import replace
import inspect
import os
from time import perf_counter
import numpy as np
from source_model_revamp.fbis.collision_parameters import collision_parameter_metadata, energy_J_from_keV
from source_model_revamp.fbis.modal.density_weighting import EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, build_eq42_density_profile, compare_eq42_density_profiles, eq42_density_weighting_model, relax_eq42_density_profile
from source_model_revamp.fbis.modal.source_projection import NO_CONFINED_FBIS_SOURCE, audit_attenuated_modal_source_projection, modal_source_projection_audit_metadata
from source_model_revamp.fbis.modal.species_state import FastIonSpeciesRequest, FastIonSpeciesState, FastIonSystemState
from source_model_revamp.fbis.modal.solver import solve_fast_ion_system
from source_model_revamp.fbis.modal.eq59.cross_species import build_external_fast_field_collision_states, external_fast_field_rosenbluth_coefficients
from source_model_revamp.fbis.modal.eq59.operator import _assemble_eq59_mode_tridiagonal
from source_model_revamp.fbis.species import DEUTERON, TRITON
from source_model_revamp.fbis.velocity_space_grid import gyrotropic_velocity_cell_volumes
from source_model_revamp.integration.evaluation import OperatingPointEvaluationMode
from source_model_revamp.integration.modal_stage.basis_density import _build_modal_basis_bundle_for_eq42_profile, _collision_state_for_density, _modal_numerics, build_reusable_modal_basis
from source_model_revamp.integration.modal_stage.density_closure import _compare_eq42_bases, _eq42_metadata, _seed_reference_species_on_confined_grid, _relative_array_change, _relative_scalar_change
from source_model_revamp.integration.modal_stage.full_device import _full_device_fast_ion_population
from source_model_revamp.integration.modal_stage.metadata_adapter import KineticMetadataContract, ordered_active_failure_reasons
from source_model_revamp.integration.modal_stage.types import OperatingPointDensityState, _ReusableModalBasisBundle
from source_model_revamp.integration.pipeline_common import distribution_roundoff
from source_model_revamp.integration.pipeline_types import BeamEnsembleResult, GeometryStageResult, KineticStageResult
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.plasma_closure import uses_startup_seed_only
from source_model_revamp.runtime_profile import finalize_runtime_profile, new_runtime_profile, record_runtime

_EPS_DENSITY_M3 = 1.0e6
_SPECIES_BY_ID = {DEUTERON.species_id: DEUTERON, TRITON.species_id: TRITON}
_EQ70_INTERPOLATION_AUDIT_ENV = "SOURCE_MODEL_EQ70_INTERPOLATION_AUDIT"

def _eq70_confined_interpolation_audit_enabled() -> bool:
    """Return whether the optional confined Eq 70 interpolation audit is enabled by the environment"""
    value = os.environ.get(_EQ70_INTERPOLATION_AUDIT_ENV, "").strip().lower()
    return value in {"1", "true", "yes", "on"}

def _coupled_source_projection_contract(source_audits: dict[str, object], tolerance: float) -> tuple[bool, float, dict[str, float], dict[str, bool]]:
    """Check modal source projection rate conservation for every active fast ion species and return the worst relative error"""
    errors: dict[str, float] = {}
    passed: dict[str, bool] = {}
    for species_id in sorted(source_audits):
        audit = source_audits[species_id]
        confined = float(audit.confined_physical_birth_rate_s)
        projected = float(audit.modal_projected_confined_rate_s)
        error = 0.0 if confined <= 0.0 else (projected - confined) / confined
        errors[species_id] = error
        passed[species_id] = bool(confined <= 0.0 or (projected > 0.0 and abs(error) <= tolerance))
    worst = max(errors.values(), key=abs, default=0.0)
  
    return bool(passed and all(passed.values())), float(worst), errors, passed

def _coupled_geometry_contract(config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult) -> tuple[bool, bool, float]:
    """Check kinetic domain, beam pitch, plasma boundary, and left to right mirror symmetry requirements for the coupled stage"""
    mirror_throat_inside_domain = geometry.metadata.get("fitted_mirror_throat_inside_kinetic_domain") is True
    beam_pitch_selection_valid = beam.metadata.get("beam_pitch_selection_valid") is True
    boundary_selection_valid = bool(geometry.plasma_boundaries is not None and geometry.plasma_boundaries.valid and geometry.B_T_function is not None and geometry.full_device_z_centers_m is not None)
    left_mirror_ratio = geometry.metadata.get("left_fitted_mirror_ratio")
    right_mirror_ratio = geometry.metadata.get("right_fitted_mirror_ratio")
    if left_mirror_ratio is None or right_mirror_ratio is None:
        symmetry_error = 0.0
    else:
        symmetry_error = abs(float(left_mirror_ratio) - float(right_mirror_ratio)) / max(abs(float(left_mirror_ratio)), abs(float(right_mirror_ratio)), np.finfo(float).tiny)
    symmetric = bool(symmetry_error <= float(config.kinetic_electrostatic.modal_loss_convention_relative_tolerance))
  
    return bool(mirror_throat_inside_domain and beam_pitch_selection_valid and boundary_selection_valid), symmetric, float(symmetry_error)

def _coupled_density_roles_separated(config: SourceModelRunConfig, modal, density_state: OperatingPointDensityState) -> tuple[bool, float | None, float | None]:
    """
    Check that Eq 68, Eq 70, and collision density roles remain distinct and consistent
    
    Eq 68 must use the solved electron midplane density, collisions must use the solved collision density, and Eq 70 must derive rather than prescribe its midplane density
    """
    tolerance = float(config.kinetic_electrostatic.modal_density_closure_relative_tolerance)
    eq68_midplane_density = modal.metadata.get("electron_midplane_density_used_by_Eq68_m3")
    collision_density = modal.metadata.get("electron_collision_density_m3")
    midpoint_error = None if eq68_midplane_density is None else abs(float(eq68_midplane_density) - density_state.electron_midplane_density_m3) / max(abs(float(eq68_midplane_density)), abs(density_state.electron_midplane_density_m3), 1.0)
    collision_error = None if collision_density is None else abs(float(collision_density) - density_state.electron_collision_density_m3) / max(abs(float(collision_density)), abs(density_state.electron_collision_density_m3), 1.0)
    eq70_derived = modal.metadata.get("electron_midplane_density_prescribed_to_Eq70_m3") is None
    separated = bool(midpoint_error is not None and midpoint_error <= tolerance and collision_error is not None and collision_error <= tolerance and eq70_derived)
   
    return separated, midpoint_error, collision_error

def _positive_relaxed(current: float, target: float, relaxation: float) -> float:
    """Return logarithmic relaxation between positive current and target density values"""
    current_value = max(float(current), _EPS_DENSITY_M3)
    target_value = max(float(target), _EPS_DENSITY_M3)
   
    return float(np.exp((1.0 - relaxation) * np.log(current_value) + relaxation * np.log(target_value)))

def _eq59_iteration_fixed_point_converged(modal_result: object) -> bool:
    """Return whether every active Eq 59 fixed point, residual, inventory, effective energy, and tail indicator check passed for one species result"""
    metadata = modal_result.metadata
    if metadata.get("modal_rosenbluth_convergence_applicable") is not True:
        return True
   
    return bool(
        metadata.get("modal_rosenbluth_fixed_point_iteration_converged") is True
        and metadata.get("modal_rosenbluth_nonlinear_residual_converged") is True
        and metadata.get("modal_rosenbluth_inventory_converged") is True
        and metadata.get("modal_rosenbluth_effective_energy_converged") is True
        and metadata.get("modal_rosenbluth_single_grid_tail_indicator_passed") is True
    )

def _density_iteration_relative_tolerance(config: SourceModelRunConfig) -> tuple[float, bool]:
    """Return the scalar density iteration tolerance and whether it falls back to the final closure tolerance"""
    configured = config.kinetic_electrostatic.modal_density_iteration_relative_tolerance
    if configured is None:
        return float(config.kinetic_electrostatic.modal_density_closure_relative_tolerance), True
   
    return float(configured), False

def _metadata_first(metadata: dict[str, object], *keys: str) -> object | None:
    """Return the first nonempty metadata value from an ordered sequence of keys"""
    for key in keys:
        value = metadata.get(key)
        if value is not None:
            return value
   
    return None

def _float_or_none(value: object | None) -> float | None:
    """Return one finite float or None when the value is absent or invalid"""
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
   
    return result if np.isfinite(result) else None

def _object_float_first(value: object | None, *names: str) -> float | None:
    """Return the first finite float found on an object from an ordered sequence of attribute names"""
    if value is None:
        return None
    for name in names:
        candidate = getattr(value, name, None)
        result = _float_or_none(candidate)
        if result is not None:
            return result
   
    return None

def _apply_eq59_tridiagonal_audit(lower: np.ndarray, diagonal: np.ndarray, upper: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Apply one tridiagonal Eq 59 matrix to a speed vector without modifying the input arrays"""
    result = np.asarray(diagonal, dtype=float) * np.asarray(values, dtype=float)
    result[1:] += np.asarray(lower, dtype=float) * values[:-1]
    result[:-1] += np.asarray(upper, dtype=float) * values[1:]
  
    return result

def _audit_relative_difference(value: float | None, reference: float | None) -> float | None:
    """Return the symmetric relative difference between two finite scalar audit values"""
    if value is None or reference is None:
        return None
    scale = max(abs(float(value)), abs(float(reference)), np.finfo(float).tiny)
   
    return abs(float(value) - float(reference)) / scale

def _eq59_cx_discrete_particle_budget(*, modal: object, volume_m3: float) -> dict[str, object]:
    """
    Audit the discrete Eq 59 source, pitch loss, charge exchange sink, and particle residual for one solved species
    
    The audit compares exact speed shell measure and matrix center width forms
    When the Eq 59 matrix signature exposes one charge exchange frequency argument, it also reassembles the matrix with and without that sink to isolate its discrete contribution
    """
    base: dict[str, object] = {"available": False, "failure_reason": None, "matrix_reassembly_available": False, "matrix_reassembly_failure_reason": None}
    try:
        volume = float(volume_m3)
        if not np.isfinite(volume) or volume <= 0.0:
            raise ValueError("Eq 59 particle budget requires a positive finite confined volume")
        warm = modal.eq59_warm_start_state
        collision = modal.eq59_collision_operator_state
        if warm is None or collision is None:
            raise ValueError("Eq 59 particle budget requires final warm start and collision operator states")
        mode_eta_raw = getattr(warm, "mode_eta_integrals", None)
        source_coefficients_raw = getattr(warm, "source_coefficients_by_component_j", None)
        source_speeds_raw = getattr(warm, "source_speeds_m_s", None)
        eigenvalues_raw = getattr(warm, "active_eigenvalues_by_speed", None)
        cx_frequency_raw = getattr(warm, "charge_exchange_loss_frequency_s", None)
        if any(value is None for value in (mode_eta_raw, source_coefficients_raw, source_speeds_raw, eigenvalues_raw, cx_frequency_raw)):
            raise ValueError("Eq 59 particle budget requires source, eigenvalue, eta measure, and charge exchange frequency state")

        distribution = np.asarray(modal.modal_distribution_f_j_v, dtype=float)
        speed = np.asarray(modal.speed_grid.centers_m_s, dtype=float)
        widths = np.asarray(modal.speed_grid.widths_m_s, dtype=float)
        shell = np.asarray(modal.speed_grid.shell_volumes_m3_s3, dtype=float)
        mode_eta = np.asarray(mode_eta_raw, dtype=float)
        source_coefficients = np.asarray(source_coefficients_raw, dtype=float)
        source_speeds = np.asarray(source_speeds_raw, dtype=float)
        eigenvalues = np.asarray(eigenvalues_raw, dtype=float)
        cx_frequency = np.asarray(cx_frequency_raw, dtype=float)
        pitch = np.asarray(collision.ion_pitch_scattering_velocity_cubed_m3_s3, dtype=float)
        tau_s = float(collision.spitzer_slowing_down_time_s)

        if distribution.ndim != 2 or distribution.shape[1] != speed.size:
            raise ValueError("Eq 59 particle budget modal distribution shape is invalid")
        mode_count, speed_count = distribution.shape
        if mode_eta.shape != (mode_count,) or np.any(~np.isfinite(mode_eta)):
            raise ValueError("Eq 59 particle budget eta integrals do not match the modal distribution")
        if source_coefficients.ndim != 2 or source_coefficients.shape[1] != mode_count:
            raise ValueError("Eq 59 particle budget source coefficients do not match the modal distribution")
        if source_speeds.shape != (source_coefficients.shape[0],):
            raise ValueError("Eq 59 particle budget source speeds do not match the source coefficients")
        if eigenvalues.shape == (mode_count,):
            eigenvalues = np.broadcast_to(eigenvalues[:, None], distribution.shape).copy()
        if eigenvalues.shape != distribution.shape:
            raise ValueError("Eq 59 particle budget eigenvalues do not match the modal distribution")
        if any(array.shape != (speed_count,) for array in (widths, shell, cx_frequency, pitch)):
            raise ValueError("Eq 59 particle budget speed arrays do not match the production speed grid")
        if np.any(~np.isfinite(distribution)) or np.any(~np.isfinite(source_coefficients)) or np.any(~np.isfinite(source_speeds)):
            raise ValueError("Eq 59 particle budget source or distribution state is nonfinite")
        if np.any(~np.isfinite(eigenvalues)) or np.any(eigenvalues <= 0.0):
            raise ValueError("Eq 59 particle budget eigenvalues must be positive and finite")
        if np.any(~np.isfinite(speed)) or np.any(speed <= 0.0) or np.any(widths <= 0.0) or np.any(shell <= 0.0):
            raise ValueError("Eq 59 particle budget speed grid must be positive and finite")
        if np.any(~np.isfinite(cx_frequency)) or np.any(cx_frequency < 0.0):
            raise ValueError("Eq 59 particle budget charge exchange frequency must be finite and nonnegative")
        if np.any(~np.isfinite(pitch)) or np.any(pitch < 0.0):
            raise ValueError("Eq 59 particle budget pitch scattering state must be finite and nonnegative")
        if not np.isfinite(tau_s) or tau_s <= 0.0:
            raise ValueError("Eq 59 particle budget slowing time must be positive and finite")

        projected_source_factor = 4.0 * np.pi * volume
        shell_measure_factor = volume
        matrix_factor = projected_source_factor / tau_s
        eta_integrated = mode_eta @ distribution
        projected_source = float(projected_source_factor * np.sum(source_coefficients * mode_eta[None, :]))
        inventory_from_shell = float(shell_measure_factor * np.sum(shell * eta_integrated))
        inventory_reference = _float_or_none(getattr(modal, "inventory_particles", None))
        shell_pitch_sink = float(shell_measure_factor * np.sum(distribution * eigenvalues * mode_eta[:, None] * pitch[None, :] * shell[None, :] / (tau_s * speed[None, :] ** 3)))
        matrix_pitch_sink = float(matrix_factor * np.sum(distribution * eigenvalues * mode_eta[:, None] * pitch[None, :] * widths[None, :] / speed[None, :]))
        realized_shell_cx_sink = float(shell_measure_factor * np.sum(eta_integrated * cx_frequency * shell))
        center_width_cx_sink = float(projected_source_factor * np.sum(eta_integrated * cx_frequency * speed**2 * widths))
        metadata_projected_source = _float_or_none(_metadata_first(modal.metadata, "modal_projected_source_rate_s", "modal_projected_source_particle_rate_from_coefficients_s"))
        metadata_eq61 = _float_or_none(_metadata_first(modal.metadata, "modal_two_end_eq61_loss_rate_s", "modal_confined_two_end_eq61_particle_loss_rate_s"))
        metadata_cx_reference = _float_or_none(modal.metadata.get("modal_charge_exchange_target_sink_rate_s"))

        base.update({
            "available": True,
            "phase_space_measure": "speed_shell_volumes_already_include_4_pi_v_squared_dv_folded_even_parallel_signs",
            "source_and_matrix_measure_factor": "4_pi_times_confined_volume",
            "shell_measure_factor": "confined_volume_only",
            "speed_cell_count": int(speed_count),
            "mode_count": int(mode_count),
            "component_count": int(source_coefficients.shape[0]),
            "spitzer_slowing_down_time_s": tau_s,
            "charge_exchange_loss_frequency_nonzero_cell_count": int(np.count_nonzero(cx_frequency > 0.0)),
            "charge_exchange_loss_frequency_max_s": float(np.max(cx_frequency)) if cx_frequency.size else 0.0,
            "eta_integrated_distribution_minimum": float(np.min(eta_integrated)) if eta_integrated.size else 0.0,
            "inventory_from_exact_speed_shell_measure_particles": inventory_from_shell,
            "modal_result_inventory_particles": inventory_reference,
            "inventory_relative_difference": _audit_relative_difference(inventory_from_shell, inventory_reference),
            "projected_source_from_final_eq59_coefficients_rate_s": projected_source,
            "metadata_projected_source_rate_s": metadata_projected_source,
            "projected_source_metadata_relative_difference": _audit_relative_difference(projected_source, metadata_projected_source),
            "shell_measure_pitch_sink_rate_s": shell_pitch_sink,
            "matrix_center_width_pitch_sink_rate_s": matrix_pitch_sink,
            "eq61_two_end_reference_rate_s": metadata_eq61,
            "shell_pitch_to_eq61_relative_difference": _audit_relative_difference(shell_pitch_sink, metadata_eq61),
            "matrix_pitch_to_eq61_relative_difference": _audit_relative_difference(matrix_pitch_sink, metadata_eq61),
            "shell_to_matrix_pitch_relative_difference": _audit_relative_difference(shell_pitch_sink, matrix_pitch_sink),
            "realized_shell_measure_charge_exchange_sink_rate_s": realized_shell_cx_sink,
            "center_width_charge_exchange_sink_rate_s": center_width_cx_sink,
            "charge_exchange_reference_rate_s": metadata_cx_reference,
            "realized_shell_cx_to_reference_relative_difference": _audit_relative_difference(realized_shell_cx_sink, metadata_cx_reference),
            "center_width_cx_to_reference_relative_difference": _audit_relative_difference(center_width_cx_sink, metadata_cx_reference),
            "shell_to_center_width_cx_relative_difference": _audit_relative_difference(realized_shell_cx_sink, center_width_cx_sink),
            "projected_shell_particle_residual_rate_s": float(projected_source - shell_pitch_sink - realized_shell_cx_sink),
            "projected_shell_particle_residual_relative": float((projected_source - shell_pitch_sink - realized_shell_cx_sink) / max(abs(projected_source), np.finfo(float).tiny)),
            "frozen_operator_speed_refinement_snapshot": {
                "speed_faces_m_s": np.asarray(modal.speed_grid.faces_m_s, dtype=float).copy(),
                "modal_distribution_f_j_v": distribution.copy(),
                "mode_eta_integrals": mode_eta.copy(),
                "source_coefficients_by_component_j": source_coefficients.copy(),
                "source_speeds_m_s": source_speeds.copy(),
                "active_eigenvalues_by_speed": eigenvalues.copy(),
                "charge_exchange_loss_frequency_s": cx_frequency.copy(),
                "ion_drag_velocity_cubed_m3_s3": np.asarray(collision.ion_drag_velocity_cubed_m3_s3, dtype=float).copy(),
                "ion_energy_diffusion_velocity_fourth_m4_s4": np.asarray(collision.ion_energy_diffusion_velocity_fourth_m4_s4, dtype=float).copy(),
                "ion_pitch_scattering_velocity_cubed_m3_s3": pitch.copy(),
                "spitzer_slowing_down_time_s": tau_s,
                "volume_m3": volume,
                "physical_confined_source_rate_s": _float_or_none(_metadata_first(modal.metadata, "modal_source_birth_rate_s", "modal_physical_birth_source_rate_s")),
                "projected_source_rate_s": projected_source,
                "eq61_two_end_reference_rate_s": metadata_eq61,
                "charge_exchange_reference_rate_s": metadata_cx_reference,
                "interpolation_model": "piecewise_linear_center_values_with_constant_endpoint_extension",
                "refinement_model": "subdivide_each_production_speed_cell_uniformly_while_retaining_all_original_faces",
                "operator_model": getattr(collision, "operator_model", None),
            },
        })

        try:
            signature = inspect.signature(_assemble_eq59_mode_tridiagonal)
            parameter_names = tuple(signature.parameters)
            base["matrix_assembly_parameter_names"] = parameter_names
            cx_parameter_candidates = tuple(
                name for name in parameter_names
                if ("charge_exchange" in name.lower() or "cx" in name.lower()) and "frequency" in name.lower()
            )
            base["matrix_charge_exchange_frequency_parameter_candidates"] = cx_parameter_candidates
            if len(cx_parameter_candidates) != 1:
                base["matrix_reassembly_failure_reason"] = "exact charge exchange frequency parameter was not uniquely identifiable in the Eq 59 matrix assembler"
                return base

            cx_parameter = cx_parameter_candidates[0]
            zero_cx = np.zeros_like(cx_frequency)
            full_lhs = np.zeros_like(distribution)
            zero_cx_lhs = np.zeros_like(distribution)
            full_rhs = np.zeros_like(distribution)
            zero_rhs = np.zeros_like(distribution)
            for mode_index in range(mode_count):
                common = {
                    "speed_grid": modal.speed_grid,
                    "source_coefficients_by_component": source_coefficients[:, mode_index],
                    "source_speeds_m_s": source_speeds,
                    "eigenvalue": eigenvalues[mode_index],
                    "collision_operator_state": collision,
                }
                full_arguments = dict(common)
                full_arguments[cx_parameter] = cx_frequency
                zero_arguments = dict(common)
                zero_arguments[cx_parameter] = zero_cx
                lower, diagonal, upper, rhs = _assemble_eq59_mode_tridiagonal(**full_arguments)
                zero_lower, zero_diagonal, zero_upper, rhs_zero = _assemble_eq59_mode_tridiagonal(**zero_arguments)
                full_lhs[mode_index] = _apply_eq59_tridiagonal_audit(lower, diagonal, upper, distribution[mode_index])
                zero_cx_lhs[mode_index] = _apply_eq59_tridiagonal_audit(zero_lower, zero_diagonal, zero_upper, distribution[mode_index])
                full_rhs[mode_index] = np.asarray(rhs, dtype=float)
                zero_rhs[mode_index] = np.asarray(rhs_zero, dtype=float)

            if np.any(~np.isfinite(full_lhs)) or np.any(~np.isfinite(zero_cx_lhs)) or np.any(~np.isfinite(full_rhs)):
                raise ValueError("Eq 59 matrix particle budget reassembly produced nonfinite values")
            rhs_difference = float(np.max(np.abs(full_rhs - zero_rhs)))
            rhs_scale = max(float(np.max(np.abs(full_rhs))) if full_rhs.size else 0.0, np.finfo(float).tiny)
            rhs_relative_difference = rhs_difference / rhs_scale
            weighted_full_lhs = float(np.sum(mode_eta[:, None] * full_lhs))
            weighted_rhs = float(np.sum(mode_eta[:, None] * full_rhs))
            weighted_cx_lhs = float(np.sum(mode_eta[:, None] * (full_lhs - zero_cx_lhs)))
            matrix_source = float(-matrix_factor * weighted_rhs)
            matrix_cx_sink = float(-matrix_factor * weighted_cx_lhs)
            full_lhs_rate = float(matrix_factor * weighted_full_lhs)
            speed_operator_signed_rate = float(full_lhs_rate + matrix_pitch_sink + matrix_cx_sink)
            matrix_residual_rate = float(matrix_factor * np.sum(mode_eta[:, None] * (full_lhs - full_rhs)))
            decomposed_residual_rate = float(matrix_source + speed_operator_signed_rate - matrix_pitch_sink - matrix_cx_sink)
            base.update({
                "matrix_reassembly_available": True,
                "matrix_charge_exchange_frequency_parameter": cx_parameter,
                "matrix_rhs_zero_cx_maximum_relative_difference": rhs_relative_difference,
                "matrix_source_rate_s": matrix_source,
                "matrix_source_to_projected_coefficients_relative_difference": _audit_relative_difference(matrix_source, projected_source),
                "matrix_charge_exchange_sink_rate_s": matrix_cx_sink,
                "matrix_cx_to_realized_shell_relative_difference": _audit_relative_difference(matrix_cx_sink, realized_shell_cx_sink),
                "matrix_cx_to_reference_relative_difference": _audit_relative_difference(matrix_cx_sink, metadata_cx_reference),
                "matrix_speed_operator_signed_particle_rate_s": speed_operator_signed_rate,
                "matrix_full_lhs_rate_s": full_lhs_rate,
                "matrix_particle_residual_rate_s": matrix_residual_rate,
                "matrix_particle_residual_relative_to_source": float(matrix_residual_rate / max(abs(matrix_source), np.finfo(float).tiny)),
                "matrix_decomposed_particle_residual_rate_s": decomposed_residual_rate,
                "matrix_residual_decomposition_relative_difference": _audit_relative_difference(matrix_residual_rate, decomposed_residual_rate),
            })
        except Exception as exc:
            base["matrix_reassembly_failure_reason"] = f"{type(exc).__name__}: {exc}"
            return base
        return base
    except Exception as exc:
        base["failure_reason"] = f"{type(exc).__name__}: {exc}"
      
        return base

def _source_to_loss_particle_authority_chain(*, beam: BeamEnsembleResult, system: FastIonSystemState, source_audits: dict[str, object], confined_ids: tuple[str, ...], volume_m3: float) -> dict[str, dict[str, object]]:
    """
    Build species resolved source to loss particle accounting from beam deposition through the modal sinks
    
    The record compares attenuated births, prompt magnetic loss, confined projected source, Eq 61 throat loss, Eq 63 directed throat crossing, charge exchange sink, and final discrete Eq 59 particle budget diagnostics
    """
    charge_exchange_sinks = dict(beam.charge_exchange_sink_by_target_species or {})
    records: dict[str, dict[str, object]] = {}
    for species_id in confined_ids:
        state = system.state_for(species_id)
        modal = state.modal_result
        if modal is None:
            continue
        audit = source_audits[species_id]
        attenuated_source = beam.attenuated_source_by_species[species_id]
        metadata = modal.metadata
        prompt_rate = float(audit.prompt_magnetic_loss_birth_rate_s)
        confined_rate = float(audit.confined_physical_birth_rate_s)
        audit_total_rate = confined_rate + prompt_rate
        eq61_left = _float_or_none(_metadata_first(metadata, "modal_eq61_left_one_end_particle_loss_rate_s", "left_ion_loss_rate_s"))
        eq61_right = _float_or_none(_metadata_first(metadata, "modal_eq61_right_one_end_particle_loss_rate_s", "right_ion_loss_rate_s"))
        eq61_two_end = _float_or_none(_metadata_first(metadata, "modal_two_end_eq61_loss_rate_s", "modal_confined_two_end_eq61_particle_loss_rate_s"))
        if eq61_two_end is None and eq61_left is not None and eq61_right is not None:
            eq61_two_end = eq61_left + eq61_right
        eq63_left = _float_or_none(metadata.get("modal_eq63_left_throat_crossing_rate_s"))
        eq63_right = _float_or_none(metadata.get("modal_eq63_right_throat_crossing_rate_s"))
        eq63_two_end = None if eq63_left is None or eq63_right is None else eq63_left + eq63_right
        sink_state = charge_exchange_sinks.get(species_id)
        warm_state = modal.eq59_warm_start_state
        warm_frequency = None if warm_state is None else getattr(warm_state, "charge_exchange_loss_frequency_s", None)
        if warm_frequency is None:
            frequency_available = False
            frequency_nonzero_count = 0
            frequency_max = None
        else:
            frequency = np.asarray(warm_frequency, dtype=float)
            frequency_available = bool(frequency.ndim == 1 and frequency.size > 0 and np.all(np.isfinite(frequency)) and np.all(frequency >= 0.0))
            frequency_nonzero_count = int(np.count_nonzero(frequency > 0.0)) if frequency_available else 0
            frequency_max = float(np.max(frequency)) if frequency_available and frequency.size else None
        records[species_id] = {
            "schema_version": 2,
            "source_state": str(audit.source_state),
            "system_request_total_physical_birth_rate_s": float(state.source_particle_rate_s),
            "beam_attenuated_total_fast_birth_rate_s": _object_float_first(attenuated_source, "total_birth_rate_s"),
            "beam_attenuated_total_ionization_rate_s": _object_float_first(attenuated_source, "total_ionization_rate_s"),
            "beam_attenuated_total_charge_exchange_birth_rate_s": _object_float_first(attenuated_source, "total_charge_exchange_rate_s"),
            "beam_attenuated_total_net_fueling_rate_s": _object_float_first(attenuated_source, "total_net_fueling_rate_s"),
            "source_audit_total_physical_birth_rate_s": float(audit_total_rate),
            "source_audit_confined_physical_birth_rate_s": confined_rate,
            "source_audit_prompt_magnetic_loss_birth_rate_s": prompt_rate,
            "source_audit_modal_projected_confined_rate_s": float(audit.modal_projected_confined_rate_s),
            "eq59_physical_confined_source_rate_s": _float_or_none(_metadata_first(metadata, "modal_source_birth_rate_s", "modal_physical_birth_source_rate_s")),
            "eq59_total_physical_birth_rate_s": _float_or_none(_metadata_first(metadata, "modal_total_physical_birth_source_rate_s", "modal_total_beam_birth_rate_s")),
            "eq59_projected_source_rate_s": _float_or_none(_metadata_first(metadata, "modal_projected_source_rate_s", "modal_projected_source_particle_rate_from_coefficients_s")),
            "eq59_eigenvalue_sink_rate_s": _float_or_none(_metadata_first(metadata, "modal_eigenvalue_sink_particle_rate_s", "modal_eigenvalue_sink_rate_s")),
            "eq61_left_one_end_rate_s": eq61_left,
            "eq61_right_one_end_rate_s": eq61_right,
            "eq61_two_end_rate_s": eq61_two_end,
            "eq63_left_throat_rate_s": eq63_left,
            "eq63_right_throat_rate_s": eq63_right,
            "eq63_two_end_throat_rate_s": eq63_two_end,
            "prompt_magnetic_loss_rate_s": prompt_rate,
            "modal_total_device_particle_loss_rate_s": _float_or_none(metadata.get("modal_total_device_particle_loss_rate_s")),
            "modal_result_ion_particle_loss_rate_s": float(modal.ion_particle_loss_rate_s),
            "modal_particle_balance_relative_error": _float_or_none(_metadata_first(metadata, "modal_particle_balance_relative_error", "modal_particle_balance_relative_error_loss_minus_source_over_source")),
            "modal_particle_balance_scope": metadata.get("modal_particle_balance_scope"),
            "modal_projected_source_to_birth_ratio": _float_or_none(_metadata_first(metadata, "modal_projected_source_to_birth_ratio", "modal_projected_source_to_birth_particle_ratio")),
            "modal_eigenvalue_sink_to_birth_ratio": _float_or_none(_metadata_first(metadata, "modal_eigenvalue_sink_to_birth_ratio", "modal_eigenvalue_sink_to_source_particle_ratio")),
            "modal_two_end_eq61_to_birth_ratio": _float_or_none(metadata.get("modal_two_end_eq61_to_birth_ratio")),
            "modal_two_end_eq61_to_eigenvalue_ratio": _float_or_none(metadata.get("modal_two_end_eq61_to_eigenvalue_ratio")),
            "modal_source_projection_relative_tolerance": _float_or_none(metadata.get("modal_source_projection_relative_tolerance")),
            "modal_loss_convention_relative_tolerance": _float_or_none(metadata.get("modal_loss_convention_relative_tolerance")),
            "modal_loss_convention_validation_passed": metadata.get("modal_loss_convention_validation_passed"),
            "modal_loss_convention_failure_reason": metadata.get("modal_loss_convention_failure_reason"),
            "modal_eq63_left_rate_relative_error": _float_or_none(metadata.get("modal_eq63_left_rate_relative_error")),
            "modal_eq63_right_rate_relative_error": _float_or_none(metadata.get("modal_eq63_right_rate_relative_error")),
            "charge_exchange_target_sink_supplied_to_fast_ion_request": sink_state is not None,
            "charge_exchange_target_sink_request_state_type": None if sink_state is None else type(sink_state).__name__,
            "charge_exchange_target_sink_request_rate_s": _object_float_first(sink_state, "total_sink_rate_s", "total_particle_sink_rate_s", "target_sink_rate_s", "total_rate_s", "particle_rate_s"),
            "eq59_charge_exchange_target_sink_active": metadata.get("modal_charge_exchange_target_sink_active") is True,
            "eq59_charge_exchange_target_sink_rate_s": _float_or_none(metadata.get("modal_charge_exchange_target_sink_rate_s")),
            "eq59_charge_exchange_target_sink_energy_W": _float_or_none(metadata.get("modal_charge_exchange_target_sink_energy_W")),
            "eq59_charge_exchange_target_sink_model": metadata.get("modal_charge_exchange_target_sink_model"),
            "eq59_charge_exchange_target_sink_reference_identity_relative_error": _float_or_none(metadata.get("modal_charge_exchange_target_sink_reference_identity_relative_error")),
            "eq59_charge_exchange_target_sink_pitch_resolution": metadata.get("modal_charge_exchange_target_sink_pitch_resolution"),
            "eq59_charge_exchange_loss_frequency_available": frequency_available,
            "eq59_charge_exchange_loss_frequency_nonzero_cell_count": frequency_nonzero_count,
            "eq59_charge_exchange_loss_frequency_max_s": frequency_max,
            "eq59_cx_discrete_particle_budget": _eq59_cx_discrete_particle_budget(modal=modal, volume_m3=volume_m3),
        }
   
    return records

def _collision_operator_seed_from_system(system: FastIonSystemState, species_ids: tuple[str, ...]) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Extract finite Eq 59 drag, energy diffusion, and pitch scattering arrays from a solved system for continuation"""
    seeds: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for species_id in species_ids:
        modal = system.state_for(species_id).modal_result
        operator = modal.eq59_collision_operator_state
        if operator is None:
            continue
        arrays = (
            np.asarray(operator.ion_drag_velocity_cubed_m3_s3, dtype=float),
            np.asarray(operator.ion_energy_diffusion_velocity_fourth_m4_s4, dtype=float),
            np.asarray(operator.ion_pitch_scattering_velocity_cubed_m3_s3, dtype=float),
        )
        if any(np.any(~np.isfinite(value)) for value in arrays):
            raise ValueError("collision operator continuation seed must be finite")
        seeds[species_id] = tuple(value.copy() for value in arrays)
   
    return seeds

def _collision_operator_continuation_warm_start_from_system(system: FastIonSystemState, species_ids: tuple[str, ...]) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """Extract the collision operator arrays and collision model identity required for scalar continuation across evaluations"""
    operator_seed = _collision_operator_seed_from_system(system, species_ids)
    warm_start: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    for species_id, arrays in operator_seed.items():
        speed = np.asarray(system.state_for(species_id).speed_grid.centers_m_s, dtype=float)
        if np.any(~np.isfinite(speed)):
            raise ValueError("collision operator continuation speed grid must be finite")
        warm_start[species_id] = (speed.copy(), *(value.copy() for value in arrays))
   
    return warm_start

def _compatible_scalar_continuation_seeds(*, confined_ids: tuple[str, ...], beam: BeamEnsembleResult, cross_species_required: bool, collision_operator_seed_by_species: dict[str, tuple[object, object, object, object]] | None, external_fast_field_seed_by_test_species: dict[str, tuple[object, ...]] | None, wall_barrier_energy_seed_J: float | None) -> tuple[dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]], dict[str, tuple[object, ...]], float | None, dict[str, object]]:
    """
    Validate optional collision operator, fast cross species field, and wall barrier continuation seeds
    
    Seeds are accepted only when species, speed grids, pair identities, array shapes, and required components match the current source state
    """
    reasons: list[str] = []
    operator_seed: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    requested_operator = dict(collision_operator_seed_by_species or {})
    for species_id in confined_ids:
        arrays = requested_operator.get(species_id)
        if arrays is None:
            continue
        if len(arrays) != 4:
            reasons.append(f"collision_operator_shape_{species_id}")
            continue
        current_speed = np.asarray(beam.speed_grid_by_species[species_id].centers_m_s, dtype=float)
        stored_speed = np.asarray(arrays[0], dtype=float)
        copied = tuple(np.asarray(value, dtype=float).copy() for value in arrays[1:])
        if stored_speed.shape != current_speed.shape or not np.array_equal(stored_speed, current_speed):
            reasons.append(f"collision_operator_grid_{species_id}")
            continue
        if any(value.shape != current_speed.shape for value in copied):
            reasons.append(f"collision_operator_shape_{species_id}")
            continue
        if np.any(~np.isfinite(stored_speed)) or any(np.any(~np.isfinite(value)) for value in copied):
            reasons.append(f"collision_operator_nonfinite_{species_id}")
            continue
        operator_seed[species_id] = copied

    external_seed: dict[str, tuple[object, ...]] = {}
    requested_external = {key: tuple(values) for key, values in dict(external_fast_field_seed_by_test_species or {}).items()}
    if requested_external and cross_species_required:
        expected_pairs = {(test_species_id, field_species_id) for test_species_id in confined_ids for field_species_id in confined_ids if field_species_id != test_species_id}
        actual_pairs: set[tuple[str, str]] = set()
        valid_external = True
        for test_species_id, states in requested_external.items():
            if test_species_id not in confined_ids:
                valid_external = False
                reasons.append("external_fast_field_test_species")
                break
            for state in states:
                try:
                    state_test_species_id = str(state.test_species.species_id)
                    field_species_id = str(state.field_species.species_id)
                    field_grid = state.field_speed_grid
                    field_values = np.asarray(state.field_first_mode_distribution, dtype=float)
                except Exception:
                    valid_external = False
                    reasons.append("external_fast_field_structure")
                    break
                if state_test_species_id != test_species_id or field_species_id not in confined_ids or field_species_id == test_species_id:
                    valid_external = False
                    reasons.append("external_fast_field_pair_identity")
                    break
                if str(getattr(state, "model", "")) != "reduced_isotropic_first_mode_Rosenbluth_cross_species":
                    valid_external = False
                    reasons.append("external_fast_field_model")
                    break
                if not np.isfinite(float(state.field_physical_density_m3)) or float(state.field_physical_density_m3) <= 0.0 or not np.isfinite(float(state.coulomb_log_value)) or float(state.coulomb_log_value) <= 0.0:
                    valid_external = False
                    reasons.append("external_fast_field_scalar_state")
                    break
                current_field_grid = beam.speed_grid_by_species[field_species_id]
                if not np.array_equal(np.asarray(field_grid.centers_m_s, dtype=float), np.asarray(current_field_grid.centers_m_s, dtype=float)) or not np.array_equal(np.asarray(field_grid.faces_m_s, dtype=float), np.asarray(current_field_grid.faces_m_s, dtype=float)):
                    valid_external = False
                    reasons.append("external_fast_field_grid")
                    break
                if field_values.shape != np.asarray(current_field_grid.centers_m_s, dtype=float).shape or np.any(~np.isfinite(field_values)):
                    valid_external = False
                    reasons.append("external_fast_field_distribution")
                    break
                actual_pairs.add((test_species_id, field_species_id))
            if not valid_external:
                break
        if valid_external and actual_pairs == expected_pairs:
            external_seed = requested_external
        elif valid_external:
            reasons.append("external_fast_field_pair_set")

    barrier_seed = None
    if wall_barrier_energy_seed_J is not None:
        barrier_value = float(wall_barrier_energy_seed_J)
        if np.isfinite(barrier_value) and barrier_value >= 0.0:
            barrier_seed = barrier_value
        else:
            reasons.append("wall_barrier_nonfinite")

    expected_operator_count = len(confined_ids)
    expected_external_pair_count = len(confined_ids) * max(len(confined_ids) - 1, 0)
    actual_external_pair_count = sum(len(states) for states in external_seed.values())
    requested = bool(requested_operator or requested_external or wall_barrier_energy_seed_J is not None)
    complete = bool(len(operator_seed) == expected_operator_count and barrier_seed is not None and (not cross_species_required or actual_external_pair_count == expected_external_pair_count))
    diagnostics = {
        "requested": requested,
        "active": bool(operator_seed or external_seed or barrier_seed is not None),
        "complete": complete,
        "operator_species_count": len(operator_seed),
        "external_pair_count": actual_external_pair_count,
        "barrier_active": barrier_seed is not None,
        "rejection_reasons": tuple(dict.fromkeys(reasons)),
    }
   
    return operator_seed, external_seed, barrier_seed, diagnostics

def _external_fast_field_pair_map(states_by_test_species: dict[str, tuple[object, ...]]) -> dict[str, object]:
    """Index reduced fast cross species collision states by their stable pair identifier"""
    pairs: dict[str, object] = {}
    for test_species_id, states in states_by_test_species.items():
        for state in states:
            if state.test_species.species_id != test_species_id:
                raise ValueError("external fast field test species mapping is inconsistent")
            pair_id = str(state.pair_id)
            if pair_id in pairs:
                raise ValueError("external fast field pair identifiers must be unique")
            pairs[pair_id] = state
  
    return pairs

def _external_fast_field_fixed_point_change(*, current_by_test_species: dict[str, tuple[object, ...]], target_by_test_species: dict[str, tuple[object, ...]], test_speed_grids_by_species: dict[str, object]) -> tuple[dict[str, dict[str, float]], float, bool]:
    """
    Compare two reduced fast cross species collision field states on the current test species speed grids
    
    The comparison includes field density, Coulomb logarithm, and `h_tilde`, `g_tilde_1`, and `g_tilde_2` Rosenbluth coefficients
    """
    current_pairs = _external_fast_field_pair_map(current_by_test_species)
    target_pairs = _external_fast_field_pair_map(target_by_test_species)
    pair_sets_match = set(current_pairs) == set(target_pairs)
    pair_changes: dict[str, dict[str, float]] = {}
    maximum_change = 0.0
    for pair_id in sorted(set(current_pairs) | set(target_pairs)):
        current = current_pairs.get(pair_id)
        target = target_pairs.get(pair_id)
        if current is None or target is None:
            pair_changes[pair_id] = {"maximum_relative_change": float("inf")}
            maximum_change = float("inf")
            continue
        if current.test_species != target.test_species or current.field_species != target.field_species or current.model != target.model:
            raise ValueError("external fast field pair identity changed during the density fixed point")
        test_species_id = current.test_species.species_id
        test_speed_grid = test_speed_grids_by_species.get(test_species_id)
        if test_speed_grid is None:
            raise ValueError("external fast field convergence requires the test species speed grid")
        current_coefficients = external_fast_field_rosenbluth_coefficients(current, test_speed_grid)
        target_coefficients = external_fast_field_rosenbluth_coefficients(target, test_speed_grid)
        density_change = _relative_scalar_change(current.field_physical_density_m3, target.field_physical_density_m3, _EPS_DENSITY_M3)
        coulomb_log_change = _relative_scalar_change(current.coulomb_log_value, target.coulomb_log_value)
        h_tilde_change = _relative_array_change(current_coefficients.h_tilde, target_coefficients.h_tilde)
        g_tilde_1_change = _relative_array_change(current_coefficients.g_tilde_1, target_coefficients.g_tilde_1)
        g_tilde_2_change = _relative_array_change(current_coefficients.g_tilde_2, target_coefficients.g_tilde_2)
        pair_maximum = float(max(density_change, coulomb_log_change, h_tilde_change, g_tilde_1_change, g_tilde_2_change))
        pair_changes[pair_id] = {
            "field_density_relative_change": float(density_change),
            "coulomb_log_relative_change": float(coulomb_log_change),
            "h_tilde_relative_change": float(h_tilde_change),
            "g_tilde_1_relative_change": float(g_tilde_1_change),
            "g_tilde_2_relative_change": float(g_tilde_2_change),
            "maximum_relative_change": pair_maximum,
        }
        maximum_change = max(maximum_change, pair_maximum)
    if not pair_sets_match:
        maximum_change = float("inf")
  
    return pair_changes, float(maximum_change), bool(pair_sets_match)

def _system_density_state(*, config: SourceModelRunConfig, geometry: GeometryStageResult, system: FastIonSystemState, electron_collision_density_m3: float, history: list[dict[str, object]]) -> OperatingPointDensityState:
    """
    Build the typed confined density state from the shared Eq 70 solution
    
    The state retains fast D, fast T, total positive charge, and electron profiles on cells and exact electrostatic nodes together with midplane and volume average densities
    """
    profile = system.shared_electrostatic_profile
    if profile is None:
        raise ValueError("coupled operating point density state requires a shared Eq 70 profile")
    volumes = np.asarray(geometry.cell_volumes_m3, dtype=float)
    fast_by_species = dict(profile.fast_ion_density_by_species_m3 or {})
    if not fast_by_species and len(system.active_species_ids) == 1:
        fast_by_species[system.active_species_ids[0]] = np.asarray(profile.fast_ion_density_m3, dtype=float)
    fast_D_cells = np.asarray(fast_by_species.get(DEUTERON.species_id, np.zeros_like(volumes)), dtype=float)
    fast_T_cells = np.asarray(fast_by_species.get(TRITON.species_id, np.zeros_like(volumes)), dtype=float)
    electron_cells = np.asarray(profile.electron_density_m3, dtype=float)
    if not (fast_D_cells.shape == fast_T_cells.shape == electron_cells.shape == volumes.shape):
        raise ValueError("coupled operating point cell densities must share the confined grid")
    total_positive_cells = np.asarray(profile.ion_density_m3, dtype=float)
    zeta_nodes = np.asarray(profile.electrostatic_node_zeta, dtype=float)
    node_fast_by_species = dict(profile.electrostatic_node_fast_ion_density_by_species_m3 or {})
    if not node_fast_by_species and len(system.active_species_ids) == 1:
        node_fast_by_species[system.active_species_ids[0]] = np.asarray(profile.electrostatic_node_fast_ion_density_m3, dtype=float)
    fast_D_nodes = np.asarray(node_fast_by_species.get(DEUTERON.species_id, np.zeros_like(zeta_nodes)), dtype=float)
    fast_T_nodes = np.asarray(node_fast_by_species.get(TRITON.species_id, np.zeros_like(zeta_nodes)), dtype=float)
    electron_nodes = np.asarray(profile.electrostatic_node_electron_density_m3, dtype=float)
    if not (fast_D_nodes.shape == fast_T_nodes.shape == electron_nodes.shape == zeta_nodes.shape):
        raise ValueError("coupled operating point node densities are unavailable or inconsistent")
    total_positive_nodes = np.asarray(profile.electrostatic_node_ion_density_m3, dtype=float)
    midplane_index = int(profile.electrostatic_midplane_index)
    fast_D_midplane = float(fast_D_nodes[midplane_index])
    fast_T_midplane = float(fast_T_nodes[midplane_index])
    electron_midplane = float(electron_nodes[midplane_index])
    total_positive_midplane = float(total_positive_nodes[midplane_index])
    exact_midplane_error = abs(electron_midplane - total_positive_midplane) / max(abs(electron_midplane), abs(total_positive_midplane), _EPS_DENSITY_M3)
    volume = float(np.sum(volumes))
    fast_D_average = float(np.sum(fast_D_cells * volumes) / volume)
    fast_T_average = float(np.sum(fast_T_cells * volumes) / volume)
    electron_average = float(np.sum(electron_cells * volumes) / volume)
    positive_inventory = float(np.sum(total_positive_cells * volumes))
    electron_inventory = float(np.sum(electron_cells * volumes))
    inventory_error = abs(electron_inventory - positive_inventory) / max(abs(electron_inventory), abs(positive_inventory), 1.0)
    support = np.asarray(profile.density_support_mask, dtype=bool) if profile.density_support_mask is not None else total_positive_cells > 0.0
   
    return OperatingPointDensityState(zeta_cells=np.asarray(geometry.zeta_centers, dtype=float), zeta_nodes=zeta_nodes, cell_volumes_m3=volumes, fast_deuterium_cell_density_m3=fast_D_cells, fast_deuterium_node_density_m3=fast_D_nodes, total_positive_charge_cell_density_m3=total_positive_cells, total_positive_charge_node_density_m3=total_positive_nodes, electron_cell_density_m3=electron_cells, electron_node_density_m3=electron_nodes, fast_deuterium_midplane_density_m3=fast_D_midplane, electron_midplane_density_m3=electron_midplane, fast_deuterium_confined_volume_average_density_m3=fast_D_average, electron_confined_volume_average_density_m3=electron_average, electron_collision_density_m3=float(electron_collision_density_m3), electron_parent_maxwellian_n0_m3=float(profile.electron_parent_maxwellian_n0_m3), profile_support_mask=np.ones(total_positive_cells.shape, dtype=bool), eq70_residual_support_mask=support, profile_scope='confined_throat_to_throat', electron_profile_scope='confined_throat_to_throat', expander_electron_profile_available=False, scalar_convergence_history=tuple((dict(item) for item in history)), exact_midplane_quasineutrality_relative_error=float(exact_midplane_error), confined_inventory_relative_error=float(inventory_error), fast_tritium_cell_density_m3=fast_T_cells, fast_tritium_node_density_m3=fast_T_nodes, fast_tritium_midplane_density_m3=fast_T_midplane, fast_tritium_confined_volume_average_density_m3=fast_T_average)

def _system_solver_arguments(*, config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, electron_temperature_J: float, ion_temperature_J: float, electron_midplane_density_m3: float, electron_collision_density_m3: float, numerics, modal_basis: _ReusableModalBasisBundle, eq70_initial_potential_energy_J: np.ndarray | None, run_numerical_convergence_assessments: bool, fixed_boundary_reconstruction_policy: str, evaluation_mode: OperatingPointEvaluationMode) -> dict[str, object]:
    """Assemble geometry, basis, electrostatic, loss, and numerical arguments shared by each solve_fast_ion_system call"""
    background_positive = np.zeros_like(np.asarray(geometry.zeta_centers, dtype=float))
    background_left = 0.0
    background_right = 0.0
    background_midplane = 0.0
   
    return {'lambda_grid': beam.lambda_grid, 'pitch_grid': beam.pitch_grid, 'volume_m3': geometry.volume_m3, 'midplane_area_m2': geometry.midplane_area_m2, 'half_length_m': geometry.half_length_m, 'mirror_ratio': geometry.mirror_ratio, 'B_tilde_function': geometry.B_tilde_function, 'zeta_faces': geometry.zeta_edges, 'B_tilde_midpoints': geometry.B_tilde_midpoints, 'electron_temperature_J': float(electron_temperature_J), 'electron_midplane_density_m3': float(electron_midplane_density_m3), 'electron_collision_density_m3': float(electron_collision_density_m3), 'eq70_prescribed_electron_midplane_density_m3': None, 'ion_temperature_J': float(ion_temperature_J), 'cell_volumes_m3': geometry.cell_volumes_m3, 'numerics': numerics, 'auxiliary_input_power_W': 0.0, 'global_loss_rate_model': config.kinetic_electrostatic.ion_loss_closure_model, 'electrostatic_feedback_model': config.kinetic_electrostatic.electrostatic_feedback_model, 'ion_loss_closure_model': config.kinetic_electrostatic.ion_loss_closure_model, 'lost_ion_parallel_temperature_model': config.kinetic_electrostatic.lost_ion_parallel_temperature_model, 'background_positive_charge_density_m3': background_positive, 'background_midplane_positive_charge_density_m3': background_midplane, 'background_left_throat_positive_charge_density_m3': background_left, 'background_right_throat_positive_charge_density_m3': background_right, 'basis': modal_basis.basis, 'basis_convergence_assessment': modal_basis.convergence_assessment, '_published_reference_lambda1': modal_basis.published_reference_lambda1, '_eq70_initial_potential_energy_J': eq70_initial_potential_energy_J, '_run_numerical_convergence_assessments': bool(run_numerical_convergence_assessments), '_fixed_boundary_reconstruction_policy': fixed_boundary_reconstruction_policy, '_evaluation_mode': evaluation_mode}

def _solve_system_scalar_closure(*, config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, source_audits: dict[str, object], electron_temperature_J: float, ion_temperature_J: float, numerics, modal_basis: _ReusableModalBasisBundle, initial_eq59_warm_start_states: dict[str, object] | None, initial_fast_midpoint_density_m3: dict[str, float] | None, initial_fast_average_density_m3: dict[str, float] | None, initial_electron_midpoint_density_m3: float | None, initial_electron_confined_average_density_m3: float | None, initial_electron_collision_density_m3: float | None, initial_eq70_potential_energy_J: np.ndarray | None, run_numerical_convergence_assessments: bool, fixed_boundary_reconstruction_policy: str, evaluation_mode: OperatingPointEvaluationMode, ion_loss_rates_s_by_component_override: dict[str, float] | None=None, initial_external_fast_field_collision_states_by_test_species: dict[str, tuple[object, ...]] | None=None, initial_previous_barrier_energy_J: float | None=None, initial_previous_collision_operators_by_species: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] | None=None, require_same_source_fixed_point_confirmation: bool=False, runtime_context_label: str='scalar_density_closure', runtime_profile: dict[str, object] | None=None) -> tuple[FastIonSystemState, OperatingPointDensityState, dict[str, object], bool, bool, list[dict[str, object]], float]:
    """
    Solve the coupled scalar density fixed point for the current Eq 42 basis
    
    The loop updates fast D and fast T densities, electron midplane and collision densities, Eq 59 collision operators, shared wall barrier, Eq 70 profile, and reduced fast cross species collision fields
    Optional final assessment reruns the same state with numerical convergence checks enabled
    """
    closure_model = str(config.kinetic_electrostatic.modal_density_closure_model).strip().lower()
    if closure_model not in {"nbi_supported_stationary"}:
        raise ValueError("coupled fast D and fast T operation requires a supported quasineutral closure")
    source_ids = tuple(sorted(beam.attenuated_source_by_species))
    confined_ids = tuple(key for key in source_ids if source_audits[key].source_state != NO_CONFINED_FBIS_SOURCE)
    charge_exchange_sinks = dict(beam.charge_exchange_sink_by_target_species or {})
    missing_sink_species = tuple(sorted(set(charge_exchange_sinks) - set(source_ids)))
    if missing_sink_species:
        raise ValueError(f"charge exchange target sinks require kinetic source states for species {missing_sink_species}")
    prompt_sink_species = tuple(sorted(set(charge_exchange_sinks) - set(confined_ids)))
    if prompt_sink_species:
        raise ValueError(f"charge exchange target sinks require confined kinetic states for species {prompt_sink_species}")
    requests = {species_id: FastIonSpeciesRequest(species=_SPECIES_BY_ID[species_id], speed_grid=beam.speed_grid_by_species[species_id], attenuated_source=beam.attenuated_source_by_species[species_id], charge_exchange_sink_state=charge_exchange_sinks.get(species_id)) for species_id in source_ids}
    prompt_rates = {key: float(source_audits[key].prompt_magnetic_loss_birth_rate_s) for key in source_ids if source_audits[key].source_state == NO_CONFINED_FBIS_SOURCE}
    prompt_powers = {key: float(sum(item.prompt_magnetic_loss_birth_rate_s * component.spatial_source.energy_J for item, component in zip(source_audits[key].component_audits, beam.attenuated_source_by_species[key].component_sources, strict=True))) for key in prompt_rates}
    volumes = np.asarray(geometry.cell_volumes_m3, dtype=float)
    configured_D_cells, configured_T_cells = _seed_reference_species_on_confined_grid(geometry)
    volume = float(np.sum(volumes))
    configured_D_average = float(np.sum(configured_D_cells * volumes) / volume)
    configured_T_average = float(np.sum(configured_T_cells * volumes) / volume)
    ion_charges = np.ones(2, dtype=float)
    ion_masses = np.asarray([DEUTERON.mass_kg, TRITON.mass_kg], dtype=float)
    seed = _EPS_DENSITY_M3
    default_midpoint = {DEUTERON.species_id: float(config.plasma_closure.background_deuterium_midplane_density_m3), TRITON.species_id: float(config.plasma_closure.background_tritium_midplane_density_m3)}
    default_average = {DEUTERON.species_id: configured_D_average, TRITON.species_id: configured_T_average}
    fast_midpoint = {key: max(float((initial_fast_midpoint_density_m3 or {}).get(key, default_midpoint.get(key, seed))), _EPS_DENSITY_M3) for key in confined_ids}
    fast_average = {key: max(float((initial_fast_average_density_m3 or {}).get(key, default_average.get(key, seed))), _EPS_DENSITY_M3) for key in confined_ids}
    for values in (fast_midpoint, fast_average):
        if any(not np.isfinite(value) or value < 0.0 for value in values.values()):
            raise ValueError("initial fast species densities must be finite and nonnegative")
    background_midplane = 0.0
    electron_midpoint = float(sum(fast_midpoint.values())) if initial_electron_midpoint_density_m3 is None else float(initial_electron_midpoint_density_m3)
    electron_average = float(sum(fast_average.values())) if initial_electron_confined_average_density_m3 is None else float(initial_electron_confined_average_density_m3)
    electron_collision = electron_average if initial_electron_collision_density_m3 is None else float(initial_electron_collision_density_m3)
    if any(not np.isfinite(value) or value <= 0.0 for value in (electron_midpoint, electron_average, electron_collision)):
        raise ValueError("initial electron densities must be finite and positive")
    warm_starts = dict(initial_eq59_warm_start_states or {})
    potential = initial_eq70_potential_energy_J
    limit = max(int(config.kinetic_electrostatic.modal_density_closure_iterations), 1)
    cross_species_required = len(confined_ids) > 1
    if cross_species_required and limit < 2:
        raise ValueError("NBI supported D and T cross collisions require at least two density closure iterations")
    qualification_tolerance = float(config.kinetic_electrostatic.modal_density_closure_relative_tolerance)
    iteration_tolerance, iteration_tolerance_uses_qualification_fallback = _density_iteration_relative_tolerance(config)
    relaxation = float(config.kinetic_electrostatic.modal_density_closure_relaxation)
    if initial_previous_barrier_energy_J is None:
        previous_barrier = None
    else:
        previous_barrier = float(initial_previous_barrier_energy_J)
        if not np.isfinite(previous_barrier) or previous_barrier < 0.0:
            raise ValueError("initial wall barrier seed must be finite and nonnegative")
    previous_operators: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for species_id, arrays in (initial_previous_collision_operators_by_species or {}).items():
        if len(arrays) != 3:
            raise ValueError("initial collision operator seed requires drag diffusion and pitch arrays")
        copied = tuple(np.asarray(value, dtype=float).copy() for value in arrays)
        if any(np.any(~np.isfinite(value)) for value in copied):
            raise ValueError("initial collision operator seed must be finite")
        previous_operators[species_id] = copied
    history: list[dict[str, object]] = []
    system = None
    density_state = None
    collision_states = {}
    relative_error = float("inf")
    converged = False
    iteration_fixed_point_converged = False
    external_fast_fields_by_test_species: dict[str, tuple[object, ...]] = {key: tuple(values) for key, values in (initial_external_fast_field_collision_states_by_test_species or {}).items()}
    assessment_requested = bool(run_numerical_convergence_assessments)
    last_solver_arguments: dict[str, object] | None = None
    last_warm_starts: dict[str, object] = {}
    last_collision_states: dict[str, object] = {}
    last_electron_collision = float(electron_collision)
    # Iterate densities, collision operators, cross species fields, barrier, and Eq 70 state together
    for iteration in range(1, limit + 1):
        collision_scope = 'nbi_supported_fast_ion_self_background'
        field_ion_averages = np.asarray([fast_average.get(DEUTERON.species_id, 0.0), fast_average.get(TRITON.species_id, 0.0)], dtype=float)
        collision_runtime_started = perf_counter()
        collision_states = {species_id: _collision_state_for_density(fast_ion_species=_SPECIES_BY_ID[species_id], electron_density_m3=electron_collision, ion_densities_m3=field_ion_averages, ion_charge_numbers=ion_charges, ion_masses_kg=ion_masses, electron_temperature_J=electron_temperature_J, ion_temperature_J=ion_temperature_J, fast_ion_density_m3=fast_average[species_id], collision_scope=f"{collision_scope}_{species_id}") for species_id in confined_ids}
        record_runtime(runtime_profile, "collision_state_build", perf_counter() - collision_runtime_started, nested=True)
        arguments = _system_solver_arguments(config=config, geometry=geometry, beam=beam, electron_temperature_J=electron_temperature_J, ion_temperature_J=ion_temperature_J, electron_midplane_density_m3=electron_midpoint, electron_collision_density_m3=electron_collision, numerics=numerics, modal_basis=modal_basis, eq70_initial_potential_energy_J=potential, run_numerical_convergence_assessments=False, fixed_boundary_reconstruction_policy=fixed_boundary_reconstruction_policy, evaluation_mode=evaluation_mode)
        arguments["_external_fast_field_collision_states_by_test_species"] = external_fast_fields_by_test_species
        arguments["_runtime_profile_accumulator"] = runtime_profile
        arguments["_runtime_eq70_confined_interpolation_audit_enabled"] = _eq70_confined_interpolation_audit_enabled()
        arguments["_runtime_eq70_solve_context"] = f"{runtime_context_label}:scalar_iteration_{iteration}"
        cross_species_active_this_iteration = bool(not cross_species_required or all(external_fast_fields_by_test_species.get(species_id) for species_id in confined_ids))
        last_solver_arguments = dict(arguments)
        last_warm_starts = dict(warm_starts)
        last_collision_states = dict(collision_states)
        last_electron_collision = float(electron_collision)
        system_runtime_started = perf_counter()
        system = solve_fast_ion_system(requests=requests, collision_states_by_species=collision_states, eq59_warm_start_states_by_species=warm_starts, prompt_only_loss_rates_s_by_species=prompt_rates, prompt_only_midplane_power_W_by_species=prompt_powers, ion_loss_rates_s_by_component_override=ion_loss_rates_s_by_component_override, **arguments)
        record_runtime(runtime_profile, "fast_ion_system_solve", perf_counter() - system_runtime_started, nested=True)
        if system.shared_electrostatic_profile is None or system.shared_current_balance is None:
            raise RuntimeError("coupled fast ion solve did not produce the shared electrostatic state")
        density_runtime_started = perf_counter()
        density_state = _system_density_state(config=config, geometry=geometry, system=system, electron_collision_density_m3=electron_collision, history=history)
        record_runtime(runtime_profile, "operating_density_state_build", perf_counter() - density_runtime_started, nested=True)
        cross_runtime_started = perf_counter()
        next_external_fast_fields = build_external_fast_field_collision_states(system=system, collision_states_by_species=collision_states) if cross_species_required else {}
        record_runtime(runtime_profile, "fast_cross_species_collision_build", perf_counter() - cross_runtime_started, nested=True)
        if cross_species_required:
            cross_species_pair_changes, cross_species_maximum_change, cross_species_pair_sets_match = _external_fast_field_fixed_point_change(
                current_by_test_species=external_fast_fields_by_test_species,
                target_by_test_species=next_external_fast_fields,
                test_speed_grids_by_species={species_id: requests[species_id].speed_grid for species_id in confined_ids},
            )
        else:
            cross_species_pair_changes = {}
            cross_species_maximum_change = 0.0
            cross_species_pair_sets_match = True
        species_changes = {}
        active_changes = []
        new_operators = {}
        for species_id in confined_ids:
            state = system.state_for(species_id)
            modal = state.modal_result
            midpoint_target = density_state.fast_deuterium_midplane_density_m3 if species_id == DEUTERON.species_id else density_state.fast_tritium_midplane_density_m3
            average_target = density_state.fast_deuterium_confined_volume_average_density_m3 if species_id == DEUTERON.species_id else density_state.fast_tritium_confined_volume_average_density_m3
            midpoint_change = _relative_scalar_change(fast_midpoint[species_id], midpoint_target, _EPS_DENSITY_M3)
            average_change = _relative_scalar_change(fast_average[species_id], average_target, _EPS_DENSITY_M3)
            operator = modal.eq59_collision_operator_state
            if operator is None:
                operator_changes = (float("inf"), float("inf"), float("inf")) if species_id not in previous_operators else (0.0, 0.0, 0.0)
            else:
                arrays = (np.asarray(operator.ion_drag_velocity_cubed_m3_s3, dtype=float), np.asarray(operator.ion_energy_diffusion_velocity_fourth_m4_s4, dtype=float), np.asarray(operator.ion_pitch_scattering_velocity_cubed_m3_s3, dtype=float))
                previous = previous_operators.get(species_id)
                operator_changes = (float("inf"), float("inf"), float("inf")) if previous is None else tuple(_relative_array_change(previous[index], arrays[index]) for index in range(3))
                new_operators[species_id] = tuple(value.copy() for value in arrays)
            species_changes[species_id] = {"midpoint_relative_change": midpoint_change, "volume_average_relative_change": average_change, "ion_drag_relative_change": operator_changes[0], "energy_diffusion_relative_change": operator_changes[1], "pitch_scattering_relative_change": operator_changes[2], "Eq59_converged": bool(modal.metadata.get("modal_rosenbluth_converged", False)), "Eq61_Eq63_left_relative_error": modal.metadata.get("modal_eq63_left_rate_relative_error"), "Eq61_Eq63_right_relative_error": modal.metadata.get("modal_eq63_right_rate_relative_error")}
            active_changes.extend((midpoint_change, average_change, *operator_changes))
        electron_midpoint_change = _relative_scalar_change(electron_midpoint, density_state.electron_midplane_density_m3, _EPS_DENSITY_M3)
        electron_average_change = _relative_scalar_change(electron_average, density_state.electron_confined_volume_average_density_m3, _EPS_DENSITY_M3)
        electron_collision_change = _relative_scalar_change(electron_collision, density_state.electron_confined_volume_average_density_m3, _EPS_DENSITY_M3)
        barrier = float(system.shared_current_balance.electron_balance.barrier_energy_J)
        barrier_change = float("inf") if previous_barrier is None else _relative_scalar_change(previous_barrier, barrier)
        active_changes.extend((electron_midpoint_change, electron_average_change, electron_collision_change, barrier_change))
        if cross_species_required:
            active_changes.append(cross_species_maximum_change)
        relative_error = float(max(active_changes))
        identity_tolerance = max(128.0 * np.finfo(float).eps, min(qualification_tolerance, 1.0e-10))
        inventory_tolerance = float(numerics.phi_z_relative_tolerance)
        identity_valid = density_state.exact_midplane_quasineutrality_relative_error <= identity_tolerance
        inventory_valid = density_state.confined_inventory_relative_error <= inventory_tolerance
        species_iteration_converged = {key: _eq59_iteration_fixed_point_converged(system.state_for(key).modal_result) for key in confined_ids}
        species_iteration_solved = bool(all(species_iteration_converged.values()))
        species_solved = all(bool(system.state_for(key).modal_result.metadata.get("modal_rosenbluth_converged", False)) for key in confined_ids)
        cross_species_iteration_fixed_point_converged = bool(not cross_species_required or (cross_species_active_this_iteration and cross_species_pair_sets_match and cross_species_maximum_change <= iteration_tolerance))
        cross_species_qualification_fixed_point_converged = bool(not cross_species_required or (cross_species_active_this_iteration and cross_species_pair_sets_match and cross_species_maximum_change <= qualification_tolerance))
        iteration_fixed_point_candidate = bool(relative_error <= iteration_tolerance and identity_valid and inventory_valid and species_iteration_solved and cross_species_iteration_fixed_point_converged)
        qualification_fixed_point_candidate = bool(relative_error <= qualification_tolerance and identity_valid and inventory_valid and cross_species_qualification_fixed_point_converged)
        same_source_confirmation_satisfied = bool(not require_same_source_fixed_point_confirmation or iteration >= 2)
        iteration_fixed_point_converged = bool(iteration_fixed_point_candidate and same_source_confirmation_satisfied)
        qualification_fixed_point_converged = bool(qualification_fixed_point_candidate and same_source_confirmation_satisfied)
        converged = bool(qualification_fixed_point_converged and system.shared_electrostatic_profile.converged and species_solved)
        history.append({"iteration": iteration, "status": "evaluated", "species_changes": species_changes, "species_Eq59_iteration_fixed_point_converged": species_iteration_converged, "electron_midplane_relative_change": electron_midpoint_change, "electron_volume_average_relative_change": electron_average_change, "electron_collision_density_relative_change": electron_collision_change, "wall_barrier_relative_change": barrier_change, "fast_cross_species_pair_changes": cross_species_pair_changes, "fast_cross_species_maximum_relative_change": float(cross_species_maximum_change), "fast_cross_species_pair_sets_match": bool(cross_species_pair_sets_match), "fast_cross_species_iteration_fixed_point_converged": bool(cross_species_iteration_fixed_point_converged), "fast_cross_species_qualification_fixed_point_converged": bool(cross_species_qualification_fixed_point_converged), "joint_relative_change": relative_error, "iteration_relative_tolerance": float(iteration_tolerance), "qualification_relative_tolerance": float(qualification_tolerance), "iteration_tolerance_uses_qualification_fallback": bool(iteration_tolerance_uses_qualification_fallback), "same_source_fixed_point_confirmation_required": bool(require_same_source_fixed_point_confirmation), "same_source_fixed_point_confirmation_satisfied": bool(same_source_confirmation_satisfied), "iteration_fixed_point_candidate_before_same_source_confirmation": bool(iteration_fixed_point_candidate), "qualification_fixed_point_candidate_before_same_source_confirmation": bool(qualification_fixed_point_candidate), "iteration_fixed_point_converged": bool(iteration_fixed_point_converged), "qualification_fixed_point_converged": bool(qualification_fixed_point_converged), "final_qualification_converged": bool(converged), "Eq70_profile_converged": bool(system.shared_electrostatic_profile.converged), "exact_midplane_quasineutrality_relative_error": density_state.exact_midplane_quasineutrality_relative_error, "confined_inventory_relative_error": density_state.confined_inventory_relative_error, "hard_density_identities_valid": bool(identity_valid), "confined_inventory_consistent": bool(inventory_valid), "fast_cross_species_collisions_required": bool(cross_species_required), "fast_cross_species_collisions_active": bool(cross_species_active_this_iteration), "fast_cross_species_pair_ids": sorted(state.pair_id for states in external_fast_fields_by_test_species.values() for state in states)})
        density_state = replace(density_state, scalar_convergence_history=tuple(dict(item) for item in history))
        if iteration_fixed_point_converged:
            break
        external_fast_fields_by_test_species = next_external_fast_fields
        for species_id in confined_ids:
            fast_midpoint[species_id] = _positive_relaxed(fast_midpoint[species_id], density_state.fast_deuterium_midplane_density_m3 if species_id == DEUTERON.species_id else density_state.fast_tritium_midplane_density_m3, relaxation)
            fast_average[species_id] = _positive_relaxed(fast_average[species_id], density_state.fast_deuterium_confined_volume_average_density_m3 if species_id == DEUTERON.species_id else density_state.fast_tritium_confined_volume_average_density_m3, relaxation)
            warm_starts[species_id] = system.state_for(species_id).modal_result.eq59_warm_start_state
        electron_midpoint = _positive_relaxed(electron_midpoint, density_state.electron_midplane_density_m3, relaxation)
        electron_average = _positive_relaxed(electron_average, density_state.electron_confined_volume_average_density_m3, relaxation)
        electron_collision = _positive_relaxed(electron_collision, density_state.electron_confined_volume_average_density_m3, relaxation)
        potential = np.asarray(system.shared_electrostatic_profile.potential_energy_J, dtype=float)
        previous_barrier = barrier
        previous_operators = new_operators
    if system is None or density_state is None:
        raise RuntimeError("coupled density closure did not evaluate a system state")
    if assessment_requested:
        if last_solver_arguments is None:
            raise RuntimeError("final numerical assessment requires an evaluated scalar closure state")
        assessment_arguments = dict(last_solver_arguments)
        assessment_arguments["_run_numerical_convergence_assessments"] = True
        assessment_arguments["_runtime_eq70_solve_context"] = f"{runtime_context_label}:final_numerical_assessment_confirmation"
        assessment_runtime_started = perf_counter()
        system = solve_fast_ion_system(requests=requests, collision_states_by_species=last_collision_states, eq59_warm_start_states_by_species=last_warm_starts, prompt_only_loss_rates_s_by_species=prompt_rates, prompt_only_midplane_power_W_by_species=prompt_powers, ion_loss_rates_s_by_component_override=ion_loss_rates_s_by_component_override, **assessment_arguments)
        assessment_runtime_s = perf_counter() - assessment_runtime_started
        record_runtime(runtime_profile, "fast_ion_system_solve", assessment_runtime_s, nested=True)
        record_runtime(runtime_profile, "final_numerical_assessment_confirmation", assessment_runtime_s, nested=True)
        if system.shared_electrostatic_profile is None or system.shared_current_balance is None:
            raise RuntimeError("final numerical assessment did not produce the shared electrostatic state")
        if cross_species_required:
            assessment_input_external_fast_fields = dict(assessment_arguments.get("_external_fast_field_collision_states_by_test_species", {}) or {})
            assessment_output_external_fast_fields = build_external_fast_field_collision_states(system=system, collision_states_by_species=last_collision_states)
            assessment_cross_species_pair_changes, assessment_cross_species_maximum_change, assessment_cross_species_pair_sets_match = _external_fast_field_fixed_point_change(
                current_by_test_species=assessment_input_external_fast_fields,
                target_by_test_species=assessment_output_external_fast_fields,
                test_speed_grids_by_species={species_id: requests[species_id].speed_grid for species_id in confined_ids},
            )
            assessment_cross_species_fixed_point_converged = bool(assessment_cross_species_pair_sets_match and assessment_cross_species_maximum_change <= iteration_tolerance)
            assessment_cross_species_qualification_fixed_point_converged = bool(assessment_cross_species_pair_sets_match and assessment_cross_species_maximum_change <= qualification_tolerance)
        else:
            assessment_cross_species_pair_changes = {}
            assessment_cross_species_maximum_change = 0.0
            assessment_cross_species_pair_sets_match = True
            assessment_cross_species_fixed_point_converged = True
            assessment_cross_species_qualification_fixed_point_converged = True
        if history:
            history[-1] = {
                **history[-1],
                "final_assessment_confirmation_cross_species_pair_changes": assessment_cross_species_pair_changes,
                "final_assessment_confirmation_cross_species_maximum_relative_change": float(assessment_cross_species_maximum_change),
                "final_assessment_confirmation_cross_species_pair_sets_match": bool(assessment_cross_species_pair_sets_match),
                "final_assessment_confirmation_cross_species_fixed_point_converged": bool(assessment_cross_species_fixed_point_converged),
                "final_assessment_confirmation_cross_species_qualification_fixed_point_converged": bool(assessment_cross_species_qualification_fixed_point_converged),
            }
        density_state = _system_density_state(config=config, geometry=geometry, system=system, electron_collision_density_m3=last_electron_collision, history=history)
        density_state = replace(density_state, scalar_convergence_history=tuple(dict(item) for item in history))
        identity_tolerance = max(128.0 * np.finfo(float).eps, min(qualification_tolerance, 1.0e-10))
        identity_valid = density_state.exact_midplane_quasineutrality_relative_error <= identity_tolerance
        inventory_valid = density_state.confined_inventory_relative_error <= float(numerics.phi_z_relative_tolerance)
        species_solved = all(bool(system.state_for(key).modal_result.metadata.get("modal_rosenbluth_converged", False)) for key in confined_ids)
        qualification_fixed_point_converged = bool(relative_error <= qualification_tolerance and identity_valid and inventory_valid and assessment_cross_species_qualification_fixed_point_converged)
        converged = bool(qualification_fixed_point_converged and system.shared_electrostatic_profile.converged and species_solved)
        if history:
            history[-1] = {
                **history[-1],
                "qualification_fixed_point_converged": bool(qualification_fixed_point_converged),
                "final_qualification_converged": bool(converged),
                "Eq70_profile_converged": bool(system.shared_electrostatic_profile.converged),
            }
            density_state = replace(density_state, scalar_convergence_history=tuple(dict(item) for item in history))
  
    return system, density_state, collision_states, converged, iteration_fixed_point_converged, history, relative_error

def _eq42_profile_from_system_density(*, config: SourceModelRunConfig, geometry: GeometryStageResult, density_state: OperatingPointDensityState, numerics) -> object:
    """Build the self consistent Eq 42 density profile from the solved total positive charge density state"""
    density_cells = np.asarray(density_state.fast_deuterium_cell_density_m3, dtype=float) + np.asarray(density_state.fast_tritium_cell_density_m3, dtype=float)
    density_nodes = np.asarray(density_state.fast_deuterium_node_density_m3, dtype=float) + np.asarray(density_state.fast_tritium_node_density_m3, dtype=float)
    source = 'operating_point_density_state_confined_kinetic_D_plus_kinetic_T'
   
    return build_eq42_density_profile(model=EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, source=source, zeta_cells=np.asarray(geometry.zeta_centers, dtype=float), density_cells_m3=density_cells, cell_volumes_m3=np.asarray(geometry.cell_volumes_m3, dtype=float), zeta_nodes=np.asarray(density_state.zeta_nodes, dtype=float), density_nodes_m3=density_nodes, symmetry_tolerance=numerics.eq42_density_symmetry_tolerance)

def _multispecies_prompt_only_kinetic_stage(*, config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, source_audits: dict[str, object], modal_basis: _ReusableModalBasisBundle, electron_temperature_J: float, ion_temperature_J: float, collision_state_recompute_policy: str | None, collision_state_is_self_consistent_electron_temperature: bool) -> KineticStageResult:
    """Return species indexed prompt loss bookkeeping and exact zero confined arrays when no active source enters the modal domain"""
    source_ids = tuple(sorted(source_audits))
    primary_id = DEUTERON.species_id if DEUTERON.species_id in source_ids else source_ids[0]
    primary_grid = beam.speed_grid_by_species[primary_id]
    n_z = geometry.z_centers_m.size
    n_speed = primary_grid.centers_m_s.size
    n_lambda = beam.lambda_grid.centers.size
    n_pitch = beam.pitch_grid.centers.size
    zero_v_lambda = np.zeros((n_speed, n_lambda), dtype=float)
    zero_z_v_lambda = np.zeros((n_z, n_speed, n_lambda), dtype=float)
    zero_z_v_pitch = np.zeros((n_z, n_speed, n_pitch), dtype=float)
    zero_density = np.zeros(n_z, dtype=float)
    species_states = {}
    prompt_rates = {}
    prompt_powers = {}
    for species_id in source_ids:
        audit = source_audits[species_id]
        source = beam.attenuated_source_by_species[species_id]
        prompt_rate = float(audit.prompt_magnetic_loss_birth_rate_s)
        prompt_power = float(sum(item.prompt_magnetic_loss_birth_rate_s * component.spatial_source.energy_J for item, component in zip(audit.component_audits, source.component_sources, strict=True)))
        prompt_rates[species_id] = prompt_rate
        prompt_powers[species_id] = prompt_power
        species_states[species_id] = FastIonSpeciesState(species=_SPECIES_BY_ID[species_id], speed_grid=beam.speed_grid_by_species[species_id], active=False, source_particle_rate_s=float(source.total_birth_rate_s), modal_result=None, prompt_only_particle_loss_rate_s=prompt_rate, prompt_only_midplane_kinetic_power_loss_W=prompt_power, status="prompt_only_no_confined_source")
    system = FastIonSystemState(species_states=species_states, shared_basis=modal_basis.basis, status="no_confined_fbis_source", metadata={"fast_ion_operating_point_model": "no_confined_fbis_source", "source_fast_species": source_ids, "active_fast_species": ()})
    metadata = {'kinetic_model': 'egedal_modal_fbis', 'kinetic_source_state': NO_CONFINED_FBIS_SOURCE, 'collision_state_recompute_policy': str(collision_state_recompute_policy or 'computed_once_for_prompt_only_background_reference'), 'collision_state_is_self_consistent_electron_temperature': bool(collision_state_is_self_consistent_electron_temperature), 'fast_ion_operating_point_model': 'no_confined_fbis_source', 'source_fast_species': source_ids, 'active_fast_species': (), 'species_prompt_loss_birth_rate_s': prompt_rates, 'species_prompt_loss_birth_power_W': prompt_powers, 'modal_prompt_loss_birth_rate_s': float(sum(prompt_rates.values())), 'modal_prompt_loss_birth_power_W': float(sum(prompt_powers.values())), 'modal_ion_particle_loss_rate_s': float(sum(prompt_rates.values())), 'modal_ion_midplane_kinetic_power_loss_W': float(sum(prompt_powers.values())), 'modal_ion_wall_power_loss_W': None, 'modal_electron_wall_power_loss_W': None, 'modal_end_loss_power_terms_available': False, 'modal_end_loss_power_unavailability_reason': 'electrostatic_barrier_and_electron_current_balance_not_solved_for_prompt_only_source', 'kinetic_convergence_status': 'not_applicable_no_confined_fbis_source', 'kinetic_convergence_failure_reason': 'no_confined_fbis_source_operating_point_unavailable', 'kinetic_eq70_profile_present': False, 'kinetic_eq70_profile_applicable': False, 'kinetic_current_balance_applicable': False, 'kinetic_confined_arrays_are_exact_zero': True, 'kinetic_prompt_only_source_is_hard_invalid': False, 'fast_D_fast_T_cross_collisions_available': False, 'global_current_balance_claimed': False, 'physical_electron_energy_equation_active': False, 'fast_T_fusion_channels_available': True, 'species_source_projection_audits': {key: modal_source_projection_audit_metadata(source_audits[key]) for key in source_ids}}
    if collision_state_is_self_consistent_electron_temperature:
        metadata["collision_state_temperature_closure"] = "self_consistent_electron_temperature_trial_background_reference"
        metadata["plasma_temperature_closure_status"] = "self_consistent_electron_temperature_power_balance"
 
    return KineticStageResult(speed_grid=primary_grid, lambda_grid=beam.lambda_grid, pitch_grid=beam.pitch_grid, final_distribution_v_lambda=zero_v_lambda, metadata=metadata, local_speed_grid=primary_grid, local_distribution_z_v_lambda=zero_z_v_lambda, local_distribution_z_v_pitch=zero_z_v_pitch, local_density_m3=zero_density, full_device_local_distribution_z_v_pitch=None, eq59_warm_start_state=None, fast_ion_system_state=system, speed_grid_by_species={key: beam.speed_grid_by_species[key] for key in source_ids}, final_distribution_v_lambda_by_species={key: np.zeros((beam.speed_grid_by_species[key].centers_m_s.size, n_lambda), dtype=float) for key in source_ids}, local_speed_grid_by_species={key: beam.speed_grid_by_species[key] for key in source_ids}, local_distribution_z_v_lambda_by_species={key: np.zeros((n_z, beam.speed_grid_by_species[key].centers_m_s.size, n_lambda), dtype=float) for key in source_ids}, local_distribution_z_v_pitch_by_species={key: np.zeros((n_z, beam.speed_grid_by_species[key].centers_m_s.size, n_pitch), dtype=float) for key in source_ids}, local_density_m3_by_species={key: zero_density.copy() for key in source_ids}, full_device_local_distribution_z_v_pitch_by_species=None, eq59_warm_start_state_by_species=None)

def build_multispecies_modal_kinetic_stage(config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, *, electron_temperature_J: float | None=None, collision_state_recompute_policy: str | None=None, collision_state_is_self_consistent_electron_temperature: bool=False, modal_basis: _ReusableModalBasisBundle | None=None, eq59_warm_start_state_by_species: dict[str, object] | None=None, eq70_initial_potential_energy_J: np.ndarray | None=None, initial_fast_midpoint_density_m3_by_species: dict[str, float] | None=None, initial_fast_average_density_m3_by_species: dict[str, float] | None=None, initial_electron_midpoint_density_m3: float | None=None, initial_electron_confined_average_density_m3: float | None=None, initial_electron_collision_density_m3: float | None=None, initial_eq42_density_profile: object | None=None, ion_loss_rates_s_by_component_override: dict[str, float] | None=None, scalar_collision_operator_warm_start_by_species: dict[str, tuple[object, object, object, object]] | None=None, scalar_external_fast_field_warm_start_by_test_species: dict[str, tuple[object, ...]] | None=None, scalar_wall_barrier_energy_warm_start_J: float | None=None, evaluation_mode: OperatingPointEvaluationMode=OperatingPointEvaluationMode.FIXED_TEMPERATURE_FINAL) -> KineticStageResult:
    """
    Build the coupled fast D and fast T modal kinetic stage
    
    The stage audits each species source, iterates the scalar density and cross species collision fixed point, couples the solved density shape back into Eq 42, performs the requested final numerical assessment, reconstructs full device fast ion populations, and returns species indexed warm states and convergence metadata
    """
    kinetic_runtime_started = perf_counter()
    runtime_profile = new_runtime_profile()
    evaluation_mode = OperatingPointEvaluationMode.parse(evaluation_mode)
    active_ids = tuple(sorted(beam.attenuated_source_by_species))
    if not active_ids or any(key not in _SPECIES_BY_ID for key in active_ids):
        raise ValueError("coupled fast ion stage supports only deuterium and tritium")
    T_e_J = float(energy_J_from_keV(config.plasma_closure.electron_temperature_initial_guess_keV) if electron_temperature_J is None else electron_temperature_J)
    T_i_J = float(energy_J_from_keV(config.plasma_closure.background_ion_temperature_keV))
    setup_runtime_started = perf_counter()
    numerics = _modal_numerics(config)
    record_runtime(runtime_profile, "kinetic_setup", perf_counter() - setup_runtime_started)
    basis_runtime_started = perf_counter()
    bundle = build_reusable_modal_basis(config, geometry) if modal_basis is None else modal_basis
    record_runtime(runtime_profile, "initial_modal_basis", perf_counter() - basis_runtime_started)
    source_runtime_started = perf_counter()
    source_audits = {key: audit_attenuated_modal_source_projection(attenuated_source=beam.attenuated_source_by_species[key], basis=bundle.basis, volume_m3=geometry.volume_m3) for key in active_ids}
    record_runtime(runtime_profile, "source_projection_audit", perf_counter() - source_runtime_started)
    confined_ids = tuple(key for key in active_ids if source_audits[key].source_state != NO_CONFINED_FBIS_SOURCE)
    if not confined_ids:
        prompt_runtime_started = perf_counter()
        prompt_result = _multispecies_prompt_only_kinetic_stage(config=config, geometry=geometry, beam=beam, source_audits=source_audits, modal_basis=bundle, electron_temperature_J=T_e_J, ion_temperature_J=T_i_J, collision_state_recompute_policy=collision_state_recompute_policy, collision_state_is_self_consistent_electron_temperature=collision_state_is_self_consistent_electron_temperature)
        record_runtime(runtime_profile, "kinetic_prompt_only", perf_counter() - prompt_runtime_started)
        runtime_kinetic_profile = finalize_runtime_profile(runtime_profile, total_s=perf_counter() - kinetic_runtime_started)
        return replace(prompt_result, metadata={**prompt_result.metadata, "runtime_kinetic_profile": runtime_kinetic_profile})
    model = eq42_density_weighting_model(config.kinetic_electrostatic.modal_eq42_density_weighting_model)
    current_bundle = bundle
    if initial_eq42_density_profile is not None and model == EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY:
        eq42_initial_basis_runtime_started = perf_counter()
        current_bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=initial_eq42_density_profile, assess_grid_convergence=False, published_reference_lambda1=bundle.published_reference_lambda1, basis_cache=bundle.basis_cache)
        record_runtime(runtime_profile, "eq42_initial_basis", perf_counter() - eq42_initial_basis_runtime_started)
    outer_history = []
    system = None
    density_state = None
    collision_states = {}
    scalar_converged = False
    scalar_iteration_fixed_point_converged = False
    scalar_history = []
    scalar_error = float("inf")
    final_profile = current_bundle.eq42_density_profile
    final_basis_change = None
    final_profile_change = None
    outer_converged = model != EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY
    outer_iteration_fixed_point_converged = model != EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY
    outer_limit = max(int(numerics.eq42_density_basis_iterations), 1) if model == EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY else 1
    warm_starts = dict(eq59_warm_start_state_by_species or {})
    fast_midpoint = dict(initial_fast_midpoint_density_m3_by_species or {})
    fast_average = dict(initial_fast_average_density_m3_by_species or {})
    potential = eq70_initial_potential_energy_J
    electron_midpoint = initial_electron_midpoint_density_m3
    electron_average = initial_electron_confined_average_density_m3
    electron_collision = initial_electron_collision_density_m3
    cross_species_required = len(confined_ids) > 1
    scalar_operator_seed, scalar_external_seed, scalar_barrier_seed, interevaluation_seed_diagnostics = _compatible_scalar_continuation_seeds(
        confined_ids=confined_ids,
        beam=beam,
        cross_species_required=cross_species_required,
        collision_operator_seed_by_species=scalar_collision_operator_warm_start_by_species,
        external_fast_field_seed_by_test_species=scalar_external_fast_field_warm_start_by_test_species,
        wall_barrier_energy_seed_J=scalar_wall_barrier_energy_warm_start_J,
    )
    outer_continuation_refresh_count = 0
    outer_scalar_iteration_counts: list[int] = []
    interevaluation_same_source_confirmation_required = bool(interevaluation_seed_diagnostics["complete"])
    interevaluation_first_outer_scalar_history: tuple[dict[str, object], ...] = ()
    # Wrap the scalar closure in the Eq 42 density shape and physical basis fixed point
    for outer_iteration in range(1, outer_limit + 1):
        scalar_runtime_started = perf_counter()
        require_same_source_confirmation = bool(outer_iteration == 1 and interevaluation_same_source_confirmation_required)
        system, density_state, collision_states, scalar_converged, scalar_iteration_fixed_point_converged, scalar_history, scalar_error = _solve_system_scalar_closure(config=config, geometry=geometry, beam=beam, source_audits=source_audits, electron_temperature_J=T_e_J, ion_temperature_J=T_i_J, numerics=numerics, modal_basis=current_bundle, initial_eq59_warm_start_states=warm_starts, initial_fast_midpoint_density_m3=fast_midpoint, initial_fast_average_density_m3=fast_average, initial_electron_midpoint_density_m3=electron_midpoint, initial_electron_confined_average_density_m3=electron_average, initial_electron_collision_density_m3=electron_collision, initial_eq70_potential_energy_J=potential, run_numerical_convergence_assessments=False, fixed_boundary_reconstruction_policy='required', evaluation_mode=evaluation_mode, ion_loss_rates_s_by_component_override=ion_loss_rates_s_by_component_override, initial_external_fast_field_collision_states_by_test_species=scalar_external_seed, initial_previous_barrier_energy_J=scalar_barrier_seed, initial_previous_collision_operators_by_species=scalar_operator_seed, require_same_source_fixed_point_confirmation=require_same_source_confirmation, runtime_context_label=f'eq42_outer_{outer_iteration}', runtime_profile=runtime_profile)
        record_runtime(runtime_profile, "scalar_density_closure", perf_counter() - scalar_runtime_started)
        outer_scalar_iteration_counts.append(len(scalar_history))
        if outer_iteration == 1:
            interevaluation_first_outer_scalar_history = tuple(dict(item) for item in scalar_history)
        if model != EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY:
            break
        eq42_profile_runtime_started = perf_counter()
        target_profile = _eq42_profile_from_system_density(config=config, geometry=geometry, density_state=density_state, numerics=numerics)
        record_runtime(runtime_profile, "eq42_target_profile", perf_counter() - eq42_profile_runtime_started)
        eq42_target_basis_runtime_started = perf_counter()
        target_bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=target_profile, assess_grid_convergence=False, published_reference_lambda1=current_bundle.published_reference_lambda1, basis_cache=current_bundle.basis_cache)
        record_runtime(runtime_profile, "eq42_target_basis", perf_counter() - eq42_target_basis_runtime_started)
        eq42_compare_runtime_started = perf_counter()
        profile_change = compare_eq42_density_profiles(current_bundle.eq42_density_profile, target_profile)
        basis_change = _compare_eq42_bases(current_bundle.basis, target_bundle.basis)
        record_runtime(runtime_profile, "eq42_convergence_comparison", perf_counter() - eq42_compare_runtime_started)
        final_profile = target_profile
        final_profile_change = float(profile_change.volume_weighted_l2_relative_change)
        final_basis_change = basis_change
        shape_basis_converged = bool(final_profile_change <= numerics.eq42_density_shape_relative_tolerance and basis_change.maximum_eigenvalue_relative_change <= numerics.eq42_basis_eigenvalue_relative_tolerance and basis_change.minimum_eigenfunction_overlap >= numerics.eq42_basis_eigenfunction_overlap_tolerance and target_profile.symmetric_within_tolerance)
        outer_iteration_fixed_point_converged = bool(scalar_iteration_fixed_point_converged and shape_basis_converged)
        outer_converged = bool(scalar_converged and system.shared_electrostatic_profile.converged and shape_basis_converged)
        outer_history.append({"iteration": outer_iteration, "scalar_density_iteration_fixed_point_converged": bool(scalar_iteration_fixed_point_converged), "scalar_density_closure_converged": bool(scalar_converged), "Eq70_profile_converged": bool(system.shared_electrostatic_profile.converged), "density_shape_volume_weighted_l2_relative_change": final_profile_change, "basis_maximum_eigenvalue_relative_change": basis_change.maximum_eigenvalue_relative_change, "basis_minimum_eigenfunction_overlap": basis_change.minimum_eigenfunction_overlap, "shape_and_basis_fixed_point_converged": bool(shape_basis_converged), "outer_iteration_fixed_point_converged": bool(outer_iteration_fixed_point_converged), "outer_iteration_converged": outer_converged})
        if outer_iteration_fixed_point_converged or outer_iteration == outer_limit:
            break
        warm_starts = {key: system.state_for(key).modal_result.eq59_warm_start_state for key in confined_ids}
        fast_midpoint = {DEUTERON.species_id: density_state.fast_deuterium_midplane_density_m3, TRITON.species_id: density_state.fast_tritium_midplane_density_m3}
        fast_average = {DEUTERON.species_id: density_state.fast_deuterium_confined_volume_average_density_m3, TRITON.species_id: density_state.fast_tritium_confined_volume_average_density_m3}
        electron_midpoint = density_state.electron_midplane_density_m3
        electron_average = density_state.electron_confined_volume_average_density_m3
        electron_collision = density_state.electron_collision_density_m3
        potential = np.asarray(system.shared_electrostatic_profile.potential_energy_J, dtype=float)
        scalar_operator_seed = _collision_operator_seed_from_system(system, confined_ids)
        scalar_external_seed = build_external_fast_field_collision_states(system=system, collision_states_by_species=collision_states) if cross_species_required else {}
        scalar_barrier_seed = float(system.shared_current_balance.electron_balance.barrier_energy_J)
        outer_continuation_refresh_count += 1
        eq42_relaxed_basis_runtime_started = perf_counter()
        relaxed_profile = relax_eq42_density_profile(current=current_bundle.eq42_density_profile, target=target_profile, relaxation=numerics.eq42_density_shape_relaxation, model=EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, source=f"relaxed_coupled_Eq42_profile_iteration_{outer_iteration}")
        current_bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=relaxed_profile, assess_grid_convergence=False, published_reference_lambda1=current_bundle.published_reference_lambda1, basis_cache=current_bundle.basis_cache)
        record_runtime(runtime_profile, "eq42_relaxed_basis", perf_counter() - eq42_relaxed_basis_runtime_started)
    if system is None or density_state is None:
        raise RuntimeError("coupled fast ion stage did not produce a system state")
    assessment_profile = final_profile if final_profile is not None else current_bundle.eq42_density_profile
    assessment_basis_runtime_started = perf_counter()
    assessed_bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=assessment_profile, assess_grid_convergence=evaluation_mode.runs_final_qualification, published_reference_lambda1=current_bundle.published_reference_lambda1, basis_cache=current_bundle.basis_cache, bypass_cache=evaluation_mode.runs_final_qualification)
    record_runtime(runtime_profile, "eq42_assessment_basis", perf_counter() - assessment_basis_runtime_started)
    warm_starts = {key: system.state_for(key).modal_result.eq59_warm_start_state for key in confined_ids}
    fast_midpoint = {DEUTERON.species_id: density_state.fast_deuterium_midplane_density_m3, TRITON.species_id: density_state.fast_tritium_midplane_density_m3}
    fast_average = {DEUTERON.species_id: density_state.fast_deuterium_confined_volume_average_density_m3, TRITON.species_id: density_state.fast_tritium_confined_volume_average_density_m3}
    final_assessment_cross_species_seed = {}
    final_assessment_operator_seed = _collision_operator_seed_from_system(system, confined_ids)
    final_assessment_barrier_seed = float(system.shared_current_balance.electron_balance.barrier_energy_J)
    final_scalar_runtime_started = perf_counter()
    system, density_state, collision_states, scalar_converged, scalar_iteration_fixed_point_converged, scalar_history, scalar_error = _solve_system_scalar_closure(config=config, geometry=geometry, beam=beam, source_audits=source_audits, electron_temperature_J=T_e_J, ion_temperature_J=T_i_J, numerics=numerics, modal_basis=assessed_bundle, initial_eq59_warm_start_states=warm_starts, initial_fast_midpoint_density_m3=fast_midpoint, initial_fast_average_density_m3=fast_average, initial_electron_midpoint_density_m3=density_state.electron_midplane_density_m3, initial_electron_confined_average_density_m3=density_state.electron_confined_volume_average_density_m3, initial_electron_collision_density_m3=density_state.electron_collision_density_m3, initial_eq70_potential_energy_J=np.asarray(system.shared_electrostatic_profile.potential_energy_J, dtype=float), run_numerical_convergence_assessments=evaluation_mode.runs_final_qualification, fixed_boundary_reconstruction_policy='best_effort_diagnostic' if evaluation_mode.runs_final_qualification else 'required', evaluation_mode=evaluation_mode, ion_loss_rates_s_by_component_override=ion_loss_rates_s_by_component_override, initial_external_fast_field_collision_states_by_test_species=final_assessment_cross_species_seed, initial_previous_barrier_energy_J=final_assessment_barrier_seed, initial_previous_collision_operators_by_species=final_assessment_operator_seed, require_same_source_fixed_point_confirmation=False, runtime_context_label='final_scalar_closure', runtime_profile=runtime_profile)
    record_runtime(runtime_profile, "final_scalar_density_closure", perf_counter() - final_scalar_runtime_started)
    if model == EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY:
        # Reconcile the selected basis against the density state returned by the final scalar solve
        final_eq42_profile_runtime_started = perf_counter()
        final_assessed_profile = _eq42_profile_from_system_density(config=config, geometry=geometry, density_state=density_state, numerics=numerics)
        record_runtime(runtime_profile, "eq42_final_assessment_profile", perf_counter() - final_eq42_profile_runtime_started)
        final_eq42_basis_runtime_started = perf_counter()
        final_assessed_bundle = _build_modal_basis_bundle_for_eq42_profile(config=config, geometry=geometry, profile=final_assessed_profile, assess_grid_convergence=False, published_reference_lambda1=assessed_bundle.published_reference_lambda1, basis_cache=assessed_bundle.basis_cache)
        record_runtime(runtime_profile, "eq42_final_assessment_basis", perf_counter() - final_eq42_basis_runtime_started)
        final_eq42_compare_runtime_started = perf_counter()
        final_assessed_profile_change = compare_eq42_density_profiles(assessment_profile, final_assessed_profile)
        final_assessed_basis_change = _compare_eq42_bases(assessed_bundle.basis, final_assessed_bundle.basis)
        record_runtime(runtime_profile, "eq42_final_assessment_comparison", perf_counter() - final_eq42_compare_runtime_started)
        final_profile_change = float(final_assessed_profile_change.volume_weighted_l2_relative_change)
        final_basis_change = final_assessed_basis_change
        final_shape_basis_converged = bool(
            final_profile_change <= numerics.eq42_density_shape_relative_tolerance
            and final_basis_change.maximum_eigenvalue_relative_change <= numerics.eq42_basis_eigenvalue_relative_tolerance
            and final_basis_change.minimum_eigenfunction_overlap >= numerics.eq42_basis_eigenfunction_overlap_tolerance
            and final_assessed_profile.symmetric_within_tolerance
        )
        outer_iteration_fixed_point_converged = bool(scalar_iteration_fixed_point_converged and final_shape_basis_converged)
        outer_converged = bool(scalar_converged and system.shared_electrostatic_profile.converged and final_shape_basis_converged)
        if outer_history:
            previous_outer = dict(outer_history[-1])
            outer_history[-1] = {
                **previous_outer,
                "pre_final_assessment_scalar_density_closure_converged": bool(previous_outer.get("scalar_density_closure_converged", False)),
                "pre_final_assessment_outer_iteration_converged": bool(previous_outer.get("outer_iteration_converged", False)),
                "scalar_density_iteration_fixed_point_converged": bool(scalar_iteration_fixed_point_converged),
                "scalar_density_closure_converged": bool(scalar_converged),
                "Eq70_profile_converged": bool(system.shared_electrostatic_profile.converged),
                "density_shape_volume_weighted_l2_relative_change": final_profile_change,
                "basis_maximum_eigenvalue_relative_change": final_basis_change.maximum_eigenvalue_relative_change,
                "basis_minimum_eigenfunction_overlap": final_basis_change.minimum_eigenfunction_overlap,
                "shape_and_basis_fixed_point_converged": bool(final_shape_basis_converged),
                "outer_iteration_fixed_point_converged": bool(outer_iteration_fixed_point_converged),
                "outer_iteration_converged": bool(outer_converged),
                "final_assessment_reconciled": True,
            }
    if outer_converged:
        eq42_failure_reason = None
    elif not scalar_converged:
        eq42_failure_reason = "final_scalar_density_closure_not_converged"
    elif not system.shared_electrostatic_profile.converged:
        eq42_failure_reason = "final_Eq70_profile_not_converged"
    else:
        eq42_failure_reason = "final_Eq42_density_shape_or_basis_not_converged"
    eq42_metadata = _eq42_metadata(model=model, profile=assessment_profile, numerics=numerics, converged=outer_converged, iterations=len(outer_history), history=outer_history, quasineutral_profile_converged=bool(system.shared_electrostatic_profile.converged), final_profile_change=final_profile_change, final_basis_change=final_basis_change, outer_iteration_status="converged" if outer_converged else "not_converged", outer_iteration_failure_reason=eq42_failure_reason, last_valid_iteration=len(outer_history), candidate_rejection_count=0)
    full_device_by_species = {}
    updated_states = dict(system.species_states)
    full_device_metadata = {}
    for species_id in confined_ids:
        state = system.state_for(species_id)
        modal = replace(state.modal_result, metadata={**state.modal_result.metadata, **eq42_metadata})
        full_device_runtime_started = perf_counter()
        full_device, full_metadata = _full_device_fast_ion_population(geometry=geometry, beam=beam, modal=modal, collision_state=state.collision_state, quadrature_order=numerics.local_velocity_quadrature_order)
        record_runtime(runtime_profile, f"full_device_mapping_{species_id}", perf_counter() - full_device_runtime_started)
        if full_device is not None:
            full_device_by_species[species_id] = full_device
        full_device_metadata[species_id] = full_metadata
        updated_states[species_id] = replace(state, modal_result=modal, full_device_local_distribution_z_v_pitch=full_device)
    system = replace(system, species_states=updated_states, shared_basis=assessed_bundle.basis, metadata={**system.metadata, **eq42_metadata})
    postprocessing_runtime_started = perf_counter()
    primary_id = DEUTERON.species_id if DEUTERON.species_id in confined_ids else confined_ids[0]
    primary = system.state_for(primary_id)
    modal = primary.modal_result
    rounded, adjusted, raw_min = distribution_roundoff(np.asarray(modal.distribution_v_lambda, dtype=float), relative_tolerance=float(config.kinetic_electrostatic.negative_distribution_roundoff_tolerance), name=f"modal fast-{primary_id} distribution")
    speed_by_species = {key: system.state_for(key).speed_grid for key in confined_ids}
    distribution_by_species = {key: np.asarray(system.state_for(key).modal_result.distribution_v_lambda, dtype=float) for key in confined_ids}
    local_speed_by_species = {key: system.state_for(key).modal_result.local_reconstruction.local_speed_grid for key in confined_ids}
    local_lambda_by_species = {key: np.asarray(system.state_for(key).modal_result.local_reconstruction.local_distribution_z_v_lambda, dtype=float) for key in confined_ids}
    local_pitch_by_species = {key: np.asarray(system.state_for(key).modal_result.local_reconstruction.local_distribution_z_v_pitch, dtype=float) for key in confined_ids}
    local_density_by_species = {key: np.asarray(system.state_for(key).modal_result.local_reconstruction.local_density_m3, dtype=float) for key in confined_ids}
    full_device_density_by_species = {key: np.sum(np.maximum(np.asarray(full_device_by_species[key], dtype=float), 0.0) * gyrotropic_velocity_cell_volumes(local_speed_by_species[key], beam.pitch_grid)[None, :, :], axis=(1, 2)) for key in confined_ids if key in full_device_by_species}
    warm_by_species = {key: system.state_for(key).modal_result.eq59_warm_start_state for key in confined_ids}
    scalar_continuation_operator_warm_start = _collision_operator_continuation_warm_start_from_system(system, confined_ids)
    scalar_continuation_external_warm_start = build_external_fast_field_collision_states(system=system, collision_states_by_species=collision_states) if cross_species_required else {}
    scalar_continuation_barrier_warm_start = float(system.shared_current_balance.electron_balance.barrier_energy_J)
    current_scale = max(abs(float(system.shared_current_balance.ion_current_loss.total_ion_current_A)), 1.0e-300)
    current_error = abs(float(system.shared_current_balance.current_residual_A)) / current_scale
    current_tolerance = float(getattr(config.kinetic_electrostatic, "modal_current_balance_relative_tolerance", numerics.loss_convention_relative_tolerance))
    shared_profile = system.shared_electrostatic_profile
    shared_eq70_production_state = {
        "converged": bool(shared_profile.converged),
        "failure_reason": shared_profile.failure_reason,
        "electron_parent_maxwellian_n0_m3": float(shared_profile.electron_parent_maxwellian_n0_m3),
        "electron_midplane_density_m3": float(shared_profile.electron_midplane_density_m3),
        "electron_collision_density_m3": float(shared_profile.electron_collision_density_m3),
        "wall_barrier_energy_J": float(system.shared_current_balance.electron_balance.barrier_energy_J),
        "left_throat_potential_drop_energy_J": float(shared_profile.throat_potential_left_energy_J),
        "right_throat_potential_drop_energy_J": float(shared_profile.throat_potential_right_energy_J),
        "exact_midplane_potential_drop_energy_J": float(shared_profile.exact_midplane_potential_energy_J),
        "floating_nonnegative_reference_active": bool(shared_profile.floating_nonnegative_reference_active),
        "floating_zero_reference_kind": shared_profile.floating_zero_reference_kind,
        "floating_zero_reference_zeta": shared_profile.floating_zero_reference_zeta,
        "floating_zero_reference_B_tilde": shared_profile.floating_zero_reference_B_tilde,
        "floating_scalar_closure_converged": bool(shared_profile.floating_scalar_closure_converged),
        "eq68_floating_reference_iteration_fixed_point_converged": None if modal.metadata.get("modal_eq68_eq70_floating_reference_iteration_converged") is None else bool(modal.metadata.get("modal_eq68_eq70_floating_reference_iteration_converged")),
        "eq68_eq70_full_coupling_converged": None if modal.metadata.get("modal_eq68_eq70_floating_reference_coupling_converged") is None else bool(modal.metadata.get("modal_eq68_eq70_floating_reference_coupling_converged")),
        "direct_exact_midplane_root_valid": bool(shared_profile.direct_exact_midplane_root_valid),
        "direct_throat_roots_valid": bool(shared_profile.direct_throat_roots_valid),
        "direct_central_roots_valid": bool(shared_profile.direct_central_roots_valid),
        "direct_central_no_root_cell_count": int(shared_profile.direct_central_no_root_cell_count),
        "maximum_relative_quasineutrality_error": float(shared_profile.max_relative_quasineutrality_error),
        "maximum_absolute_density_residual_normalized_to_reference": float(shared_profile.maximum_absolute_density_residual_normalized_to_reference),
        "volume_integrated_absolute_particle_mismatch_fraction": float(shared_profile.volume_integrated_absolute_particle_mismatch_fraction),
    }
    source_tolerance = float(config.kinetic_electrostatic.modal_source_projection_relative_tolerance)
    source_projection_converged, source_projection_error, source_projection_errors, species_source_projection_passed = _coupled_source_projection_contract(source_audits, source_tolerance)
    geometry_boundary_valid, symmetric_end_model_applicable, end_symmetry_relative_error = _coupled_geometry_contract(config, geometry, beam)
    density_roles_separated, midpoint_role_relative_error, collision_role_relative_error = _coupled_density_roles_separated(config, modal, density_state)
    density_identity_tolerance = max(128.0 * np.finfo(float).eps, min(float(config.kinetic_electrostatic.modal_density_closure_relative_tolerance), 1.0e-10))
    density_inventory_tolerance = float(numerics.phi_z_relative_tolerance)
    exact_midplane_quasineutrality_valid = bool(density_state.exact_midplane_quasineutrality_relative_error <= density_identity_tolerance)
    confined_inventory_consistent = bool(density_state.confined_inventory_relative_error <= density_inventory_tolerance)
    species_eq59 = {key: bool(system.state_for(key).modal_result.metadata.get("modal_rosenbluth_converged", False)) for key in confined_ids}
    species_local_reconstruction_convergence = {}
    for key in confined_ids:
        reconstruction = system.state_for(key).modal_result.local_reconstruction
        diagnostics = {} if reconstruction is None else dict(reconstruction.local_speed_grid_diagnostics or {})
        species_local_reconstruction_convergence[key] = {
            "assessed": diagnostics.get("local_grid_refinement_assessed") is True,
            "converged": diagnostics.get("local_grid_refinement_converged") is True,
            "failure_reason": diagnostics.get("local_grid_refinement_failure_reason"),
            "relative_tolerance": diagnostics.get("local_grid_refinement_relative_tolerance"),
            "density_profile_relative_change": diagnostics.get("local_grid_final_density_relative_change"),
            "energy_density_profile_relative_change": diagnostics.get("local_grid_final_energy_density_relative_change"),
            "particle_inventory_relative_change": diagnostics.get("local_grid_final_particle_inventory_relative_change"),
            "energy_inventory_relative_change": diagnostics.get("local_grid_final_energy_inventory_relative_change"),
            "domain_sufficient": diagnostics.get("local_speed_domain_sufficient"),
            "overflow_particle_fraction": diagnostics.get("local_speed_overflow_particle_fraction"),
            "overflow_energy_fraction": diagnostics.get("local_speed_overflow_energy_fraction"),
            "tail_clipped": diagnostics.get("local_speed_tail_clipped"),
        }
    local_reconstruction_required = bool(evaluation_mode.runs_final_qualification)
    local_reconstruction_converged = bool(
        not local_reconstruction_required
        or (species_local_reconstruction_convergence and all(item["assessed"] is True and item["converged"] is True for item in species_local_reconstruction_convergence.values())))
    species_loss = {key: {"Eq61_left_rate_s": float(system.state_for(key).modal_result.metadata["modal_eq61_left_one_end_particle_loss_rate_s"]), "Eq63_left_rate_s": system.state_for(key).modal_result.metadata.get("modal_eq63_left_throat_crossing_rate_s"), "Eq63_left_relative_error": system.state_for(key).modal_result.metadata.get("modal_eq63_left_rate_relative_error"), "Eq63_right_relative_error": system.state_for(key).modal_result.metadata.get("modal_eq63_right_rate_relative_error")} for key in confined_ids}
    species_eq61_boundary_reconstruction_diagnostics = {key: system.state_for(key).modal_result.metadata.get("modal_eq61_boundary_reconstruction_diagnostics") for key in confined_ids}
    particle_authority_chain = _source_to_loss_particle_authority_chain(beam=beam, system=system, source_audits=source_audits, confined_ids=confined_ids, volume_m3=geometry.volume_m3)
    source_metadata = {key: modal_source_projection_audit_metadata(source_audits[key]) for key in active_ids}
    collision_policy = str(collision_state_recompute_policy or "rebuilt_from_each_outer_closure_state_and_each_fast_self_density_iterate")
    collision_metadata_by_species = {}
    for key in confined_ids:
        item = collision_parameter_metadata(system.state_for(key).collision_state)
        item["collision_state_recompute_policy"] = collision_policy
        item["collision_state_is_self_consistent_electron_temperature"] = bool(collision_state_is_self_consistent_electron_temperature)
        if collision_state_is_self_consistent_electron_temperature:
            item["collision_state_temperature_closure"] = "self_consistent_electron_temperature_power_balance"
        collision_metadata_by_species[key] = item
    full_device_includes_confined = bool(confined_ids and all(full_device_metadata[key].get("full_device_fast_ion_distribution_includes_confined") is True for key in confined_ids))
    full_device_includes_directed_lost = bool(confined_ids and all(full_device_metadata[key].get("full_device_fast_ion_distribution_includes_directed_lost") is True for key in confined_ids))
    basis_converged = bool(not evaluation_mode.runs_final_qualification or modal.metadata.get("modal_basis_convergence_assessed") is True and modal.metadata.get("modal_basis_converged") is True)
    kinetic_converged = bool(scalar_converged and outer_converged and system.shared_electrostatic_profile.converged and all(species_eq59.values()) and local_reconstruction_converged and basis_converged and source_projection_converged and geometry_boundary_valid and symmetric_end_model_applicable and density_roles_separated and exact_midplane_quasineutrality_valid and confined_inventory_consistent and current_error <= current_tolerance)
    failure_states: list[tuple[str, bool]] = [
        ("geometry_or_plasma_boundary_invalid", not geometry_boundary_valid),
        ("asymmetric_branch_kinetics_not_implemented", not symmetric_end_model_applicable),
        ("modal_source_projection_rate_not_converged", not source_projection_converged),
        (str(eq42_metadata.get("modal_eq42_outer_iteration_failure_reason") or "modal_eq42_density_basis_iteration_not_converged"), not outer_converged),
        (str(modal.metadata.get("modal_basis_convergence_failure_reason") or "modal_basis_grid_convergence_failed"), not basis_converged),
    ]
    for species_id in sorted(species_eq59):
        species_reason = system.state_for(species_id).modal_result.metadata.get("modal_rosenbluth_convergence_failure_reason")
        failure_states.append((f"fast_{species_id}_eq59:{species_reason or 'modal_eq59_nonlinear_or_speed_convergence_failed'}", not species_eq59[species_id]))
    if local_reconstruction_required:
        for species_id in sorted(species_local_reconstruction_convergence):
            item = species_local_reconstruction_convergence[species_id]
            species_passed = bool(item["assessed"] is True and item["converged"] is True)
            if item["assessed"] is not True:
                species_reason = "local_physical_speed_grid_convergence_not_assessed"
            else:
                species_reason = str(item.get("failure_reason") or "local_physical_speed_grid_convergence_failed")
            failure_states.append((f"fast_{species_id}_local_reconstruction:{species_reason}", not species_passed))
    failure_states.extend((
        (str(modal.metadata.get("modal_phi_z_failure_reason") or "modal_eq70_profile_invalid"), not system.shared_electrostatic_profile.converged),
        ("modal_collision_density_closure_not_converged", not scalar_converged),
        ("exact_midplane_quasineutrality_identity_failed", not exact_midplane_quasineutrality_valid),
        ("electron_ion_confined_inventory_identity_failed", not confined_inventory_consistent),
        ("electron_density_roles_not_consistent", not density_roles_separated),
        ("modal_current_balance_not_converged", current_error > current_tolerance),
    ))
    convergence_failures = ordered_active_failure_reasons(failure_states)
    minimum_distribution = min(float(np.min(value)) if value.size else 0.0 for value in distribution_by_species.values())
    metadata_contract = KineticMetadataContract(exact_midplane_quasineutrality_check_passed=exact_midplane_quasineutrality_valid, modal_source_projection_rate_converged=source_projection_converged, modal_source_projection_relative_error=source_projection_error, modal_source_projection_relative_tolerance=source_tolerance, kinetic_min_distribution_value=minimum_distribution, kinetic_geometry_boundary_check_passed=geometry_boundary_valid, kinetic_symmetric_end_model_applicable=symmetric_end_model_applicable, density_roles_separated=density_roles_separated, electron_ion_confined_inventory_check_passed=confined_inventory_consistent, kinetic_convergence_status="converged" if kinetic_converged else "not_converged", kinetic_convergence_failure_reasons=convergence_failures)
    record_runtime(runtime_profile, "kinetic_postprocessing", perf_counter() - postprocessing_runtime_started)
    runtime_kinetic_profile = finalize_runtime_profile(runtime_profile, total_s=perf_counter() - kinetic_runtime_started)
    final_scalar_history_record = scalar_history[-1] if scalar_history else {}
    final_cross_species_fixed_point_converged = final_scalar_history_record.get("fast_cross_species_iteration_fixed_point_converged")
    final_cross_species_maximum_relative_change = final_scalar_history_record.get("fast_cross_species_maximum_relative_change")
    final_cross_species_pair_changes = final_scalar_history_record.get("fast_cross_species_pair_changes", {})
    final_assessment_confirmation_cross_species_fixed_point_converged = final_scalar_history_record.get("final_assessment_confirmation_cross_species_fixed_point_converged")
    final_assessment_confirmation_cross_species_maximum_relative_change = final_scalar_history_record.get("final_assessment_confirmation_cross_species_maximum_relative_change")
    final_assessment_confirmation_cross_species_pair_changes = final_scalar_history_record.get("final_assessment_confirmation_cross_species_pair_changes", {})
    metadata = {**modal.metadata, 'runtime_kinetic_profile': runtime_kinetic_profile, 'runtime_interevaluation_same_source_confirmation_required': bool(interevaluation_same_source_confirmation_required), 'runtime_interevaluation_same_source_confirmation_performed': bool(not interevaluation_same_source_confirmation_required or len(interevaluation_first_outer_scalar_history) >= 2), 'runtime_interevaluation_first_outer_scalar_iterations': int(len(interevaluation_first_outer_scalar_history)), 'runtime_interevaluation_first_outer_first_iteration_fixed_point_candidate': None if not interevaluation_first_outer_scalar_history else bool(interevaluation_first_outer_scalar_history[0].get('iteration_fixed_point_candidate_before_same_source_confirmation', False)), 'runtime_interevaluation_first_outer_first_iteration_confirmation_satisfied': None if not interevaluation_first_outer_scalar_history else bool(interevaluation_first_outer_scalar_history[0].get('same_source_fixed_point_confirmation_satisfied', False)), 'runtime_interevaluation_first_outer_final_iteration_confirmation_satisfied': None if not interevaluation_first_outer_scalar_history else bool(interevaluation_first_outer_scalar_history[-1].get('same_source_fixed_point_confirmation_satisfied', False)), 'runtime_interevaluation_scalar_continuation_seed_requested': bool(interevaluation_seed_diagnostics['requested']), 'runtime_interevaluation_scalar_continuation_seed_active': bool(interevaluation_seed_diagnostics['active']), 'runtime_interevaluation_scalar_continuation_seed_complete': bool(interevaluation_seed_diagnostics['complete']), 'runtime_interevaluation_scalar_continuation_operator_species_count': int(interevaluation_seed_diagnostics['operator_species_count']), 'runtime_interevaluation_scalar_continuation_external_pair_count': int(interevaluation_seed_diagnostics['external_pair_count']), 'runtime_interevaluation_scalar_continuation_barrier_active': bool(interevaluation_seed_diagnostics['barrier_active']), 'runtime_interevaluation_scalar_continuation_rejection_reasons': tuple(interevaluation_seed_diagnostics['rejection_reasons']), 'runtime_eq42_outer_scalar_continuation_refresh_count': int(outer_continuation_refresh_count), 'runtime_eq42_outer_scalar_iteration_counts': tuple((int(value) for value in outer_scalar_iteration_counts)), 'runtime_scalar_continuation_output_operator_species_count': int(len(scalar_continuation_operator_warm_start)), 'runtime_scalar_continuation_output_external_pair_count': int(sum((len(states) for states in scalar_continuation_external_warm_start.values()))), 'runtime_scalar_continuation_output_barrier_active': bool(np.isfinite(scalar_continuation_barrier_warm_start)), 'kinetic_model': 'egedal_modal_fbis', 'fast_ion_operating_point_model': 'coupled_fast_D_fast_T_shared_electrostatic_state', 'active_fast_species': confined_ids, 'source_fast_species': active_ids, 'operating_point_evaluation_mode': evaluation_mode.value, 'state_finite': bool(all((np.all(np.isfinite(value)) for value in distribution_by_species.values()))), 'state_physically_evaluable': bool(all((np.all(value >= 0.0) for value in distribution_by_species.values()))), 'kinetic_convergence_status': 'converged' if kinetic_converged else 'not_converged', 'kinetic_numerical_convergence_passed': kinetic_converged, 'kinetic_all_inner_steps_converged': kinetic_converged, 'runtime_density_iteration_fixed_point_converged': bool(scalar_iteration_fixed_point_converged), 'runtime_density_scalar_iterations': len(scalar_history), 'runtime_density_iteration_relative_tolerance': float(_density_iteration_relative_tolerance(config)[0]), 'runtime_density_iteration_tolerance_uses_qualification_fallback': bool(_density_iteration_relative_tolerance(config)[1]), 'runtime_eq42_iteration_fixed_point_converged': bool(outer_iteration_fixed_point_converged), 'runtime_eq42_outer_iterations': len(outer_history), 'runtime_final_assessment_continuation_seeded': bool(final_assessment_cross_species_seed or final_assessment_operator_seed or np.isfinite(final_assessment_barrier_seed)), 'runtime_final_assessment_cross_species_seed_active': bool(final_assessment_cross_species_seed), 'runtime_final_assessment_cross_species_fixed_point_converged': final_cross_species_fixed_point_converged, 'runtime_final_assessment_cross_species_maximum_relative_change': final_cross_species_maximum_relative_change, 'runtime_final_assessment_cross_species_pair_changes': final_cross_species_pair_changes, 'runtime_final_assessment_confirmation_cross_species_fixed_point_converged': final_assessment_confirmation_cross_species_fixed_point_converged, 'runtime_final_assessment_confirmation_cross_species_maximum_relative_change': final_assessment_confirmation_cross_species_maximum_relative_change, 'runtime_final_assessment_confirmation_cross_species_pair_changes': final_assessment_confirmation_cross_species_pair_changes, 'runtime_final_assessment_operator_seed_species_count': int(len(final_assessment_operator_seed)), 'runtime_final_assessment_barrier_seed_active': bool(np.isfinite(final_assessment_barrier_seed)), 'runtime_final_assessment_scalar_iterations': len(scalar_history), 'runtime_final_assessment_numerical_confirmation_requested': bool(evaluation_mode.runs_final_qualification), 'runtime_final_assessment_numerical_confirmation_performed': bool(evaluation_mode.runs_final_qualification), 'modal_density_closure_model': config.kinetic_electrostatic.modal_density_closure_model, 'modal_density_closure_applicable': True, 'modal_density_closure_converged': bool(scalar_converged), 'modal_density_closure_iterations': len(scalar_history), 'modal_density_closure_relative_error': float(scalar_error), 'modal_density_closure_relative_tolerance': float(config.kinetic_electrostatic.modal_density_closure_relative_tolerance), 'modal_density_iteration_relative_tolerance': float(_density_iteration_relative_tolerance(config)[0]), 'modal_density_closure_history': scalar_history, 'modal_eq42_outer_iteration_converged': bool(outer_converged), 'modal_phi_z_converged': bool(system.shared_electrostatic_profile.converged), 'kinetic_eq70_profile_present': True, 'kinetic_eq70_profile_applicable': True, 'kinetic_current_balance_applicable': True, 'modal_current_balance_relative_error': current_error, 'kinetic_current_balance_relative_tolerance': current_tolerance, 'modal_ion_wall_power_loss_W': system.total_ion_wall_power_W, 'modal_electron_wall_power_loss_W': system.electron_wall_power_loss_W, 'modal_end_loss_power_terms_available': system.total_ion_wall_power_W is not None and system.electron_wall_power_loss_W is not None, 'electron_temperature_power_balance_scope': 'fast_D_fast_T_electron_FBIS_subsystem', 'global_plasma_power_balance_available': False, 'shared_eq42_basis': True, 'shared_eq68_wall_barrier': True, 'shared_eq70_potential': True, 'shared_eq70_production_state': shared_eq70_production_state, 'fast_D_fast_T_cross_collisions_available': False, 'global_current_balance_claimed': False, 'physical_electron_energy_equation_active': False, 'fast_T_fusion_channels_available': True, 'species_eq59_converged': species_eq59, 'kinetic_local_speed_refinement_applicable': local_reconstruction_required, 'species_local_reconstruction_convergence': species_local_reconstruction_convergence, 'all_species_local_reconstruction_converged': local_reconstruction_converged, 'species_loss_reconstruction': species_loss, 'species_eq61_boundary_reconstruction_diagnostics': species_eq61_boundary_reconstruction_diagnostics, 'nbi_particle_authority_chain_by_species': particle_authority_chain, 'species_source_projection_audits': source_metadata, 'collision_state_recompute_policy': collision_policy, 'collision_state_is_self_consistent_electron_temperature': bool(collision_state_is_self_consistent_electron_temperature), 'species_collision_metadata': collision_metadata_by_species, 'species_full_device_population_metadata': full_device_metadata, 'full_device_fast_ion_distribution_includes_confined': full_device_includes_confined, 'full_device_fast_ion_distribution_includes_directed_lost': full_device_includes_directed_lost, 'fast_ion_species_mass_kg_by_species': {key: float(_SPECIES_BY_ID[key].mass_kg) for key in confined_ids}, 'fast_ion_species_charge_number_by_species': {key: int(_SPECIES_BY_ID[key].charge_number) for key in confined_ids}, 'fast_ion_speed_grid_m_s_by_species': {key: np.asarray(speed_by_species[key].centers_m_s, dtype=float) for key in confined_ids}, 'fast_ion_speed_grid_faces_m_s_by_species': {key: np.asarray(speed_by_species[key].faces_m_s, dtype=float) for key in confined_ids}, 'fast_ion_distribution_v_lambda_by_species': distribution_by_species, 'fast_ion_local_speed_grid_m_s_by_species': {key: np.asarray(local_speed_by_species[key].centers_m_s, dtype=float) for key in confined_ids}, 'fast_ion_local_speed_grid_faces_m_s_by_species': {key: np.asarray(local_speed_by_species[key].faces_m_s, dtype=float) for key in confined_ids}, 'fast_ion_local_distribution_z_v_lambda_by_species': local_lambda_by_species, 'fast_ion_local_distribution_z_v_pitch_by_species': local_pitch_by_species, 'fast_ion_local_density_m3_by_species': local_density_by_species, 'fast_ion_full_device_distribution_z_v_pitch_by_species': full_device_by_species, 'fast_ion_full_device_density_m3_by_species': full_device_density_by_species, 'kinetic_raw_final_distribution_min': float(raw_min), 'kinetic_roundoff_negative_cells_zeroed_for_fusion': int(adjusted), **density_state.as_metadata(), **eq42_metadata, **system.metadata}
    metadata.update({"density_identity_relative_tolerance": density_identity_tolerance, "density_inventory_relative_tolerance": density_inventory_tolerance, "electron_midplane_density_role_relative_error": midpoint_role_relative_error, "electron_collision_density_role_relative_error": collision_role_relative_error, "kinetic_left_right_mirror_ratio_relative_difference": end_symmetry_relative_error, "kinetic_beam_pitch_selection_valid": beam.metadata.get("beam_pitch_selection_valid") is True, "kinetic_beam_injection_geometry_consistent": beam.metadata.get("beam_injection_geometry_consistent") is True, "species_source_projection_relative_error": source_projection_errors, "species_source_projection_passed": species_source_projection_passed, "eq68_root_residual_passed": bool(current_error <= current_tolerance), "eq68_root_relative_residual": current_error, **metadata_contract.as_metadata()})
    if assessed_bundle.basis_cache is not None:
        metadata.update(assessed_bundle.basis_cache.as_metadata())
    if collision_state_is_self_consistent_electron_temperature:
        metadata["collision_state_temperature_closure"] = "self_consistent_electron_temperature_power_balance"
        metadata["plasma_temperature_closure_status"] = "self_consistent_electron_temperature_power_balance"
   
    return KineticStageResult(
        speed_grid=primary.speed_grid, 
        lambda_grid=beam.lambda_grid, 
        pitch_grid=beam.pitch_grid, 
        final_distribution_v_lambda=rounded, 
        metadata=metadata, 
        local_speed_grid=modal.local_reconstruction.local_speed_grid, 
        local_distribution_z_v_lambda=local_lambda_by_species[primary_id], 
        local_distribution_z_v_pitch=local_pitch_by_species[primary_id], 
        local_density_m3=local_density_by_species[primary_id], 
        full_device_local_distribution_z_v_pitch=full_device_by_species.get(primary_id), 
        eq59_warm_start_state=warm_by_species[primary_id], 
        operating_point_density_state=density_state, 
        modal_basis_warm_start=assessed_bundle, 
        eq70_potential_energy_warm_start_J=np.asarray(system.shared_electrostatic_profile.potential_energy_J, dtype=float), 
        eq42_density_profile_warm_start=assessment_profile, 
        fast_ion_system_state=system, speed_grid_by_species=speed_by_species, 
        final_distribution_v_lambda_by_species=distribution_by_species, 
        local_speed_grid_by_species=local_speed_by_species, 
        local_distribution_z_v_lambda_by_species=local_lambda_by_species, 
        local_distribution_z_v_pitch_by_species=local_pitch_by_species, 
        local_density_m3_by_species=local_density_by_species, 
        full_device_local_distribution_z_v_pitch_by_species=full_device_by_species or None,
        eq59_warm_start_state_by_species=warm_by_species, 
        scalar_collision_operator_warm_start_by_species=scalar_continuation_operator_warm_start, 
        scalar_external_fast_field_warm_start_by_test_species=scalar_continuation_external_warm_start, 
        scalar_wall_barrier_energy_warm_start_J=scalar_continuation_barrier_warm_start
    )
