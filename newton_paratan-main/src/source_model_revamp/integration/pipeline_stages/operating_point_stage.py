"""
Fixed electron temperature beam and plasma operating point coupling

The stage iterates beam attenuation and the modal kinetic state until solved density profiles, beam birth rate and power, axial source shape, and component deposition are mutually consistent
"""
from __future__ import annotations
from collections.abc import Callable, Mapping
from dataclasses import dataclass, replace
from typing import Any
from time import perf_counter
import numpy as np
from source_model_revamp.fbis.collision_parameters import energy_J_from_keV
from source_model_revamp.integration.modal_stage.types import OperatingPointDensityState
from source_model_revamp.integration.config.operating_point import BeamDensityCouplingConfig
from source_model_revamp.integration.pipeline_stages.beam_stage import build_beam_stage
from source_model_revamp.integration.pipeline_stages.kinetic_stage import build_kinetic_stage
from source_model_revamp.integration.pipeline_types import BeamDensityCouplingIteration, BeamEnsembleResult, GeometryStageResult, KineticStageResult, OperatingPointStageResult, OperatingPointWarmStartState
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.plasma_closure import uses_startup_seed_only
from source_model_revamp.integration.nbi_supported_balance import StationaryParticleBalanceState, build_stationary_particle_balance
from source_model_revamp.integration.assessment import kinetic_numerical_convergence
from source_model_revamp.integration.evaluation import OperatingPointEvaluationMode, OperatingPointProgressCallback, OperatingPointProgressEvent, emit_progress
from source_model_revamp.runtime_profile import compact_runtime_summary, finalize_runtime_profile, new_runtime_profile, record_runtime

class OperatingPointWarmStartIncompatibleError(ValueError):
    """Raised when a same run operating point warm state does not match the requested geometry or grids"""

@dataclass(frozen=True)
class _BeamDiagnostics:
    """Species resolved beam totals, axial birth vectors, component deposited powers, support masks, and physical cell volumes used by the fixed point"""
    total_birth_rate_s: float
    deposited_power_W: float
    birth_profile_m3_s: np.ndarray
    birth_profile_support_mask: np.ndarray
    birth_profile_cell_volumes_m3: np.ndarray
    component_deposited_power_W: np.ndarray
    component_support_mask: np.ndarray
    total_birth_rate_s_by_species: Mapping[str, float]
    deposited_power_W_by_species: Mapping[str, float]

@dataclass(frozen=True)
class BeamSourceProfileChangeMetrics:
    """Support aware pointwise, weighted L2, and absolute reference change metrics for one beam source vector"""
    pointwise_relative_change: float
    weighted_L2_relative_change: float
    absolute_reference_change: float
    support_mask: np.ndarray
    support_count: int
    floor: float
    reference: float
    weights: np.ndarray

@dataclass(frozen=True)
class DensityProfileChangeMetrics:
    """Support aware pointwise, volume L2, and absolute reference density change metrics"""
    pointwise_relative_change: float
    volume_L2_relative_change: float
    absolute_reference_change: float
    support_mask: np.ndarray
    support_count: int
    density_floor_m3: float
    reference_density_m3: float
    physical_cell_volumes_m3: np.ndarray

@dataclass(frozen=True)
class _OperatingPointSolveStatus:
    """Independent finite state, kinetic convergence, and end loss availability status for one operating point evaluation"""
    state_finite: bool
    state_physically_evaluable: bool
    kinetic_numerical_convergence_passed: bool
    end_loss_power_terms_available: bool
    convergence_check_names: tuple[str, ...]

def _relative_change(current: np.ndarray | float, previous: np.ndarray | float) -> float:
    """Return the largest symmetric relative change with a finite zero scale"""
    a = np.asarray(current, dtype=float)
    b = np.asarray(previous, dtype=float)
   
    if a.shape != b.shape:
        raise ValueError("coupled quantities must have matching shapes")
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        raise ValueError("coupled quantities must be finite")
  
    scale = max(float(np.max(np.abs(a))) if a.size else 0.0, float(np.max(np.abs(b))) if b.size else 0.0, 1.0)

    return float(np.max(np.abs(a - b)) / scale) if a.size else 0.0

def _relative_mapping_change(current: Mapping[str, float], previous: Mapping[str, float]) -> float:
    """Return the largest relative change across matching scalar mapping entries"""
    keys = tuple(sorted(set(current) | set(previous)))
    if not keys:
        return 0.0
  
    return float(max(_relative_change(float(current.get(key, 0.0)), float(previous.get(key, 0.0))) for key in keys))

def _supported_beam_vector_change_metrics( current: np.ndarray, previous: np.ndarray, *, current_support_mask: np.ndarray, previous_support_mask: np.ndarray, weights: np.ndarray, floor_absolute: float, floor_reference_fraction: float, absolute_reference_model: str) -> BeamSourceProfileChangeMetrics:
    """
    Evaluate beam vector convergence only on a shared physical support
    
    The metric combines a floor protected pointwise change, weighted L2 change, and an absolute reference change
    """
    a = np.asarray(current, dtype=float)
    b = np.asarray(previous, dtype=float)
    current_support = np.asarray(current_support_mask, dtype=bool)
    previous_support = np.asarray(previous_support_mask, dtype=bool)
    metric_weights = np.asarray(weights, dtype=float)
   
    if not ( a.ndim == 1 and a.shape == b.shape == current_support.shape == previous_support.shape == metric_weights.shape):
        raise ValueError("beam source metric inputs must be matching one dimensional arrays")
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        raise ValueError("beam source metric values must be finite")
    if np.any(a < 0.0) or np.any(b < 0.0):
        raise ValueError("beam source metric values must be nonnegative")
    if np.any(~np.isfinite(metric_weights)) or np.any(metric_weights <= 0.0):
        raise ValueError("beam source metric weights must be finite and positive")
   
    mismatch = current_support ^ previous_support
 
    if np.any(mismatch):
        indices = np.flatnonzero(mismatch).tolist()
        raise ValueError(f"beam source metric support changed at indices {indices}")
    
    support = current_support
    support_count = int(np.count_nonzero(support))
  
    if support_count == 0:
        raise ValueError("beam source metric support is empty")
    if absolute_reference_model == "maximum_supported_amplitude":
        reference = max( float(np.max(np.abs(a[support]))), float(np.max(np.abs(b[support]))))
    elif absolute_reference_model == "total_supported_physical_value":
        reference = max( float(np.sum(np.abs(a[support]))), float(np.sum(np.abs(b[support]))))
    else:
        raise ValueError("unsupported beam source absolute reference model")   
    floor = max(float(floor_absolute), float(floor_reference_fraction) * reference)
   
    if not np.isfinite(floor) or floor < 0.0:
        raise ValueError("beam source metric floor must be finite and nonnegative")
    difference = np.abs(a[support] - b[support])
  
    if reference == 0.0 and np.all(difference == 0.0):
        pointwise = 0.0
        weighted_L2 = 0.0
        absolute_reference = 0.0
    else:
        local_scale = np.maximum( np.maximum(np.abs(a[support]), np.abs(b[support])), floor)
        pointwise = float(np.max(difference / local_scale))
        supported_weights = metric_weights[support]
        numerator = float(np.sum(supported_weights * difference**2))
        denominator = float(np.sum(supported_weights * local_scale**2))
        weighted_L2 = float(np.sqrt(numerator / denominator))
        absolute_reference = ( float(np.max(difference) / reference) if reference > 0.0 else float("inf"))

    return BeamSourceProfileChangeMetrics(
        pointwise_relative_change=pointwise,
        weighted_L2_relative_change=weighted_L2,
        absolute_reference_change=absolute_reference,
        support_mask=support.copy(),
        support_count=support_count,
        floor=floor,
        reference=reference,
        weights=metric_weights.copy(),
    )

def beam_birth_profile_change_metrics( current: np.ndarray, previous: np.ndarray, *, current_support_mask: np.ndarray, previous_support_mask: np.ndarray, physical_cell_volumes_m3: np.ndarray, floor_absolute_m3_s: float, floor_reference_fraction: float) -> BeamSourceProfileChangeMetrics:
    """Evaluate axial beam birth profile convergence with physical cell volume weighting"""
    return _supported_beam_vector_change_metrics(
        current,
        previous,
        current_support_mask=current_support_mask,
        previous_support_mask=previous_support_mask,
        weights=physical_cell_volumes_m3,
        floor_absolute=floor_absolute_m3_s,
        floor_reference_fraction=floor_reference_fraction,
        absolute_reference_model="maximum_supported_amplitude",
    )

def beam_component_deposition_change_metrics( current: np.ndarray, previous: np.ndarray, *, current_support_mask: np.ndarray, previous_support_mask: np.ndarray, floor_absolute_W: float, floor_reference_fraction: float) -> BeamSourceProfileChangeMetrics:
    """Evaluate beam component deposited power convergence on the shared component support"""
    values = np.asarray(current, dtype=float)
    return _supported_beam_vector_change_metrics(
        values,
        previous,
        current_support_mask=current_support_mask,
        previous_support_mask=previous_support_mask,
        weights=np.ones(values.shape, dtype=float),
        floor_absolute=floor_absolute_W,
        floor_reference_fraction=floor_reference_fraction,
        absolute_reference_model="total_supported_physical_value",
    )

def density_profile_change_metrics( current: np.ndarray, previous: np.ndarray, *, current_support_mask: np.ndarray, previous_support_mask: np.ndarray, physical_cell_volumes_m3: np.ndarray, density_floor_absolute_m3: float, density_floor_reference_fraction: float, physical_zero_mask: np.ndarray | None = None) -> DensityProfileChangeMetrics:
    """
    Evaluate confined density profile convergence on the shared physical support
    
    An optional physical zero mask removes cells that are structurally outside the represented density domain
    """
    a = np.asarray(current, dtype=float)
    b = np.asarray(previous, dtype=float)
    current_support = np.asarray(current_support_mask, dtype=bool)
    previous_support = np.asarray(previous_support_mask, dtype=bool)
    volumes = np.asarray(physical_cell_volumes_m3, dtype=float)
  
    if not ( a.shape == b.shape == current_support.shape == previous_support.shape == volumes.shape):
        raise ValueError("density comparison inputs must have matching shapes")
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        raise ValueError("density comparison values must be finite")
    if np.any(a < 0.0) or np.any(b < 0.0):
        raise ValueError("density comparison values must be nonnegative")
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("density comparison volumes must be finite and positive")
    support_mismatch = previous_support ^ current_support
   
    if np.any(support_mismatch):
        indices = np.flatnonzero(support_mismatch).tolist()
        raise ValueError( "density comparison support is inconsistent on physical cells " f"{indices}")
    if physical_zero_mask is None:
        physical_zero = np.zeros(a.shape, dtype=bool)
    else:
        physical_zero = np.asarray(physical_zero_mask, dtype=bool)
        if physical_zero.shape != a.shape:
            raise ValueError("physical zero mask must match the density profile")
        if np.any(physical_zero & ~current_support):
            raise ValueError("physical zero classification must lie inside profile support")
   
    support = current_support & ~physical_zero
    support_count = int(np.count_nonzero(support))
  
    if support_count == 0:
        if not np.any(previous_support):
            raise ValueError("density comparison has no physically represented support")
        reference = float(np.max(np.abs(b[previous_support])))
    else:
        reference = max( float(np.max(np.abs(a[support]))), float(np.max(np.abs(b[support]))),)
  
    floor = max( float(density_floor_absolute_m3), float(density_floor_reference_fraction) * reference,)
  
    if not np.isfinite(floor) or floor < 0.0:
        raise ValueError("density comparison floor must be finite and nonnegative")
    zero_reference_profile = bool( support_count > 0 and reference == 0.0 and np.all(a[support] == 0.0) and np.all(b[support] == 0.0))
   
    if support_count == 0 or zero_reference_profile:
        pointwise = 0.0
        volume_L2 = 0.0
        absolute_reference = 0.0
    else:
        difference = np.abs(a[support] - b[support])
        local_scale = np.maximum( np.maximum(np.abs(a[support]), np.abs(b[support])), floor)
        pointwise = float(np.max(difference / local_scale))
        support_volumes = volumes[support]
        numerator = float(np.sum(support_volumes * difference**2))
        denominator = float(np.sum(support_volumes * local_scale**2))
        volume_L2 = float(np.sqrt(numerator / denominator))
        absolute_reference = float(np.max(difference) / reference)

    return DensityProfileChangeMetrics(
        pointwise_relative_change=pointwise,
        volume_L2_relative_change=volume_L2,
        absolute_reference_change=absolute_reference,
        support_mask=support.copy(),
        support_count=support_count,
        density_floor_m3=floor,
        reference_density_m3=reference,
        physical_cell_volumes_m3=volumes.copy(),
    )

def _initial_background_target( config: SourceModelRunConfig, geometry: GeometryStageResult) -> tuple[np.ndarray, np.ndarray, bool, str]:
    """
    Build the electron stopping target used to initialize the beam and density fixed point
    
    Startup seed closure uses the configured positive charge seed only as an initialization target
    """
    z = np.asarray( getattr(geometry, "z_centers_m", getattr(geometry, "zeta_centers", ())), dtype=float)
    background = getattr(geometry, "background_profiles", None)
   
    if background is None:
        density = np.full( z.shape, float( config.plasma_closure.background_deuterium_midplane_density_m3 + config.plasma_closure.background_tritium_midplane_density_m3), dtype=float)
        support = np.ones(z.shape, dtype=bool)
    else:
        density = np.asarray(background.positive_charge_density_m3_at(z), dtype=float)
        support = np.asarray(background.support_domain.contains(z), dtype=bool)
  
    if density.shape != z.shape or support.shape != z.shape:
        raise ValueError("startup target must match the geometry axial grid")
    if np.any(~np.isfinite(density)) or np.any(density < 0.0):
        raise ValueError("startup target must be finite and nonnegative")
   
    density = np.where(support, density, 0.0)
   
    if uses_startup_seed_only(config.plasma_closure):
        return density, support, False, "magnetically_mapped_startup_D_plus_T_quasineutral"
   
    guess_supplied = bool(getattr(config.plasma_closure, "electron_density_initial_guess_supplied", False))
    maintained_midpoint = float( config.plasma_closure.background_deuterium_midplane_density_m3 + config.plasma_closure.background_tritium_midplane_density_m3)
   
    if guess_supplied:
        guess = float(config.plasma_closure.electron_density_initial_guess_m3)
        if not np.isfinite(guess) or guess <= 0.0:
            raise ValueError("electron density initial guess must be finite and positive")
        if maintained_midpoint > 0.0:
            density = np.where(support, density * guess / maintained_midpoint, 0.0)
        else:
            density = np.where(support, guess, 0.0)
        guess_source = "plasma_closure_electrons_density_initial_guess_m3"
    else:
        fast_seed = max(1.0e6, 1.0e-12 * max(maintained_midpoint, 1.0e6))
        if maintained_midpoint > 0.0:
            density = np.where( support, density * (maintained_midpoint + fast_seed) / maintained_midpoint, 0.0)
        else:
            density = np.where(support, fast_seed, 0.0)
        guess_source = "startup_positive_charge_plus_small_kinetic_seed"

    return density, support, guess_supplied, guess_source

def _warm_start_target(geometry: GeometryStageResult, warm_start_state: OperatingPointWarmStartState) -> tuple[np.ndarray, np.ndarray]:
    """Validate a same run warm state against the current confined geometry and return its target density profile and support"""
    target = np.asarray(warm_start_state.target_density_profile_m3, dtype=float)
    defined = np.asarray(warm_start_state.target_density_defined_mask, dtype=bool)
    expected_shape = np.asarray(geometry.zeta_centers).shape
   
    if target.shape != expected_shape or defined.shape != expected_shape:
        raise OperatingPointWarmStartIncompatibleError("operating point warm start target must match the geometry axial grid")
    if np.any(~np.isfinite(target)) or np.any(target < 0.0):
        raise OperatingPointWarmStartIncompatibleError("operating point warm start target must be finite and nonnegative")

    return target.copy(), defined.copy()

def _beam_diagnostics(beam: BeamEnsembleResult) -> _BeamDiagnostics:
    """Aggregate beam observables by species and return the vectors required by the fixed point convergence tests"""
    total_rate = float(beam.total_fast_birth_rate_s)
    total_power = float(beam.total_deposited_birth_power_W)
    rates_by_species = {str(key): float(value) for key, value in beam.total_fast_birth_rate_by_species.items()}
    powers_by_species: dict[str, float] = {}
    profiles_by_species: dict[str, np.ndarray] = {}
    supports_by_species: dict[str, np.ndarray] = {}
    component_power: list[float] = []
    component_support: list[bool] = []
    volumes = None
  
    for beam_id in sorted(beam.beam_results_by_id):
        result = beam.beam_results_by_id[beam_id]
        species_id = str(result.projectile_species)
        powers_by_species[species_id] = powers_by_species.get(species_id, 0.0) + float(result.total_deposited_birth_power_W)
        profile = np.asarray(result.metadata.get("beam_source_axial_physical_birth_rate_density_m3_s"), dtype=float)
        support = np.asarray(result.metadata.get("beam_birth_profile_support_mask", np.ones(profile.shape, dtype=bool)), dtype=bool)
        result_volumes = np.asarray(result.metadata.get("beam_birth_profile_physical_cell_volumes_m3", np.ones(profile.shape, dtype=float)), dtype=float)
     
        if profile.ndim != 1 or support.shape != profile.shape or result_volumes.shape != profile.shape:
            raise ValueError("beam species diagnostics must share one axial profile grid")
        if volumes is None:
            volumes = result_volumes
        elif not np.array_equal(volumes, result_volumes):
            raise ValueError("beam birth profile physical cell volumes changed between beams")
       
        profiles_by_species[species_id] = profiles_by_species.get(species_id, np.zeros_like(profile)) + profile
        supports_by_species[species_id] = supports_by_species.get(species_id, np.zeros_like(support)) | support
        values = np.asarray(result.metadata.get("beam_component_deposited_birth_power_W", [float(component.deposited_birth_power_W) for component in result.attenuated_source.component_sources]), dtype=float)
        configured = np.asarray(result.metadata.get("beam_component_support_mask", np.ones(values.shape, dtype=bool)), dtype=bool)
        component_power.extend(float(value) for value in values)
        component_support.extend(bool(value) for value in configured)
   
    species_ids = tuple(sorted(rates_by_species))
  
    if not species_ids or volumes is None:
        raise ValueError("beam fixed point diagnostics require at least one species source")
  
    birth_profile = np.concatenate([profiles_by_species[key] for key in species_ids])
    support = np.concatenate([supports_by_species[key] for key in species_ids])
    profile_volumes = np.concatenate([volumes for _ in species_ids])
    component_power_array = np.asarray(component_power, dtype=float)
    component_support_array = np.asarray(component_support, dtype=bool)
   
    if (
        not np.isfinite(total_rate)
        or total_rate < 0.0
        or not np.isfinite(total_power)
        or total_power < 0.0
        or any(not np.isfinite(value) or value < 0.0 for value in rates_by_species.values())
        or any(not np.isfinite(value) or value < 0.0 for value in powers_by_species.values())
        or np.any(~np.isfinite(birth_profile))
        or np.any(birth_profile < 0.0)
        or np.any(~np.isfinite(profile_volumes))
        or np.any(profile_volumes <= 0.0)
        or component_power_array.ndim != 1
        or component_power_array.size == 0
        or component_support_array.shape != component_power_array.shape
        or np.any(~np.isfinite(component_power_array))
        or np.any(component_power_array < 0.0)
        or not np.any(support)
        or not np.any(component_support_array)):
        raise ValueError("beam fixed point diagnostics must be finite and nonnegative")
   
    return _BeamDiagnostics(total_birth_rate_s=total_rate, deposited_power_W=total_power, birth_profile_m3_s=birth_profile, birth_profile_support_mask=support, birth_profile_cell_volumes_m3=profile_volumes, component_deposited_power_W=component_power_array, component_support_mask=component_support_array, total_birth_rate_s_by_species=rates_by_species, deposited_power_W_by_species=powers_by_species)

def _beam_source_change_metrics(current: _BeamDiagnostics, previous: _BeamDiagnostics, controls: BeamDensityCouplingConfig) -> tuple[BeamSourceProfileChangeMetrics, BeamSourceProfileChangeMetrics]:
    """Evaluate axial birth profile and component deposited power changes between consecutive beam states"""
    birth_metrics = beam_birth_profile_change_metrics(
        current.birth_profile_m3_s,
        previous.birth_profile_m3_s,
        current_support_mask=current.birth_profile_support_mask,
        previous_support_mask=previous.birth_profile_support_mask,
        physical_cell_volumes_m3=current.birth_profile_cell_volumes_m3,
        floor_absolute_m3_s=controls.birth_profile_floor_absolute_m3_s,
        floor_reference_fraction=controls.birth_profile_floor_reference_fraction,
    )
   
    if not np.array_equal( current.birth_profile_cell_volumes_m3, previous.birth_profile_cell_volumes_m3):
        raise ValueError("beam birth profile physical cell volumes changed")
   
    component_metrics = beam_component_deposition_change_metrics(
        current.component_deposited_power_W,
        previous.component_deposited_power_W,
        current_support_mask=current.component_support_mask,
        previous_support_mask=previous.component_support_mask,
        floor_absolute_W=controls.component_power_floor_absolute_W,
        floor_reference_fraction=controls.component_power_floor_reference_fraction,
    )

    return birth_metrics, component_metrics

def _solved_density_state( kinetic: KineticStageResult, geometry: GeometryStageResult) -> OperatingPointDensityState:
    """Return the typed solved confined density state from the kinetic result and validate its geometry"""
    state = kinetic.operating_point_density_state
  
    if state is None:
        raise ValueError("kinetic stage did not return an operating point density state")
    
    electron = np.asarray(state.electron_cell_density_m3, dtype=float)
    support = np.asarray(state.profile_support_mask, dtype=bool)
    expected_shape = np.asarray(geometry.zeta_centers).shape
   
    if electron.shape != expected_shape or support.shape != expected_shape:
        raise ValueError("solved electron profile must match the geometry axial grid")
    if np.any(~np.isfinite(electron)) or np.any(electron < 0.0):
        raise ValueError("solved electron profile must be finite and nonnegative")
    if state.electron_profile_scope != "confined_throat_to_throat":
        raise ValueError("solved electron profile has an unsupported scope")

    return state

def _kinetic_species_density_change_metrics(current: KineticStageResult, previous: KineticStageResult | None, *, support_mask: np.ndarray, physical_cell_volumes_m3: np.ndarray, controls: BeamDensityCouplingConfig) -> tuple[DensityProfileChangeMetrics | None, dict[str, float] | None, dict[str, dict[str, float]] | None]:
    """Return aggregate and species resolved fast ion density profile changes between kinetic states"""
    current_map = dict(current.local_density_m3_by_species or {})
    previous_map = {} if previous is None else dict(previous.local_density_m3_by_species or {})
   
    if not current_map or not previous_map:
        return None, None, None
  
    support = np.asarray(support_mask, dtype=bool)
    volumes = np.asarray(physical_cell_volumes_m3, dtype=float)
    metrics_by_species: dict[str, DensityProfileChangeMetrics] = {}
  
    for species_id in sorted(set(current_map) | set(previous_map)):
        current_density = np.asarray(current_map.get(species_id, np.zeros(support.shape, dtype=float)), dtype=float)
        previous_density = np.asarray(previous_map.get(species_id, np.zeros(support.shape, dtype=float)), dtype=float)
        if current_density.shape != support.shape or previous_density.shape != support.shape:
            raise ValueError("kinetic species density profile must match the geometry axial grid")
       
        metrics_by_species[species_id] = density_profile_change_metrics(
            current_density,
            previous_density,
            current_support_mask=support,
            previous_support_mask=support,
            physical_cell_volumes_m3=volumes,
            density_floor_absolute_m3=controls.density_floor_absolute_m3,
            density_floor_reference_fraction=controls.density_floor_reference_fraction,
        )
    pointwise = max(value.pointwise_relative_change for value in metrics_by_species.values())
    volume_L2 = max(value.volume_L2_relative_change for value in metrics_by_species.values())
    absolute = max(value.absolute_reference_change for value in metrics_by_species.values())
    representative = next(iter(metrics_by_species.values()))
    aggregate = DensityProfileChangeMetrics(
        pointwise_relative_change=pointwise,
        volume_L2_relative_change=volume_L2,
        absolute_reference_change=absolute,
        support_mask=representative.support_mask.copy(),
        support_count=representative.support_count,
        density_floor_m3=max(value.density_floor_m3 for value in metrics_by_species.values()),
        reference_density_m3=max(value.reference_density_m3 for value in metrics_by_species.values()),
        physical_cell_volumes_m3=representative.physical_cell_volumes_m3.copy(),
    )
    volume_L2_by_species = {key: value.volume_L2_relative_change for key, value in metrics_by_species.items()}
    full_metrics_by_species = {key: {
            "pointwise_relative_change": value.pointwise_relative_change,
            "volume_L2_relative_change": value.volume_L2_relative_change,
            "absolute_reference_change": value.absolute_reference_change,
            "density_floor_m3": value.density_floor_m3,
            "reference_density_m3": value.reference_density_m3,
        } for key, value in metrics_by_species.items()
    }
   
    return aggregate, volume_L2_by_species, full_metrics_by_species


def _kinetic_check_status( kinetic: KineticStageResult, density_state: OperatingPointDensityState) -> _OperatingPointSolveStatus:
    """Return whether one kinetic state is finite, evaluable, numerically converged, and has usable end loss terms"""
    metadata = kinetic.metadata
    arrays = (np.asarray(kinetic.final_distribution_v_lambda, dtype=float), np.asarray(density_state.fast_deuterium_cell_density_m3, dtype=float), np.asarray(density_state.fast_tritium_cell_density_m3, dtype=float), np.asarray(density_state.electron_cell_density_m3, dtype=float), np.asarray(density_state.total_positive_charge_cell_density_m3, dtype=float), np.asarray(density_state.cell_volumes_m3, dtype=float))
    scalars = np.asarray([density_state.fast_deuterium_midplane_density_m3, density_state.fast_tritium_midplane_density_m3, density_state.electron_midplane_density_m3, density_state.fast_deuterium_confined_volume_average_density_m3, density_state.fast_tritium_confined_volume_average_density_m3, density_state.electron_confined_volume_average_density_m3, density_state.electron_collision_density_m3, density_state.electron_parent_maxwellian_n0_m3], dtype=float,)
    potential = kinetic.eq70_potential_energy_warm_start_J
    state_finite = bool( all(np.all(np.isfinite(values)) for values in arrays) and np.all(np.isfinite(scalars)) and (potential is None or np.all(np.isfinite(np.asarray(potential, dtype=float)))) and metadata.get("kinetic_eq59_state_finite", True) is True)
    state_physically_evaluable = bool( state_finite and all(np.all(values >= 0.0) for values in arrays[1:-1]) and np.all(scalars >= 0.0) and np.all(arrays[-1] > 0.0) and density_state.electron_profile_scope == "confined_throat_to_throat" and not density_state.expander_electron_profile_available and metadata.get("kinetic_hard_model_applicability_passed", True) is True)
    numerical, numerical_records = kinetic_numerical_convergence(metadata)

    return _OperatingPointSolveStatus(
        state_finite=state_finite,
        state_physically_evaluable=state_physically_evaluable,
        kinetic_numerical_convergence_passed=numerical,
        end_loss_power_terms_available=metadata.get("modal_end_loss_power_terms_available") is True,
        convergence_check_names=numerical_records,
    )

def _all_coupling_metrics_pass(controls: BeamDensityCouplingConfig, *, profile_pointwise_change: float, profile_volume_L2_change: float, profile_absolute_reference_change: float, species_density_metrics: DensityProfileChangeMetrics | None, birth_rate_change: float | None, deposited_power_change: float | None, birth_profile_metrics: BeamSourceProfileChangeMetrics | None, component_metrics: BeamSourceProfileChangeMetrics | None) -> bool:
    """Return whether every required beam and density coupling metric satisfies its configured tolerance"""
    if species_density_metrics is None or birth_rate_change is None or deposited_power_change is None or birth_profile_metrics is None or component_metrics is None:
        return False

    return bool(
        profile_pointwise_change <= controls.electron_profile_relative_tolerance
        and profile_volume_L2_change <= controls.profile_volume_L2_tolerance
        and profile_absolute_reference_change
        <= controls.profile_absolute_reference_tolerance
        and species_density_metrics.pointwise_relative_change <= controls.electron_profile_relative_tolerance
        and species_density_metrics.volume_L2_relative_change <= controls.profile_volume_L2_tolerance
        and species_density_metrics.absolute_reference_change <= controls.profile_absolute_reference_tolerance
        and float(birth_rate_change) <= controls.total_birth_rate_relative_tolerance
        and float(deposited_power_change) <= controls.deposited_power_relative_tolerance
        and birth_profile_metrics.pointwise_relative_change
        <= controls.axial_birth_profile_relative_tolerance
        and birth_profile_metrics.weighted_L2_relative_change
        <= controls.axial_birth_profile_relative_tolerance
        and birth_profile_metrics.absolute_reference_change
        <= controls.axial_birth_profile_relative_tolerance
        and component_metrics.pointwise_relative_change
        <= controls.energy_component_deposition_relative_tolerance
        and component_metrics.weighted_L2_relative_change
        <= controls.energy_component_deposition_relative_tolerance
        and component_metrics.absolute_reference_change
        <= controls.energy_component_deposition_relative_tolerance
    )

def _gate_record(value: float | None, tolerance: float) -> dict[str, object]:
    """Return one convergence gate value, tolerance, availability, and pass state"""
    numeric = None if value is None else float(value)
    limit = float(tolerance)
    passed = bool(numeric is not None and np.isfinite(numeric) and numeric <= limit)
    ratio = None if numeric is None else (0.0 if limit == 0.0 and numeric == 0.0 else float("inf") if limit == 0.0 else numeric / limit)
   
    return {"value": numeric, "tolerance": limit, "ratio": ratio, "passed": passed}

def _final_consistency_metric_gates(controls: BeamDensityCouplingConfig, *, profile_metrics: DensityProfileChangeMetrics, species_density_metrics: DensityProfileChangeMetrics | None, species_density_change_metrics_by_species: Mapping[str, Mapping[str, float]] | None, birth_rate_change: float | None, deposited_power_change: float | None, birth_profile_metrics: BeamSourceProfileChangeMetrics | None, component_metrics: BeamSourceProfileChangeMetrics | None) -> dict[str, object]:
    """Build the final beam and density consistency gate records after the candidate fixed point is rebuilt"""
    gates = {
        "electron_density_pointwise": _gate_record(profile_metrics.pointwise_relative_change, controls.electron_profile_relative_tolerance),
        "electron_density_volume_L2": _gate_record(profile_metrics.volume_L2_relative_change, controls.profile_volume_L2_tolerance),
        "electron_density_absolute_reference": _gate_record(profile_metrics.absolute_reference_change, controls.profile_absolute_reference_tolerance),
        "species_density_pointwise": _gate_record(None if species_density_metrics is None else species_density_metrics.pointwise_relative_change, controls.electron_profile_relative_tolerance),
        "species_density_volume_L2": _gate_record(None if species_density_metrics is None else species_density_metrics.volume_L2_relative_change, controls.profile_volume_L2_tolerance),
        "species_density_absolute_reference": _gate_record(None if species_density_metrics is None else species_density_metrics.absolute_reference_change, controls.profile_absolute_reference_tolerance),
        "beam_birth_rate": _gate_record(birth_rate_change, controls.total_birth_rate_relative_tolerance),
        "beam_deposited_power": _gate_record(deposited_power_change, controls.deposited_power_relative_tolerance),
        "beam_birth_profile_pointwise": _gate_record(None if birth_profile_metrics is None else birth_profile_metrics.pointwise_relative_change, controls.axial_birth_profile_relative_tolerance),
        "beam_birth_profile_volume_L2": _gate_record(None if birth_profile_metrics is None else birth_profile_metrics.weighted_L2_relative_change, controls.axial_birth_profile_relative_tolerance),
        "beam_birth_profile_absolute_reference": _gate_record(None if birth_profile_metrics is None else birth_profile_metrics.absolute_reference_change, controls.axial_birth_profile_relative_tolerance),
        "beam_component_pointwise": _gate_record(None if component_metrics is None else component_metrics.pointwise_relative_change, controls.energy_component_deposition_relative_tolerance),
        "beam_component_L2": _gate_record(None if component_metrics is None else component_metrics.weighted_L2_relative_change, controls.energy_component_deposition_relative_tolerance),
        "beam_component_absolute_reference": _gate_record(None if component_metrics is None else component_metrics.absolute_reference_change, controls.energy_component_deposition_relative_tolerance),
    }
    species = {}
   
    for species_id, metrics in (species_density_change_metrics_by_species or {}).items():
        species[str(species_id)] = {
            "pointwise": _gate_record(float(metrics["pointwise_relative_change"]), controls.electron_profile_relative_tolerance),
            "volume_L2": _gate_record(float(metrics["volume_L2_relative_change"]), controls.profile_volume_L2_tolerance),
            "absolute_reference": _gate_record(float(metrics["absolute_reference_change"]), controls.profile_absolute_reference_tolerance),
        }
   
    return {"aggregate": gates, "by_species": species}

def _progress_gate_failure_text(gates: Mapping[str, object]) -> str:
    """Return compact text for beam and density coupling gates that have not passed"""
    labels = {
        "electron_density_pointwise": "electron profile pointwise",
        "electron_density_volume_L2": "electron profile volume L2",
        "electron_density_absolute_reference": "electron profile absolute reference",
        "species_density_pointwise": "species density pointwise",
        "species_density_volume_L2": "species density volume L2",
        "species_density_absolute_reference": "species density absolute reference",
        "beam_birth_rate": "beam birth rate",
        "beam_deposited_power": "beam deposited power",
        "beam_birth_profile_pointwise": "axial birth profile pointwise",
        "beam_birth_profile_volume_L2": "axial birth profile volume L2",
        "beam_birth_profile_absolute_reference": "axial birth profile absolute reference",
        "beam_component_pointwise": "component deposition pointwise",
        "beam_component_L2": "component deposition L2",
        "beam_component_absolute_reference": "component deposition absolute reference",
    }
    failures: list[str] = []
    aggregate = gates.get("aggregate", {})
  
    if not isinstance(aggregate, Mapping):
        return "beam convergence gates unavailable"
  
    for key, record in aggregate.items():
        if not isinstance(record, Mapping) or record.get("passed") is True:
            continue
       
        label = labels.get(str(key), str(key))
        value = record.get("value")
        tolerance = record.get("tolerance")
      
        if value is None:
            failures.append(f"{label}=unavailable")
        else:
            failures.append(f"{label}={float(value):.6g}>{float(tolerance):.6g}")
  
    return "all beam coupling gates passed" if not failures else "not passing: " + "; ".join(failures)

def _candidate_beam_progress_failure_text(controls: BeamDensityCouplingConfig, *, birth_rate_change: float, deposited_power_change: float, birth_profile_metrics: BeamSourceProfileChangeMetrics, component_metrics: BeamSourceProfileChangeMetrics) -> str:
    """Return compact text for failed candidate beam confirmation gates"""
    records = {
        "beam birth rate": _gate_record(birth_rate_change, controls.total_birth_rate_relative_tolerance),
        "beam deposited power": _gate_record(deposited_power_change, controls.deposited_power_relative_tolerance),
        "axial birth profile pointwise": _gate_record(birth_profile_metrics.pointwise_relative_change, controls.axial_birth_profile_relative_tolerance),
        "axial birth profile volume L2": _gate_record(birth_profile_metrics.weighted_L2_relative_change, controls.axial_birth_profile_relative_tolerance),
        "axial birth profile absolute reference": _gate_record(birth_profile_metrics.absolute_reference_change, controls.axial_birth_profile_relative_tolerance),
        "component deposition pointwise": _gate_record(component_metrics.pointwise_relative_change, controls.energy_component_deposition_relative_tolerance),
        "component deposition L2": _gate_record(component_metrics.weighted_L2_relative_change, controls.energy_component_deposition_relative_tolerance),
        "component deposition absolute reference": _gate_record(component_metrics.absolute_reference_change, controls.energy_component_deposition_relative_tolerance),
    }
    failures = [
        f"{label}={float(record['value']):.6g}>{float(record['tolerance']):.6g}"
        for label, record in records.items()
        if record["passed"] is not True
    ]
  
    return "candidate beam confirmation passed" if not failures else "candidate beam not passing: " + "; ".join(failures)

def build_fixed_temperature_operating_point(config: SourceModelRunConfig, geometry: GeometryStageResult, *, electron_temperature_J: float, modal_basis: Any | None = None, warm_start_state: OperatingPointWarmStartState | None = None, collision_state_recompute_policy: str = "computed_for_fixed_electron_temperature_operating_point", collision_state_is_self_consistent_electron_temperature: bool = False, controls: BeamDensityCouplingConfig | None = None, beam_builder: Callable[..., BeamEnsembleResult] = build_beam_stage, kinetic_builder: Callable[..., KineticStageResult] = build_kinetic_stage, evaluation_mode: OperatingPointEvaluationMode = OperatingPointEvaluationMode.FIXED_TEMPERATURE_FINAL, progress_callback: OperatingPointProgressCallback | None = None, temperature_trial_index: int | None = None, warm_start_source_temperature_keV: float | None = None, clock: Callable[[], float] = perf_counter) -> OperatingPointStageResult:
    """
    Solve the beam attenuation and plasma state fixed point at one electron temperature
    
    Each iteration builds beams from the current stopping target, solves the modal kinetic state, evaluates density and beam changes, and updates warm states
    A candidate fixed point is rebuilt with its own solved density target before final consistency and stationary particle balance checks are accepted
    """
    operating_runtime_started = clock()
    runtime_profile = new_runtime_profile()
    kinetic_runtime_evaluations: list[dict[str, object]] = []
    evaluation_mode = OperatingPointEvaluationMode.parse(evaluation_mode)
    temperature = float(electron_temperature_J)
   
    if not np.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("electron_temperature_J must be finite and positive")
  
    active_controls = controls or config.beam_density_coupling
    benchmark_closure = str(config.kinetic_electrostatic.modal_density_closure_model).strip().lower() == "egedal_beam_plasma_quasineutral"
    iteration_evaluation_mode = (OperatingPointEvaluationMode.TEMPERATURE_TRIAL if evaluation_mode.runs_final_qualification and not benchmark_closure else evaluation_mode)
    startup_seed_only = uses_startup_seed_only(config.plasma_closure)
   
    if warm_start_state is None:
        target_density, target_defined, initial_guess_used, initial_target_source = _initial_background_target(config, geometry)
        density_warm_start = None
        eq59_warm_start_state = None
        eq59_warm_start_state_by_species = None
        eq70_warm_start = None
        eq42_density_profile_warm_start = None
        scalar_collision_operator_warm_start_by_species = None
        scalar_external_fast_field_warm_start_by_test_species = None
        scalar_wall_barrier_energy_warm_start_J = None
        operating_point_warm_start_used = False
        kinetic_target_state = None
    else:
        target_density, target_defined = _warm_start_target(geometry, warm_start_state)
        density_warm_start = warm_start_state.density_state
        eq59_warm_start_state = warm_start_state.eq59_warm_start_state
        eq59_warm_start_state_by_species = warm_start_state.eq59_warm_start_state_by_species
        eq70_warm_start = warm_start_state.eq70_potential_energy_warm_start_J
        eq42_density_profile_warm_start = warm_start_state.eq42_density_profile_warm_start
        scalar_collision_operator_warm_start_by_species = warm_start_state.scalar_collision_operator_warm_start_by_species
        scalar_external_fast_field_warm_start_by_test_species = warm_start_state.scalar_external_fast_field_warm_start_by_test_species
        scalar_wall_barrier_energy_warm_start_J = warm_start_state.scalar_wall_barrier_energy_warm_start_J
       
        if warm_start_state.modal_basis_warm_start is not None:
            modal_basis = warm_start_state.modal_basis_warm_start
      
        initial_guess_used = False
        initial_target_source = "same_run_operating_point_warm_start"
        operating_point_warm_start_used = True
        kinetic_target_state = warm_start_state.kinetic_target_state

    def density_warm_values(state: OperatingPointDensityState | None) -> tuple[float | None, float | None, float | None, float | None, float | None, float | None, float | None, float | None]:
        """Extract scalar D, T, and electron density warm values from one solved OperatingPointDensityState"""
        if state is None:
            volumes = np.asarray(geometry.cell_volumes_m3, dtype=float)
            average = float(np.sum(target_density[target_defined] * volumes[target_defined]) / np.sum(volumes[target_defined]))
            startup_charge = float(config.plasma_closure.background_deuterium_midplane_density_m3 + config.plasma_closure.background_tritium_midplane_density_m3)
            midpoint = startup_charge if startup_seed_only else float(config.plasma_closure.electron_density_initial_guess_m3) if initial_guess_used else startup_charge + max(1.0e6, 1.0e-12 * max(startup_charge, 1.0e6))
           
            return None, None, None, None, midpoint, average, None, average
      
        return (
            state.fast_deuterium_midplane_density_m3,
            state.fast_deuterium_confined_volume_average_density_m3,
            state.fast_tritium_midplane_density_m3,
            state.fast_tritium_confined_volume_average_density_m3,
            state.electron_midplane_density_m3,
            state.electron_confined_volume_average_density_m3,
            state.electron_parent_maxwellian_n0_m3,
            state.electron_collision_density_m3,
        )

    fast_midpoint_warm_start, fast_average_warm_start, fast_T_midpoint_warm_start, fast_T_average_warm_start, electron_midpoint_warm_start, electron_average_warm_start, electron_parent_n0_warm_start, electron_collision_density_warm_start = density_warm_values(density_warm_start)
    history: list[BeamDensityCouplingIteration] = []
    previous_beam_diagnostics: _BeamDiagnostics | None = None
    last_beam: BeamEnsembleResult | None = None
    last_kinetic: KineticStageResult | None = None
    last_density_state: OperatingPointDensityState | None = None
    final_beam: BeamEnsembleResult | None = None
    final_birth_rate_change: float | None = None
    final_deposited_power_change: float | None = None
    final_birth_profile_metrics: BeamSourceProfileChangeMetrics | None = None
    final_component_metrics: BeamSourceProfileChangeMetrics | None = None
    final_species_density_metrics: DensityProfileChangeMetrics | None = None
    final_species_density_changes: dict[str, float] | None = None
    final_species_density_change_metrics_by_species: dict[str, dict[str, float]] | None = None
    particle_balance_state: StationaryParticleBalanceState | None = None
    profile_metrics: DensityProfileChangeMetrics | None = None
    solve_status: _OperatingPointSolveStatus | None = None
    failure_reason: str | None = None
    converged = False
    benchmark_bypass = False
    candidate_fixed_point_passed = False
    final_consistency_recompute_passed = False
    final_consistency_overall_qualification_passed = False
    beam_density_termination_reason: str | None = None
    final_consistency_attempt_history: list[dict[str, object]] = []
    latest_final_consistency_record: dict[str, object] | None = None
   
    # Iterate beam attenuation and the solved density state until both source and density metrics close
    for iteration in range(1, active_controls.max_iterations + 1):
        previous_kinetic_target_state = kinetic_target_state
        iteration_started = clock()
      
        emit_progress(progress_callback, OperatingPointProgressEvent(
                phase="beam_density_iteration_started",
                evaluation_mode=evaluation_mode,
                temperature_trial_index=temperature_trial_index,
                temperature_keV=float(temperature / energy_J_from_keV(1.0)),
                warm_start_source_temperature_keV=warm_start_source_temperature_keV,
                beam_iteration=iteration,
                beam_iteration_limit=active_controls.max_iterations,
            ),
        )
      
        beam_runtime_started = clock()
        beam = beam_builder(config, geometry, electron_temperature_J=temperature) if benchmark_closure else beam_builder(config, geometry, electron_temperature_J=temperature, target_density_profile_m3=target_density, target_density_defined_mask=target_defined, kinetic_target_state=kinetic_target_state)
       
        record_runtime(runtime_profile, "beam_build", clock() - beam_runtime_started)
       
        target_physical_support = np.asarray(beam.metadata.get("beam_target_density_profile_support_mask", target_defined), dtype=bool)
       
        if target_physical_support.shape != target_density.shape or not np.any(target_physical_support):
            raise ValueError("beam target physical support is invalid")
       
        beam_diagnostics = _beam_diagnostics(beam)
        kinetic_started = clock()
        kinetic = kinetic_builder(
            config,
            geometry,
            beam,
            electron_temperature_J=temperature,
            collision_state_recompute_policy=collision_state_recompute_policy,
            collision_state_is_self_consistent_electron_temperature=collision_state_is_self_consistent_electron_temperature,
            modal_basis=modal_basis,
            eq59_warm_start_state=eq59_warm_start_state,
            eq59_warm_start_state_by_species=eq59_warm_start_state_by_species,
            eq70_initial_potential_energy_J=eq70_warm_start,
            initial_representative_fast_density_m3=fast_average_warm_start,
            initial_fast_D_midpoint_density_m3=fast_midpoint_warm_start,
            initial_fast_D_confined_average_density_m3=fast_average_warm_start,
            initial_fast_T_midpoint_density_m3=fast_T_midpoint_warm_start,
            initial_fast_T_confined_average_density_m3=fast_T_average_warm_start,
            initial_electron_midpoint_density_m3=electron_midpoint_warm_start,
            initial_electron_confined_average_density_m3=electron_average_warm_start,
            initial_electron_parent_n0_m3=electron_parent_n0_warm_start,
            initial_electron_collision_density_m3=electron_collision_density_warm_start,
            initial_eq42_density_profile=eq42_density_profile_warm_start,
            scalar_collision_operator_warm_start_by_species=None if scalar_collision_operator_warm_start_by_species is None else dict(scalar_collision_operator_warm_start_by_species),
            scalar_external_fast_field_warm_start_by_test_species=None if scalar_external_fast_field_warm_start_by_test_species is None else {key: tuple(values) for key, values in scalar_external_fast_field_warm_start_by_test_species.items()},
            scalar_wall_barrier_energy_warm_start_J=scalar_wall_barrier_energy_warm_start_J,
            evaluation_mode=iteration_evaluation_mode,
        )
        kinetic_elapsed_s = max(0.0, float(clock() - kinetic_started))
       
        record_runtime(runtime_profile, "kinetic_build", kinetic_elapsed_s)
      
        kinetic_profile = kinetic.metadata.get("runtime_kinetic_profile")
        kinetic_runtime_evaluations.append({
            "kind": "beam_iteration",
            "beam_iteration": iteration,
            "evaluation_mode": iteration_evaluation_mode.value,
            "runtime_s": kinetic_elapsed_s,
            "profile": kinetic_profile,
        })
       
        emit_progress(progress_callback, OperatingPointProgressEvent(
                phase="kinetic_stage_completed",
                evaluation_mode=evaluation_mode,
                temperature_trial_index=temperature_trial_index,
                temperature_keV=float(temperature / energy_J_from_keV(1.0)),
                beam_iteration=iteration,
                elapsed_stage_s=kinetic_elapsed_s,
                message=compact_runtime_summary(kinetic_profile),
            ),
        )
       
        if kinetic.operating_point_density_state is None:
            source_state = str(kinetic.metadata.get("kinetic_source_state", "")).strip().lower()
            last_beam = beam
            last_kinetic = kinetic
           
            if benchmark_closure:
                benchmark_bypass = True
                converged = True
            elif source_state == "no_confined_fbis_source":
                history.append(BeamDensityCouplingIteration(iteration=iteration, density_pointwise_relative_change=float("inf"), density_volume_L2_relative_change=float("inf"), density_absolute_reference_change=float("inf"), species_density_pointwise_relative_change=None, species_density_volume_L2_relative_change=None, species_density_absolute_reference_change=None, species_density_relative_change_by_species=None, total_fast_birth_rate_relative_change=None, deposited_beam_power_relative_change=None, birth_profile_pointwise_relative_change=None, birth_profile_volume_L2_relative_change=None, birth_profile_absolute_reference_change=None, component_pointwise_relative_change=None, component_L2_relative_change=None, component_absolute_reference_change=None, state_finite=True, state_physically_evaluable=False, kinetic_numerical_convergence_passed=False, end_loss_power_terms_available=False))
                failure_reason = "no_confined_fbis_source_operating_point_unavailable"
            else:
                raise ValueError("kinetic stage did not return an operating point density state")
            break
       
        density_state = _solved_density_state(kinetic, geometry)
        solved_electron = np.asarray(density_state.electron_cell_density_m3, dtype=float)
        solved_scope = np.asarray(density_state.profile_support_mask, dtype=bool)
      
        if np.any(target_physical_support & ~solved_scope):
            raise ValueError("solved electron profile is undefined on required beam target cells")
       
        solved_electron = np.where(solved_scope, solved_electron, 0.0)
        profile_metrics = density_profile_change_metrics(
            solved_electron,
            target_density,
            current_support_mask=target_physical_support,
            previous_support_mask=target_physical_support,
            physical_cell_volumes_m3=np.asarray(geometry.cell_volumes_m3, dtype=float),
            density_floor_absolute_m3=active_controls.density_floor_absolute_m3,
            density_floor_reference_fraction=active_controls.density_floor_reference_fraction,
        )
        species_density_metrics, species_density_changes, species_density_change_metrics_by_species = _kinetic_species_density_change_metrics(
            kinetic,
            previous_kinetic_target_state,
            support_mask=target_physical_support,
            physical_cell_volumes_m3=np.asarray(geometry.cell_volumes_m3, dtype=float),
            controls=active_controls,
        )
      
        if previous_beam_diagnostics is None:
            birth_rate_change = None
            deposited_power_change = None
            birth_profile_metrics = None
            component_metrics = None
        else:
            birth_rate_change = _relative_mapping_change(beam_diagnostics.total_birth_rate_s_by_species, previous_beam_diagnostics.total_birth_rate_s_by_species)
            deposited_power_change = _relative_mapping_change(beam_diagnostics.deposited_power_W_by_species, previous_beam_diagnostics.deposited_power_W_by_species)
            birth_profile_metrics, component_metrics = _beam_source_change_metrics(beam_diagnostics, previous_beam_diagnostics, active_controls)
       
        solve_status = _kinetic_check_status(kinetic, density_state)
        if not solve_status.state_finite:
            raise ValueError("operating point state is nonfinite")
        if not solve_status.state_physically_evaluable:
            raise ValueError("operating point state is not physically evaluable")
      
        coupling_metrics_passed = _all_coupling_metrics_pass(
            active_controls,
            profile_pointwise_change=profile_metrics.pointwise_relative_change,
            profile_volume_L2_change=profile_metrics.volume_L2_relative_change,
            profile_absolute_reference_change=profile_metrics.absolute_reference_change,
            species_density_metrics=species_density_metrics,
            birth_rate_change=birth_rate_change,
            deposited_power_change=deposited_power_change,
            birth_profile_metrics=birth_profile_metrics,
            component_metrics=component_metrics,
        )
        candidate_beam = None
        candidate_passed = False
        final_birth_rate_change = None
        final_deposited_power_change = None
        final_birth_profile_metrics = None
        final_component_metrics = None
      
        if coupling_metrics_passed:
            candidate_beam_runtime_started = clock()
            candidate_beam = beam_builder(config, geometry, electron_temperature_J=temperature, operating_point_density_state=density_state, kinetic_target_state=kinetic)
            record_runtime(runtime_profile, "candidate_beam_build", clock() - candidate_beam_runtime_started)
            candidate_diagnostics = _beam_diagnostics(candidate_beam)
            final_birth_rate_change = _relative_mapping_change(candidate_diagnostics.total_birth_rate_s_by_species, beam_diagnostics.total_birth_rate_s_by_species)
            final_deposited_power_change = _relative_mapping_change(candidate_diagnostics.deposited_power_W_by_species, beam_diagnostics.deposited_power_W_by_species)
            final_birth_profile_metrics, final_component_metrics = _beam_source_change_metrics(candidate_diagnostics, beam_diagnostics, active_controls)
            candidate_passed = _all_coupling_metrics_pass(
                active_controls,
                profile_pointwise_change=0.0,
                profile_volume_L2_change=0.0,
                profile_absolute_reference_change=0.0,
                species_density_metrics=species_density_metrics,
                birth_rate_change=final_birth_rate_change,
                deposited_power_change=final_deposited_power_change,
                birth_profile_metrics=final_birth_profile_metrics,
                component_metrics=final_component_metrics,
            )
        history.append(
            BeamDensityCouplingIteration(
                iteration=iteration,
                density_pointwise_relative_change=profile_metrics.pointwise_relative_change,
                density_volume_L2_relative_change=profile_metrics.volume_L2_relative_change,
                density_absolute_reference_change=profile_metrics.absolute_reference_change,
                species_density_pointwise_relative_change=None if species_density_metrics is None else species_density_metrics.pointwise_relative_change,
                species_density_volume_L2_relative_change=None if species_density_metrics is None else species_density_metrics.volume_L2_relative_change,
                species_density_absolute_reference_change=None if species_density_metrics is None else species_density_metrics.absolute_reference_change,
                species_density_relative_change_by_species=species_density_changes,
                species_density_change_metrics_by_species=species_density_change_metrics_by_species,
                total_fast_birth_rate_relative_change=birth_rate_change,
                deposited_beam_power_relative_change=deposited_power_change,
                birth_profile_pointwise_relative_change=None if birth_profile_metrics is None else birth_profile_metrics.pointwise_relative_change,
                birth_profile_volume_L2_relative_change=None if birth_profile_metrics is None else birth_profile_metrics.weighted_L2_relative_change,
                birth_profile_absolute_reference_change=None if birth_profile_metrics is None else birth_profile_metrics.absolute_reference_change,
                component_pointwise_relative_change=None if component_metrics is None else component_metrics.pointwise_relative_change,
                component_L2_relative_change=None if component_metrics is None else component_metrics.weighted_L2_relative_change,
                component_absolute_reference_change=None if component_metrics is None else component_metrics.absolute_reference_change,
                state_finite=solve_status.state_finite,
                state_physically_evaluable=solve_status.state_physically_evaluable,
                kinetic_numerical_convergence_passed=solve_status.kinetic_numerical_convergence_passed,
                end_loss_power_terms_available=solve_status.end_loss_power_terms_available,
                candidate_fixed_point_passed=candidate_passed,
            )
        )
        iteration_gate_records = _final_consistency_metric_gates(
            active_controls,
            profile_metrics=profile_metrics,
            species_density_metrics=species_density_metrics,
            species_density_change_metrics_by_species=species_density_change_metrics_by_species,
            birth_rate_change=birth_rate_change,
            deposited_power_change=deposited_power_change,
            birth_profile_metrics=birth_profile_metrics,
            component_metrics=component_metrics,
        )
        progress_message = _progress_gate_failure_text(iteration_gate_records)
       
        if coupling_metrics_passed:
            if candidate_passed:
                progress_message = progress_message + "; candidate beam confirmation passed"
            elif (
                final_birth_rate_change is not None
                and final_deposited_power_change is not None
                and final_birth_profile_metrics is not None
                and final_component_metrics is not None
            ):
                progress_message = progress_message + "; " + _candidate_beam_progress_failure_text(
                    active_controls,
                    birth_rate_change=final_birth_rate_change,
                    deposited_power_change=final_deposited_power_change,
                    birth_profile_metrics=final_birth_profile_metrics,
                    component_metrics=final_component_metrics,
                )
       
        emit_progress(progress_callback, OperatingPointProgressEvent(
                phase="beam_density_iteration_completed",
                evaluation_mode=evaluation_mode,
                temperature_trial_index=temperature_trial_index,
                temperature_keV=float(temperature / energy_J_from_keV(1.0)),
                warm_start_source_temperature_keV=warm_start_source_temperature_keV,
                beam_iteration=iteration,
                beam_iteration_limit=active_controls.max_iterations,
                elapsed_stage_s=max(0.0, float(clock() - iteration_started)),
                beam_profile_change=profile_metrics.volume_L2_relative_change,
                message=progress_message,
            ),
        )
       
        last_beam = beam
        last_kinetic = kinetic
        last_density_state = density_state
        kinetic_target_state = kinetic
        eq59_warm_start_state = kinetic.eq59_warm_start_state
        eq59_warm_start_state_by_species = kinetic.eq59_warm_start_state_by_species
      
        if kinetic.modal_basis_warm_start is not None:
            modal_basis = kinetic.modal_basis_warm_start
       
        eq42_density_profile_warm_start = kinetic.eq42_density_profile_warm_start
        eq70_warm_start = kinetic.eq70_potential_energy_warm_start_J
        scalar_collision_operator_warm_start_by_species = kinetic.scalar_collision_operator_warm_start_by_species
        scalar_external_fast_field_warm_start_by_test_species = kinetic.scalar_external_fast_field_warm_start_by_test_species
        scalar_wall_barrier_energy_warm_start_J = kinetic.scalar_wall_barrier_energy_warm_start_J
        fast_midpoint_warm_start, fast_average_warm_start, fast_T_midpoint_warm_start, fast_T_average_warm_start, electron_midpoint_warm_start, electron_average_warm_start, electron_parent_n0_warm_start, electron_collision_density_warm_start = density_warm_values(density_state)
       
        # Final consistency recomputes the kinetic state from the confirmed candidate beam
        if candidate_passed and candidate_beam is not None:
            candidate_fixed_point_passed = True
            final_kinetic_started = clock()
            final_kinetic = kinetic_builder(
                config,
                geometry,
                candidate_beam,
                electron_temperature_J=temperature,
                collision_state_recompute_policy=collision_state_recompute_policy,
                collision_state_is_self_consistent_electron_temperature=collision_state_is_self_consistent_electron_temperature,
                modal_basis=modal_basis,
                eq59_warm_start_state=eq59_warm_start_state,
                eq59_warm_start_state_by_species=eq59_warm_start_state_by_species,
                eq70_initial_potential_energy_J=eq70_warm_start,
                initial_representative_fast_density_m3=fast_average_warm_start,
                initial_fast_D_midpoint_density_m3=fast_midpoint_warm_start,
                initial_fast_D_confined_average_density_m3=fast_average_warm_start,
                initial_fast_T_midpoint_density_m3=fast_T_midpoint_warm_start,
                initial_fast_T_confined_average_density_m3=fast_T_average_warm_start,
                initial_electron_midpoint_density_m3=electron_midpoint_warm_start,
                initial_electron_confined_average_density_m3=electron_average_warm_start,
                initial_electron_parent_n0_m3=electron_parent_n0_warm_start,
                initial_electron_collision_density_m3=electron_collision_density_warm_start,
                initial_eq42_density_profile=eq42_density_profile_warm_start,
                scalar_collision_operator_warm_start_by_species=None if scalar_collision_operator_warm_start_by_species is None else dict(scalar_collision_operator_warm_start_by_species),
                scalar_external_fast_field_warm_start_by_test_species=None if scalar_external_fast_field_warm_start_by_test_species is None else {key: tuple(values) for key, values in scalar_external_fast_field_warm_start_by_test_species.items()},
                scalar_wall_barrier_energy_warm_start_J=scalar_wall_barrier_energy_warm_start_J,
                evaluation_mode=evaluation_mode,
            )
            final_kinetic_elapsed_s = max(0.0, float(clock() - final_kinetic_started))
            record_runtime(runtime_profile, "final_consistency_kinetic_build", final_kinetic_elapsed_s)
            final_kinetic_profile = final_kinetic.metadata.get("runtime_kinetic_profile")
            kinetic_runtime_evaluations.append({
                "kind": "final_consistency",
                "beam_iteration": iteration,
                "evaluation_mode": evaluation_mode.value,
                "runtime_s": final_kinetic_elapsed_s,
                "profile": final_kinetic_profile,
            })
            emit_progress(progress_callback, OperatingPointProgressEvent(
                    phase="final_beam_kinetic_consistency_completed",
                    evaluation_mode=evaluation_mode,
                    temperature_trial_index=temperature_trial_index,
                    temperature_keV=float(temperature / energy_J_from_keV(1.0)),
                    beam_iteration=iteration,
                    elapsed_stage_s=final_kinetic_elapsed_s,
                    message=compact_runtime_summary(final_kinetic_profile),
                ),
            )
            final_density_state = _solved_density_state(final_kinetic, geometry)
            final_solve_status = _kinetic_check_status(final_kinetic, final_density_state)
          
            if not final_solve_status.state_finite or not final_solve_status.state_physically_evaluable:
                raise ValueError("final consistency operating point state is invalid")
            
            candidate_target = np.asarray(candidate_beam.metadata.get("beam_target_density_profile_m3", density_state.electron_cell_density_m3), dtype=float)
            candidate_support = np.asarray(candidate_beam.metadata.get("beam_target_density_profile_support_mask", density_state.profile_support_mask), dtype=bool)
            final_solved_electron = np.asarray(final_density_state.electron_cell_density_m3, dtype=float)
            final_solved_support = np.asarray(final_density_state.profile_support_mask, dtype=bool)
          
            if np.any(candidate_support & ~final_solved_support):
                raise ValueError("final consistency electron profile is undefined on beam support")
            
            final_profile_metrics = density_profile_change_metrics(
                np.where(final_solved_support, final_solved_electron, 0.0),
                candidate_target,
                current_support_mask=candidate_support,
                previous_support_mask=candidate_support,
                physical_cell_volumes_m3=np.asarray(geometry.cell_volumes_m3, dtype=float),
                density_floor_absolute_m3=active_controls.density_floor_absolute_m3,
                density_floor_reference_fraction=active_controls.density_floor_reference_fraction,
            )
            final_species_density_metrics, final_species_density_changes, final_species_density_change_metrics_by_species = _kinetic_species_density_change_metrics(
                final_kinetic,
                kinetic,
                support_mask=candidate_support,
                physical_cell_volumes_m3=np.asarray(geometry.cell_volumes_m3, dtype=float),
                controls=active_controls,
            )
            particle_balance_state = build_stationary_particle_balance(
                beam=candidate_beam,
                kinetic=final_kinetic,
                relative_tolerance=active_controls.total_birth_rate_relative_tolerance,
            ) if startup_seed_only else None
            implied_beam_runtime_started = clock()
            implied_beam = beam_builder(config, geometry, electron_temperature_J=temperature, operating_point_density_state=final_density_state, kinetic_target_state=final_kinetic)
            record_runtime(runtime_profile, "final_implied_beam_build", clock() - implied_beam_runtime_started)
            implied_diagnostics = _beam_diagnostics(implied_beam)
            candidate_diagnostics = _beam_diagnostics(candidate_beam)
            final_birth_rate_change = _relative_change(implied_diagnostics.total_birth_rate_s, candidate_diagnostics.total_birth_rate_s)
            final_deposited_power_change = _relative_change(implied_diagnostics.deposited_power_W, candidate_diagnostics.deposited_power_W)
            final_birth_profile_metrics, final_component_metrics = _beam_source_change_metrics(implied_diagnostics, candidate_diagnostics, active_controls)
            final_consistency_beam_metrics_passed = _all_coupling_metrics_pass(
                active_controls,
                profile_pointwise_change=final_profile_metrics.pointwise_relative_change,
                profile_volume_L2_change=final_profile_metrics.volume_L2_relative_change,
                profile_absolute_reference_change=final_profile_metrics.absolute_reference_change,
                species_density_metrics=final_species_density_metrics,
                birth_rate_change=final_birth_rate_change,
                deposited_power_change=final_deposited_power_change,
                birth_profile_metrics=final_birth_profile_metrics,
                component_metrics=final_component_metrics,
            )
            final_consistency_particle_balance_required = bool(startup_seed_only)
            final_consistency_particle_balance_passed = bool(
                not startup_seed_only
                or particle_balance_state is not None and particle_balance_state.converged
            )
            final_consistency_recompute_passed = bool(final_consistency_beam_metrics_passed)
            final_consistency_overall_qualification_passed = bool(
                final_consistency_beam_metrics_passed
                and final_consistency_particle_balance_passed
            )
            latest_final_consistency_record = {
                "iteration": iteration,
                "beam_metrics_passed": final_consistency_beam_metrics_passed,
                "stationary_particle_balance_required": final_consistency_particle_balance_required,
                "stationary_particle_balance_passed": final_consistency_particle_balance_passed,
                "overall_passed": final_consistency_overall_qualification_passed,
                "kinetic_state_finite": final_solve_status.state_finite,
                "kinetic_state_physically_evaluable": final_solve_status.state_physically_evaluable,
                "kinetic_numerical_convergence_passed": final_solve_status.kinetic_numerical_convergence_passed,
                "metric_gates": _final_consistency_metric_gates(
                    active_controls,
                    profile_metrics=final_profile_metrics,
                    species_density_metrics=final_species_density_metrics,
                    species_density_change_metrics_by_species=final_species_density_change_metrics_by_species,
                    birth_rate_change=final_birth_rate_change,
                    deposited_power_change=final_deposited_power_change,
                    birth_profile_metrics=final_birth_profile_metrics,
                    component_metrics=final_component_metrics,
                ),
                "stationary_particle_balance": None if particle_balance_state is None else particle_balance_state.as_metadata(),
            }
            final_consistency_attempt_history.append(latest_final_consistency_record)
            history[-1] = replace(history[-1], final_consistency_recompute_passed=final_consistency_recompute_passed)
            final_beam = candidate_beam
            last_beam = candidate_beam
            last_kinetic = final_kinetic
            last_density_state = final_density_state
            solve_status = final_solve_status
            profile_metrics = final_profile_metrics
            if final_consistency_recompute_passed:
                converged = True
                beam_density_termination_reason = "same_source_final_beam_density_metrics_converged"
                failure_reason = None
                break
            failure_reason = "final_beam_kinetic_consistency_recompute_not_converged"
            previous_beam_diagnostics = candidate_diagnostics
            target_density = np.where(candidate_support, (1.0 - active_controls.target_density_relaxation) * candidate_target + active_controls.target_density_relaxation * final_solved_electron, 0.0)
            target_defined = candidate_support
            density_warm_start = final_density_state
            kinetic_target_state = final_kinetic
            eq59_warm_start_state = final_kinetic.eq59_warm_start_state
            eq59_warm_start_state_by_species = final_kinetic.eq59_warm_start_state_by_species
            if final_kinetic.modal_basis_warm_start is not None:
                modal_basis = final_kinetic.modal_basis_warm_start
            eq42_density_profile_warm_start = final_kinetic.eq42_density_profile_warm_start
            eq70_warm_start = final_kinetic.eq70_potential_energy_warm_start_J
            scalar_collision_operator_warm_start_by_species = final_kinetic.scalar_collision_operator_warm_start_by_species
            scalar_external_fast_field_warm_start_by_test_species = final_kinetic.scalar_external_fast_field_warm_start_by_test_species
            scalar_wall_barrier_energy_warm_start_J = final_kinetic.scalar_wall_barrier_energy_warm_start_J
            fast_midpoint_warm_start, fast_average_warm_start, fast_T_midpoint_warm_start, fast_T_average_warm_start, electron_midpoint_warm_start, electron_average_warm_start, electron_parent_n0_warm_start, electron_collision_density_warm_start = density_warm_values(final_density_state)
            continue
        previous_beam_diagnostics = beam_diagnostics
        next_target_support = solved_scope if startup_seed_only else target_physical_support
        target_density = np.where(next_target_support, (1.0 - active_controls.target_density_relaxation) * target_density + active_controls.target_density_relaxation * solved_electron, 0.0)
        target_defined = next_target_support

    if last_beam is None or last_kinetic is None:
        raise RuntimeError("operating point loop produced no paired beam and kinetic state")
    terminal_failure_qualification_performed = False
    if (
        evaluation_mode.runs_final_qualification
        and not converged
        and last_density_state is not None
        and str(last_kinetic.metadata.get("operating_point_evaluation_mode", ""))
        != evaluation_mode.value
    ):
        terminal_qualification_started = clock()
        terminal_kinetic = kinetic_builder(
            config,
            geometry,
            last_beam,
            electron_temperature_J=temperature,
            collision_state_recompute_policy=collision_state_recompute_policy,
            collision_state_is_self_consistent_electron_temperature=collision_state_is_self_consistent_electron_temperature,
            modal_basis=modal_basis,
            eq59_warm_start_state=last_kinetic.eq59_warm_start_state,
            eq59_warm_start_state_by_species=last_kinetic.eq59_warm_start_state_by_species,
            eq70_initial_potential_energy_J=last_kinetic.eq70_potential_energy_warm_start_J,
            initial_representative_fast_density_m3=last_density_state.fast_deuterium_confined_volume_average_density_m3,
            initial_fast_D_midpoint_density_m3=last_density_state.fast_deuterium_midplane_density_m3,
            initial_fast_D_confined_average_density_m3=last_density_state.fast_deuterium_confined_volume_average_density_m3,
            initial_fast_T_midpoint_density_m3=last_density_state.fast_tritium_midplane_density_m3,
            initial_fast_T_confined_average_density_m3=last_density_state.fast_tritium_confined_volume_average_density_m3,
            initial_electron_midpoint_density_m3=last_density_state.electron_midplane_density_m3,
            initial_electron_confined_average_density_m3=last_density_state.electron_confined_volume_average_density_m3,
            initial_electron_parent_n0_m3=last_density_state.electron_parent_maxwellian_n0_m3,
            initial_electron_collision_density_m3=last_density_state.electron_collision_density_m3,
            initial_eq42_density_profile=last_kinetic.eq42_density_profile_warm_start,
            scalar_collision_operator_warm_start_by_species=None if last_kinetic.scalar_collision_operator_warm_start_by_species is None else dict(last_kinetic.scalar_collision_operator_warm_start_by_species),
            scalar_external_fast_field_warm_start_by_test_species=None if last_kinetic.scalar_external_fast_field_warm_start_by_test_species is None else {key: tuple(values) for key, values in last_kinetic.scalar_external_fast_field_warm_start_by_test_species.items()},
            scalar_wall_barrier_energy_warm_start_J=last_kinetic.scalar_wall_barrier_energy_warm_start_J,
            evaluation_mode=evaluation_mode,
        )
        terminal_qualification_runtime_s = max(0.0, float(clock() - terminal_qualification_started))
        record_runtime(runtime_profile, "terminal_failure_qualification_kinetic_build", terminal_qualification_runtime_s)
        terminal_profile = terminal_kinetic.metadata.get("runtime_kinetic_profile")
        kinetic_runtime_evaluations.append({
            "kind": "terminal_failure_qualification",
            "beam_iteration": len(history),
            "evaluation_mode": evaluation_mode.value,
            "runtime_s": terminal_qualification_runtime_s,
            "profile": terminal_profile,
        })
        terminal_density_state = _solved_density_state(terminal_kinetic, geometry)
        terminal_solve_status = _kinetic_check_status(terminal_kinetic, terminal_density_state)
        if not terminal_solve_status.state_finite or not terminal_solve_status.state_physically_evaluable:
            raise ValueError("terminal qualification operating point state is invalid")
        last_kinetic = terminal_kinetic
        last_density_state = terminal_density_state
        solve_status = terminal_solve_status
        terminal_failure_qualification_performed = True
    if last_density_state is None:
        kinetic_converged, numerical_records = kinetic_numerical_convergence(last_kinetic.metadata)
        operating_point_converged = bool(benchmark_bypass and kinetic_converged)
        if benchmark_bypass and not operating_point_converged:
            failure_reason = "required_kinetic_numerical_convergence_not_passed"
        runtime_fixed_temperature_operating_point_profile = finalize_runtime_profile(runtime_profile, total_s=clock() - operating_runtime_started)
        metadata = {'runtime_fixed_temperature_operating_point_profile': runtime_fixed_temperature_operating_point_profile, 'runtime_kinetic_evaluations': tuple(kinetic_runtime_evaluations), 'runtime_beam_intermediate_kinetic_evaluation_mode': iteration_evaluation_mode.value, 'runtime_beam_intermediate_final_qualification_deferred': bool(iteration_evaluation_mode is not evaluation_mode), 'runtime_terminal_failure_qualification_performed': bool(terminal_failure_qualification_performed), 'beam_density_fixed_point_converged': bool(benchmark_bypass), 'beam_density_coupling_iteration_count': len(history), 'beam_density_coupling_iteration_limit': active_controls.max_iterations, 'beam_density_coupling_reached_iteration_limit': bool(len(history) >= active_controls.max_iterations), 'beam_density_coupling_exhausted_iteration_limit': bool(len(history) >= active_controls.max_iterations and (not converged)), 'beam_density_coupling_history_export_requested': bool(config.output.write_iteration_histories), 'operating_point_numerically_converged': operating_point_converged, 'operating_point_failure_reason': failure_reason, 'operating_point_density_state_available': False, 'state_finite': bool(benchmark_bypass), 'state_physically_evaluable': bool(benchmark_bypass), 'kinetic_numerical_convergence_passed': kinetic_converged, 'kinetic_required_numerical_check_names': list(numerical_records), 'end_loss_power_terms_available': last_kinetic.metadata.get('modal_end_loss_power_terms_available') is True, 'beam_target_density_profile_valid': last_beam.metadata.get('beam_target_density_profile_valid') is True, 'beam_density_coupling_numerics_settings': active_controls.as_metadata(), 'operating_point_evaluation_mode': evaluation_mode.value, 'electron_temperature_power_balance_scope': 'fast_D_fast_T_electron_FBIS_subsystem' if last_kinetic.fast_ion_system_state is not None and 'tritium' in last_kinetic.fast_ion_system_state.active_species_ids else 'fast_D_electron_FBIS_subsystem', 'global_plasma_power_balance_available': False, 'startup_seed_profile_role': 'initialization_only' if startup_seed_only else 'prescribed_reference_population', 'startup_seed_replenished': not startup_seed_only, 'startup_seed_final_density_authority': False if startup_seed_only else True, 'startup_seed_final_fusion_authority': False if startup_seed_only else True, 'startup_seed_initial_target_source': initial_target_source, 'startup_seed_heavy_particle_target_transition': 'solved_kinetic_D_T_after_first_iteration' if startup_seed_only else 'not_applicable'}
        if config.output.write_iteration_histories:
            metadata["beam_density_coupling_history"] = [record.as_metadata() for record in history]
        last_kinetic.metadata.update(metadata)
        warm_start = OperatingPointWarmStartState(
            target_density_profile_m3=np.asarray(last_beam.metadata.get("beam_target_density_profile_m3", target_density), dtype=float).copy(),
            target_density_defined_mask=np.asarray(last_beam.metadata.get("beam_target_density_profile_support_mask", target_defined), dtype=bool).copy(),
            kinetic_target_state=last_kinetic,
            eq59_warm_start_state=last_kinetic.eq59_warm_start_state,
            eq59_warm_start_state_by_species=last_kinetic.eq59_warm_start_state_by_species,
            modal_basis_warm_start=last_kinetic.modal_basis_warm_start,
            eq70_potential_energy_warm_start_J=last_kinetic.eq70_potential_energy_warm_start_J,
            eq42_density_profile_warm_start=last_kinetic.eq42_density_profile_warm_start,
            scalar_collision_operator_warm_start_by_species=last_kinetic.scalar_collision_operator_warm_start_by_species,
            scalar_external_fast_field_warm_start_by_test_species=last_kinetic.scalar_external_fast_field_warm_start_by_test_species,
            scalar_wall_barrier_energy_warm_start_J=last_kinetic.scalar_wall_barrier_energy_warm_start_J,
        )
        result = OperatingPointStageResult(beam=last_beam, kinetic=last_kinetic, density_state=None, beam_density_coupling_history=tuple(history), converged=operating_point_converged, failure_reason=failure_reason, warm_start_state=warm_start, metadata=metadata, state_finite=bool(benchmark_bypass), state_physically_evaluable=bool(benchmark_bypass), beam_density_fixed_point_converged=bool(benchmark_bypass), kinetic_numerical_convergence_passed=kinetic_converged, end_loss_power_terms_available=metadata["end_loss_power_terms_available"], electron_temperature_assignment_valid=True)
        return result

    if final_beam is None:
        final_beam = last_beam
        if history:
            final_birth_rate_change = history[-1].total_fast_birth_rate_relative_change
            final_deposited_power_change = history[-1].deposited_beam_power_relative_change
            final_birth_profile_metrics = birth_profile_metrics
            final_component_metrics = component_metrics
    if not converged and failure_reason is None:
        failure_reason = "beam_density_coupling_not_converged"
    if solve_status is None or profile_metrics is None:
        raise RuntimeError("operating point convergence state is unavailable")
    operating_point_converged = bool(converged and solve_status.kinetic_numerical_convergence_passed)
    if converged and not operating_point_converged:
        failure_reason = "required_kinetic_numerical_convergence_not_passed"
    if startup_seed_only and particle_balance_state is None:
        particle_balance_state = build_stationary_particle_balance(beam=final_beam, kinetic=last_kinetic, relative_tolerance=active_controls.total_birth_rate_relative_tolerance)
    particle_balance_metadata = particle_balance_state.as_metadata() if particle_balance_state is not None else {
        "nbi_supported_particle_balance_model": "not_applicable_reference_mode",
        "nbi_supported_particle_balance_relative_tolerance": None,
        "nbi_supported_particle_balance_converged": None,
        "nbi_supported_particle_balance_by_species": {},
    }
    runtime_fixed_temperature_operating_point_profile = finalize_runtime_profile(runtime_profile, total_s=clock() - operating_runtime_started)
    metadata = {**last_density_state.as_metadata(), 'runtime_fixed_temperature_operating_point_profile': runtime_fixed_temperature_operating_point_profile, 'runtime_kinetic_evaluations': tuple(kinetic_runtime_evaluations), 'runtime_beam_intermediate_kinetic_evaluation_mode': iteration_evaluation_mode.value, 'runtime_beam_intermediate_final_qualification_deferred': bool(iteration_evaluation_mode is not evaluation_mode), 'runtime_terminal_failure_qualification_performed': bool(terminal_failure_qualification_performed), 'beam_target_density_profile_valid': final_beam.metadata.get('beam_target_density_profile_valid') is True, 'beam_density_fixed_point_converged': converged, 'beam_density_coupling_iteration_count': len(history), 'beam_density_coupling_iteration_limit': active_controls.max_iterations, 'beam_density_coupling_reached_iteration_limit': bool(len(history) >= active_controls.max_iterations), 'beam_density_coupling_exhausted_iteration_limit': bool(len(history) >= active_controls.max_iterations and (not converged)), 'beam_density_coupling_history_export_requested': bool(config.output.write_iteration_histories), 'density_pointwise_relative_change': profile_metrics.pointwise_relative_change, 'density_volume_L2_relative_change': profile_metrics.volume_L2_relative_change, 'density_absolute_reference_change': profile_metrics.absolute_reference_change, 'species_density_pointwise_relative_change': None if final_species_density_metrics is None else final_species_density_metrics.pointwise_relative_change, 'species_density_volume_L2_relative_change': None if final_species_density_metrics is None else final_species_density_metrics.volume_L2_relative_change, 'species_density_absolute_reference_change': None if final_species_density_metrics is None else final_species_density_metrics.absolute_reference_change, 'species_density_relative_change_by_species': final_species_density_changes, 'species_density_change_metrics_by_species': final_species_density_change_metrics_by_species, 'beam_birth_rate_relative_change': final_birth_rate_change, 'beam_deposited_power_relative_change': final_deposited_power_change, 'beam_birth_profile_pointwise_relative_change': None if final_birth_profile_metrics is None else final_birth_profile_metrics.pointwise_relative_change, 'beam_birth_profile_volume_L2_relative_change': None if final_birth_profile_metrics is None else final_birth_profile_metrics.weighted_L2_relative_change, 'beam_birth_profile_absolute_reference_change': None if final_birth_profile_metrics is None else final_birth_profile_metrics.absolute_reference_change, 'beam_component_pointwise_relative_change': None if final_component_metrics is None else final_component_metrics.pointwise_relative_change, 'beam_component_L2_relative_change': None if final_component_metrics is None else final_component_metrics.weighted_L2_relative_change, 'beam_component_absolute_reference_change': None if final_component_metrics is None else final_component_metrics.absolute_reference_change, 'beam_axial_birth_profile_converged': bool(converged and final_birth_profile_metrics is not None and (final_birth_profile_metrics.pointwise_relative_change <= active_controls.axial_birth_profile_relative_tolerance) and (final_birth_profile_metrics.weighted_L2_relative_change <= active_controls.axial_birth_profile_relative_tolerance) and (final_birth_profile_metrics.absolute_reference_change <= active_controls.axial_birth_profile_relative_tolerance)), 'beam_energy_component_deposition_converged': bool(converged and final_component_metrics is not None and (final_component_metrics.pointwise_relative_change <= active_controls.energy_component_deposition_relative_tolerance) and (final_component_metrics.weighted_L2_relative_change <= active_controls.energy_component_deposition_relative_tolerance) and (final_component_metrics.absolute_reference_change <= active_controls.energy_component_deposition_relative_tolerance)), 'candidate_fixed_point_passed': candidate_fixed_point_passed, 'final_consistency_recompute_passed': final_consistency_recompute_passed, 'beam_density_final_consistency_attempt_count': len(final_consistency_attempt_history), 'beam_density_final_consistency_numerical_recompute_passed': final_consistency_recompute_passed, 'beam_density_final_consistency_stationary_particle_balance_in_stopping_condition': False, 'beam_density_final_consistency_overall_qualification_passed': final_consistency_overall_qualification_passed, 'beam_density_numerical_termination_reason': beam_density_termination_reason, 'beam_density_final_consistency_recompute_performed': latest_final_consistency_record is not None, 'beam_density_final_consistency_beam_metrics_passed': None if latest_final_consistency_record is None else latest_final_consistency_record['beam_metrics_passed'], 'beam_density_final_consistency_stationary_particle_balance_required': None if latest_final_consistency_record is None else latest_final_consistency_record['stationary_particle_balance_required'], 'beam_density_final_consistency_stationary_particle_balance_passed': None if latest_final_consistency_record is None else latest_final_consistency_record['stationary_particle_balance_passed'], 'beam_density_final_consistency_overall_passed': None if latest_final_consistency_record is None else latest_final_consistency_record['overall_passed'], 'beam_density_final_consistency_metric_gates': None if latest_final_consistency_record is None else latest_final_consistency_record['metric_gates'], 'beam_density_final_consistency_stationary_particle_balance': None if latest_final_consistency_record is None else latest_final_consistency_record['stationary_particle_balance'], 'operating_point_numerically_converged': operating_point_converged, 'operating_point_failure_reason': failure_reason, 'operating_point_density_state_available': True, 'operating_point_electron_temperature_J': temperature, 'state_finite': solve_status.state_finite, 'state_physically_evaluable': solve_status.state_physically_evaluable, 'kinetic_numerical_convergence_passed': solve_status.kinetic_numerical_convergence_passed, 'kinetic_required_numerical_check_names': list(solve_status.convergence_check_names), 'end_loss_power_terms_available': solve_status.end_loss_power_terms_available, 'operating_point_warm_start_used': operating_point_warm_start_used, 'electron_density_initial_guess_used': initial_guess_used, 'beam_density_coupling_numerics_settings': active_controls.as_metadata(), 'operating_point_evaluation_mode': evaluation_mode.value, 'electron_temperature_power_balance_scope': 'fast_D_fast_T_electron_FBIS_subsystem' if last_kinetic.fast_ion_system_state is not None and 'tritium' in last_kinetic.fast_ion_system_state.active_species_ids else 'fast_D_electron_FBIS_subsystem', 'global_plasma_power_balance_available': False, 'startup_seed_profile_role': 'initialization_only' if startup_seed_only else 'prescribed_reference_population', 'startup_seed_replenished': not startup_seed_only, 'startup_seed_final_density_authority': False if startup_seed_only else True, 'startup_seed_final_fusion_authority': False if startup_seed_only else True, 'startup_seed_initial_target_source': initial_target_source, 'startup_seed_heavy_particle_target_transition': 'solved_kinetic_D_T_after_first_iteration' if startup_seed_only else 'not_applicable', 'startup_seed_final_authority_fraction': 0.0 if startup_seed_only else 1.0, 'nbi_supported_temperature_closure_baseline': 'self_consistent_electron_energy' if startup_seed_only else None, 'nbi_supported_fixed_temperature_role': 'benchmark_or_sensitivity_only' if startup_seed_only and config.power_balance.electron_temperature_mode == 'fixed_closure' else None, **particle_balance_metadata}
    if config.output.write_iteration_histories:
        metadata["beam_density_coupling_history"] = [record.as_metadata() for record in history]
        metadata["beam_density_final_consistency_attempt_history"] = list(final_consistency_attempt_history)
    last_kinetic.metadata.update(metadata)
    warm_start = OperatingPointWarmStartState(
        target_density_profile_m3=np.asarray(last_density_state.electron_cell_density_m3, dtype=float).copy(),
        target_density_defined_mask=np.asarray(last_density_state.profile_support_mask, dtype=bool).copy(),
        kinetic_target_state=last_kinetic,
        eq59_warm_start_state=last_kinetic.eq59_warm_start_state,
        eq59_warm_start_state_by_species=last_kinetic.eq59_warm_start_state_by_species,
        density_state=last_density_state,
        modal_basis_warm_start=last_kinetic.modal_basis_warm_start,
        eq70_potential_energy_warm_start_J=last_kinetic.eq70_potential_energy_warm_start_J,
        eq42_density_profile_warm_start=last_kinetic.eq42_density_profile_warm_start,
        scalar_collision_operator_warm_start_by_species=last_kinetic.scalar_collision_operator_warm_start_by_species,
        scalar_external_fast_field_warm_start_by_test_species=last_kinetic.scalar_external_fast_field_warm_start_by_test_species,
        scalar_wall_barrier_energy_warm_start_J=last_kinetic.scalar_wall_barrier_energy_warm_start_J,
    )
    result = OperatingPointStageResult(
        beam=final_beam,
        kinetic=last_kinetic,
        density_state=last_density_state,
        beam_density_coupling_history=tuple(history),
        converged=operating_point_converged,
        failure_reason=failure_reason,
        warm_start_state=warm_start,
        metadata=metadata,
        state_finite=solve_status.state_finite,
        state_physically_evaluable=solve_status.state_physically_evaluable,
        beam_density_fixed_point_converged=converged,
        kinetic_numerical_convergence_passed=solve_status.kinetic_numerical_convergence_passed,
        end_loss_power_terms_available=solve_status.end_loss_power_terms_available,
        electron_temperature_assignment_valid=True,
    )
    return result

__all__ = [
    "BeamSourceProfileChangeMetrics",
    "DensityProfileChangeMetrics",
    "beam_birth_profile_change_metrics",
    "beam_component_deposition_change_metrics",
    "build_fixed_temperature_operating_point",
    "density_profile_change_metrics",
]
