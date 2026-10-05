"""Typed plotting view for the species resolved source model metadata"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any
import numpy as np
from source_model_revamp.fbis.species import DEUTERON, TRITON, ion_species
from source_model_revamp.plotting.common import MissingMetadata

PLOTTING_METADATA_CONTRACT_VERSION = 1
PLOTTING_SPECIES_ORDER = (DEUTERON.species_id, TRITON.species_id)
CURRENT_COMPONENT_ORDER = ('fast_deuterium', 'fast_tritium', 'prompt_deuterium', 'prompt_tritium')

def _readonly_array(value: Any, name: str, *, ndim: int | None = None, dtype: Any = float, required: bool = False, nonnegative: bool = False, allow_nan: bool = False) -> np.ndarray:
    if value is None:
        if required:
            raise MissingMetadata(f"missing {name}")
        result = np.asarray([], dtype=dtype)
    else:
        try:
            result = np.asarray(value, dtype=dtype)
        except Exception as exc:
            raise MissingMetadata(f"{name} could not be converted to an array") from exc
    if ndim is not None and result.size and result.ndim != ndim:
        raise MissingMetadata(f"{name} must be {ndim}D, got shape {result.shape}")
    if required and result.size == 0:
        raise MissingMetadata(f"{name} must be nonempty")
    if result.size and np.issubdtype(result.dtype, np.number):
        if allow_nan and np.any(np.isinf(result)):
            raise MissingMetadata(f"{name} must not contain infinite values")
        if not allow_nan and not np.all(np.isfinite(result)):
            raise MissingMetadata(f"{name} must contain only finite values")
    if nonnegative and result.size and np.any(result < 0):
        raise MissingMetadata(f"{name} must be nonnegative")
    result = np.array(result, copy=True)
    result.setflags(write=False)
    return result

def _mapping(value: Any, name: str, *, required: bool = False) -> Mapping[str, Any]:
    if value is None:
        if required:
            raise MissingMetadata(f"missing {name}")
        return MappingProxyType({})
    if not isinstance(value, Mapping):
        raise MissingMetadata(f"{name} must be a mapping")
    if required and not value:
        raise MissingMetadata(f"{name} must be nonempty")
    return MappingProxyType({str(key): item for key, item in value.items()})

def _frozen_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _frozen_value(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_frozen_value(item) for item in value)
    if isinstance(value, np.ndarray):
        return _readonly_array(value, "contract array")
    return value

def _optional_float(value: Any, name: str) -> float | None:
    if value is None:
        return None
    try:
        result = float(value)
    except Exception as exc:
        raise MissingMetadata(f"{name} must be numeric") from exc
    if not np.isfinite(result):
        raise MissingMetadata(f"{name} must be finite")
    return result

def _required_float(value: Any, name: str) -> float:
    result = _optional_float(value, name)
    if result is None:
        raise MissingMetadata(f"missing {name}")
    return result

def _required_nonnegative_float(value: Any, name: str) -> float:
    result = _required_float(value, name)
    if result < 0.0:
        raise MissingMetadata(f"{name} must be nonnegative")
    return result

def _required_positive_float(value: Any, name: str) -> float:
    result = _required_float(value, name)
    if result <= 0.0:
        raise MissingMetadata(f"{name} must be positive")
    return result

def _optional_nonnegative_float(value: Any, name: str) -> float | None:
    result = _optional_float(value, name)
    if result is not None and result < 0.0:
        raise MissingMetadata(f"{name} must be nonnegative")
    return result

def _string_tuple(value: Any, name: str) -> tuple[str, ...]:
    if value is None:
        return ()
    if not isinstance(value, (list, tuple)):
        raise MissingMetadata(f"{name} must be a list or tuple")
    return tuple(str(item) for item in value)

def _required_string(value: Any, name: str) -> str:
    result = str(value or "").strip()
    if not result:
        raise MissingMetadata(f"missing {name}")
    return result

def _ordered_ids(values: Mapping[str, Any], preferred: tuple[str, ...]) -> tuple[str, ...]:
    ordered = [value for value in preferred if value in values]
    ordered.extend(sorted(value for value in values if value not in preferred))
    return tuple(ordered)

def _validate_profile(profile: np.ndarray, grid: "AxialGridRecord", name: str) -> None:
    if profile.size and profile.shape != grid.centers_m.shape:
        raise MissingMetadata(f"{name} must match {grid.grid_id} centers")

def _validate_centers_faces(centers: np.ndarray, faces: np.ndarray, name: str) -> None:
    if faces.size != centers.size + 1:
        raise MissingMetadata(f"{name} faces must contain one more value than centers")
    if np.any(np.diff(faces) <= 0.0) or (centers.size > 1 and np.any(np.diff(centers) <= 0.0)):
        raise MissingMetadata(f"{name} centers and faces must be strictly increasing")
    if np.any(centers <= faces[:-1]) or np.any(centers >= faces[1:]):
        raise MissingMetadata(f"{name} centers must lie inside their cells")

def _require_matching_keys(reference: Mapping[str, Any], named_mappings: tuple[tuple[str, Mapping[str, Any]], ...]) -> None:
    expected = set(reference)
    for name, mapping in named_mappings:
        if set(mapping) != expected:
            raise MissingMetadata(f"{name} keys must match the reference species keys")

@dataclass(frozen=True)
class AxialGridRecord:
    grid_id: str
    scope: str
    edges_m: np.ndarray
    centers_m: np.ndarray
    cell_volumes_m3: np.ndarray

@dataclass(frozen=True)
class BeamPlotRecord:
    beam_id: str
    species_id: str
    enabled: bool
    power_W: float
    energy_keV: float
    fast_birth_rate_s: float
    deposited_birth_power_W: float
    shine_through_power_W: float
    beam_radius_m: float
    configured_injection_angle_deg: float | None
    effective_injection_angle_deg: float | None
    beamline_axis_angle_deg: float | None
    injection_geometry_consistent: bool | None
    component_ids: tuple[str, ...]
    component_energies_J: np.ndarray
    component_configured_power_fraction: np.ndarray
    component_deposited_birth_power_W: np.ndarray
    start_m: np.ndarray
    end_m: np.ndarray
    direction_unit: np.ndarray
    path_center_points_m: np.ndarray
    beam_radius_profile_m: np.ndarray
    physical_path_lengths_m: np.ndarray
    effective_path_lengths_m: np.ndarray
    radial_overlap_fractions: np.ndarray
    plasma_radius_profile_m: np.ndarray
    birth_pitch_angle_profile_deg: np.ndarray
    birth_lambda_profile: np.ndarray
    axial_grid_id: str
    axial_birth_rate_density_m3_s: np.ndarray
    axial_birth_rate_s: np.ndarray
    axial_birth_profile_support_mask: np.ndarray
    axial_birth_power_density_W_m3: np.ndarray
    axial_birth_power_W: np.ndarray
    source_z_v_lambda_per_s: np.ndarray
    volume_averaged_source_v_lambda_per_s: np.ndarray

@dataclass(frozen=True)
class BeamSpeciesProfileRecord:
    species_id: str
    axial_grid_id: str
    axial_birth_rate_density_m3_s: np.ndarray
    axial_birth_rate_s: np.ndarray
    axial_birth_power_density_W_m3: np.ndarray
    axial_birth_power_W: np.ndarray

@dataclass(frozen=True)
class FastSpeciesPlotRecord:
    species_id: str
    symbol: str
    mass_kg: float
    charge_number: float
    confined_grid_id: str
    full_device_grid_id: str
    invariant_speed_centers_m_s: np.ndarray
    invariant_speed_faces_m_s: np.ndarray
    lambda_centers: np.ndarray
    lambda_faces: np.ndarray
    distribution_v_lambda: np.ndarray
    local_speed_centers_m_s: np.ndarray
    local_speed_faces_m_s: np.ndarray
    pitch_centers: np.ndarray
    pitch_faces: np.ndarray
    local_distribution_z_v_lambda: np.ndarray
    local_distribution_z_v_pitch: np.ndarray
    local_density_m3: np.ndarray
    full_device_distribution_z_v_pitch: np.ndarray
    full_device_density_m3: np.ndarray


@dataclass(frozen=True)
class ElectrostaticPlotRecord:
    grid_id: str
    potential_relative_to_midplane_V: np.ndarray
    electron_density_m3: np.ndarray
    ion_density_m3: np.ndarray
    wall_barrier_energy_J: float | None
    normalized_wall_barrier: float | None
    converged: bool | None

@dataclass(frozen=True)
class FusionComponentPlotRecord:
    label: str
    kind: str
    reaction: str
    reactant_a_population_id: str
    reactant_b_population_id: str
    reactant_a_species_id: str
    reactant_b_species_id: str
    axial_grid_id: str
    axial_neutron_rate_s: np.ndarray
    axial_neutron_rate_density_m3_s: np.ndarray
    total_reaction_rate_s: float
    total_neutron_rate_s: float

@dataclass(frozen=True)
class CurrentComponentPlotRecord:
    component_id: str
    population_kind: str
    species_id: str
    charge_number: float
    total_particle_loss_rate_s: float
    left_particle_loss_rate_s: float | None
    right_particle_loss_rate_s: float | None
    side_resolved: bool
    current_A: float

@dataclass(frozen=True)
class CurrentBalancePlotRecord:
    components: tuple[CurrentComponentPlotRecord, ...]
    scope: str
    converged: bool
    target_ion_current_A: float | None
    electron_current_A: float | None
    current_residual_A: float | None
    current_relative_residual: float | None
    wall_barrier_energy_J: float | None
    wall_barrier_energy_keV: float | None
    normalized_wall_barrier: float | None
    left_ion_current_A: float | None
    right_ion_current_A: float | None
    end_asymmetry_relative_error: float | None
    end_asymmetry_relative_tolerance: float | None
    end_asymmetry_applicable: bool
    end_asymmetry_check_passed: bool
    component_identity_relative_error: float | None
    component_identity_check_passed: bool
    side_sum_identity_applicable: bool
    side_sum_identity_relative_error: float | None
    side_sum_identity_check_passed: bool
    represented_confined_plasma_balance_claimed: bool
    full_device_wall_balance_claimed: bool

@dataclass(frozen=True)
class PlottingModelView:
    contract_version: int
    contract: Mapping[str, Any]
    grids: Mapping[str, AxialGridRecord]
    beams: tuple[BeamPlotRecord, ...]
    beam_profiles_by_species: Mapping[str, BeamSpeciesProfileRecord]
    fast_species: tuple[FastSpeciesPlotRecord, ...]
    electrostatic: ElectrostaticPlotRecord | None
    fusion_components: tuple[FusionComponentPlotRecord, ...]
    current_balance: CurrentBalancePlotRecord | None

    def grid(self, grid_id: str) -> AxialGridRecord:
        try:
            return self.grids[grid_id]
        except KeyError as exc:
            raise MissingMetadata(f"unknown plotting grid {grid_id!r}") from exc

    def beam(self, beam_id: str) -> BeamPlotRecord:
        for record in self.beams:
            if record.beam_id == beam_id:
                return record
        raise MissingMetadata(f"unknown plotting beam {beam_id!r}")

    def fast(self, species_id: str) -> FastSpeciesPlotRecord:
        for record in self.fast_species:
            if record.species_id == species_id:
                return record
        raise MissingMetadata(f"fast plotting species {species_id!r} is unavailable")



    @property
    def beam_total_birth_rate_density_m3_s(self) -> np.ndarray:
        if not self.beam_profiles_by_species:
            return _readonly_array(None, "beam total birth rate density")
        return _readonly_array(np.sum(tuple(record.axial_birth_rate_density_m3_s for record in self.beam_profiles_by_species.values()), axis=0), "beam total birth rate density", ndim=1, nonnegative=True)

    @property
    def beam_total_birth_power_density_W_m3(self) -> np.ndarray:
        if not self.beam_profiles_by_species:
            return _readonly_array(None, "beam total birth power density")
        return _readonly_array(np.sum(tuple(record.axial_birth_power_density_W_m3 for record in self.beam_profiles_by_species.values()), axis=0), "beam total birth power density", ndim=1, nonnegative=True)

    @property
    def fast_confined_total_density_m3(self) -> np.ndarray:
        if not self.fast_species:
            return _readonly_array(None, "fast confined total density")
        return _readonly_array(np.sum(tuple(record.local_density_m3 for record in self.fast_species), axis=0), "fast confined total density", ndim=1, nonnegative=True)

    @property
    def fast_full_device_total_density_m3(self) -> np.ndarray:
        available = tuple(record.full_device_density_m3 for record in self.fast_species if record.full_device_density_m3.size)
        if not available:
            return _readonly_array(None, "fast full device total density")
        if len(available) != len(self.fast_species):
            raise MissingMetadata("full device fast density is unavailable for one or more active species")
        return _readonly_array(np.sum(available, axis=0), "fast full device total density", ndim=1, nonnegative=True)

    @property
    def fusion_total_axial_neutron_rate_s(self) -> np.ndarray:
        if not self.fusion_components:
            return _readonly_array(None, "fusion total axial neutron rate")
        return _readonly_array(np.sum(tuple(record.axial_neutron_rate_s for record in self.fusion_components), axis=0), "fusion total axial neutron rate", ndim=1, nonnegative=True)

    @property
    def fusion_total_axial_neutron_rate_density_m3_s(self) -> np.ndarray:
        if not self.fusion_components:
            return _readonly_array(None, "fusion total axial neutron rate density")
        return _readonly_array(np.sum(tuple(record.axial_neutron_rate_density_m3_s for record in self.fusion_components), axis=0), "fusion total axial neutron rate density", ndim=1, nonnegative=True)

def _build_grid(metadata: Mapping[str, Any], grid_id: str, *, scope: str, edges_key: str | None, centers_key: str, volumes_key: str | None) -> AxialGridRecord:
    centers = _readonly_array(metadata.get(centers_key), centers_key, ndim=1, required=True)
    edges = _readonly_array(None if edges_key is None else metadata.get(edges_key), str(edges_key), ndim=1)
    volumes = _readonly_array(None if volumes_key is None else metadata.get(volumes_key), str(volumes_key), ndim=1, nonnegative=True)
    if edges.size:
        _validate_centers_faces(centers, edges, grid_id)
    elif centers.size > 1 and np.any(np.diff(centers) <= 0.0):
        raise MissingMetadata(f"{centers_key} must be strictly increasing")
    if volumes.size and volumes.size != centers.size:
        raise MissingMetadata(f"{volumes_key} must match {centers_key}")
    return AxialGridRecord(grid_id=grid_id, scope=scope, edges_m=edges, centers_m=centers, cell_volumes_m3=volumes)

def _build_grids(metadata: Mapping[str, Any]) -> Mapping[str, AxialGridRecord]:
    grids = {
        "confined_kinetic": _build_grid(metadata, "confined_kinetic", scope="confined_throat_to_throat", edges_key="coordinate_edges_m", centers_key="coordinate_centers_m", volumes_key="cell_volumes_m3"),
        "full_device_population": _build_grid(metadata, "full_device_population", scope="full_device", edges_key="full_device_z_edges_m", centers_key="full_device_z_centers_m", volumes_key="full_device_cell_volumes_m3"),
        "magnetic_visual": _build_grid(metadata, "magnetic_visual", scope="full_device_visual", edges_key=None, centers_key="magnetic_field_visual_coordinate_m", volumes_key=None),
    }
    return MappingProxyType(grids)

def _build_beam_record(beam_id: str, raw: Mapping[str, Any], confined_grid: AxialGridRecord) -> BeamPlotRecord:
    species_id = str(raw.get("projectile_species", "")).strip()
    if not species_id:
        raise MissingMetadata(f"beam {beam_id!r} is missing projectile_species")
    ion_species(species_id)
    path_centers = _readonly_array(raw.get("path_center_points_m"), f"beam_deposition_by_id.{beam_id}.path_center_points_m", ndim=2, required=True, allow_nan=True)
    if path_centers.shape[1:] != (3,):
        raise MissingMetadata(f"beam {beam_id!r} path_center_points_m must have shape N by 3")
    if path_centers.shape[0] != confined_grid.centers_m.size:
        raise MissingMetadata(f"beam {beam_id!r} path_center_points_m must match the confined grid")
    start = _readonly_array(raw.get("start_m"), f"beam_deposition_by_id.{beam_id}.start_m", ndim=1, required=True)
    end = _readonly_array(raw.get("end_m"), f"beam_deposition_by_id.{beam_id}.end_m", ndim=1, required=True)
    direction = _readonly_array(raw.get("direction_unit"), f"beam_deposition_by_id.{beam_id}.direction_unit", ndim=1, required=True)
    if start.shape != (3,) or end.shape != (3,) or direction.shape != (3,):
        raise MissingMetadata(f"beam {beam_id!r} start, end, and direction must be length 3 vectors")
    if not np.isclose(float(np.linalg.norm(direction)), 1.0, rtol=0.0, atol=1.0e-12):
        raise MissingMetadata(f"beam {beam_id!r} direction must be a unit vector")
    profile_names = (
        "beam_radius_profile_m",
        "physical_path_lengths_m",
        "effective_path_lengths_m",
        "radial_overlap_fractions",
        "plasma_radius_profile_m",
        "birth_pitch_angle_profile_deg",
        "birth_lambda_profile",
        "axial_birth_rate_density_m3_s",
        "axial_birth_rate_s",
        "axial_birth_power_density_W_m3",
        "axial_birth_power_W",
    )
    profiles = {name: _readonly_array(raw.get(name), f"beam_deposition_by_id.{beam_id}.{name}", ndim=1, required=True, nonnegative=True) for name in profile_names}
    for name, profile in profiles.items():
        _validate_profile(profile, confined_grid, f"beam_deposition_by_id.{beam_id}.{name}")
    finite_path_rows = np.all(np.isfinite(path_centers), axis=1)
    nan_path_rows = np.all(np.isnan(path_centers), axis=1)
    if not np.all(finite_path_rows | nan_path_rows):
        raise MissingMetadata(f"beam {beam_id!r} path center rows must be fully finite or fully NaN")
    if not np.array_equal(finite_path_rows, profiles["physical_path_lengths_m"] > 0.0):
        raise MissingMetadata(f"beam {beam_id!r} finite path center rows must match positive physical path lengths")
    support = _readonly_array(raw.get("axial_birth_profile_support_mask"), f"beam_deposition_by_id.{beam_id}.axial_birth_profile_support_mask", ndim=1, dtype=bool, required=True)
    _validate_profile(support, confined_grid, f"beam_deposition_by_id.{beam_id}.axial_birth_profile_support_mask")
    component_ids = _string_tuple(raw.get("component_ids"), f"beam_deposition_by_id.{beam_id}.component_ids")
    if len(set(component_ids)) != len(component_ids):
        raise MissingMetadata(f"beam {beam_id!r} component ids contain duplicates")
    component_energies = _readonly_array(raw.get("component_energies_J"), f"beam_deposition_by_id.{beam_id}.component_energies_J", ndim=1, required=True, nonnegative=True)
    component_fractions = _readonly_array(raw.get("component_configured_power_fraction"), f"beam_deposition_by_id.{beam_id}.component_configured_power_fraction", ndim=1, required=True, nonnegative=True)
    component_power = _readonly_array(raw.get("component_deposited_birth_power_W"), f"beam_deposition_by_id.{beam_id}.component_deposited_birth_power_W", ndim=1, required=True, nonnegative=True)
    if not (len(component_ids) == component_energies.size == component_fractions.size == component_power.size):
        raise MissingMetadata(f"beam {beam_id!r} component arrays must have matching lengths")
    source = _readonly_array(raw.get("source_z_v_lambda_per_s"), f"beam_deposition_by_id.{beam_id}.source_z_v_lambda_per_s", ndim=3, required=True, nonnegative=True)
    if source.shape[0] != confined_grid.centers_m.size:
        raise MissingMetadata(f"beam {beam_id!r} source_z_v_lambda_per_s must match the confined grid")
    volume_source = _readonly_array(raw.get("volume_averaged_source_v_lambda_per_s"), f"beam_deposition_by_id.{beam_id}.volume_averaged_source_v_lambda_per_s", ndim=2, required=True, nonnegative=True)
    if source.shape[1:] != volume_source.shape:
        raise MissingMetadata(f"beam {beam_id!r} source arrays must use matching velocity grids")
    return BeamPlotRecord(
        beam_id=beam_id,
        species_id=species_id,
        enabled=bool(raw.get("enabled", True)),
        power_W=_required_nonnegative_float(raw.get("power_W"), f"beam_deposition_by_id.{beam_id}.power_W"),
        energy_keV=_required_nonnegative_float(raw.get("energy_keV"), f"beam_deposition_by_id.{beam_id}.energy_keV"),
        fast_birth_rate_s=_required_nonnegative_float(raw.get("fast_birth_rate_s"), f"beam_deposition_by_id.{beam_id}.fast_birth_rate_s"),
        deposited_birth_power_W=_required_nonnegative_float(raw.get("deposited_birth_power_W"), f"beam_deposition_by_id.{beam_id}.deposited_birth_power_W"),
        shine_through_power_W=_required_nonnegative_float(raw.get("shine_through_power_W"), f"beam_deposition_by_id.{beam_id}.shine_through_power_W"),
        beam_radius_m=_required_nonnegative_float(raw.get("beam_radius_m"), f"beam_deposition_by_id.{beam_id}.beam_radius_m"),
        configured_injection_angle_deg=_optional_float(raw.get("configured_injection_angle_deg"), f"beam_deposition_by_id.{beam_id}.configured_injection_angle_deg"),
        effective_injection_angle_deg=_optional_float(raw.get("effective_injection_angle_deg"), f"beam_deposition_by_id.{beam_id}.effective_injection_angle_deg"),
        beamline_axis_angle_deg=_optional_float(raw.get("beamline_axis_angle_deg"), f"beam_deposition_by_id.{beam_id}.beamline_axis_angle_deg"),
        injection_geometry_consistent=None if raw.get("injection_geometry_consistent") is None else bool(raw.get("injection_geometry_consistent")),
        component_ids=component_ids,
        component_energies_J=component_energies,
        component_configured_power_fraction=component_fractions,
        component_deposited_birth_power_W=component_power,
        start_m=start,
        end_m=end,
        direction_unit=direction,
        path_center_points_m=path_centers,
        beam_radius_profile_m=profiles["beam_radius_profile_m"],
        physical_path_lengths_m=profiles["physical_path_lengths_m"],
        effective_path_lengths_m=profiles["effective_path_lengths_m"],
        radial_overlap_fractions=profiles["radial_overlap_fractions"],
        plasma_radius_profile_m=profiles["plasma_radius_profile_m"],
        birth_pitch_angle_profile_deg=profiles["birth_pitch_angle_profile_deg"],
        birth_lambda_profile=profiles["birth_lambda_profile"],
        axial_grid_id=confined_grid.grid_id,
        axial_birth_rate_density_m3_s=profiles["axial_birth_rate_density_m3_s"],
        axial_birth_rate_s=profiles["axial_birth_rate_s"],
        axial_birth_profile_support_mask=support,
        axial_birth_power_density_W_m3=profiles["axial_birth_power_density_W_m3"],
        axial_birth_power_W=profiles["axial_birth_power_W"],
        source_z_v_lambda_per_s=source,
        volume_averaged_source_v_lambda_per_s=volume_source,
    )

def _build_beams(metadata: Mapping[str, Any], contract: Mapping[str, Any], confined_grid: AxialGridRecord) -> tuple[BeamPlotRecord, ...]:
    raw_records = _mapping(metadata.get("beam_deposition_by_id"), "beam_deposition_by_id", required=True)
    contract_order = _string_tuple(contract.get("canonical_beam_order"), "plotting_metadata_contract.canonical_beam_order")
    if not contract_order:
        raise MissingMetadata("plotting metadata contract canonical beam order must be nonempty")
    if len(set(contract_order)) != len(contract_order):
        raise MissingMetadata("plotting metadata contract canonical beam order contains duplicates")
    if set(contract_order) != set(raw_records):
        raise MissingMetadata("plotting metadata contract canonical beam order does not match the beam records")
    return tuple(_build_beam_record(beam_id, _mapping(raw_records[beam_id], f"beam_deposition_by_id.{beam_id}", required=True), confined_grid) for beam_id in contract_order)

def _build_beam_species_profiles(metadata: Mapping[str, Any], confined_grid: AxialGridRecord) -> Mapping[str, BeamSpeciesProfileRecord]:
    rate_density_map = _mapping(metadata.get("beam_axial_birth_rate_density_m3_s_by_species"), "beam_axial_birth_rate_density_m3_s_by_species", required=True)
    rate_map = _mapping(metadata.get("beam_axial_birth_rate_s_by_species"), "beam_axial_birth_rate_s_by_species", required=True)
    power_density_map = _mapping(metadata.get("beam_axial_birth_power_density_W_m3_by_species"), "beam_axial_birth_power_density_W_m3_by_species", required=True)
    power_map = _mapping(metadata.get("beam_axial_birth_power_W_by_species"), "beam_axial_birth_power_W_by_species", required=True)
    _require_matching_keys(rate_density_map, (("beam_axial_birth_rate_s_by_species", rate_map), ("beam_axial_birth_power_density_W_m3_by_species", power_density_map), ("beam_axial_birth_power_W_by_species", power_map)))
    species_ids = _ordered_ids(rate_density_map, PLOTTING_SPECIES_ORDER)
    result = {}
    for species_id in species_ids:
        ion_species(species_id)
        rate_density = _readonly_array(rate_density_map.get(species_id), f"beam_axial_birth_rate_density_m3_s_by_species.{species_id}", ndim=1, required=True, nonnegative=True)
        rate = _readonly_array(rate_map.get(species_id), f"beam_axial_birth_rate_s_by_species.{species_id}", ndim=1, required=True, nonnegative=True)
        power_density = _readonly_array(power_density_map.get(species_id), f"beam_axial_birth_power_density_W_m3_by_species.{species_id}", ndim=1, required=True, nonnegative=True)
        power = _readonly_array(power_map.get(species_id), f"beam_axial_birth_power_W_by_species.{species_id}", ndim=1, required=True, nonnegative=True)
        for name, profile in (("rate_density", rate_density), ("rate", rate), ("power_density", power_density), ("power", power)):
            _validate_profile(profile, confined_grid, f"beam species {species_id} {name}")
        result[species_id] = BeamSpeciesProfileRecord(species_id=species_id, axial_grid_id=confined_grid.grid_id, axial_birth_rate_density_m3_s=rate_density, axial_birth_rate_s=rate, axial_birth_power_density_W_m3=power_density, axial_birth_power_W=power)
    return MappingProxyType(result)

def _build_fast_species(metadata: Mapping[str, Any], confined_grid: AxialGridRecord, full_grid: AxialGridRecord) -> tuple[FastSpeciesPlotRecord, ...]:
    density_map = _mapping(metadata.get("fast_ion_local_density_m3_by_species"), "fast_ion_local_density_m3_by_species")
    if not density_map:
        return ()
    mass_map = _mapping(metadata.get("fast_ion_species_mass_kg_by_species"), "fast_ion_species_mass_kg_by_species", required=True)
    charge_map = _mapping(metadata.get("fast_ion_species_charge_number_by_species"), "fast_ion_species_charge_number_by_species", required=True)
    invariant_speed_map = _mapping(metadata.get("fast_ion_speed_grid_m_s_by_species"), "fast_ion_speed_grid_m_s_by_species", required=True)
    invariant_faces_map = _mapping(metadata.get("fast_ion_speed_grid_faces_m_s_by_species"), "fast_ion_speed_grid_faces_m_s_by_species", required=True)
    distribution_map = _mapping(metadata.get("fast_ion_distribution_v_lambda_by_species"), "fast_ion_distribution_v_lambda_by_species", required=True)
    local_speed_map = _mapping(metadata.get("fast_ion_local_speed_grid_m_s_by_species"), "fast_ion_local_speed_grid_m_s_by_species", required=True)
    local_faces_map = _mapping(metadata.get("fast_ion_local_speed_grid_faces_m_s_by_species"), "fast_ion_local_speed_grid_faces_m_s_by_species", required=True)
    local_lambda_map = _mapping(metadata.get("fast_ion_local_distribution_z_v_lambda_by_species"), "fast_ion_local_distribution_z_v_lambda_by_species", required=True)
    local_pitch_map = _mapping(metadata.get("fast_ion_local_distribution_z_v_pitch_by_species"), "fast_ion_local_distribution_z_v_pitch_by_species", required=True)
    full_distribution_map = _mapping(metadata.get("fast_ion_full_device_distribution_z_v_pitch_by_species"), "fast_ion_full_device_distribution_z_v_pitch_by_species")
    full_density_map = _mapping(metadata.get("fast_ion_full_device_density_m3_by_species"), "fast_ion_full_device_density_m3_by_species")
    _require_matching_keys(density_map, (("fast_ion_species_mass_kg_by_species", mass_map), ("fast_ion_species_charge_number_by_species", charge_map), ("fast_ion_speed_grid_m_s_by_species", invariant_speed_map), ("fast_ion_speed_grid_faces_m_s_by_species", invariant_faces_map), ("fast_ion_distribution_v_lambda_by_species", distribution_map), ("fast_ion_local_speed_grid_m_s_by_species", local_speed_map), ("fast_ion_local_speed_grid_faces_m_s_by_species", local_faces_map), ("fast_ion_local_distribution_z_v_lambda_by_species", local_lambda_map), ("fast_ion_local_distribution_z_v_pitch_by_species", local_pitch_map)))
    if set(full_distribution_map) != set(full_density_map):
        raise MissingMetadata("full device fast distribution and density species keys must match")
    if not set(full_distribution_map).issubset(density_map) or not set(full_density_map).issubset(density_map):
        raise MissingMetadata("full device fast species keys must be a subset of active fast species")
    lambda_centers = _readonly_array(metadata.get("lambda_centers"), "lambda_centers", ndim=1, required=True)
    lambda_faces = _readonly_array(metadata.get("lambda_faces"), "lambda_faces", ndim=1, required=True)
    pitch_centers_value = metadata.get("pitch_centers")
    pitch_faces_value = metadata.get("pitch_faces") if metadata.get("pitch_faces") is not None else metadata.get("pitch_edges")
    if pitch_centers_value is None and pitch_faces_value is None:
        pitch_counts = set()
        for mapping_name, mapping in (("fast_ion_local_distribution_z_v_pitch_by_species", local_pitch_map), ("fast_ion_full_device_distribution_z_v_pitch_by_species", full_distribution_map)):
            for species_id, values in mapping.items():
                array = np.asarray(values)
                if array.ndim != 3 or array.shape[2] < 1:
                    raise MissingMetadata(f"{mapping_name}.{species_id} cannot define the missing pitch grid")
                pitch_counts.add(int(array.shape[2]))
        if len(pitch_counts) != 1:
            raise MissingMetadata("missing pitch grid cannot be reconstructed consistently from the represented fast ion distributions")
        pitch_count = pitch_counts.pop()
        pitch_faces_value = np.linspace(-1.0, 1.0, pitch_count + 1)
        pitch_centers_value = 0.5 * (pitch_faces_value[:-1] + pitch_faces_value[1:])
    elif pitch_centers_value is None or pitch_faces_value is None:
        raise MissingMetadata("pitch_centers and pitch_faces must be supplied together")
    pitch_centers = _readonly_array(pitch_centers_value, "pitch_centers", ndim=1, required=True)
    pitch_faces = _readonly_array(pitch_faces_value, "pitch_faces", ndim=1, required=True)
    _validate_centers_faces(lambda_centers, lambda_faces, "lambda grid")
    _validate_centers_faces(pitch_centers, pitch_faces, "pitch grid")
    result = []
    for species_id in _ordered_ids(density_map, PLOTTING_SPECIES_ORDER):
        species = ion_species(species_id)
        invariant_speed = _readonly_array(invariant_speed_map.get(species_id), f"fast_ion_speed_grid_m_s_by_species.{species_id}", ndim=1, required=True, nonnegative=True)
        invariant_faces = _readonly_array(invariant_faces_map.get(species_id), f"fast_ion_speed_grid_faces_m_s_by_species.{species_id}", ndim=1, required=True, nonnegative=True)
        distribution = _readonly_array(distribution_map.get(species_id), f"fast_ion_distribution_v_lambda_by_species.{species_id}", ndim=2, required=True)
        local_speed = _readonly_array(local_speed_map.get(species_id), f"fast_ion_local_speed_grid_m_s_by_species.{species_id}", ndim=1, required=True, nonnegative=True)
        local_faces = _readonly_array(local_faces_map.get(species_id), f"fast_ion_local_speed_grid_faces_m_s_by_species.{species_id}", ndim=1, required=True, nonnegative=True)
        local_lambda = _readonly_array(local_lambda_map.get(species_id), f"fast_ion_local_distribution_z_v_lambda_by_species.{species_id}", ndim=3, required=True)
        local_pitch = _readonly_array(local_pitch_map.get(species_id), f"fast_ion_local_distribution_z_v_pitch_by_species.{species_id}", ndim=3, required=True)
        local_density = _readonly_array(density_map.get(species_id), f"fast_ion_local_density_m3_by_species.{species_id}", ndim=1, required=True, nonnegative=True)
        full_distribution = _readonly_array(full_distribution_map.get(species_id), f"fast_ion_full_device_distribution_z_v_pitch_by_species.{species_id}", ndim=3)
        full_density = _readonly_array(full_density_map.get(species_id), f"fast_ion_full_device_density_m3_by_species.{species_id}", ndim=1, nonnegative=True)
        _validate_centers_faces(invariant_speed, invariant_faces, f"fast species {species_id} invariant speed grid")
        _validate_centers_faces(local_speed, local_faces, f"fast species {species_id} local speed grid")
        if distribution.shape != (invariant_speed.size, lambda_centers.size):
            raise MissingMetadata(f"fast species {species_id} invariant distribution shape is inconsistent")
        if local_lambda.shape != (confined_grid.centers_m.size, local_speed.size, lambda_centers.size):
            raise MissingMetadata(f"fast species {species_id} local lambda distribution shape is inconsistent")
        if local_pitch.shape != (confined_grid.centers_m.size, local_speed.size, pitch_centers.size):
            raise MissingMetadata(f"fast species {species_id} local pitch distribution shape is inconsistent")
        _validate_profile(local_density, confined_grid, f"fast species {species_id} local density")
        if full_distribution.size and full_distribution.shape != (full_grid.centers_m.size, local_speed.size, pitch_centers.size):
            raise MissingMetadata(f"fast species {species_id} full device distribution shape is inconsistent")
        if full_density.size:
            _validate_profile(full_density, full_grid, f"fast species {species_id} full device density")
        result.append(FastSpeciesPlotRecord(
            species_id=species_id,
            symbol=species.symbol,
            mass_kg=_required_positive_float(mass_map.get(species_id), f"fast_ion_species_mass_kg_by_species.{species_id}"),
            charge_number=_required_positive_float(charge_map.get(species_id), f"fast_ion_species_charge_number_by_species.{species_id}"),
            confined_grid_id=confined_grid.grid_id,
            full_device_grid_id=full_grid.grid_id,
            invariant_speed_centers_m_s=invariant_speed,
            invariant_speed_faces_m_s=invariant_faces,
            lambda_centers=lambda_centers,
            lambda_faces=lambda_faces,
            distribution_v_lambda=distribution,
            local_speed_centers_m_s=local_speed,
            local_speed_faces_m_s=local_faces,
            pitch_centers=pitch_centers,
            pitch_faces=pitch_faces,
            local_distribution_z_v_lambda=local_lambda,
            local_distribution_z_v_pitch=local_pitch,
            local_density_m3=local_density,
            full_device_distribution_z_v_pitch=full_distribution,
            full_device_density_m3=full_density,
        ))
    return tuple(result)



def _build_electrostatic(metadata: Mapping[str, Any], confined_grid: AxialGridRecord) -> ElectrostaticPlotRecord | None:
    potential = _readonly_array(metadata.get("modal_phi_z_profile_V"), "modal_phi_z_profile_V", ndim=1)
    electron = _readonly_array(metadata.get("modal_phi_z_electron_density_profile_m3"), "modal_phi_z_electron_density_profile_m3", ndim=1, nonnegative=True)
    ion = _readonly_array(metadata.get("modal_phi_z_ion_density_profile_m3"), "modal_phi_z_ion_density_profile_m3", ndim=1, nonnegative=True)
    if potential.size == electron.size == ion.size == 0:
        return None
    if not (potential.size and electron.size and ion.size):
        return None
    for name, profile in (("modal_phi_z_profile_V", potential), ("modal_phi_z_electron_density_profile_m3", electron), ("modal_phi_z_ion_density_profile_m3", ion)):
        _validate_profile(profile, confined_grid, name)
    barrier = metadata.get("total_current_balance_wall_barrier_energy_J")
    if barrier is None:
        barrier = metadata.get("modal_wall_barrier_energy_J")
    normalized = metadata.get("total_current_balance_normalized_wall_barrier")
    if normalized is None:
        normalized = metadata.get("modal_wall_barrier_over_Te")
    return ElectrostaticPlotRecord(grid_id=confined_grid.grid_id, potential_relative_to_midplane_V=potential, electron_density_m3=electron, ion_density_m3=ion, wall_barrier_energy_J=_optional_nonnegative_float(barrier, "wall_barrier_energy_J"), normalized_wall_barrier=_optional_nonnegative_float(normalized, "normalized_wall_barrier"), converged=None if metadata.get("modal_phi_z_converged") is None else bool(metadata.get("modal_phi_z_converged")))

def _build_fusion_components(metadata: Mapping[str, Any], full_grid: AxialGridRecord) -> tuple[FusionComponentPlotRecord, ...]:
    pairs = _mapping(metadata.get("fusion_component_population_pairs"), "fusion_component_population_pairs")
    if not pairs:
        return ()
    kinds = _mapping(metadata.get("fusion_component_kind_by_label"), "fusion_component_kind_by_label", required=True)
    reactions = _mapping(metadata.get("fusion_component_reaction_by_label"), "fusion_component_reaction_by_label", required=True)
    rates = _mapping(metadata.get("fusion_component_rates_s"), "fusion_component_rates_s", required=True)
    neutron_rates = _mapping(metadata.get("fusion_component_neutron_rates_s"), "fusion_component_neutron_rates_s", required=True)
    axial_rates = _mapping(metadata.get("fusion_channel_axial_neutron_rate_s"), "fusion_channel_axial_neutron_rate_s", required=True)
    axial_density = _mapping(metadata.get("fusion_channel_axial_neutron_rate_density_m3_s"), "fusion_channel_axial_neutron_rate_density_m3_s", required=True)
    population_species = _mapping(metadata.get("fusion_population_species_by_id"), "fusion_population_species_by_id", required=True)
    _require_matching_keys(pairs, (("fusion_component_kind_by_label", kinds), ("fusion_component_reaction_by_label", reactions), ("fusion_component_rates_s", rates), ("fusion_component_neutron_rates_s", neutron_rates), ("fusion_channel_axial_neutron_rate_s", axial_rates), ("fusion_channel_axial_neutron_rate_density_m3_s", axial_density)))
    grid_id = str(metadata.get("fusion_component_axial_grid_id", full_grid.grid_id))
    if grid_id != full_grid.grid_id:
        raise MissingMetadata(f"unsupported fusion component axial grid {grid_id!r}")
    component_order = _string_tuple(metadata.get("active_fusion_channels"), "active_fusion_channels")
    if not component_order:
        component_order = tuple(sorted(pairs))
    if len(set(component_order)) != len(component_order) or set(component_order) != set(pairs):
        raise MissingMetadata("active fusion channel order does not match the fusion component records")
    result = []
    for label in component_order:
        population_ids = _string_tuple(pairs[label], f"fusion_component_population_pairs.{label}")
        if len(population_ids) != 2:
            raise MissingMetadata(f"fusion component {label!r} must contain two population ids")
        if population_ids[0] not in population_species or population_ids[1] not in population_species:
            raise MissingMetadata(f"fusion component {label!r} references an unknown population")
        reactant_a_species_id = _required_string(population_species[population_ids[0]], f"fusion_population_species_by_id.{population_ids[0]}")
        reactant_b_species_id = _required_string(population_species[population_ids[1]], f"fusion_population_species_by_id.{population_ids[1]}")
        ion_species(reactant_a_species_id)
        ion_species(reactant_b_species_id)
        rate_profile = _readonly_array(axial_rates.get(label), f"fusion_channel_axial_neutron_rate_s.{label}", ndim=1, required=True, nonnegative=True)
        density_profile = _readonly_array(axial_density.get(label), f"fusion_channel_axial_neutron_rate_density_m3_s.{label}", ndim=1, required=True, nonnegative=True)
        _validate_profile(rate_profile, full_grid, f"fusion_channel_axial_neutron_rate_s.{label}")
        _validate_profile(density_profile, full_grid, f"fusion_channel_axial_neutron_rate_density_m3_s.{label}")
        result.append(FusionComponentPlotRecord(
            label=label,
            kind=_required_string(kinds.get(label), f"fusion_component_kind_by_label.{label}"),
            reaction=_required_string(reactions.get(label), f"fusion_component_reaction_by_label.{label}"),
            reactant_a_population_id=population_ids[0],
            reactant_b_population_id=population_ids[1],
            reactant_a_species_id=reactant_a_species_id,
            reactant_b_species_id=reactant_b_species_id,
            axial_grid_id=grid_id,
            axial_neutron_rate_s=rate_profile,
            axial_neutron_rate_density_m3_s=density_profile,
            total_reaction_rate_s=_required_nonnegative_float(rates.get(label), f"fusion_component_rates_s.{label}"),
            total_neutron_rate_s=_required_nonnegative_float(neutron_rates.get(label), f"fusion_component_neutron_rates_s.{label}"),
        ))
    return tuple(result)

def _build_current_balance(metadata: Mapping[str, Any]) -> CurrentBalancePlotRecord | None:
    records = _mapping(metadata.get("total_current_balance_components"), "total_current_balance_components")
    if not records:
        return None
    currents = _mapping(metadata.get("total_current_balance_component_currents_A"), "total_current_balance_component_currents_A", required=True)
    if set(currents) != set(records):
        raise MissingMetadata("total current component current keys must match the component records")
    order = _ordered_ids(records, CURRENT_COMPONENT_ORDER)
    components = []
    for component_id in order:
        raw = _mapping(records[component_id], f"total_current_balance_components.{component_id}", required=True)
        species_id = _required_string(raw.get("species_id"), f"total_current_balance_components.{component_id}.species_id")
        ion_species(species_id)
        left_rate = _optional_nonnegative_float(raw.get("left_particle_loss_rate_s"), f"total_current_balance_components.{component_id}.left_particle_loss_rate_s")
        right_rate = _optional_nonnegative_float(raw.get("right_particle_loss_rate_s"), f"total_current_balance_components.{component_id}.right_particle_loss_rate_s")
        side_resolved = bool(raw.get("side_resolved", False))
        if side_resolved and (left_rate is None or right_rate is None):
            raise MissingMetadata(f"side resolved current component {component_id!r} requires left and right rates")
        components.append(CurrentComponentPlotRecord(
            component_id=component_id,
            population_kind=_required_string(raw.get("population_kind"), f"total_current_balance_components.{component_id}.population_kind"),
            species_id=species_id,
            charge_number=_required_positive_float(raw.get("charge_number"), f"total_current_balance_components.{component_id}.charge_number"),
            total_particle_loss_rate_s=_required_nonnegative_float(raw.get("total_particle_loss_rate_s"), f"total_current_balance_components.{component_id}.total_particle_loss_rate_s"),
            left_particle_loss_rate_s=left_rate,
            right_particle_loss_rate_s=right_rate,
            side_resolved=side_resolved,
            current_A=_required_nonnegative_float(currents.get(component_id), f"total_current_balance_component_currents_A.{component_id}"),
        ))
    return CurrentBalancePlotRecord(
        components=tuple(components),
        scope=_required_string(metadata.get("total_current_balance_scope"), "total_current_balance_scope"),
        converged=bool(metadata.get("total_current_balance_converged", False)),
        target_ion_current_A=_optional_nonnegative_float(metadata.get("total_current_balance_target_ion_current_A"), "total_current_balance_target_ion_current_A"),
        electron_current_A=_optional_nonnegative_float(metadata.get("total_current_balance_electron_current_A"), "total_current_balance_electron_current_A"),
        current_residual_A=_optional_float(metadata.get("total_current_balance_current_residual_A"), "total_current_balance_current_residual_A"),
        current_relative_residual=_optional_nonnegative_float(metadata.get("total_current_balance_current_relative_residual"), "total_current_balance_current_relative_residual"),
        wall_barrier_energy_J=_optional_nonnegative_float(metadata.get("total_current_balance_wall_barrier_energy_J"), "total_current_balance_wall_barrier_energy_J"),
        wall_barrier_energy_keV=_optional_nonnegative_float(metadata.get("total_current_balance_wall_barrier_energy_keV"), "total_current_balance_wall_barrier_energy_keV"),
        normalized_wall_barrier=_optional_nonnegative_float(metadata.get("total_current_balance_normalized_wall_barrier"), "total_current_balance_normalized_wall_barrier"),
        left_ion_current_A=_optional_nonnegative_float(metadata.get("total_current_balance_side_resolved_left_ion_current_A"), "total_current_balance_side_resolved_left_ion_current_A"),
        right_ion_current_A=_optional_nonnegative_float(metadata.get("total_current_balance_side_resolved_right_ion_current_A"), "total_current_balance_side_resolved_right_ion_current_A"),
        end_asymmetry_relative_error=_optional_nonnegative_float(metadata.get("total_current_balance_end_asymmetry_relative_error"), "total_current_balance_end_asymmetry_relative_error"),
        end_asymmetry_relative_tolerance=_optional_nonnegative_float(metadata.get("total_current_balance_end_asymmetry_relative_tolerance"), "total_current_balance_end_asymmetry_relative_tolerance"),
        end_asymmetry_applicable=bool(metadata.get("total_current_balance_end_asymmetry_applicable", False)),
        end_asymmetry_check_passed=bool(metadata.get("total_current_balance_end_asymmetry_check_passed", False)),
        component_identity_relative_error=_optional_nonnegative_float(metadata.get("total_current_balance_component_current_identity_relative_error"), "total_current_balance_component_current_identity_relative_error"),
        component_identity_check_passed=bool(metadata.get("total_current_balance_component_current_identity_check_passed", False)),
        side_sum_identity_relative_error=_optional_nonnegative_float(metadata.get("total_current_balance_side_sum_identity_relative_error"), "total_current_balance_side_sum_identity_relative_error"),
        side_sum_identity_applicable=bool(metadata.get("total_current_balance_side_sum_identity_applicable", False)),
        side_sum_identity_check_passed=bool(metadata.get("total_current_balance_side_sum_identity_check_passed", False)),
        represented_confined_plasma_balance_claimed=bool(metadata.get("global_current_balance_claimed", False)),
        full_device_wall_balance_claimed=False,
    )

def build_model_view(metadata: Mapping[str, Any]) -> PlottingModelView:
    """Build the validated immutable plotting view"""
    version_value = metadata.get("plotting_metadata_contract_version")
    try:
        version = int(version_value)
    except Exception as exc:
        raise MissingMetadata("plotting metadata contract is missing, rerun the source model with Plotting Pass 1") from exc
    if version != PLOTTING_METADATA_CONTRACT_VERSION:
        raise MissingMetadata(f"unsupported plotting metadata contract version {version}")
    raw_contract = _mapping(metadata.get("plotting_metadata_contract"), "plotting_metadata_contract", required=True)
    try:
        contract_version = int(raw_contract.get("version"))
    except Exception as exc:
        raise MissingMetadata("plotting metadata contract version is invalid") from exc
    if contract_version != version:
        raise MissingMetadata("plotting metadata contract versions do not agree")
    contract_species_order = _string_tuple(raw_contract.get("canonical_species_order"), "plotting_metadata_contract.canonical_species_order")
    if contract_species_order != PLOTTING_SPECIES_ORDER:
        raise MissingMetadata("plotting metadata contract species order is unsupported")
    contract = _frozen_value(raw_contract)
    grids = _build_grids(metadata)
    confined_grid = grids["confined_kinetic"]
    full_grid = grids["full_device_population"]
    beams = _build_beams(metadata, contract, confined_grid)
    beam_profiles = _build_beam_species_profiles(metadata, confined_grid)
    fast_species = _build_fast_species(metadata, confined_grid, full_grid)
    active_fast_species = _string_tuple(raw_contract.get("active_fast_species"), "plotting_metadata_contract.active_fast_species")
    if set(active_fast_species) != {record.species_id for record in fast_species}:
        raise MissingMetadata("active fast species do not match the species resolved plotting records")
    if {beam.species_id for beam in beams} != set(beam_profiles):
        raise MissingMetadata("beam species profiles do not match the represented beam species")
    electrostatic = _build_electrostatic(metadata, confined_grid)
    fusion_components = _build_fusion_components(metadata, full_grid)
    current_balance = _build_current_balance(metadata)
    return PlottingModelView(contract_version=version, contract=contract, grids=grids, beams=beams, beam_profiles_by_species=beam_profiles, fast_species=fast_species, electrostatic=electrostatic, fusion_components=fusion_components, current_balance=current_balance)

__all__ = ['PLOTTING_METADATA_CONTRACT_VERSION', 'PLOTTING_SPECIES_ORDER', 'CURRENT_COMPONENT_ORDER', 'AxialGridRecord', 'BeamPlotRecord', 'BeamSpeciesProfileRecord', 'FastSpeciesPlotRecord', 'ElectrostaticPlotRecord', 'FusionComponentPlotRecord', 'CurrentComponentPlotRecord', 'CurrentBalancePlotRecord', 'PlottingModelView', 'build_model_view']
