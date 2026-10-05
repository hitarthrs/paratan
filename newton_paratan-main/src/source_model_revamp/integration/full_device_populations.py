"""Full device population grid assembly without extrapolating confined electron or fast ion states beyond their represented support"""
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Mapping
import numpy as np
from source_model_revamp.integration.pipeline_types import GeometryStageResult, KineticStageResult

@dataclass(frozen=True)
class FullDevicePopulationGrid:
    """
    Full device geometry and represented population arrays
    
    Fast distributions use shape `(n_z, n_speed, n_pitch)`
    The solved electron density is finite only on `electron_profile_support_mask` and is NaN outside the confined solved scope
    """
    z_edges_m: np.ndarray
    z_centers_m: np.ndarray
    B_tilde_centers: np.ndarray
    cell_volumes_m3: np.ndarray
    radial_outer_radius_m_by_z: np.ndarray
    radial_outer_radius_m_by_edge: np.ndarray
    radial_inner_radius_m_by_z: np.ndarray
    source_connected_radial_rho_max_edges: np.ndarray
    radial_source_profile_model: str | None
    fast_distribution_z_v_pitch: np.ndarray | None
    fast_distribution_z_v_pitch_by_species: Mapping[str, np.ndarray]
    solved_confined_electron_density_m3: np.ndarray | None
    electron_profile_support_mask: np.ndarray
    electron_profile_scope: str | None
    electron_outside_solved_scope_unavailable: bool
    includes_lost_fast_ions: bool
    includes_lost_fast_ions_by_species: Mapping[str, bool]

    @property
    def confined_fast_distribution_z_v_pitch(self) -> np.ndarray | None:
        """Return the legacy deuterium fast distribution alias"""
        return self.fast_distribution_z_v_pitch

    @property
    def expander_electron_profile_available(self) -> bool:
        """Return False because this mapping does not construct an expander electron population"""
        return False

    @property
    def electron_density_outside_solved_scope_available(self) -> bool:
        """Return False because the confined electron solution is not extrapolated outside its solved scope"""
        return False

def conservative_confined_distribution_on_grid(*, source_edges_m: np.ndarray, source_cell_volumes_m3: np.ndarray, source_distribution_z_v_pitch: np.ndarray, target_edges_m: np.ndarray, target_cell_volumes_m3: np.ndarray) -> np.ndarray:
    """
    Remap axial cell averaged values by physical overlap volume
    
    For each target cell the source cell averages are weighted by their overlap volume and divided by the target cell volume
    """
    source_edges = np.asarray(source_edges_m, dtype=float)
    source_volumes = np.asarray(source_cell_volumes_m3, dtype=float)
    local = np.asarray(source_distribution_z_v_pitch, dtype=float)
    target_edges = np.asarray(target_edges_m, dtype=float)
    target_volumes = np.asarray(target_cell_volumes_m3, dtype=float)
    if source_edges.size != source_volumes.size + 1 or local.shape[0] != source_volumes.size:
        raise ValueError("source distribution and geometry arrays have inconsistent shapes")
    if target_edges.size != target_volumes.size + 1:
        raise ValueError("target geometry arrays have inconsistent shapes")
    if np.any(np.diff(source_edges) <= 0.0) or np.any(np.diff(target_edges) <= 0.0):
        raise ValueError("source and target grid edges must be strictly increasing")
    if np.any(~np.isfinite(source_volumes)) or np.any(source_volumes <= 0.0):
        raise ValueError("source cell volumes must be positive and finite")
    if np.any(~np.isfinite(target_volumes)) or np.any(target_volumes <= 0.0):
        raise ValueError("target cell volumes must be positive and finite")
    if np.any(~np.isfinite(local)):
        raise ValueError("source distribution must be finite")
    scale = max(float(np.max(np.abs(local))) if local.size else 0.0, np.finfo(float).tiny)
    if float(np.min(local)) < -1.0e-12 * scale:
        raise ValueError("source distribution has material negative values")
    remapped = np.zeros((target_volumes.size, *local.shape[1:]), dtype=float)
    source_widths = np.diff(source_edges)
    for target_index, (target_left, target_right) in enumerate(zip(target_edges[:-1], target_edges[1:], strict=True)):
        overlap_left = np.maximum(source_edges[:-1], target_left)
        overlap_right = np.minimum(source_edges[1:], target_right)
        overlap_width = np.maximum(overlap_right - overlap_left, 0.0)
        # Scale each source cell volume by the fraction of its axial width inside the target cell
        overlap_volume = overlap_width * source_volumes / source_widths
        if np.any(overlap_volume > 0.0):
            remapped[target_index] = np.tensordot(overlap_volume, local, axes=(0, 0)) / target_volumes[target_index]

    return remapped

def _solved_electron_profile_on_population_grid(geometry: GeometryStageResult, kinetic: KineticStageResult, *, target_edges_m: np.ndarray, target_volumes_m3: np.ndarray) -> tuple[np.ndarray | None, np.ndarray, str | None]:
    """
    Map the solved confined electron density onto the requested population grid without extending its physical scope
    
    Cells outside the complete mapped support are returned as NaN and marked False in the support mask
    """

    target_edges = np.asarray(target_edges_m, dtype=float)
    target_volumes = np.asarray(target_volumes_m3, dtype=float)
    target_size = target_volumes.size
    state = getattr(kinetic, "operating_point_density_state", None)
    if state is None:
        return None, np.zeros(target_size, dtype=bool), None

    if bool(state.expander_electron_profile_available):
        raise ValueError("the operating point must not claim an expander electron profile")
    scope = str(state.electron_profile_scope).strip()
    if scope != "confined_throat_to_throat":
        raise ValueError("the solved electron profile scope must be confined_throat_to_throat")

    source_edges = np.asarray(geometry.z_edges_m, dtype=float)
    source_volumes = np.asarray(geometry.cell_volumes_m3, dtype=float)
    source_density = np.asarray(state.electron_cell_density_m3, dtype=float)
    source_support = np.asarray(state.profile_support_mask)
    state_zeta = np.asarray(state.zeta_cells, dtype=float)
    geometry_zeta = np.asarray(geometry.zeta_centers, dtype=float)
    if source_density.shape != source_volumes.shape:
        raise ValueError("the solved electron profile must match the confined geometry grid")
    if source_support.shape != source_volumes.shape:
        raise ValueError("the electron profile support mask must match the confined geometry grid")
    if not np.issubdtype(source_support.dtype, np.bool_):
        raise ValueError("the electron profile support mask must contain booleans")
    if (
        state_zeta.shape != geometry_zeta.shape
        or np.any(~np.isfinite(state_zeta))
        or not np.allclose(state_zeta, geometry_zeta, rtol=0.0, atol=1.0e-12)
    ):
        raise ValueError("the operating point density state must match the confined geometry coordinates")
    if np.any(~np.isfinite(source_density)) or np.any(source_density < 0.0):
        raise ValueError("the solved electron profile must be finite and nonnegative on its confined grid")

    if np.array_equal(source_edges, target_edges) and np.array_equal(source_volumes, target_volumes):
        mapped_density = source_density.copy()
        mapped_support = source_support.astype(bool, copy=True)
    else:
        mapped_density = conservative_confined_distribution_on_grid(
            source_edges_m=source_edges,
            source_cell_volumes_m3=source_volumes,
            source_distribution_z_v_pitch=source_density,
            target_edges_m=target_edges,
            target_cell_volumes_m3=target_volumes,
        )
        mapped_scope_weight = conservative_confined_distribution_on_grid(
            source_edges_m=source_edges,
            source_cell_volumes_m3=source_volumes,
            source_distribution_z_v_pitch=np.ones(source_volumes.shape, dtype=float),
            target_edges_m=target_edges,
            target_cell_volumes_m3=target_volumes,
        )
        mapped_support_weight = conservative_confined_distribution_on_grid(
            source_edges_m=source_edges,
            source_cell_volumes_m3=source_volumes,
            source_distribution_z_v_pitch=source_support.astype(float),
            target_edges_m=target_edges,
            target_cell_volumes_m3=target_volumes,
        )
        edge_scale = max(float(np.max(np.abs(target_edges))), float(np.max(np.abs(source_edges))), 1.0)
        edge_tolerance = 64.0 * np.finfo(float).eps * edge_scale
        inside_confined_scope = (
            (target_edges[:-1] >= source_edges[0] - edge_tolerance)
            & (target_edges[1:] <= source_edges[-1] + edge_tolerance)
        )
        mapped_support = (
            inside_confined_scope
            & (mapped_scope_weight > 0.0)
            & np.isclose(mapped_support_weight, mapped_scope_weight, rtol=1.0e-12, atol=1.0e-15)
        )

    # Keep the unsolved electron region explicit rather than extrapolating the confined closure
    electron = np.full(target_size, np.nan, dtype=float)
    electron[mapped_support] = mapped_density[mapped_support]

    return electron, mapped_support, scope

def build_full_device_population_grid(geometry: GeometryStageResult, kinetic: KineticStageResult) -> FullDevicePopulationGrid:
    """
    Build the common full device geometry and represented fast and electron populations
    
    Explicit full device fast ion distributions are preferred when present
    Otherwise confined local distributions are conservatively remapped and remain zero outside their confined source support
    """
    has_full_grid = all(value is not None for value in (geometry.full_device_z_edges_m, geometry.full_device_z_centers_m, geometry.full_device_B_tilde_centers, geometry.full_device_cell_volumes_m3))
    if has_full_grid:
        edges = np.asarray(geometry.full_device_z_edges_m, dtype=float)
        centers = np.asarray(geometry.full_device_z_centers_m, dtype=float)
        B_tilde = np.asarray(geometry.full_device_B_tilde_centers, dtype=float)
        volumes = np.asarray(geometry.full_device_cell_volumes_m3, dtype=float)
    else:
        edges = np.asarray(geometry.z_edges_m, dtype=float)
        centers = np.asarray(geometry.z_centers_m, dtype=float)
        B_tilde = np.asarray(geometry.B_tilde_centers, dtype=float)
        volumes = np.asarray(geometry.cell_volumes_m3, dtype=float)
    if edges.size != centers.size + 1 or B_tilde.shape != centers.shape or volumes.shape != centers.shape:
        raise ValueError("full device geometry arrays have inconsistent shapes")
    fast_distribution_by_species: dict[str, np.ndarray] = {}
    lost_population_by_species: dict[str, bool] = {}
    full_device_mapping = getattr(kinetic, "full_device_local_distribution_z_v_pitch_by_species", None)
    local_mapping = getattr(kinetic, "local_distribution_z_v_pitch_by_species", None)
    full_device_metadata = kinetic.metadata.get("species_full_device_population_metadata", {})
    if full_device_mapping is not None:
        for species_id, values in full_device_mapping.items():
            distribution = np.asarray(values, dtype=float)
            if distribution.shape[0] != centers.size:
                raise ValueError(f"full device fast {species_id} distribution does not match the full device grid")
            if np.any(~np.isfinite(distribution)):
                raise ValueError(f"full device fast {species_id} distribution must be finite")
            scale = max(float(np.max(np.abs(distribution))) if distribution.size else 0.0, np.finfo(float).tiny)
            if float(np.min(distribution)) < -1.0e-12 * scale:
                raise ValueError(f"full device fast {species_id} distribution has material negative values")
            fast_distribution_by_species[str(species_id)] = distribution
            metadata = full_device_metadata.get(species_id, {}) if isinstance(full_device_metadata, Mapping) else {}
            flag = metadata.get("full_device_fast_ion_distribution_includes_directed_lost", kinetic.metadata.get("full_device_fast_ion_distribution_includes_directed_lost", False))
            if not isinstance(flag, (bool, np.bool_)):
                raise ValueError("full device lost ion population metadata must be a boolean")
            lost_population_by_species[str(species_id)] = bool(flag)
    elif local_mapping is not None:
        for species_id, values in local_mapping.items():
            local = np.asarray(values, dtype=float)
            if local.shape[0] != geometry.z_centers_m.size:
                raise ValueError(f"confined local fast {species_id} distribution does not match the confined geometry grid")
            fast_distribution_by_species[str(species_id)] = conservative_confined_distribution_on_grid(source_edges_m=np.asarray(geometry.z_edges_m, dtype=float), source_cell_volumes_m3=np.asarray(geometry.cell_volumes_m3, dtype=float), source_distribution_z_v_pitch=local, target_edges_m=edges, target_cell_volumes_m3=volumes)
            lost_population_by_species[str(species_id)] = False
    fast_distribution = fast_distribution_by_species.get("deuterium")
    if fast_distribution is None:
        full_device_distribution = getattr(kinetic, "full_device_local_distribution_z_v_pitch", None)
        if full_device_distribution is not None:
            fast_distribution = np.asarray(full_device_distribution, dtype=float)
            fast_distribution_by_species["deuterium"] = fast_distribution
            lost_population_by_species["deuterium"] = bool(kinetic.metadata.get("full_device_fast_ion_distribution_includes_directed_lost", False))
        elif kinetic.local_distribution_z_v_pitch is not None:
            fast_distribution = conservative_confined_distribution_on_grid(source_edges_m=np.asarray(geometry.z_edges_m, dtype=float), source_cell_volumes_m3=np.asarray(geometry.cell_volumes_m3, dtype=float), source_distribution_z_v_pitch=np.asarray(kinetic.local_distribution_z_v_pitch, dtype=float), target_edges_m=edges, target_cell_volumes_m3=volumes)
            fast_distribution_by_species["deuterium"] = fast_distribution
            lost_population_by_species["deuterium"] = False
    includes_lost_fast_ions = any(lost_population_by_species.values())
    electron, electron_support, electron_scope = _solved_electron_profile_on_population_grid(geometry, kinetic, target_edges_m=edges, target_volumes_m3=volumes)
    radial_outer = float(geometry.plasma_radius_m) / np.sqrt(B_tilde)
    if geometry.B_T_function is None:
        radial_outer_edges = np.interp(edges, centers, radial_outer, left=radial_outer[0], right=radial_outer[-1])
    else:
        B_edges = np.asarray(geometry.B_T_function(edges), dtype=float) / float(geometry.B0_T)
        if B_edges.shape != edges.shape or np.any(~np.isfinite(B_edges)) or np.any(B_edges <= 0.0):
            raise ValueError("full device edge magnetic field must be positive and finite")
        radial_outer_edges = float(geometry.plasma_radius_m) / np.sqrt(B_edges)
    radial_inner = np.zeros_like(radial_outer)
    radial_rho_raw = kinetic.metadata.get("expander_source_connected_radial_rho_limit_edges_full_device")
    if radial_rho_raw is None:
        radial_rho_edges = np.ones(edges.size, dtype=float)
        radial_profile_model = None
    else:
        radial_rho_edges = np.asarray(radial_rho_raw, dtype=float)
        if radial_rho_edges.shape != edges.shape:
            raise ValueError("full device source connected radial support must match the full device edges")
        if np.any(~np.isfinite(radial_rho_edges)) or np.any(radial_rho_edges < 0.0) or np.any(radial_rho_edges > 1.0):
            raise ValueError("full device source connected radial support must lie inside zero to one")
        radial_profile_model = str(kinetic.metadata.get("expander_radial_profile_model", "")).strip() or None

    return FullDevicePopulationGrid(z_edges_m=edges, z_centers_m=centers, B_tilde_centers=B_tilde, cell_volumes_m3=volumes, radial_outer_radius_m_by_z=radial_outer, radial_outer_radius_m_by_edge=radial_outer_edges, radial_inner_radius_m_by_z=radial_inner, source_connected_radial_rho_max_edges=radial_rho_edges, radial_source_profile_model=radial_profile_model, fast_distribution_z_v_pitch=fast_distribution, fast_distribution_z_v_pitch_by_species=fast_distribution_by_species, solved_confined_electron_density_m3=electron, electron_profile_support_mask=electron_support, electron_profile_scope=electron_scope, electron_outside_solved_scope_unavailable=True, includes_lost_fast_ions=includes_lost_fast_ions, includes_lost_fast_ions_by_species=lost_population_by_species)
