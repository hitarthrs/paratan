"""Shared metadata and plotting helpers"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping
import numpy as np

METADATA_FILENAME = "source_model_metadata.json"
EV_TO_J = 1.602176634e-19
KEV_TO_J = 1.0e3 * EV_TO_J
MEV_TO_J = 1.0e6 * EV_TO_J
DEFAULT_LOCAL_B_OVER_B0 = (1.1, 3.0, 6.0, 15.5)

class MissingMetadata(ValueError):
    """Raised when a plot cannot be made from the available metadata"""

@dataclass(frozen=True)
class PlotSpec:
    name: str
    group: str
    filename: str
    function: Callable[[Mapping[str, Any], Path, argparse.Namespace], None]
    description: str

def _configure_matplotlib():
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    return plt

def _as_array(metadata: Mapping[str, Any], key: str, *, default: Any = None, ndim: int | None = None) -> np.ndarray:
    value = metadata.get(key, default)
    if value is None:
        arr = np.asarray([], dtype=float)
    else:
        arr = np.asarray(value, dtype=float)
    if ndim is not None and arr.size and arr.ndim != ndim:
        raise MissingMetadata(f"{key} must be {ndim}D, got shape {arr.shape}")
    return arr

def _scalar(metadata: Mapping[str, Any], key: str, default: float | None = None) -> float | None:
    value = metadata.get(key, default)
    if value is None:
        return None
    try:
        value = float(value)
    except Exception as exc:
        raise MissingMetadata(f"{key} must be numeric") from exc
    if not np.isfinite(value):
        return None
    return value

def _first_scalar(metadata: Mapping[str, Any], keys: Iterable[str]) -> float | None:
    for key in keys:
        value = _scalar(metadata, key)
        if value is not None:
            return value
    return None

def _centers_from_edges(edges: np.ndarray) -> np.ndarray:
    edges = np.asarray(edges, dtype=float)
    if edges.size < 2:
        return np.asarray([], dtype=float)
    return 0.5 * (edges[:-1] + edges[1:])

def _cm(values: Any) -> np.ndarray:
    return 100.0 * np.asarray(values, dtype=float)

def _m_to_cm(value: float | int) -> float:
    return 100.0 * float(value)

def _safe_divide(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    return np.divide(num, den, out=np.zeros_like(num, dtype=float), where=den != 0.0)

def _rate_density(rates: np.ndarray, volumes: np.ndarray) -> np.ndarray:
    if rates.size == 0 or volumes.size != rates.size:
        return np.asarray([], dtype=float)
    return _safe_divide(rates, volumes)

def _positive_floor(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    positive = values[np.isfinite(values) & (values > 0.0)]
    if positive.size == 0:
        return np.finfo(float).tiny
    return max(float(np.min(positive)) * 0.1, np.finfo(float).tiny)

def _log10_positive(values: np.ndarray) -> np.ndarray:
    floor = _positive_floor(values)
    return np.log10(np.maximum(np.asarray(values, dtype=float), floor))

def _z_edges_centers_m(metadata: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    edges = _as_array(metadata, "coordinate_edges_m", ndim=1)
    centers = _as_array(metadata, "coordinate_centers_m", ndim=1)
    if edges.size == 0:
        edges = _as_array(metadata, "z_edges_m", ndim=1)
    if edges.size == 0:
        edges_cm = _as_array(metadata, "z_edges_cm", ndim=1)
        if edges_cm.size:
            edges = edges_cm / 100.0
    if centers.size == 0 and edges.size >= 2:
        centers = _centers_from_edges(edges)
    if centers.size == 0:
        raise MissingMetadata("missing coordinate_centers_m or coordinate_edges_m")
    return edges, centers

def _matching_edges( metadata: Mapping[str, Any], keys: Iterable[str], bin_count: int) -> np.ndarray:
    for key in keys:
        edges = _as_array(metadata, key, ndim=1)
        if edges.size == int(bin_count) + 1:
            return edges
    return np.asarray([], dtype=float)

def _volumes_for_length( metadata: Mapping[str, Any], length: int, *, prefer_full_device: bool = False) -> np.ndarray:
    keys = ( ("full_device_cell_volumes_m3", "cell_volumes_m3") if prefer_full_device else ("cell_volumes_m3", "full_device_cell_volumes_m3"))
    for key in keys:
        volumes = _as_array(metadata, key, ndim=1)
        if volumes.size == int(length):
            return volumes
        
    return np.asarray([], dtype=float)

def _field_arrays(metadata: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z = _as_array(metadata, "magnetic_field_visual_coordinate_m", ndim=1)
    B = _as_array(metadata, "magnetic_field_visual_B_T", ndim=1)
    Bt = _as_array(metadata, "magnetic_field_visual_B_tilde", ndim=1)
    if z.size and B.size == z.size:
        if Bt.size != z.size:
            B0 = _scalar(metadata, "B0_T", 1.0) or 1.0
            Bt = B / B0
        return z, B, Bt
    _, zc = _z_edges_centers_m(metadata)
    Bt = _as_array(metadata, "B_tilde_centers", ndim=1)
    if Bt.size != zc.size:
        raise MissingMetadata("missing B_tilde_centers matching coordinate centers")
    B0 = _scalar(metadata, "B0_T", 1.0) or 1.0

    return zc, B0 * Bt, Bt

def _source_rate_matrix(metadata: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    matrix = _as_array(metadata, "source_rate_matrix_s")
    if matrix.size == 0:
        matrix = _as_array(metadata, "neutron_source_rate_matrix_s")
    if matrix.size == 0:
        raise MissingMetadata("missing source_rate_matrix_s")
    if matrix.ndim == 3:
        matrix = np.sum(matrix, axis=2)
    if matrix.ndim != 2:
        raise MissingMetadata("source_rate_matrix_s must be 2D or 3D")
    edges = _matching_edges(metadata, ("full_device_z_edges_m", "neutron_z_edges_m", "openmc_source_z_edges_m", "coordinate_edges_m", "z_edges_m"), matrix.shape[0])
    e_edges = _matching_edges(metadata, ("energy_edges_J",), matrix.shape[1])
    if e_edges.size == 0:
        e_edges_eV = _as_array(metadata, "energy_edges_eV", ndim=1)
        if e_edges_eV.size == matrix.shape[1] + 1:
            e_edges = e_edges_eV * EV_TO_J
    if e_edges.size == 0:
        e_edges_MeV = _as_array(metadata, "energy_edges_MeV", ndim=1)
        if e_edges_MeV.size == matrix.shape[1] + 1:
            e_edges = e_edges_MeV * MEV_TO_J
    if edges.size != matrix.shape[0] + 1:
        raise MissingMetadata("coordinate_edges_m does not match source_rate_matrix_s")
    if e_edges.size != matrix.shape[1] + 1:
        raise MissingMetadata("energy_edges do not match source_rate_matrix_s")
    
    return matrix, edges, e_edges

def _plasma_radius_profile(metadata: Mapping[str, Any], z_centers: np.ndarray) -> np.ndarray:
    B_tilde = _as_array(metadata, "B_tilde_centers", ndim=1)
    plasma_radius = _scalar(metadata, "plasma_radius_m")
    if B_tilde.size == z_centers.size and plasma_radius is not None:
        return plasma_radius / np.sqrt(np.maximum(B_tilde, np.finfo(float).tiny))
    if plasma_radius is not None:
        return np.full(z_centers.size, plasma_radius, dtype=float)
    
    return np.asarray([], dtype=float)

def _source_inner_radius_profile(z_centers: np.ndarray) -> np.ndarray:
    return np.zeros(z_centers.size, dtype=float)

def _vessel_radius_m(metadata: Mapping[str, Any]) -> float | None:
    for key in ( "beam_attenuation_vessel_radius_m", "source_geometry_linkage_vacuum_central_radius_m", "vacuum_vessel_radius_m"):
        value = _scalar(metadata, key)
        if value is not None and value > 0.0:
            return value
        
    return None

def _fusion_population_label(population_id: str) -> str:
    labels = {'fast_deuterium': 'Fast D', 'fast_tritium': 'Fast T'}
    return labels.get(str(population_id), str(population_id).replace("_", " ").title())

def _fusion_component_label(component: Any) -> str:
    left = _fusion_population_label(component.reactant_a_population_id)
    right = _fusion_population_label(component.reactant_b_population_id)
    return f"{left} + {right}"

def _beam_path_arrays(beam: Any) -> tuple[np.ndarray, np.ndarray]:
    centers = np.asarray(beam.path_center_points_m, dtype=float)
    radii = np.asarray(beam.beam_radius_profile_m, dtype=float)
    valid = np.all(np.isfinite(centers), axis=1)
    centers = centers[valid]
    if radii.size == valid.size:
        radii = radii[valid]
    if radii.size != centers.shape[0]:
        radii = np.full(centers.shape[0], beam.beam_radius_m, dtype=float)
    points = np.vstack((beam.start_m, centers, beam.end_m))
    point_radii = np.concatenate(([beam.beam_radius_m], radii, [beam.beam_radius_m]))
    coordinate = (points - beam.start_m) @ np.asarray(beam.direction_unit, dtype=float)
    order = np.argsort(coordinate)
    points = points[order]
    point_radii = point_radii[order]
    keep = np.ones(points.shape[0], dtype=bool)
    if points.shape[0] > 1:
        keep[1:] = np.linalg.norm(np.diff(points, axis=0), axis=1) > 1.0e-12
    return points[keep], point_radii[keep]


def _beam_projected_polygon(beam: Any, coordinate_index: int) -> np.ndarray:
    centers, radii = _beam_path_arrays(beam)
    points = centers[:, [2, coordinate_index]]
    tangents = np.gradient(points, axis=0)
    norms = np.linalg.norm(tangents, axis=1)
    fallback = points[-1] - points[0]
    fallback_norm = float(np.linalg.norm(fallback))
    if fallback_norm <= 0.0:
        return np.asarray([], dtype=float)
    tangents[norms <= 0.0] = fallback
    norms = np.linalg.norm(tangents, axis=1)
    normals = np.column_stack((-tangents[:, 1], tangents[:, 0])) / norms[:, None]
    upper = points + radii[:, None] * normals
    lower = points - radii[:, None] * normals
    return np.vstack((upper, lower[::-1]))

def _select_z_by_B(metadata: Mapping[str, Any], targets: Iterable[float]) -> list[int]:
    B_tilde = _as_array(metadata, "B_tilde_centers", ndim=1)
    if B_tilde.size == 0:
        raise MissingMetadata("missing B_tilde_centers")
    out: list[int] = []
    for target in targets:
        idx = int(np.argmin(np.abs(B_tilde - float(target))))
        if idx not in out:
            out.append(idx)

    return out

def _save(fig: Any, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt = _configure_matplotlib()
    plt.close(fig)
