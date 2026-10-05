"""Build normalized axial density profiles for Egedal Eq 42 orbit averaging

Profiles use the physical flux tube cell volume average so n(z) divided by average n has unit weighted mean
Input profiles are symmetrized about the magnetic midplane before they are supplied to the orbit averaged scattering operator
"""
from __future__ import annotations
from dataclasses import dataclass
from hashlib import sha256
from typing import Any, Callable
import numpy as np
from numpy.typing import ArrayLike

EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY = "self_consistent_quasineutral_density"
EQ42_PRESCRIBED_AXIAL_DENSITY_PROFILE = "prescribed_axial_density_profile"
EQ42_UNIFORM_DENSITY_REFERENCE = "uniform_density_reference"
_PROFILE_IDENTITY_QUANTIZATION = 1.0e-12
_EQ42_ALLOWED_MODELS = {EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, EQ42_PRESCRIBED_AXIAL_DENSITY_PROFILE, EQ42_UNIFORM_DENSITY_REFERENCE}

def _smooth_floating_point_noise(values: ArrayLike) -> np.ndarray:
    """Round coordinates to fourteen decimal places before profile node operations"""
    array = np.asarray(values, dtype=float)
    return np.round(array, decimals=14)

def _point_coords_into_bytes(values: ArrayLike) -> bytes:
    """Quantize finite profile values and return deterministic bytes for identity keys"""
    array = np.asarray(values, dtype=float)
    if np.any(~np.isfinite(array)):
        raise ValueError("profile identity values must be finite")
    quantized = np.rint(array / _PROFILE_IDENTITY_QUANTIZATION).astype(np.int64)

    return quantized.tobytes()

def _merge_near_duplicate_nodes(coordinates: ArrayLike, values: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    """Merge coordinates that differ only by floating point roundoff and average their values"""
    coordinate_array = np.asarray(coordinates, dtype=float)
    value_array = np.asarray(values, dtype=float)
    if coordinate_array.ndim != 1 or value_array.shape != coordinate_array.shape:
        raise ValueError("node coordinates and values must be matching 1D arrays")
    canonical = _smooth_floating_point_noise(coordinate_array)
    unique, inverse = np.unique(canonical, return_inverse=True)
    sums = np.zeros(unique.shape, dtype=float)
    counts = np.zeros(unique.shape, dtype=int)
    np.add.at(sums, inverse, value_array)
    np.add.at(counts, inverse, 1)
    averaged = sums / counts
    if np.any(np.diff(unique) <= 0.0):
        raise ValueError("node coordinates must be strictly increasing")
    
    return unique, averaged

@dataclass(frozen=True)
class Eq42DensityProfile:
    """Even Eq 42 density profile normalized by its physical flux tube volume average
    
    Cell densities use m^-3 and cell volumes use m^3
    Positive node arrays span normalized axial coordinate ζ from zero at the midplane to one at the throat
    """
    model: str
    source: str
    zeta_cells: np.ndarray
    density_cells_m3: np.ndarray
    cell_volumes_m3: np.ndarray
    zeta_positive_nodes: np.ndarray
    density_positive_nodes_m3: np.ndarray
    density_ratio_cells: np.ndarray
    density_ratio_positive_nodes: np.ndarray
    volume_average_density_m3: float
    raw_volume_average_density_m3: float
    symmetry_relative_error: float
    symmetry_tolerance: float
    symmetric_within_tolerance: bool
    profile_identity: str

    def __post_init__(self) -> None:
        """Validate profile shapes, units by convention, monotonic coordinates, and unit weighted normalization"""
        model = eq42_density_weighting_model(self.model)
        zeta_cells = np.asarray(self.zeta_cells, dtype=float)
        density_cells = np.asarray(self.density_cells_m3, dtype=float)
        volumes = np.asarray(self.cell_volumes_m3, dtype=float)
        zeta_nodes = np.asarray(self.zeta_positive_nodes, dtype=float)
        density_nodes = np.asarray(self.density_positive_nodes_m3, dtype=float)
        ratio_cells = np.asarray(self.density_ratio_cells, dtype=float)
        ratio_nodes = np.asarray(self.density_ratio_positive_nodes, dtype=float)
        if zeta_cells.ndim != 1 or density_cells.shape != zeta_cells.shape:
            raise ValueError("Eq 42 cell coordinates and densities must be matching 1D arrays")
        if volumes.shape != zeta_cells.shape:
            raise ValueError("Eq 42 cell volumes must match the cell coordinate array")
        if zeta_nodes.ndim != 1 or density_nodes.shape != zeta_nodes.shape:
            raise ValueError("Eq 42 positive node coordinates and densities must be matching 1D arrays")
        if ratio_cells.shape != zeta_cells.shape or ratio_nodes.shape != zeta_nodes.shape:
            raise ValueError("Eq 42 density ratio arrays must match their coordinate arrays")
        if zeta_cells.size < 2 or zeta_nodes.size < 2:
            raise ValueError("Eq 42 density profiles require at least two cells and two positive nodes")
        if np.any(np.diff(zeta_cells) <= 0.0) or np.any(np.diff(zeta_nodes) <= 0.0):
            raise ValueError("Eq 42 density coordinates must be strictly increasing")
        if not np.isclose(zeta_nodes[0], 0.0, rtol=0.0, atol=1.0e-14):
            raise ValueError("Eq 42 positive nodes must start at the magnetic midplane")
        if not np.isclose(zeta_nodes[-1], 1.0, rtol=0.0, atol=1.0e-12):
            raise ValueError("Eq 42 positive nodes must end at the magnetic throat")
        for name, values in (("density_cells_m3", density_cells), ("cell_volumes_m3", volumes), ("density_positive_nodes_m3", density_nodes), ("density_ratio_cells", ratio_cells), ("density_ratio_positive_nodes", ratio_nodes),):
            if np.any(~np.isfinite(values)):
                raise ValueError(f"{name} must contain only finite values")
        if np.any(density_cells < 0.0) or np.any(density_nodes < 0.0):
            raise ValueError("Eq 42 density profiles must be nonnegative")
        if np.any(ratio_cells < 0.0) or np.any(ratio_nodes < 0.0):
            raise ValueError("Eq 42 density ratios must be nonnegative")
        if np.any(volumes <= 0.0):
            raise ValueError("Eq 42 cell volumes must be positive")
        volume_average = float(self.volume_average_density_m3)
        raw_average = float(self.raw_volume_average_density_m3)
        if not np.isfinite(volume_average) or volume_average <= 0.0:
            raise ValueError("Eq 42 volume average density must be finite and positive")
        if not np.isfinite(raw_average) or raw_average <= 0.0:
            raise ValueError("Eq 42 raw volume average density must be finite and positive")
        weighted_ratio_average = float(np.sum(ratio_cells * volumes) / np.sum(volumes))
        if not np.isclose(weighted_ratio_average, 1.0, rtol=5.0e-13, atol=5.0e-13):
            raise ValueError("Eq 42 density ratio must have unit physical volume average")
        
        object.__setattr__(self, "model", model)
        object.__setattr__(self, "zeta_cells", zeta_cells)
        object.__setattr__(self, "density_cells_m3", density_cells)
        object.__setattr__(self, "cell_volumes_m3", volumes)
        object.__setattr__(self, "zeta_positive_nodes", zeta_nodes)
        object.__setattr__(self, "density_positive_nodes_m3", density_nodes)
        object.__setattr__(self, "density_ratio_cells", ratio_cells)
        object.__setattr__(self, "density_ratio_positive_nodes", ratio_nodes)
        object.__setattr__(self, "volume_average_density_m3", volume_average)
        object.__setattr__(self, "raw_volume_average_density_m3", raw_average)

    def density_ratio_at_zeta(self, zeta: float | ArrayLike) -> float | np.ndarray:
        """Return the even Eq 42 ratio n(z)/<n> on the normalized well"""
        values = np.asarray(zeta, dtype=float)
        if np.any(~np.isfinite(values)):
            raise ValueError("zeta must contain only finite values")
        absolute = np.abs(values)
        tolerance = 128.0 * np.finfo(float).eps
        result = np.interp(np.clip(absolute, 0.0, 1.0), self.zeta_positive_nodes, self.density_ratio_positive_nodes)
        if np.ndim(zeta) == 0:
            return float(result)
        
        return np.asarray(result, dtype=float)

    @property
    def density_ratio_function(self) -> Callable[[ArrayLike], float | np.ndarray]:
        """Return an array capable callable for Eq 42 bounce quadrature"""
        def evaluate(zeta: ArrayLike) -> float | np.ndarray:
            """Evaluate the stored even density ratio at normalized axial coordinates"""
            return self.density_ratio_at_zeta(zeta)

        return evaluate

    @property
    def nonuniform_weighting_active(self) -> bool:
        """Return whether the normalized density shape differs numerically from unity"""
        return bool(not np.allclose(self.density_ratio_positive_nodes, 1.0, rtol=1.0e-12, atol=1.0e-12,))

    @property
    def cache_key(self) -> tuple[str, bytes, bytes]:
        """Return a basis cache key containing the canonical model and normalized profile shape"""
        return (self.model, _point_coords_into_bytes(self.zeta_positive_nodes), _point_coords_into_bytes(self.density_ratio_positive_nodes),)

    def as_metadata(self) -> dict[str, Any]:
        """Return serializable profile normalization and symmetry metadata"""
        return {
            "model": self.model,
            "source": self.source,
            "volume_average_density_m3": self.volume_average_density_m3,
            "raw_volume_average_density_m3": self.raw_volume_average_density_m3,
            "zeta_positive_nodes": [float(value) for value in self.zeta_positive_nodes],
            "density_ratio_positive_nodes": [float(value) for value in self.density_ratio_positive_nodes],
            "density_ratio_cells": [float(value) for value in self.density_ratio_cells],
            "nonuniform_weighting_active": self.nonuniform_weighting_active,
            "symmetry_relative_error": self.symmetry_relative_error,
            "symmetry_tolerance": self.symmetry_tolerance,
            "symmetric_within_tolerance": self.symmetric_within_tolerance,
            "profile_identity": self.profile_identity,
            "normalization_measure": "physical_flux_tube_cell_volume_average",
            "eta_mapping_density_weighted": False,
        }

@dataclass(frozen=True)
class Eq42DensityProfileChange:
    """Difference metrics between two normalized Eq 42 density shapes"""
    volume_weighted_l2_relative_change: float
    maximum_absolute_ratio_change: float
    maximum_supported_relative_change: float
    support_floor: float

@dataclass(frozen=True)
class Eq42BasisChange:
    """Eigenvalue and eigenfunction changes between two density weighted physical bases"""
    maximum_eigenvalue_relative_change: float
    minimum_eigenfunction_overlap: float
    eigenvalue_relative_changes: np.ndarray
    eigenfunction_overlaps: np.ndarray

def eq42_density_weighting_model(value: str) -> str:
    """Return the canonical Eq 42 density weighting model identifier"""
    model = str(value).strip().lower()
    aliases = {"self_consistent": EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, "quasineutral": EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, "self_consistent_quasineutral": EQ42_SELF_CONSISTENT_QUASINEUTRAL_DENSITY, "prescribed": EQ42_PRESCRIBED_AXIAL_DENSITY_PROFILE, "prescribed_profile": EQ42_PRESCRIBED_AXIAL_DENSITY_PROFILE, "uniform": EQ42_UNIFORM_DENSITY_REFERENCE, "uniform_density": EQ42_UNIFORM_DENSITY_REFERENCE}
    model = aliases.get(model, model)
    if model not in _EQ42_ALLOWED_MODELS:
        raise ValueError("modal_eq42_density_weighting_model must be one of " f"{sorted(_EQ42_ALLOWED_MODELS)}, got {value!r}")
    return model

def _symmetry_relative_error(zeta: np.ndarray, density: np.ndarray) -> tuple[np.ndarray, float]:
    """Symmetrize a profile about ζ equals zero and return its maximum relative asymmetry"""
    mirrored = np.interp(-zeta, zeta, density)
    symmetric = 0.5 * (density + mirrored)
    scale = max(float(np.max(np.abs(symmetric))), np.finfo(float).tiny)
    error = float(np.max(np.abs(density - mirrored)) / scale)

    return symmetric, error

def _positive_nodes(*, zeta_nodes: np.ndarray, density_nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the nonnegative ζ nodes and their even averaged density values"""
    coordinates = np.unique(_smooth_floating_point_noise(np.concatenate(( np.asarray([0.0, 1.0], dtype=float), np.abs(np.asarray(zeta_nodes, dtype=float)),))))
    coordinates = coordinates[(coordinates >= 0.0) & (coordinates <= 1.0)]
    positive = np.interp(coordinates, zeta_nodes, density_nodes)
    negative = np.interp(-coordinates, zeta_nodes, density_nodes)

    return coordinates, 0.5 * (positive + negative)

def build_eq42_density_profile(*, model: str, source: str, zeta_cells: ArrayLike, density_cells_m3: ArrayLike, cell_volumes_m3: ArrayLike, zeta_nodes: ArrayLike | None = None, density_nodes_m3: ArrayLike | None = None, symmetry_tolerance: float = 1.0e-2) -> Eq42DensityProfile:
    """Build an even Eq 42 profile normalized by n=N/V

    The supplied cell profile must span both sides of the magnetic midplane
    The returned physical density scale is the symmetrized flux tube volume average in m^-3
    """
    canonical_model = eq42_density_weighting_model(model)
    zeta = np.asarray(zeta_cells, dtype=float)
    density = np.asarray(density_cells_m3, dtype=float)
    volumes = np.asarray(cell_volumes_m3, dtype=float)
    if zeta.ndim != 1 or density.shape != zeta.shape or volumes.shape != zeta.shape:
        raise ValueError("Eq 42 cell coordinates, density, and volumes must be matching 1D arrays")
    if np.any(~np.isfinite(zeta)) or np.any(np.diff(zeta) <= 0.0):
        raise ValueError("Eq 42 cell coordinates must be finite and strictly increasing")
    if np.any(~np.isfinite(density)) or np.any(density < 0.0):
        raise ValueError("Eq 42 cell density must be finite and nonnegative")
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("Eq 42 cell volumes must be finite and positive")
    if not (zeta[0] < 0.0 < zeta[-1]):
        raise ValueError("Eq 42 profile must span both sides of the magnetic midplane")
    tolerance = float(symmetry_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("Eq 42 symmetry tolerance must be finite and nonnegative")
    if zeta_nodes is None:
        nodes = np.unique(_smooth_floating_point_noise(np.concatenate(([-1.0, 0.0, 1.0], zeta))))
        node_density = np.interp(nodes, zeta, density)
    else:
        if density_nodes_m3 is None:
            raise ValueError("density_nodes_m3 is required when zeta_nodes is supplied")
        nodes = np.asarray(zeta_nodes, dtype=float)
        node_density = np.asarray(density_nodes_m3, dtype=float)
        if nodes.ndim != 1 or node_density.shape != nodes.shape:
            raise ValueError("Eq 42 node coordinates and density must be matching 1D arrays")
        if np.any(~np.isfinite(nodes)):
            raise ValueError("Eq 42 node coordinates must be finite")
        coordinate_tolerance = 64.0 * np.finfo(float).eps
        if np.any(np.diff(nodes) < -coordinate_tolerance):
            raise ValueError("Eq 42 node coordinates must be increasing")
        if np.any(~np.isfinite(node_density)) or np.any(node_density < 0.0):
            raise ValueError("Eq 42 node density must be finite and nonnegative")
        if nodes[0] > -1.0 + 1.0e-12 or nodes[-1] < 1.0 - 1.0e-12:
            raise ValueError("Eq 42 node profile must include both magnetic throats")
        nodes, node_density = _merge_near_duplicate_nodes(nodes, node_density)
    # Symmetrize before normalization so the coupled basis uses the even axial density assumed by the model
    symmetric_density, cell_symmetry_error = _symmetry_relative_error(zeta, density)
    symmetric_node_density, node_symmetry_error = _symmetry_relative_error(nodes, node_density)
    symmetry_error = max(cell_symmetry_error, node_symmetry_error)
    raw_average = float(np.sum(density * volumes) / np.sum(volumes))
    volume_average = float(np.sum(symmetric_density * volumes) / np.sum(volumes))
    if not np.isfinite(volume_average) or volume_average <= 0.0:
        raise ValueError("Eq 42 profile has zero or invalid volume average density")
    ratio_cells = symmetric_density / volume_average
    positive_zeta, positive_density = _positive_nodes(zeta_nodes=nodes, density_nodes=symmetric_node_density)
    ratio_nodes = positive_density / volume_average
    weighted_ratio_average = float(np.sum(ratio_cells * volumes) / np.sum(volumes))
    ratio_cells = ratio_cells / weighted_ratio_average
    ratio_nodes = ratio_nodes / weighted_ratio_average
    volume_average = volume_average * weighted_ratio_average
    digest = sha256()
    digest.update(_point_coords_into_bytes(positive_zeta))
    digest.update(_point_coords_into_bytes(ratio_nodes))
    identity = f"{canonical_model}:source={source}:shape_sha256={digest.hexdigest()}"

    return Eq42DensityProfile(
        model=canonical_model,
        source=str(source),
        zeta_cells=zeta,
        density_cells_m3=ratio_cells * volume_average,
        cell_volumes_m3=volumes,
        zeta_positive_nodes=positive_zeta,
        density_positive_nodes_m3=ratio_nodes * volume_average,
        density_ratio_cells=ratio_cells,
        density_ratio_positive_nodes=ratio_nodes,
        volume_average_density_m3=volume_average,
        raw_volume_average_density_m3=raw_average,
        symmetry_relative_error=symmetry_error,
        symmetry_tolerance=tolerance,
        symmetric_within_tolerance=bool(symmetry_error <= tolerance),
        profile_identity=identity,
    )

def uniform_eq42_density_profile(*, zeta_cells: ArrayLike, cell_volumes_m3: ArrayLike, reference_density_m3: float, model: str = EQ42_UNIFORM_DENSITY_REFERENCE, symmetry_tolerance: float = 1.0e-2, source: str = "Egedal_first_iteration_uniform_density_ratio") -> Eq42DensityProfile:
    """Build the n(z) / <n> = 1 reference used before density feedback"""
    density = float(reference_density_m3)
    if not np.isfinite(density) or density <= 0.0:
        raise ValueError("reference_density_m3 must be finite and positive")
    zeta = np.asarray(zeta_cells, dtype=float)
    nodes = np.unique(np.concatenate(([-1.0, 0.0, 1.0], zeta)))

    return build_eq42_density_profile(model=model, source=source, zeta_cells=zeta, density_cells_m3=np.full(zeta.shape, density, dtype=float), cell_volumes_m3=cell_volumes_m3, zeta_nodes=nodes, density_nodes_m3=np.full(nodes.shape, density, dtype=float), symmetry_tolerance=symmetry_tolerance)

def relax_eq42_density_profile(*, current: Eq42DensityProfile, target: Eq42DensityProfile, relaxation: float, model: str, source: str) -> Eq42DensityProfile:
    """Relax normalized density shape and logarithmically interpolate its positive physical density scale"""
    fraction = float(relaxation)
    if not np.isfinite(fraction) or not 0.0 < fraction <= 1.0:
        raise ValueError("Eq 42 density shape relaxation must satisfy 0 < value <= 1")
    if current.zeta_cells.shape != target.zeta_cells.shape or not np.allclose(current.zeta_cells, target.zeta_cells, rtol=0.0, atol=1.0e-13):
        raise ValueError("Eq 42 density profiles must use the same cell grid for relaxation")
    if not np.allclose(current.cell_volumes_m3, target.cell_volumes_m3, rtol=1.0e-13, atol=0.0):
        raise ValueError("Eq 42 density profiles must use the same cell volumes for relaxation")
    common_nodes = np.unique(_smooth_floating_point_noise(np.concatenate((current.zeta_positive_nodes, target.zeta_positive_nodes))))
    current_nodes = np.interp(common_nodes, current.zeta_positive_nodes, current.density_ratio_positive_nodes)
    target_nodes = np.interp(common_nodes, target.zeta_positive_nodes, target.density_ratio_positive_nodes)
    relaxed_cells = ((1.0 - fraction) * current.density_ratio_cells + fraction * target.density_ratio_cells)
    relaxed_nodes = (1.0 - fraction) * current_nodes + fraction * target_nodes
    # Interpolate the positive physical density scale geometrically while relaxing the normalized shape linearly
    density_scale = float(np.exp((1.0 - fraction) * np.log(current.volume_average_density_m3) + fraction * np.log(target.volume_average_density_m3)))
    full_nodes = np.concatenate((-common_nodes[:0:-1], common_nodes))
    full_density_nodes = np.concatenate((relaxed_nodes[:0:-1], relaxed_nodes)) * density_scale

    return build_eq42_density_profile(model=model, source=source, zeta_cells=current.zeta_cells, density_cells_m3=relaxed_cells * density_scale, cell_volumes_m3=current.cell_volumes_m3, zeta_nodes=full_nodes, density_nodes_m3=full_density_nodes, symmetry_tolerance=max(current.symmetry_tolerance, target.symmetry_tolerance))

def compare_eq42_density_profiles(current: Eq42DensityProfile, target: Eq42DensityProfile, *, support_relative_floor: float = 1.0e-8) -> Eq42DensityProfileChange:
    """Compare normalized density shapes on their shared physical cell grid"""
    if current.zeta_cells.shape != target.zeta_cells.shape or not np.allclose(current.zeta_cells, target.zeta_cells, rtol=0.0, atol=1.0e-13):
        raise ValueError("Eq 42 profile comparison requires a common cell grid")
    volumes = current.cell_volumes_m3
    difference = target.density_ratio_cells - current.density_ratio_cells
    numerator = float(np.sum(volumes * difference**2))
    denominator = float(np.sum(volumes * target.density_ratio_cells**2))
    l2 = float(np.sqrt(numerator / max(denominator, np.finfo(float).tiny)))
    maximum_absolute = float(np.max(np.abs(difference)))
    scale = max(float(np.max(np.abs(target.density_ratio_cells))), np.finfo(float).tiny)
    floor = float(support_relative_floor) * scale
    supported = np.abs(target.density_ratio_cells) >= floor
    maximum_relative = (float(np.max(np.abs(difference[supported]) / np.maximum(np.abs(target.density_ratio_cells[supported]), floor))) if np.any(supported) else 0.0)

    return Eq42DensityProfileChange(volume_weighted_l2_relative_change=l2, maximum_absolute_ratio_change=maximum_absolute, maximum_supported_relative_change=maximum_relative, support_floor=floor)
