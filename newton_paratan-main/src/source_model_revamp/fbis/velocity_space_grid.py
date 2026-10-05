"""
Velocity space grids and Jacobian helpers for FBIS collision models

Distribution convention
    n = integral f(v) d^3v

For an isotropic distribution f(v), this becomes

    n = 4π * integral_0^∞ v^2 f(v) dv.

For a gyrotropic distribution f(v, ξ), where ξ = v_parallel / v,

    n = 2π * integral_0^∞ v^2 dv * integral_-1^1 f(v, ξ) dξ.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike

def _as_1d_finite_array(values: ArrayLike, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{name} must be a 1D array")
    if array.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array

def _require_strictly_increasing(values: np.ndarray, name: str) -> None:
    if values.size < 2:
        raise ValueError(f"{name} must contain at least two points")
    if np.any(np.diff(values) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")

_PHYSICAL_DISTRIBUTION_RELATIVE_TOLERANCE = 1.0e-12

def _distribution_nonnegativity_metrics(distribution: ArrayLike, *, relative_tolerance: float = _PHYSICAL_DISTRIBUTION_RELATIVE_TOLERANCE, scale_floor: float = 0.0, name: str = "distribution") -> tuple[float, float, float, float, bool]:
    array = np.asarray(distribution, dtype=float)
    if array.size == 0:
        raise ValueError(f"{name} must not be empty")
    if np.any(~np.isfinite(array)):
        raise ValueError(f"{name} contains a nonfinite value")
    scale = max(float(np.max(np.abs(array))), float(scale_floor), np.finfo(float).tiny)
    tolerance = float(relative_tolerance) * scale
    minimum = float(np.min(array))
    negative_fraction = float(np.count_nonzero(array < -tolerance) / array.size)
    return scale, tolerance, minimum, negative_fraction, bool(minimum < -tolerance)

@dataclass(frozen=True)
class SpeedGrid:
    """Cell centered speed grid with spherical velocity space shell measures"""
    centers_m_s: np.ndarray
    faces_m_s: np.ndarray
    widths_m_s: np.ndarray
    shell_volumes_m3_s3: np.ndarray
    face_areas_m2_s2: np.ndarray
    def __post_init__(self) -> None:
        centers = _as_1d_finite_array(self.centers_m_s, "speed centers")
        faces = _as_1d_finite_array(self.faces_m_s, "speed faces")
        widths = _as_1d_finite_array(self.widths_m_s, "speed widths")
        shell_volumes = _as_1d_finite_array(self.shell_volumes_m3_s3, "speed shell volumes")
        face_areas = _as_1d_finite_array(self.face_areas_m2_s2, "speed face areas")
        _require_strictly_increasing(faces, "speed faces")
        if np.any(faces < 0.0):
            raise ValueError("speed faces must be nonnegative")
        if centers.size != faces.size - 1:
            raise ValueError("speed centers must have one fewer entry than speed faces")
        if widths.shape != centers.shape:
            raise ValueError("speed widths must have the same shape as speed centers")
        if shell_volumes.shape != centers.shape:
            raise ValueError("speed shell volumes must have the same shape as speed centers")
        if face_areas.shape != faces.shape:
            raise ValueError("speed face areas must have the same shape as speed faces")
        if np.any(centers <= 0.0):
            raise ValueError("speed centers must be positive")
        if np.any(widths <= 0.0):
            raise ValueError("speed widths must be positive")
        if np.any(shell_volumes <= 0.0):
            raise ValueError("speed shell volumes must be positive")
        if np.any(face_areas < 0.0):
            raise ValueError("speed face areas must be nonnegative")
        object.__setattr__(self, "centers_m_s", centers)
        object.__setattr__(self, "faces_m_s", faces)
        object.__setattr__(self, "widths_m_s", widths)
        object.__setattr__(self, "shell_volumes_m3_s3", shell_volumes)
        object.__setattr__(self, "face_areas_m2_s2", face_areas)

@dataclass(frozen=True)
class PitchGrid:
    """Cell centered pitch grid in ξ = v_parallel / v"""
    centers: np.ndarray
    faces: np.ndarray
    widths: np.ndarray
    def __post_init__(self) -> None:
        centers = _as_1d_finite_array(self.centers, "pitch centers")
        faces = _as_1d_finite_array(self.faces, "pitch faces")
        widths = _as_1d_finite_array(self.widths, "pitch widths")
        _require_strictly_increasing(faces, "pitch faces")
        if centers.size != faces.size - 1:
            raise ValueError("pitch centers must have one fewer entry than pitch faces")
        if widths.shape != centers.shape:
            raise ValueError("pitch widths must have the same shape as pitch centers")
        if np.any(widths <= 0.0):
            raise ValueError("pitch widths must be positive")
        if np.any(faces < -1.0) or np.any(faces > 1.0):
            raise ValueError("pitch faces must lie inside [-1, 1]")
        if np.any(centers < -1.0) or np.any(centers > 1.0):
            raise ValueError("pitch centers must lie inside [-1, 1]")
        object.__setattr__(self, "centers", centers)
        object.__setattr__(self, "faces", faces)
        object.__setattr__(self, "widths", widths)

def cell_centers_from_faces(faces: ArrayLike):
    """Cell centers from adjacent cell faces"""
    face_values = _as_1d_finite_array(faces, "cell faces")
    _require_strictly_increasing(face_values, "cell faces")
 
    return 0.5 * (face_values[:-1] + face_values[1:])

def cell_widths_from_faces(faces: ArrayLike):
    """Cell widths from adjacent cell faces"""
    face_values = _as_1d_finite_array(faces, "cell faces")
    _require_strictly_increasing(face_values, "cell faces")
  
    return face_values[1:] - face_values[:-1]

def speed_shell_volumes_from_faces(faces_m_s: ArrayLike):
    """Exact spherical shell measures 4π/3*(v_hi^3 - v_lo^3)"""
    faces = _as_1d_finite_array(faces_m_s, "speed faces")
    _require_strictly_increasing(faces, "speed faces")
    if np.any(faces < 0.0):
        raise ValueError("speed faces must be nonnegative")
    
    return (4.0 * np.pi / 3.0) * (faces[1:] ** 3 - faces[:-1] ** 3)

def speed_face_areas_from_faces(faces_m_s: ArrayLike):
    """Spherical velocity space face areas 4π*v_face^2"""
    faces = _as_1d_finite_array(faces_m_s, "speed faces")
    if np.any(faces < 0.0):
        raise ValueError("speed faces must be nonnegative")
   
    return 4.0 * np.pi * faces**2

def speed_grid_from_faces(faces_m_s: ArrayLike):
    """Build a speed grid from speed cell faces"""
    faces = _as_1d_finite_array(faces_m_s, "speed faces")
    _require_strictly_increasing(faces, "speed faces")
    if np.any(faces < 0.0):
        raise ValueError("speed faces must be nonnegative")
  
    return SpeedGrid(centers_m_s=0.5 * (faces[:-1] + faces[1:]), faces_m_s=faces, widths_m_s=faces[1:] - faces[:-1], shell_volumes_m3_s3=speed_shell_volumes_from_faces(faces), face_areas_m2_s2=speed_face_areas_from_faces(faces))

def uniform_speed_grid(max_speed_m_s: float, num_cells: int, min_speed_m_s: float = 0.0):
    """Uniform cell-centered speed grid"""
    if int(num_cells) < 1:
        raise ValueError("num_cells must be at least 1")
    if not np.isfinite(max_speed_m_s) or not np.isfinite(min_speed_m_s):
        raise ValueError("speed bounds must be finite")
    if min_speed_m_s < 0.0:
        raise ValueError("min_speed_m_s must be nonnegative")
    if max_speed_m_s <= min_speed_m_s:
        raise ValueError("max_speed_m_s must be greater than min_speed_m_s")
    faces = np.linspace(float(min_speed_m_s), float(max_speed_m_s), int(num_cells) + 1)
  
    return speed_grid_from_faces(faces)

def source_aligned_stretched_speed_grid(*, max_speed_m_s: float, num_cells: int, source_speeds_m_s: ArrayLike, core_cell_fraction: float = 0.7, tail_stretch_power: float = 2.0) -> SpeedGrid:
    """Build a source aligned speed grid with a high speed tail"""
    upper = float(max_speed_m_s)
    cell_count = int(num_cells)
    if cell_count < 2:
        raise ValueError("num_cells must be at least 2")
    if not np.isfinite(upper) or upper <= 0.0:
        raise ValueError("max_speed_m_s must be positive and finite")
    fraction = float(core_cell_fraction)
    if not np.isfinite(fraction) or not 0.0 < fraction < 1.0:
        raise ValueError("core_cell_fraction must lie inside (0, 1)")
    stretch = float(tail_stretch_power)
    if not np.isfinite(stretch) or stretch < 1.0:
        raise ValueError("tail_stretch_power must be finite and at least one")
    source = np.asarray(source_speeds_m_s, dtype=float).reshape(-1)
    if source.size == 0 or np.any(~np.isfinite(source)) or np.any(source <= 0.0):
        raise ValueError("source_speeds_m_s must contain positive finite values")
    tolerance = 64.0 * np.finfo(float).eps * max(upper, float(np.max(source)), 1.0)
    source = np.sort(source[source < upper - tolerance])
    if source.size == 0:
        return uniform_speed_grid(max_speed_m_s=upper, num_cells=cell_count)
    source = np.asarray([value for index, value in enumerate(source) if index == 0 or value - source[index - 1] > tolerance], dtype=float)
    largest_source = float(source[-1])
    mandatory_core_intervals = source.size
    core_cells = int(round(fraction * cell_count))
    core_cells = max(core_cells, mandatory_core_intervals)
    core_cells = min(core_cells, cell_count - 1)
    tail_cells = cell_count - core_cells
    core_breaks = np.concatenate(([0.0], source))
    lengths = np.diff(core_breaks)
    allocation = np.ones(lengths.size, dtype=int)
    remaining = core_cells - allocation.size
    if remaining > 0:
        ideal = remaining * lengths / max(float(np.sum(lengths)), np.finfo(float).tiny)
        additional = np.floor(ideal).astype(int)
        allocation += additional
        leftover = remaining - int(np.sum(additional))
        if leftover > 0:
            order = np.argsort(-(ideal - additional))
            allocation[order[:leftover]] += 1
    faces: list[float] = [0.0]
    for left, right, count in zip(core_breaks[:-1], core_breaks[1:], allocation, strict=True):
        segment = np.linspace(float(left), float(right), int(count) + 1)
        faces.extend(float(value) for value in segment[1:])
    faces[-1] = largest_source
    if tail_cells > 0:
        coordinate = np.linspace(0.0, 1.0, tail_cells + 1)[1:]
        tail = largest_source + (upper - largest_source) * coordinate**stretch
        faces.extend(float(value) for value in tail)
    else:
        faces[-1] = upper
    face_array = np.asarray(faces, dtype=float)
    face_array[0] = 0.0
    face_array[-1] = upper
    for value in source:
        nearest = int(np.argmin(np.abs(face_array - value)))
        face_array[nearest] = value

    return speed_grid_from_faces(face_array)

def pitch_grid_from_faces(faces: ArrayLike):
    """Build a pitch grid from ξ cell faces"""
    face_values = _as_1d_finite_array(faces, "pitch faces")
    _require_strictly_increasing(face_values, "pitch faces")
    if np.any(face_values < -1.0) or np.any(face_values > 1.0):
        raise ValueError("pitch faces must lie inside [-1, 1]")

    return PitchGrid(centers=0.5 * (face_values[:-1] + face_values[1:]), faces=face_values, widths=face_values[1:] - face_values[:-1])

def uniform_pitch_grid(num_cells: int, min_xi: float = -1.0, max_xi: float = 1.0):
    """Uniform pitch grid in ξ = v_parallel/v"""
    if int(num_cells) < 1:
        raise ValueError("num_cells must be at least 1")
    if not np.isfinite(min_xi) or not np.isfinite(max_xi):
        raise ValueError("pitch bounds must be finite")
    if min_xi < -1.0 or max_xi > 1.0:
        raise ValueError("pitch bounds must lie inside [-1, 1]")
    if max_xi <= min_xi:
        raise ValueError("max_xi must be greater than min_xi")
    faces = np.linspace(float(min_xi), float(max_xi), int(num_cells) + 1)
  
    return pitch_grid_from_faces(faces)

def gyrotropic_distribution_shape(speed_grid: SpeedGrid, pitch_grid: PitchGrid) -> tuple[int, int]:
    """Expected array shape for a gyrotropic f(v, ξ) distribution"""
    return (speed_grid.centers_m_s.size, pitch_grid.centers.size)

def require_gyrotropic_distribution_shape(distribution: ArrayLike, speed_grid: SpeedGrid, pitch_grid: PitchGrid, name: str = "gyrotropic distribution", allow_negative: bool = True) -> np.ndarray:
    """Return a validated gyrotropic distribution array on a speed/pitch grid"""
    f = np.asarray(distribution, dtype=float)
    expected = gyrotropic_distribution_shape(speed_grid, pitch_grid)
    if f.shape != expected:
        raise ValueError(f"{name} must have shape {expected}, got {f.shape}")
    if not np.all(np.isfinite(f)):
        raise ValueError(f"{name} must contain only finite values")
    if not allow_negative and np.any(f < 0.0):
        raise ValueError(f"{name} must be nonnegative")
    
    return f

def gyrotropic_velocity_cell_volumes(speed_grid: SpeedGrid, pitch_grid: PitchGrid):
    """Velocity space cell volumes for f(v, ξ): 2π * v^2 * dv * dξ"""
    return 0.5 * np.outer(speed_grid.shell_volumes_m3_s3, pitch_grid.widths)

def isotropic_maxwellian_distribution(speed_m_s: ArrayLike, density_m3: float, temperature_energy_J: float, particle_mass_kg: float):
    """Maxwellian f(v) normalized by n = integral f d^3v"""
    v = np.asarray(speed_m_s, dtype=float)
    if np.any(v < 0.0) or not np.all(np.isfinite(v)):
        raise ValueError("speed_m_s must be finite and nonnegative")
    if density_m3 < 0.0 or not np.isfinite(density_m3):
        raise ValueError("density_m3 must be finite and nonnegative")
    if temperature_energy_J <= 0.0 or not np.isfinite(temperature_energy_J):
        raise ValueError("temperature_energy_J must be positive and finite")
    if particle_mass_kg <= 0.0 or not np.isfinite(particle_mass_kg):
        raise ValueError("particle_mass_kg must be positive and finite")
    prefactor = density_m3 * (particle_mass_kg / (2.0 * np.pi * temperature_energy_J)) ** 1.5

    return prefactor * np.exp(-particle_mass_kg * v**2 / (2.0 * temperature_energy_J))