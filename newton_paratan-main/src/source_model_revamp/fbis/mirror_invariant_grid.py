"""Mirror invariant grid helpers for FBIS"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid

def _as_1d_finite_array(values: ArrayLike, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{name} must be 1D")
    if array.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
 
    return array

def _validate_lambda_faces(faces: ArrayLike) -> np.ndarray:
    face_values = _as_1d_finite_array(faces, "Lambda faces")
    if face_values.size < 2:
        raise ValueError("Lambda grid requires at least two faces")
    if np.any(np.diff(face_values) <= 0.0):
        raise ValueError("Lambda faces must be strictly increasing")
    if face_values[0] < 0.0:
        raise ValueError("Lambda faces must be nonnegative")
    if face_values[-1] > 1.0:
        raise ValueError("Lambda faces must not exceed the magnetic invariant limit Lambda=1")
   
    return face_values

@dataclass(frozen=True)
class LambdaGrid:
    """Cell centered magnetic invariant Lambda grid"""
    centers: np.ndarray
    faces: np.ndarray
    widths: np.ndarray
    def __post_init__(self) -> None:
        faces = _validate_lambda_faces(self.faces)
        centers = _as_1d_finite_array(self.centers, "Lambda centers")
        widths = _as_1d_finite_array(self.widths, "Lambda widths")
        expected_cells = faces.size - 1
        if centers.size != expected_cells or widths.size != expected_cells:
            raise ValueError("Lambda grid centers and widths must have one value per cell")
        expected_centers = 0.5 * (faces[:-1] + faces[1:])
        expected_widths = faces[1:] - faces[:-1]
        if not np.allclose(centers, expected_centers, rtol=0.0, atol=1.0e-14):
            raise ValueError("Lambda centers must equal adjacent face midpoints")
        if not np.allclose(widths, expected_widths, rtol=0.0, atol=1.0e-14):
            raise ValueError("Lambda widths must equal adjacent face differences")
        if np.any(widths <= 0.0):
            raise ValueError("Lambda widths must be positive")

def lambda_grid_from_faces(faces: ArrayLike) -> LambdaGrid:
    """Build a cell centered Lambda grid from cell faces"""
    face_values = _validate_lambda_faces(faces)
    return LambdaGrid(centers=0.5 * (face_values[:-1] + face_values[1:]), faces=face_values, widths=face_values[1:] - face_values[:-1])

def uniform_lambda_grid(num_cells: int, lambda_min: float = 0.0, lambda_max: float = 1.0) -> LambdaGrid:
    """Uniform Lambda grid over a magneticn invariant interval"""
    count = int(num_cells)
    if count < 1:
        raise ValueError("num_cells must be at least 1")
    if not np.isfinite(lambda_min) or not np.isfinite(lambda_max):
        raise ValueError("Lambda bounds must be finite")
    if float(lambda_max) <= float(lambda_min):
        raise ValueError("lambda_max must be greater than lambda_min")
    faces = np.linspace(float(lambda_min), float(lambda_max), count + 1)
  
    return lambda_grid_from_faces(faces)

def mirror_distribution_shape(speed_grid: SpeedGrid, lambda_grid: LambdaGrid) -> tuple[int, int]:
    """Expected shape for a mirror distribution F(v, Lambda)"""
    return (speed_grid.centers_m_s.size, lambda_grid.centers.size)

def require_mirror_distribution_shape(distribution_v_lambda: ArrayLike, speed_grid: SpeedGrid, lambda_grid: LambdaGrid, name: str = "mirror distribution") -> np.ndarray:
    """Return a finite mirror distribution after checking its grid shape"""
    expected = mirror_distribution_shape(speed_grid, lambda_grid)
    distribution = np.asarray(distribution_v_lambda, dtype=float)
    if distribution.shape != expected:
        raise ValueError(f"{name} must have shape {expected}, got {distribution.shape}")
    if not np.all(np.isfinite(distribution)):
        raise ValueError(f"{name} must contain only finite values")
   
    return distribution