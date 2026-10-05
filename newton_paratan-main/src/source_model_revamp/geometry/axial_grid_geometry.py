"""Axial cell sizes, flux tube areas, and volume integrals

Use cell centered magnetic fields to assign volumes for beam and plasma profiles
"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike

def _finite_1d_array(name: str, values: ArrayLike) -> np.ndarray:
    """Convert input to a nonempty finite vector of floating point values"""
    array = np.asarray(values, dtype=float)
    if array.ndim != 1:
        raise ValueError(f"{name} must be 1D")
    if array.size == 0:
        raise ValueError(f"{name} must be non empty")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")

    return array

def cell_centers_from_faces(coordinate_faces_m: ArrayLike) -> np.ndarray:
    """Cell centers from strictly increasing axial cell faces"""
    faces = _finite_1d_array("coordinate_faces_m", coordinate_faces_m)
    if faces.size < 2:
        raise ValueError("coordinate_faces_m must contain at least two faces")
    if np.any(np.diff(faces) <= 0.0):
        raise ValueError("coordinate_faces_m must be strictly increasing")

    return 0.5 * (faces[:-1] + faces[1:])

def cell_widths_from_faces(coordinate_faces_m: ArrayLike) -> np.ndarray:
    """Cell widths from strictly increasing axial cell faces"""
    faces = _finite_1d_array("coordinate_faces_m", coordinate_faces_m)
    if faces.size < 2:
        raise ValueError("coordinate_faces_m must contain at least two faces")
    widths = faces[1:] - faces[:-1]
    if np.any(widths <= 0.0):
        raise ValueError("coordinate_faces_m must be strictly increasing")

    return widths

def flux_tube_area_profile_m2(B_tilde_centers: ArrayLike, reference_area_m2: float) -> np.ndarray:
    """Return flux tube areas from magnetic flux conservation

    reference_area_m2 is the area at B_tilde = 1
    """
    B_tilde = _finite_1d_array("B_tilde_centers", B_tilde_centers)
    if np.any(B_tilde <= 0.0):
        raise ValueError("B_tilde_centers must be positive")
    area0 = float(reference_area_m2)
    if not np.isfinite(area0) or area0 <= 0.0:
        raise ValueError("reference_area_m2 must be positive and finite")

    return area0 / B_tilde

def flux_tube_cell_volumes_m3(coordinate_faces_m: ArrayLike, B_tilde_centers: ArrayLike, reference_area_m2: float) -> np.ndarray:
    """Return cell volumes using the magnetic field at each cell center

    Each volume is the local flux tube area times the axial cell width
    """
    widths = cell_widths_from_faces(coordinate_faces_m)
    B_tilde = _finite_1d_array("B_tilde_centers", B_tilde_centers)
    if B_tilde.shape != widths.shape:
        raise ValueError("B_tilde_centers must have one value per axial cell")

    return flux_tube_area_profile_m2(B_tilde, reference_area_m2=reference_area_m2) * widths

def volume_integral(values: ArrayLike, cell_volumes_m3: ArrayLike) -> float:
    """Integral of a cell centered profile over reduced flux tube volume"""
    profile = _finite_1d_array("values", values)
    volumes = _finite_1d_array("cell_volumes_m3", cell_volumes_m3)
    if profile.shape != volumes.shape:
        raise ValueError("values and cell_volumes_m3 must have the same shape")
    if np.any(volumes <= 0.0):
        raise ValueError("cell_volumes_m3 must be positive")

    return float(np.sum(profile * volumes))
