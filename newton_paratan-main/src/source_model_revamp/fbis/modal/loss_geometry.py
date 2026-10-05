"""Evaluate axial geometry factors used by the modal ion loss closure"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike

def _positive_scalar(value: float, name: str) -> float:
    """Return a finite positive scalar or raise for invalid input"""
    scalar = float(value)
    if not np.isfinite(scalar) or scalar <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return scalar

def loss_geometry_factor_from_midpoint_profile(zeta_faces: ArrayLike, B_tilde_midpoints: ArrayLike, mirror_ratio: float, density_ratio_midpoints: ArrayLike | None = None) -> tuple[float, float]:
    """Return the dimensionless Eq 62 geometry factor and flux tube volume integral
    
    Inputs use normalized axial coordinate faces and midpoint B divided by B0
    An optional Eq 42 density ratio multiplies the loss geometry integrand while the volume integral remains purely magnetic
    """
    faces = np.asarray(zeta_faces, dtype=float)
    field = np.asarray(B_tilde_midpoints, dtype=float)
    ratio = _positive_scalar(mirror_ratio, "mirror_ratio")

    if faces.ndim != 1 or field.ndim != 1 or faces.size != field.size + 1:
        raise ValueError("zeta_faces must have one more value than B_tilde_midpoints")
    if np.any(np.diff(faces) <= 0.0):
        raise ValueError("zeta_faces must be strictly increasing")
    if np.any(~np.isfinite(field)) or np.any(field <= 0.0):
        raise ValueError("B_tilde_midpoints must be positive and finite")
    if density_ratio_midpoints is None:
        density_ratio = np.ones_like(field)
    else:
        density_ratio = np.asarray(density_ratio_midpoints, dtype=float)
        if density_ratio.shape != field.shape:
            raise ValueError("density_ratio_midpoints must match B_tilde_midpoints")
        if np.any(~np.isfinite(density_ratio)) or np.any(density_ratio < 0.0):
            raise ValueError("density_ratio_midpoints must be finite and nonnegative")

    widths = np.diff(faces)
    geometry_integrand = (density_ratio * (1.0 / field) * np.sqrt(np.maximum(1.0 - field / ratio, 0.0)))
    volume_integrand = 1.0 / field
    geometry_factor = float(np.sum(widths * geometry_integrand))
    volume_integral = float(np.sum(widths * volume_integrand))
    if volume_integral <= 0.0:
        raise ValueError("volume geometry integral must be positive")

    return geometry_factor, volume_integral