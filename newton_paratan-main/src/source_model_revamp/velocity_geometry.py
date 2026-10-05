"""
Shared gyrotropic velocity coordinate geometry helpers

File provides relative speed geometry for fusion and reactivity integrals
"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike

def pitch_perpendicular_fraction(pitch: ArrayLike):
    """Return sqrt(1 - ξ^2) for pitch ξ = v_parallel / v"""
    xi = np.asarray(pitch, dtype=float)
  
    return np.sqrt(np.maximum(1.0 - xi**2, 0.0))

def gyrotropic_cosine_between_velocities(evaluation_pitch: ArrayLike, source_pitch: ArrayLike, gyrophase_angle_rad: ArrayLike):
    """Cosine of the angle between two gyrotropic velocity directions"""
    xi = np.asarray(evaluation_pitch, dtype=float)
    mu = np.asarray(source_pitch, dtype=float)
    phi = np.asarray(gyrophase_angle_rad, dtype=float)

    return (xi * mu + pitch_perpendicular_fraction(xi) * pitch_perpendicular_fraction(mu) * np.cos(phi))

def relative_speed_m_s(evaluation_speed_m_s: ArrayLike, evaluation_pitch: ArrayLike, source_speed_m_s: ArrayLike, source_pitch: ArrayLike, gyrophase_angle_rad: ArrayLike):
    """Return |v - u| for gyrotropic speed/pitch/gyrophase coordinates"""
    v = np.asarray(evaluation_speed_m_s, dtype=float)
    u = np.asarray(source_speed_m_s, dtype=float)
    cos_gamma = gyrotropic_cosine_between_velocities(evaluation_pitch=evaluation_pitch, source_pitch=source_pitch, gyrophase_angle_rad=gyrophase_angle_rad)
    speed_squared = v**2 + u**2 - 2.0 * v * u * cos_gamma
    
    return np.sqrt(np.maximum(speed_squared, 0.0))