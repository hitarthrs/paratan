"""Analytic density normalized Maxwellian Rosenbluth functions used by Eq 59 ion pairs"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike
from scipy.special import erf
from source_model_revamp.fbis.modal.types import ModalRosenbluthCoefficients

_SMALL_X_THRESHOLD = 1.0e-4

def _positive_finite_scalar(value: float, name: str) -> float:
    """Validate one strictly positive finite scalar"""
    scalar = float(value)
    if not np.isfinite(scalar) or scalar <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
  
    return scalar

def _maxwellian_psi_and_functions(x: np.ndarray, thermal_speed_m_s: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return stable Ψ g_tilde_1 and g_tilde_2 for x = (v / v_t)²
    
    Small x series are used to avoid cancellation near v = 0
    """
    values = np.asarray(x, dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError("Maxwellian speed ratio must be finite and nonnegative")
    vt = _positive_finite_scalar(thermal_speed_m_s, "thermal_speed_m_s")
    psi = np.empty_like(values)
    g1 = np.empty_like(values)
    g2 = np.empty_like(values)
    small = values < _SMALL_X_THRESHOLD
    ordinary = ~small
    xs = values[small]
    sqrt_pi = np.sqrt(np.pi)
    sqrt_xs = np.sqrt(xs)
    psi[small] = (4.0 / (3.0 * sqrt_pi)) * xs * sqrt_xs - (4.0 / (5.0 * sqrt_pi)) * xs**2 * sqrt_xs + (2.0 / (7.0 * sqrt_pi)) * xs**3 * sqrt_xs - (2.0 / (27.0 * sqrt_pi)) * xs**4 * sqrt_xs
    g1[small] = (4.0 / (3.0 * sqrt_pi)) * sqrt_xs - (4.0 / (15.0 * sqrt_pi)) * xs * sqrt_xs + (2.0 / (35.0 * sqrt_pi)) * xs**2 * sqrt_xs - (2.0 / (189.0 * sqrt_pi)) * xs**3 * sqrt_xs
    g2[small] = vt * ((4.0 / (3.0 * sqrt_pi)) * xs - (4.0 / (5.0 * sqrt_pi)) * xs**2 + (2.0 / (7.0 * sqrt_pi)) * xs**3 - (2.0 / (27.0 * sqrt_pi)) * xs**4)
    xo = values[ordinary]
    root = np.sqrt(xo)
    erf_root = erf(root)
    psi_o = erf_root - 2.0 * root * np.exp(-xo) / sqrt_pi
    psi[ordinary] = psi_o
    g1[ordinary] = erf_root - psi_o / (2.0 * xo)
    g2[ordinary] = vt * psi_o / np.sqrt(xo)
    roundoff_scale = max(float(np.max(np.abs(np.concatenate((psi, g1, g2 / vt))))) if values.size else 0.0, 1.0)
    tolerance = 256.0 * np.finfo(float).eps * roundoff_scale
    for name, array in (("h_tilde", psi), ("g_tilde_1", g1), ("g_tilde_2", g2)):
        minimum = float(np.min(array)) if array.size else 0.0
        if minimum < -tolerance:
            raise ValueError(f"analytic Maxwellian {name} contains material negative values")
        array[array < 0.0] = 0.0
        if np.any(~np.isfinite(array)):
            raise ValueError(f"analytic Maxwellian {name} contains nonfinite values")
 
    return psi, g1, g2

def maxwellian_rosenbluth_coefficients(*, speed_m_s: ArrayLike, field_temperature_J: float, field_mass_kg: float) -> ModalRosenbluthCoefficients:
    """Return density normalized isotropic Maxwellian Eq 59 functions
    
    The field thermal speed is v_t = √(2 T / m) and the returned coefficient state uses unit density normalization
    """
    speed = np.asarray(speed_m_s, dtype=float)
    if speed.ndim != 1 or speed.size < 1 or np.any(~np.isfinite(speed)) or np.any(speed < 0.0):
        raise ValueError("speed_m_s must be a finite nonnegative one dimensional array")
    temperature = _positive_finite_scalar(field_temperature_J, "field_temperature_J")
    mass = _positive_finite_scalar(field_mass_kg, "field_mass_kg")
    thermal_speed = np.sqrt(2.0 * temperature / mass)
    x = (speed / thermal_speed) ** 2
    h_tilde, g1, g2 = _maxwellian_psi_and_functions(x, thermal_speed)
  
    return ModalRosenbluthCoefficients(h_tilde=h_tilde, g_tilde_1=g1, g_tilde_2=g2, density_normalization_m3=1.0)

__all__ = [
    "maxwellian_rosenbluth_coefficients",
]
