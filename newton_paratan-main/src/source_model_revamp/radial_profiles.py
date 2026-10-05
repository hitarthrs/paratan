"""Prescribed normalized magnetic flux radial profiles"""
from __future__ import annotations
import numpy as np

KOTELNIKOV_PARABOLIC_FLUX_K1 = "kotelnikov_parabolic_flux_k1"
SUPPORTED_RADIAL_PROFILE_MODELS = frozenset({KOTELNIKOV_PARABOLIC_FLUX_K1})

def canonical_radial_profile_model(value: str | None) -> str:
    model = KOTELNIKOV_PARABOLIC_FLUX_K1 if value is None else str(value).strip().lower()
    if model not in SUPPORTED_RADIAL_PROFILE_MODELS:
        allowed = ", ".join(sorted(SUPPORTED_RADIAL_PROFILE_MODELS))
        raise ValueError(f"radial profile model must be one of {allowed}")
 
    return model

def _rho_array(value: np.ndarray | float) -> np.ndarray:
    rho = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(rho)):
        raise ValueError("normalized radial coordinate must be finite")

    return np.clip(rho, 0.0, 1.0)

def radial_source_shape(value: np.ndarray | float, model: str = KOTELNIKOV_PARABOLIC_FLUX_K1) -> np.ndarray:
    """Return the area normalized Kotelnikov k equals one source shape"""
    canonical_radial_profile_model(model)
    rho = _rho_array(value)
   
    return 2.0 * (1.0 - rho**2)

def radial_probability_cdf(value: np.ndarray | float, model: str = KOTELNIKOV_PARABOLIC_FLUX_K1) -> np.ndarray:
    """Return the cumulative source probability inside normalized radius ρ"""
    canonical_radial_profile_model(model)
    rho = _rho_array(value)
  
    return 2.0 * rho**2 - rho**4

def radial_probability_inverse_cdf(value: np.ndarray | float, model: str = KOTELNIKOV_PARABOLIC_FLUX_K1) -> np.ndarray:
    """Invert the normalized radial source probability CDF"""
    canonical_radial_profile_model(model)
    probability = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(probability)) or np.any(probability < 0.0) or np.any(probability > 1.0):
        raise ValueError("radial profile probability must lie inside zero to one")
  
    return np.sqrt(np.maximum(0.0, 1.0 - np.sqrt(np.maximum(0.0, 1.0 - probability))))

def radial_conditional_inverse_cdf(value: np.ndarray | float, rho_max: np.ndarray | float, model: str = KOTELNIKOV_PARABOLIC_FLUX_K1) -> np.ndarray:
    """Sample the radial profile conditioned on ρ being below ρ_max"""
    canonical_radial_profile_model(model)
    probability = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(probability)) or np.any(probability < 0.0) or np.any(probability > 1.0):
        raise ValueError("conditional radial probability must lie inside zero to one")
    limit = _rho_array(rho_max)
    absolute_probability = probability * radial_probability_cdf(limit, model)
  
    return radial_probability_inverse_cdf(absolute_probability, model)

def radial_conditional_mean_rho_squared(value: np.ndarray | float, model: str = KOTELNIKOV_PARABOLIC_FLUX_K1) -> np.ndarray:
    """Return mean ρ^2 conditioned on ρ being below the supplied limit"""
    canonical_radial_profile_model(model)
    rho = _rho_array(value)
    denominator = radial_probability_cdf(rho, model)
    numerator = rho**4 - (2.0 / 3.0) * rho**6
  
    return np.divide(numerator, denominator, out=np.zeros_like(rho), where=denominator > 0.0)

def _cdf_antiderivative(rho: np.ndarray) -> np.ndarray:
    return (2.0 / 3.0) * rho**3 - (1.0 / 5.0) * rho**5

def radial_survival_average(rho_start: float, rho_end: float, segment_fraction: np.ndarray, model: str = KOTELNIKOV_PARABOLIC_FLUX_K1) -> np.ndarray:
    """Average surviving source fraction over the visited part of one path segment"""
    canonical_radial_profile_model(model)
    fraction = np.asarray(segment_fraction, dtype=float)
    if np.any(~np.isfinite(fraction)) or np.any(fraction < 0.0) or np.any(fraction > 1.0):
        raise ValueError("segment fraction must lie inside zero to one")
    start = float(np.clip(rho_start, 0.0, 1.0))
    end = float(np.clip(rho_end, 0.0, 1.0))
    stop = start + fraction * (end - start)
    delta = stop - start
    numerator = _cdf_antiderivative(stop) - _cdf_antiderivative(np.asarray(start))
    average = np.divide(numerator, delta, out=np.full_like(stop, float(radial_probability_cdf(start, model))), where=np.abs(delta) > 64.0 * np.finfo(float).eps)
  
    return np.where(fraction > 0.0, average, 0.0)  