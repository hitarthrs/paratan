"""
Isotropic Rosenbluth potentials for physical velocity distributions

Distribution convention
    n = integral f(v) d^3v = 4π integral_0^inf v^2 f(v) dv

For isotropic f, the angular integrals reduce to 1D speed integrals
    H(v) = 4π [ M2(v) / v + N1(v) ]
    G(v) = 4π [ v * M2(v) + M4(v) / (3v) + N3(v) + v^2 * N1(v) / 3 ]

with
    M2(v) = integral_0^v u^2 * f(u) du
    M4(v) = integral_0^v u^2 * f(u) du
    N1(v) = integral_v^vmax u * f(u) du
    N3(v) = integral_v^vmax u^3 * f(u) du
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import cumulative_trapezoid

@dataclass(frozen=True)
class IsotropicRosenbluthPotentials:
    """Rosenbluth potentials and radial derivatives on a speed grid"""
    speed_m_s: np.ndarray
    distribution: np.ndarray
    density_m3: float
    H: np.ndarray
    G: np.ndarray
    dH_dv: np.ndarray
    d2H_dv2: np.ndarray
    dG_dv: np.ndarray
    d2G_dv2: np.ndarray

def _validated_speed_distribution(speed_m_s: ArrayLike, distribution: ArrayLike, *, allow_negative_distribution: bool = False) -> tuple[np.ndarray, np.ndarray]:
    v = np.asarray(speed_m_s, dtype=float)
    f = np.asarray(distribution, dtype=float)
    if v.ndim != 1 or f.ndim != 1:
        raise ValueError("speed_m_s and distribution must be 1D")
    if v.shape != f.shape:
        raise ValueError("speed_m_s and distribution must have the same shape")
    if v.size < 2:
        raise ValueError("at least two speed points are required")
    if not np.all(np.isfinite(v)) or not np.all(np.isfinite(f)):
        raise ValueError("speed_m_s and distribution must contain only finite values")
    if np.any(np.diff(v) <= 0.0):
        raise ValueError("speed_m_s must be strictly increasing")
    if np.any(v < 0.0):
        raise ValueError("speed_m_s must be nonnegative")
    if not allow_negative_distribution and np.any(f < 0.0):
        raise ValueError("distribution must be nonnegative")
    
    return v, f

def cumulative_integral_from_zero(speed_m_s: ArrayLike, integrand: ArrayLike):
    """Cumulative trapezoid integral from the first speed point"""
    v, y = _validated_speed_distribution(speed_m_s, integrand, allow_negative_distribution=True)
    return cumulative_trapezoid(y, v, initial=0.0)

def cumulative_integral_to_infinity(speed_m_s: ArrayLike, integrand: ArrayLike):
    """Cumulative trapezoid integral from each grid point to the last speed point"""
    left = cumulative_integral_from_zero(speed_m_s, integrand)
    total = left[-1]
   
    return total - left

def isotropic_density_m3(speed_m_s: ArrayLike, distribution: ArrayLike):
    """n = 4π integral v^2 * f(v) dv"""
    v, f = _validated_speed_distribution(speed_m_s, distribution)
  
    return 4.0 * np.pi * np.trapezoid(v**2 * f, v)

def lower_speed_moment(speed_m_s: ArrayLike, distribution: ArrayLike, speed_power: float):
    """M_p(v) = integral_0^v u^p * f(u) du"""
    v, f = _validated_speed_distribution(speed_m_s, distribution)
   
    return cumulative_integral_from_zero(v, v**speed_power * f)

def upper_speed_moment(speed_m_s: ArrayLike, distribution: ArrayLike, speed_power: float):
    """N_p(v) = integral_v^vmax u^p * f(u) du on the finite input grid"""
    v, f = _validated_speed_distribution(speed_m_s, distribution)
  
    return cumulative_integral_to_infinity(v, v**speed_power * f)

def isotropic_rosenbluth_H(speed_m_s: ArrayLike, distribution: ArrayLike):
    """H(v) = integral f(u) / |v - u| d^3u for isotropic f"""
    v, f = _validated_speed_distribution(speed_m_s, distribution)
    M2 = lower_speed_moment(v, f, 2.0)
    N1 = upper_speed_moment(v, f, 1.0)
    H = np.empty_like(v)
    nonzero = v != 0.0
    H[nonzero] = 4.0 * np.pi * (M2[nonzero] / v[nonzero] + N1[nonzero])
    H[~nonzero] = 4.0 * np.pi * N1[~nonzero]
   
    return H

def isotropic_rosenbluth_G(speed_m_s: ArrayLike, distribution: ArrayLike):
    """G(v) = integral f(u) * |v - u| d^3u for isotropic f"""
    v, f = _validated_speed_distribution(speed_m_s, distribution)
    M2 = lower_speed_moment(v, f, 2.0)
    M4 = lower_speed_moment(v, f, 4.0)
    N1 = upper_speed_moment(v, f, 1.0)
    N3 = upper_speed_moment(v, f, 3.0)
    G = np.empty_like(v)
    nonzero = v != 0.0
    G[nonzero] = 4.0 * np.pi * (v[nonzero] * M2[nonzero] + M4[nonzero] / (3.0 * v[nonzero]) + N3[nonzero] + v[nonzero] ** 2 * N1[nonzero] / 3.0)
    G[~nonzero] = 4.0 * np.pi * N3[~nonzero]
  
    return G

def isotropic_dH_dv(speed_m_s: ArrayLike, distribution: ArrayLike):
    """dH/dv = -4π * M2(v) / v^2, with dH/dv = 0 at v = 0"""
    v, f = _validated_speed_distribution(speed_m_s, distribution)
    M2 = lower_speed_moment(v, f, 2.0)
    dH = np.zeros_like(v)
    nonzero = v != 0.0
    dH[nonzero] = -4.0 * np.pi * M2[nonzero] / v[nonzero] ** 2
  
    return dH

def isotropic_d2H_dv2(speed_m_s: ArrayLike, distribution: ArrayLike):
    """Second radial derivative of H for isotropic f"""
    v, f = _validated_speed_distribution(speed_m_s, distribution)
    M2 = lower_speed_moment(v, f, 2.0)
    d2H = np.empty_like(v)
    nonzero = v != 0.0
    d2H[nonzero] = 8.0 * np.pi * M2[nonzero] / v[nonzero] ** 3 - 4.0 * np.pi * f[nonzero]
    d2H[~nonzero] = -4.0 * np.pi * f[~nonzero] / 3.0
  
    return d2H

def isotropic_dG_dv(speed_m_s: ArrayLike, distribution: ArrayLike):
    """First radial derivative of G for isotropic f"""
    v, f = _validated_speed_distribution(speed_m_s, distribution)
    M2 = lower_speed_moment(v, f, 2.0)
    M4 = lower_speed_moment(v, f, 4.0)
    N1 = upper_speed_moment(v, f, 1.0)
    dG = np.zeros_like(v)
    nonzero = v != 0.0
    dG[nonzero] = 4.0 * np.pi * (M2[nonzero] - M4[nonzero] / (3.0 * v[nonzero] ** 2) + 2.0 * v[nonzero] * N1[nonzero] / 3.0)
 
    return dG

def isotropic_d2G_dv2(speed_m_s: ArrayLike, distribution: ArrayLike):
    """Second radial derivative of G, using ∇^2_v G = 2H"""
    v, f = _validated_speed_distribution(speed_m_s, distribution)
    H = isotropic_rosenbluth_H(v, f)
    dG = isotropic_dG_dv(v, f)
    d2G = np.empty_like(v)
    nonzero = v != 0.0
    d2G[nonzero] = 2.0 * H[nonzero] - 2.0 * dG[nonzero] / v[nonzero]
    d2G[~nonzero] = 2.0 * H[~nonzero] / 3.0
  
    return d2G

def isotropic_laplacian_from_radial_derivatives(speed_m_s: ArrayLike, first_derivative: ArrayLike, second_derivative: ArrayLike):
    """Spherical radial Laplacian F'' + 2F'/v, with the v=0 limit from arrays"""
    v = np.asarray(speed_m_s, dtype=float)
    d1 = np.asarray(first_derivative, dtype=float)
    d2 = np.asarray(second_derivative, dtype=float)
    if v.ndim != 1 or d1.shape != v.shape or d2.shape != v.shape:
        raise ValueError("derivative arrays must be one-dimensional and match speed_m_s")
    if not np.all(np.isfinite(v)) or not np.all(np.isfinite(d1)) or not np.all(np.isfinite(d2)):
        raise ValueError("speed and derivative arrays must contain only finite values")
    if v.size < 2 or np.any(np.diff(v) <= 0.0) or np.any(v < 0.0):
        raise ValueError("speed_m_s must be nonnegative and strictly increasing")
    laplacian = np.empty_like(v)
    nonzero = v != 0.0
    laplacian[nonzero] = d2[nonzero] + 2.0 * d1[nonzero] / v[nonzero]
    laplacian[~nonzero] = 3.0 * d2[~nonzero]
  
    return laplacian

def isotropic_rosenbluth_potentials(speed_m_s: ArrayLike, distribution: ArrayLike):
    """Compute isotropic H, G, radial derivatives, and density on a speed grid"""
    v, f = _validated_speed_distribution(speed_m_s, distribution)
   
    return IsotropicRosenbluthPotentials(speed_m_s=v, distribution=f, density_m3=isotropic_density_m3(v, f), H=isotropic_rosenbluth_H(v, f), G=isotropic_rosenbluth_G(v, f), dH_dv=isotropic_dH_dv(v, f), d2H_dv2=isotropic_d2H_dv2(v, f), dG_dv=isotropic_dG_dv(v, f), d2G_dv2=isotropic_d2G_dv2(v, f))