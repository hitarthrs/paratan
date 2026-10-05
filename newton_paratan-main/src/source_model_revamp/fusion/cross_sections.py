"""Fusion kinematics helpers and Bosch Hale Table VII thermal reactivity fits"""
from __future__ import annotations
from dataclasses import dataclass
from collections.abc import Callable
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fusion.reactions import CM3_TO_M3, KEV_TO_J, J_TO_KEV
from source_model_revamp.fusion.bosch_hale import bosch_hale_cross_section_m2_from_keV

@dataclass(frozen=True)
class BoschHaleReactivityCoefficients:
    """Bosch Hale Table VII Maxwellian thermal reactivity coefficients"""
    c1: float
    c2: float
    c3: float
    c4: float
    c5: float
    c6: float
    c7: float
    b_g: float
    reduced_mass_energy_keV: float

DT_N_ALPHA_REACTIVITY_COEFFICIENTS = BoschHaleReactivityCoefficients(c1=1.17302e-9, c2=1.51361e-2, c3=7.51886e-2, c4=4.60643e-3, c5=1.35000e-2, c6=-1.06750e-4, c7=1.36600e-5, b_g=34.3827, reduced_mass_energy_keV=1_124_656.0)
DD_P_T_REACTIVITY_COEFFICIENTS = BoschHaleReactivityCoefficients(c1=5.65718e-12, c2=3.41267e-3, c3=1.99167e-3, c4=0.0, c5=1.05060e-5, c6=0.0, c7=0.0, b_g=31.3970, reduced_mass_energy_keV=937_814.0)
DD_N_HE3_REACTIVITY_COEFFICIENTS = BoschHaleReactivityCoefficients(c1=5.43360e-12, c2=5.85778e-3, c3=7.68222e-3, c4=0.0, c5=-2.96400e-6, c6=0.0, c7=0.0, b_g=31.3970, reduced_mass_energy_keV=937_814.0)

def _as_energy_array(values: ArrayLike, name: str) -> np.ndarray:
    """Return a finite scalar or array energy representation"""
    array = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    
    return array

def reduced_mass_kg(mass_a_kg: float, mass_b_kg: float) -> float:
    """Reduced mass μ = m_a * m_b / (m_a + m_b)"""
    ma = float(mass_a_kg)
    mb = float(mass_b_kg)
    if not np.isfinite(ma) or not np.isfinite(mb) or ma <= 0.0 or mb <= 0.0:
        raise ValueError("masses must be positive and finite")
    
    return ma * mb / (ma + mb)

def center_of_mass_energy_J(relative_speed_m_s: ArrayLike, reduced_mass_kg: float):
    """E_cm = 0.5 * μ * g^2"""
    g = np.asarray(relative_speed_m_s, dtype=float)
    mu = float(reduced_mass_kg)
    if mu <= 0.0 or not np.isfinite(mu):
        raise ValueError("reduced_mass_kg must be positive and finite")
    
    return 0.5 * mu * g**2

def relative_speed_from_center_of_mass_energy_m_s(energy_cm_J: ArrayLike, reduced_mass_kg: float):
    """g = sqrt(2 E_cm / μ)"""
    energy = _as_energy_array(energy_cm_J, "energy_cm_J")
    mu = float(reduced_mass_kg)
    if mu <= 0.0 or not np.isfinite(mu):
        raise ValueError("reduced_mass_kg must be positive and finite")
    
    return np.sqrt(np.maximum(2.0 * energy / mu, 0.0))

def bosch_hale_reactivity_coefficients_for_reaction(reaction: str) -> BoschHaleReactivityCoefficients:
    """Return Table VII thermal reactivity coefficients for a supported reaction label or alias"""
    label = str(reaction).strip().lower().replace(" ", "")
    if label in {"dt", "dt_n", "d-t", "td", "t-d", "t(d,n)4he", "d(t,n)4he", "d(t,n)alpha", "d(t,n)a"}:
        return DT_N_ALPHA_REACTIVITY_COEFFICIENTS
    if label in {"dd", "dd_n", "d-d", "dd-neutron", "dd_neutron", "d(d,n)3he", "d(d,n)he3"}:
        return DD_N_HE3_REACTIVITY_COEFFICIENTS
    if label in {"dd_p", "dd-proton", "dd_proton", "d(d,p)t"}:
        return DD_P_T_REACTIVITY_COEFFICIENTS
    raise ValueError(f"Unsupported active Bosch-Hale thermal-reactivity reaction label {reaction!r}")

def bosch_hale_theta_keV(temperature_keV: ArrayLike, coefficients: BoschHaleReactivityCoefficients):
    """Effective temperature θ used by the Bosch Hale Table VII reactivity fit"""
    T = _as_energy_array(temperature_keV, "temperature_keV")
    if np.any(T <= 0.0):
        raise ValueError("temperature_keV must be positive")
    c = coefficients
    numerator = T * (c.c2 + T * (c.c4 + T * c.c6))
    denominator = 1.0 + T * (c.c3 + T * (c.c5 + T * c.c7))

    return T / (1.0 - numerator / denominator)

def bosch_hale_xi(temperature_keV: ArrayLike, coefficients: BoschHaleReactivityCoefficients):
    """ξ = (B_G^2 / (4θ))^(1/3)"""
    theta = bosch_hale_theta_keV(temperature_keV, coefficients)

    return (coefficients.b_g**2 / (4.0 * theta)) ** (1.0 / 3.0)

def bosch_hale_thermal_reactivity_cm3_s(temperature_keV: ArrayLike, coefficients: BoschHaleReactivityCoefficients):
    """Evaluate the Bosch Hale Table VII Maxwellian thermal reactivity in cm^3/s"""
    T = _as_energy_array(temperature_keV, "temperature_keV")
    if np.any(T <= 0.0):
        raise ValueError("temperature_keV must be positive")
    theta = bosch_hale_theta_keV(T, coefficients)
    xi = bosch_hale_xi(T, coefficients)

    return coefficients.c1 * theta * np.sqrt(xi / (coefficients.reduced_mass_energy_keV * T**3)) * np.exp(-3.0 * xi)

def bosch_hale_thermal_reactivity_m3_s(temperature_energy_J: ArrayLike, coefficients: BoschHaleReactivityCoefficients):
    """Evaluate the Bosch Hale Table VII Maxwellian thermal reactivity in m^3/s with T in joules"""
    T_J = _as_energy_array(temperature_energy_J, "temperature_energy_J")

    return bosch_hale_thermal_reactivity_cm3_s(T_J / KEV_TO_J, coefficients) * CM3_TO_M3

def bosch_hale_thermal_reactivity_for_reaction_m3_s(reaction: str, temperature_energy_J: ArrayLike):
    """Evaluate the Table VII Maxwellian thermal reactivity for a supported reaction in m^3/s"""
    return bosch_hale_thermal_reactivity_m3_s(temperature_energy_J, bosch_hale_reactivity_coefficients_for_reaction(reaction))

def bosch_hale_cross_section_m2_from_J(reaction: str, energy_cm_J: ArrayLike):
    """Evaluate a Table IV cross section with E_cm supplied in joules"""
    E_J = _as_energy_array(energy_cm_J, "energy_cm_J")

    return bosch_hale_cross_section_m2_from_keV(reaction, E_J * J_TO_KEV)

def sigma_g_from_cross_section(relative_speed_m_s: ArrayLike, reduced_mass_kg: float, cross_section_function_m2: Callable[[ArrayLike], ArrayLike]):
    """Return σ(E_cm) * g for a relative speed grid"""
    g = np.asarray(relative_speed_m_s, dtype=float)
    if not np.all(np.isfinite(g)):
        raise ValueError("relative_speed_m_s must be finite")
    energy = center_of_mass_energy_J(g, reduced_mass_kg)

    return np.asarray(cross_section_function_m2(energy), dtype=float) * g