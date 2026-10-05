"""Bosch Hale Table IV total fusion cross section fits"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike

MILLIBARN_TO_M2 = 1.0e-31
BARN_TO_M2 = 1.0e-28

@dataclass(frozen=True)
class BoschHaleCrossSectionCoefficients:
    """Bosch Hale Table IV S factor coefficients and E_cm fit interval for one reaction branch"""
    reaction: str
    B_G_sqrt_keV: float
    fit_energy_min_keV: float
    fit_energy_max_keV: float
    A1: float
    A2: float
    A3: float
    A4: float
    A5: float
    B1: float
    B2: float
    B3: float
    B4: float

DT_N_ALPHA_COEFFICIENTS = BoschHaleCrossSectionCoefficients(reaction="T(d,n)4He", B_G_sqrt_keV=34.3827, fit_energy_min_keV=0.5, fit_energy_max_keV=550.0, A1=6.927e4, A2=7.454e8, A3=2.050e6, A4=5.2002e4, A5=0.0, B1=6.38e1, B2=-9.95e-1, B3=6.981e-5, B4=1.728e-4)
DD_N_HE3_COEFFICIENTS = BoschHaleCrossSectionCoefficients(reaction="D(d,n)3He", B_G_sqrt_keV=31.3970, fit_energy_min_keV=0.5, fit_energy_max_keV=4900.0, A1=5.3701e4, A2=3.3027e2, A3=-1.2706e-1, A4=2.9327e-5, A5=-2.5151e-9, B1=0.0, B2=0.0, B3=0.0, B4=0.0)
DD_P_T_COEFFICIENTS = BoschHaleCrossSectionCoefficients(reaction="D(d,p)T", B_G_sqrt_keV=31.3970, fit_energy_min_keV=0.5, fit_energy_max_keV=5000.0, A1=5.5576e4, A2=2.1054e2, A3=-3.2638e-2, A4=1.4987e-6, A5=1.8181e-10, B1=0.0, B2=0.0, B3=0.0, B4=0.0)

def bosch_hale_s_factor_pade_keV_millibarn(E_cm_keV: ArrayLike, coefficients: BoschHaleCrossSectionCoefficients):
    """Evaluate the Table IV Pade S(E_cm) factor in keV millibarn

    Input and output follow the shape of E_cm_keV
    """
    energy = np.asarray(E_cm_keV, dtype=float)
    numerator = coefficients.A1 + energy * (coefficients.A2 + energy * (coefficients.A3 + energy * (coefficients.A4 + energy * coefficients.A5)))
    denominator = 1.0 + energy * (coefficients.B1 + energy * (coefficients.B2 + energy * (coefficients.B3 + energy * coefficients.B4)))

    return numerator / denominator

def bosch_hale_sigma_millibarn(E_cm_keV: ArrayLike, coefficients: BoschHaleCrossSectionCoefficients):
    """Evaluate σ(E_cm) = S(E_cm) / E_cm * exp(−B_G / sqrt(E_cm)) in millibarn

    Nonpositive E_cm returns zero and scalar input returns float
    """
    energy = np.asarray(E_cm_keV, dtype=float)
    scalar_input = energy.ndim == 0
    energy_work = np.atleast_1d(energy)
    sigma = np.zeros_like(energy_work, dtype=float)
    positive = energy_work > 0.0
    energy_positive = energy_work[positive]
    if energy_positive.size:
        s_factor = bosch_hale_s_factor_pade_keV_millibarn(energy_positive, coefficients)
        sigma[positive] = (s_factor / energy_positive * np.exp(-coefficients.B_G_sqrt_keV / np.sqrt(energy_positive)))
    if scalar_input:
        return float(sigma[0])
    
    return sigma

def bosch_hale_sigma_m2(E_cm_keV: ArrayLike, coefficients: BoschHaleCrossSectionCoefficients):
    """Evaluate the Table IV total fusion cross section σ(E_cm) in m^2"""
    return MILLIBARN_TO_M2 * bosch_hale_sigma_millibarn(E_cm_keV, coefficients)

def bosch_hale_sigma_DT_m2(E_cm_keV: ArrayLike):
    """Evaluate the T(d,n)4He Table IV cross section in m^2"""
    return bosch_hale_sigma_m2(E_cm_keV, DT_N_ALPHA_COEFFICIENTS)

def bosch_hale_sigma_DD_neutron_m2(E_cm_keV: ArrayLike):
    """Evaluate the D(d,n)3He Table IV cross section in m^2"""
    return bosch_hale_sigma_m2(E_cm_keV, DD_N_HE3_COEFFICIENTS)

def bosch_hale_sigma_DD_proton_m2(E_cm_keV: ArrayLike):
    """Evaluate the D(d,p)T Table IV cross section in m^2"""
    return bosch_hale_sigma_m2(E_cm_keV, DD_P_T_COEFFICIENTS)

def bosch_hale_cross_section_coefficients_for_reaction(reaction: str) -> BoschHaleCrossSectionCoefficients:
    """Return Table IV coefficients for a supported reaction label or alias"""
    label = str(reaction).strip().lower().replace(" ", "")
    if label in {"dt", "dt_n", "d-t", "td", "t-d", "t(d,n)4he", "d(t,n)4he", "d(t,n)alpha", "d(t,n)a"}:
        return DT_N_ALPHA_COEFFICIENTS
    if label in {"dd", "dd_n", "d-d", "dd-neutron", "dd_neutron", "d(d,n)3he", "d(d,n)he3"}:
        return DD_N_HE3_COEFFICIENTS
    if label in {"dd_p", "dd-proton", "dd_proton", "d(d,p)t"}:
        return DD_P_T_COEFFICIENTS
    raise ValueError("Unsupported active Bosch Hale cross section reaction label " f"{reaction!r}")

def bosch_hale_table_iv_fit_domain_keV(reaction: str,) -> tuple[float, float]:
    """Return the Table IV E_cm fit interval in keV for one reaction"""
    coefficients = bosch_hale_cross_section_coefficients_for_reaction(reaction)
    return (float(coefficients.fit_energy_min_keV), float(coefficients.fit_energy_max_keV))

def bosch_hale_cross_section_m2_from_keV(reaction: str, E_cm_keV: ArrayLike):
    """Evaluate a supported Table IV cross section from E_cm in keV"""
    coefficients = bosch_hale_cross_section_coefficients_for_reaction(reaction)
    return bosch_hale_sigma_m2(E_cm_keV, coefficients)