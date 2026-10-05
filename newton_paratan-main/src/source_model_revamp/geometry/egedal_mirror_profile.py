"""Analytic mirror field and axial gradient from Egedal 2022 Eq 56

Use the central mirror interval with abs(ζ) <= 1
The midplane is ζ = 0 and the throats are at abs(ζ) = 1
"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike

def normalized_zeta(z_phys_m: ArrayLike, half_length_m: float):
    """Return normalized axial coordinate ζ = z/L"""
    return np.asarray(z_phys_m, dtype=float) / half_length_m

def throat_field_T(B0_T: float, mirror_ratio: float) -> float:
    """Return mirror throat field B_M = R_M * B_0"""
    return float(mirror_ratio * B0_T)

def egedal_B_tilde(zeta: ArrayLike, mirror_ratio: float, transition_width: float, shape_exponent: float):
    """Return the normalized magnetic field B/B0 from Egedal Eq 56

    zeta is the axial coordinate divided by the mirror half length
    mirror_ratio sets the throat field relative to the midplane field
    transition_width and shape_exponent set the shape near each throat
    """
    zeta = np.asarray(zeta, dtype=float)
    k = float(shape_exponent)
    w = float(transition_width)
    w_k = w**k
    # The offset fixes the midplane value at one and the throat value at mirror_ratio
    exp_edge = np.exp(-1.0 / w_k)
    numerator = np.exp(-(np.abs(np.abs(zeta) - 1.0) ** k) / w_k) - exp_edge
    denominator = 1.0 - exp_edge

    return 1.0 + (float(mirror_ratio) - 1.0) * numerator / denominator

def egedal_dB_tilde_dzeta(zeta: ArrayLike, mirror_ratio: float, transition_width: float, shape_exponent: float):
    """Return the Eq 56 gradient on the central mirror interval where abs(ζ) <= 1"""
    zeta = np.asarray(zeta, dtype=float)
    k = float(shape_exponent)
    w = float(transition_width)
    w_k = w**k
    a = 1.0 - np.abs(zeta)
    exp_edge = np.exp(-1.0 / w_k)

    return ((float(mirror_ratio) - 1.0) * np.exp(-(a**k) / w_k) * k * a ** (k - 1.0) * np.sign(zeta) / (w_k * (1.0 - exp_edge)))

def B_phys_T(z_phys_m: ArrayLike, B0_T: float, half_length_m: float, mirror_ratio: float, transition_width: float, shape_exponent: float):
    """Return physical analytic field B(z) = B0 * B_tilde(z/L)"""
    zeta = normalized_zeta(z_phys_m=z_phys_m, half_length_m=half_length_m)

    return float(B0_T) * egedal_B_tilde(zeta=zeta, mirror_ratio=mirror_ratio, transition_width=transition_width, shape_exponent=shape_exponent)

def dB_phys_dz_T_per_m(z_phys_m: ArrayLike, B0_T: float, half_length_m: float, mirror_ratio: float, transition_width: float, shape_exponent: float):
    """Return dB/dz = (B0/L) * dB_tilde/dζ for the analytic profile"""
    zeta = normalized_zeta(z_phys_m=z_phys_m, half_length_m=half_length_m)

    return (float(B0_T) / float(half_length_m)) * egedal_dB_tilde_dzeta(zeta=zeta, mirror_ratio=mirror_ratio, transition_width=transition_width, shape_exponent=shape_exponent)

def throat_locations_m(half_length_m: float) -> tuple[float, float]:
    """Return the left and right throat positions relative to the midplane"""
    return (-float(half_length_m), float(half_length_m))
