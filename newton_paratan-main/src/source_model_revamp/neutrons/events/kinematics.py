"""Center of mass geometry helpers for correlated neutron event sampling"""
from __future__ import annotations
import numpy as np
from source_model_revamp.constants import SPEED_OF_LIGHT_M_S
from source_model_revamp.neutrons.kinematics import TwoBodyReactionKinematics

def pair_invariant_and_projectile_cm_direction( reaction_kinematics: TwoBodyReactionKinematics, velocity_a_m_s: np.ndarray, velocity_b_m_s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Return pair invariant `s` and the reactant a momentum direction in the pair center of mass frame
    
    Reactant velocity arrays have shape `(n, 3)` and are interpreted in the lab frame
    """
    kin = reaction_kinematics
    v_a = np.asarray(velocity_a_m_s, dtype=float)
    v_b = np.asarray(velocity_b_m_s, dtype=float)
    if v_a.ndim != 2 or v_a.shape[1] != 3 or v_b.shape != v_a.shape:
        raise ValueError("reactant velocity arrays must have shape n by 3")
    if np.any(~np.isfinite(v_a)) or np.any(~np.isfinite(v_b)):
        raise ValueError("reactant velocities must be finite")
    c = SPEED_OF_LIGHT_M_S
    speed2_a = np.einsum("ij,ij->i", v_a, v_a)
    speed2_b = np.einsum("ij,ij->i", v_b, v_b)
    beta2_a = speed2_a / c**2
    beta2_b = speed2_b / c**2
    if np.any(beta2_a >= 1.0) or np.any(beta2_b >= 1.0):
        raise ValueError("reactant speeds must be smaller than light speed")
    gamma_a = 1.0 / np.sqrt(1.0 - beta2_a)
    gamma_b = 1.0 / np.sqrt(1.0 - beta2_b)
    energy_a = gamma_a * kin.reactant_a_mass_kg * c**2
    energy_b = gamma_b * kin.reactant_b_mass_kg * c**2
    momentum_a = gamma_a[:, None] * kin.reactant_a_mass_kg * v_a
    momentum_b = gamma_b[:, None] * kin.reactant_b_mass_kg * v_b
    total_energy = energy_a + energy_b
    total_momentum = momentum_a + momentum_b
    invariant_s = total_energy**2 - np.einsum( "ij,ij->i", total_momentum, total_momentum) * c**2
    if np.any(~np.isfinite(invariant_s)) or np.any(invariant_s <= 0.0):
        raise ValueError("pair invariant mass must be positive and finite")
    beta_cm = total_momentum * c / total_energy[:, None]
    beta2_cm = np.einsum("ij,ij->i", beta_cm, beta_cm)
    if np.any(beta2_cm >= 1.0):
        raise ValueError("pair CM speed must be smaller than light speed")
    gamma_cm = 1.0 / np.sqrt(1.0 - beta2_cm)
    beta_dot_p = np.einsum("ij,ij->i", beta_cm, momentum_a)
    beta2_safe = np.where(beta2_cm > 0.0, beta2_cm, 1.0)
    coefficient = ( (gamma_cm - 1.0) * beta_dot_p / beta2_safe - gamma_cm * energy_a / c)
    coefficient = np.where(beta2_cm > 0.0, coefficient, 0.0)
    momentum_a_cm = momentum_a + coefficient[:, None] * beta_cm
    magnitude = np.linalg.norm(momentum_a_cm, axis=1)
    if np.any(magnitude <= 0.0):
        raise ValueError("projectile direction is undefined in the pair CM frame")
    
    return invariant_s, momentum_a_cm / magnitude[:, None]

def directions_about_axes( axes: np.ndarray, cosine: np.ndarray, azimuth_rad: np.ndarray) -> np.ndarray:
    """
    Build unit directions at supplied cosine and azimuth around each unit axis
    
    The returned array has shape `(n, 3)` and preserves the requested axis cosine to roundoff
    """
    axis = np.asarray(axes, dtype=float)
    mu = np.asarray(cosine, dtype=float)
    azimuth = np.asarray(azimuth_rad, dtype=float)
    if axis.ndim != 2 or axis.shape[1] != 3:
        raise ValueError("axes must have shape n by 3")
    if mu.shape != (axis.shape[0],) or azimuth.shape != mu.shape:
        raise ValueError("angular arrays must have one value per axis")
    if np.any(mu < -1.0) or np.any(mu > 1.0):
        raise ValueError("cosine values must lie inside minus one to one")
    reference = np.zeros_like(axis)
    use_z = np.abs(axis[:, 2]) < 0.9
    reference[use_z, 2] = 1.0
    reference[~use_z, 0] = 1.0
    first = np.cross(reference, axis)
    first /= np.linalg.norm(first, axis=1)[:, None]
    second = np.cross(axis, first)
    sine = np.sqrt(np.maximum(1.0 - mu**2, 0.0))
    direction = ( mu[:, None] * axis + sine[:, None] * ( np.cos(azimuth)[:, None] * first + np.sin(azimuth)[:, None] * second))
    direction /= np.linalg.norm(direction, axis=1)[:, None]

    return direction