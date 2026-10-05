"""Map arbitrary reactant pair invariant mass to the incident energy used by deuteron angular laws"""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.constants import SPEED_OF_LIGHT_M_S

ELECTRON_VOLT_J = 1.602176634e-19

def equivalent_projectile_lab_kinetic_energy_J_from_invariant_s(invariant_s_J2: float, projectile_mass_kg: float, target_mass_kg: float) -> float:
    """Return equivalent projectile laboratory kinetic energy in J from invariant `s`

    The target is stationary and `s` is supplied in `J^2` so that
    `K_lab = [s − (m_p c^2 + m_t c^2)^2] / (2 m_t c^2)`
    """
    s_value = float(invariant_s_J2)
    projectile_mass = float(projectile_mass_kg)
    target_mass = float(target_mass_kg)
    if not np.isfinite(s_value) or s_value <= 0.0:
        raise ValueError("invariant_s_J2 must be positive and finite")
    if not np.isfinite(projectile_mass) or projectile_mass <= 0.0:
        raise ValueError("projectile_mass_kg must be positive and finite")
    if not np.isfinite(target_mass) or target_mass <= 0.0:
        raise ValueError("target_mass_kg must be positive and finite")
    c2 = SPEED_OF_LIGHT_M_S**2
    projectile_rest_energy = projectile_mass * c2
    target_rest_energy = target_mass * c2
    threshold_s = (projectile_rest_energy + target_rest_energy) ** 2
    kinetic_energy = (s_value - threshold_s) / (2.0 * target_rest_energy)
    if kinetic_energy < -1.0e-12 * projectile_rest_energy:
        raise ValueError("invariant mass lies below the projectile target rest threshold")
    
    return float(max(kinetic_energy, 0.0))

def equivalent_projectile_lab_kinetic_energy_eV_from_invariant_s(invariant_s_J2: float, projectile_mass_kg: float, target_mass_kg: float) -> float:
    """Return the scalar equivalent projectile laboratory kinetic energy in eV"""
    return equivalent_projectile_lab_kinetic_energy_J_from_invariant_s(invariant_s_J2, projectile_mass_kg, target_mass_kg) / ELECTRON_VOLT_J

def equivalent_projectile_lab_kinetic_energy_J_from_invariant_s_array(invariant_s_J2: ArrayLike, projectile_mass_kg: float, target_mass_kg: float) -> np.ndarray:
    """Vectorize the exact invariant `s` to stationary target laboratory energy mapping

    The output preserves the input array shape and has units J
    """
    s_value = np.asarray(invariant_s_J2, dtype=float)
    projectile_mass = float(projectile_mass_kg)
    target_mass = float(target_mass_kg)
    if np.any(~np.isfinite(s_value)) or np.any(s_value <= 0.0):
        raise ValueError("invariant_s_J2 must be positive and finite")
    if not np.isfinite(projectile_mass) or projectile_mass <= 0.0:
        raise ValueError("projectile_mass_kg must be positive and finite")
    if not np.isfinite(target_mass) or target_mass <= 0.0:
        raise ValueError("target_mass_kg must be positive and finite")
    c2 = SPEED_OF_LIGHT_M_S**2
    projectile_rest_energy = projectile_mass * c2
    target_rest_energy = target_mass * c2
    threshold_s = (projectile_rest_energy + target_rest_energy) ** 2
    kinetic_energy = (s_value - threshold_s) / (2.0 * target_rest_energy)
    if np.any(kinetic_energy < -1.0e-12 * projectile_rest_energy):
        raise ValueError("invariant mass lies below the projectile target rest threshold")
    
    return np.maximum(kinetic_energy, 0.0)

def equivalent_projectile_lab_kinetic_energy_eV_from_invariant_s_array(invariant_s_J2: ArrayLike, projectile_mass_kg: float, target_mass_kg: float) -> np.ndarray:
    """Return the array equivalent projectile laboratory kinetic energies in eV"""
    return equivalent_projectile_lab_kinetic_energy_J_from_invariant_s_array(invariant_s_J2, projectile_mass_kg, target_mass_kg ) / ELECTRON_VOLT_J