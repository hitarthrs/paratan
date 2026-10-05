"""Electron electrostatic loss boundary and Maxwellian tail utilities

Sign convention
    Potential differences are voltages relative to the midplane
        Δϕ = ϕ_wall - ϕ_midplane

    For electrons q_e = -e, the electrostatic energy change at the wall is
        q_e * Δϕ = -e * Δϕ

    The positive energy barrier which a midplane electron must overcome is
        E_barrier = e * (ϕ_midplane - ϕ_wall) = -e * Δϕ

Thermal convention
    Temperatures are represented as energies in joules
        y = E_barrier / T_e

    The escaping parallel speed threshold is
        v_w = sqrt(2 * E_barrier / m_e)

Scalar and array inputs follow NumPy broadcasting where more than one input is combined
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from scipy.constants import elementary_charge, electron_mass
ELECTRON_CHARGE_MAGNITUDE_C = elementary_charge
ELECTRON_MASS_KG = electron_mass

@dataclass(frozen=True)
class ElectronLossBoundary:
    """Electron wall loss boundary specified by a positive barrier energy
    
    barrier_energy_J is E_barrier in joules
    normalized_barrier is y = E_barrier / T_e
    wall_potential_relative_to_midplane_V is Δϕ = ϕ_wall - ϕ_midplane
    escape_speed_m_s is v_w = sqrt(2 * E_barrier / m_e)
    """
    barrier_energy_J: float | np.ndarray
    normalized_barrier: float | np.ndarray
    wall_potential_relative_to_midplane_V: float | np.ndarray
    escape_speed_m_s: float | np.ndarray

def electron_barrier_energy_from_wall_potential_J(wall_potential_relative_to_midplane_V: ArrayLike):
    """Convert wall potential relative to the midplane into barrier energy in joules
    
    E_barrier = e * (ϕ_midplane - ϕ_wall) = -e * Δϕ
    """
    delta_phi = np.asarray(wall_potential_relative_to_midplane_V, dtype=float)

    return -ELECTRON_CHARGE_MAGNITUDE_C * delta_phi

def electron_wall_potential_from_barrier_V(barrier_energy_J: ArrayLike):
    """Convert positive electron barrier energy into wall potential relative to the midplane
    
    Δϕ_wall = -E_barrier/e for an electron confining wall
    """
    barrier = np.asarray(barrier_energy_J, dtype=float)

    return -barrier / ELECTRON_CHARGE_MAGNITUDE_C

def electron_barrier_energy_from_normalized(normalized_barrier: ArrayLike, electron_temperature_J: ArrayLike):
    """Convert normalized barrier and electron temperature into joules
    
    E_barrier = y * T_e
    """
    y = np.asarray(normalized_barrier, dtype=float)
    T_e = np.asarray(electron_temperature_J, dtype=float)

    return y * T_e

def normalized_electron_barrier(barrier_energy_J: ArrayLike, electron_temperature_J: ArrayLike):
    """Return the dimensionless electron barrier
    
    y = E_barrier / T_e
    """
    barrier = np.asarray(barrier_energy_J, dtype=float)
    T_e = np.asarray(electron_temperature_J, dtype=float)

    return barrier / T_e

def electron_escape_speed_m_s(barrier_energy_J: ArrayLike, electron_mass_kg: float = ELECTRON_MASS_KG):
    """Return the parallel speed threshold in m/s
    
    v_w = sqrt(2 * E_barrier / m_e)
    """
    barrier = np.asarray(barrier_energy_J, dtype=float)

    return np.sqrt(2.0 * barrier / electron_mass_kg)

def electron_thermal_speed_m_s(electron_temperature_J: ArrayLike, electron_mass_kg: float = ELECTRON_MASS_KG):
    """Return the electron thermal speed in m/s for temperature in joules
    
    v_te = sqrt(2 * T_e / m_e)
    """
    T_e = np.asarray(electron_temperature_J, dtype=float)

    return np.sqrt(2.0 * T_e / electron_mass_kg)

def electron_escape_speed_ratio(barrier_energy_J: ArrayLike, electron_temperature_J: ArrayLike):
    """Return the dimensionless escape to thermal speed ratio
    
    v_w / v_te = sqrt(E_barrier / T_e)
    """
    return np.sqrt(normalized_electron_barrier(barrier_energy_J=barrier_energy_J, electron_temperature_J=electron_temperature_J,))

def maxwellian_parallel_tail_probability(barrier_energy_J: ArrayLike, electron_temperature_J: ArrayLike):
    """Return the one sided 1D Maxwellian probability above the parallel threshold
    
    P(v_parallel > v_w) = 0.5 erfc(sqrt(y))
    """
    from scipy.special import erfc
    y = normalized_electron_barrier(barrier_energy_J=barrier_energy_J, electron_temperature_J=electron_temperature_J)

    return 0.5 * erfc(np.sqrt(y))

def maxwellian_flux_tail_factor(barrier_energy_J: ArrayLike, electron_temperature_J: ArrayLike):
    """Return the Maxwellian flux suppression from the barrier
    
    Flux tail factor exp(-E_barrier / T_e)
    """
    y = normalized_electron_barrier(barrier_energy_J=barrier_energy_J, electron_temperature_J=electron_temperature_J)

    return np.exp(-y)

def mean_parallel_energy_midplane_for_escaping_electrons_J(barrier_energy_J: ArrayLike, electron_temperature_J: ArrayLike):
    """Flux weighted midplane parallel kinetic energy of electrons above barrier
    
    For a 1D Maxwellian half flux with x = K_parallel/T_e and x >= y,
        the flux weighting is proportional to exp(-x)
            so <K_parallel> = E_b + T_e
    """
    barrier = np.asarray(barrier_energy_J, dtype=float)
    T_e = np.asarray(electron_temperature_J, dtype=float)

    return barrier + T_e

def mean_total_energy_midplane_for_escaping_electrons_J(barrier_energy_J: ArrayLike, electron_temperature_J: ArrayLike):
    """Flux weighted total midplane kinetic energy of escaping Maxwellian electrons
    
    The perpendicular Maxwellian energy contributes T_e,
        independent of the parallel loss threshold, so
            <K_total,mid> = E_b + 2 T_e
    """
    barrier = np.asarray(barrier_energy_J, dtype=float)
    T_e = np.asarray(electron_temperature_J, dtype=float)

    return barrier + 2.0 * T_e

def mean_parallel_energy_at_wall_for_escaping_electrons_J(electron_temperature_J: ArrayLike):
    """Return the residual flux weighted parallel energy after the barrier
    
    <K_parallel,wall> = T_e
    """
    return np.asarray(electron_temperature_J, dtype=float)

def mean_total_energy_at_wall_for_escaping_electrons_J(electron_temperature_J: ArrayLike, include_perpendicular_energy: bool = False):
    """Return the mean kinetic energy delivered to the wall after electrostatic deceleration
    
    The default retains only the residual parallel contribution T_e
    When include_perpendicular_energy is true the returned total is 2 T_e
    """
    T_e = np.asarray(electron_temperature_J, dtype=float)
    if include_perpendicular_energy:
        return 2.0 * T_e

    return T_e

def build_electron_loss_boundary(barrier_energy_J: ArrayLike, electron_temperature_J: ArrayLike):
    """Build the electron boundary state from positive barrier energy and electron temperature
    
    Inputs may be scalars or broadcast compatible arrays
    """
    barrier = np.asarray(barrier_energy_J, dtype=float)
    T_e = np.asarray(electron_temperature_J, dtype=float)

    return ElectronLossBoundary(barrier_energy_J=barrier, normalized_barrier=normalized_electron_barrier(barrier, T_e), wall_potential_relative_to_midplane_V=electron_wall_potential_from_barrier_V(barrier), escape_speed_m_s=electron_escape_speed_m_s(barrier))

def build_electron_loss_boundary_from_wall_potential(wall_potential_relative_to_midplane_V: ArrayLike, electron_temperature_J: ArrayLike):
    """Build the electron boundary state from wall potential relative to the midplane
    
    The wall potential is first converted with E_barrier = -e * Δϕ
    """
    barrier = electron_barrier_energy_from_wall_potential_J(wall_potential_relative_to_midplane_V=wall_potential_relative_to_midplane_V)

    return build_electron_loss_boundary(barrier_energy_J=barrier, electron_temperature_J=electron_temperature_J)