"""
Production beam source definitions for the modal FBIS source model

File contains Lambda space beam source states used by attenuation and conservative source placement
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import numpy as np
from numpy.typing import ArrayLike
from scipy.constants import elementary_charge
from source_model_revamp.orbits.invariants import lambda_from_pitch_angle, speed_from_energy_m_s as invariant_speed_from_energy_m_s

def _coalesce(primary: Any, alias: Any, primary_name: str, alias_name: str):
    if primary is None:
        primary = alias
    elif alias is not None:
        raise TypeError(f"Pass only one of {primary_name!r} or {alias_name!r}, not both")
    if primary is None:
        raise TypeError(f"Missing required argument {primary_name!r}")
  
    return primary

def _require_positive(name: str, value: ArrayLike) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be positive")
    
    return arr

def _require_nonnegative(name: str, value: ArrayLike) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be nonnegative")
   
    return arr

def _require_finite(name: str, value: ArrayLike) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
   
    return arr

@dataclass(frozen=True)
class BeamSourceDefinition:
    """Representation independent monoenergetic, mono pitch Λ space source"""
    energy_J: float
    speed_m_s: float
    birth_rate_s: float
    source_rate_density_m3_s: float
    pitch_angle_rad: float
    B_tilde_birth: float
    Lambda_birth: float

@dataclass(frozen=True)
class MultiEnergyBeamSourceDefinition:
    """Single beam represented by multiple energy components"""
    component_definitions: tuple[BeamSourceDefinition, ...]
    total_birth_rate_s: float
    total_source_rate_density_m3_s: float

def beam_speed_from_energy_m_s( beam_energy_J: ArrayLike | None = None, fast_ion_mass_kg: float | None = None, *, energy_J: ArrayLike | None = None, particle_mass_kg: float | None = None):
    """Return beam particle speed: v_b = sqrt(2 * E_b / m_f)"""
    beam_energy_J = _coalesce(beam_energy_J, energy_J, "beam_energy_J", "energy_J")
    fast_ion_mass_kg = _coalesce(fast_ion_mass_kg, particle_mass_kg, "fast_ion_mass_kg", "particle_mass_kg")
    energy = _require_positive("beam_energy_J", beam_energy_J)
    mass = float(_require_positive("fast_ion_mass_kg", fast_ion_mass_kg))

    return invariant_speed_from_energy_m_s(energy_J=energy, mass_kg=mass)

def particle_birth_rate_from_power_s(beam_power_W: ArrayLike | None = None, beam_energy_J: ArrayLike | None = None, *, power_W: ArrayLike | None = None, energy_J: ArrayLike | None = None):
    """Return beam particle birth rate: Ndot_b = P_b / E_b"""
    beam_power_W = _coalesce(beam_power_W, power_W, "beam_power_W", "power_W")
    beam_energy_J = _coalesce(beam_energy_J, energy_J, "beam_energy_J", "energy_J")
    power = _require_nonnegative("beam_power_W", beam_power_W)
    energy = _require_positive("beam_energy_J", beam_energy_J)

    return power / energy

def particle_birth_rate_from_current_s(beam_current_A: ArrayLike | None = None, charge_state: float | None = None, *, current_A: ArrayLike | None = None, charge_number: float | None = None):
    """Return equivalent ion birth rate: Ndot_b = I_b / (Z_f * e)"""
    beam_current_A = _coalesce(beam_current_A, current_A, "beam_current_A", "current_A")
    if charge_state is None:
        charge_state = 1.0 if charge_number is None else charge_number
    current = _require_nonnegative("beam_current_A", beam_current_A)
    charge = float(_require_positive("charge_state", charge_state))

    return current / (charge * elementary_charge)

def component_powers_from_fractions_W(total_beam_power_W: float, power_fractions: ArrayLike) -> np.ndarray:
    """Return component powers P_i = f_i * P_total"""
    total_power = float(_require_nonnegative("total_beam_power_W", total_beam_power_W))
    fractions = _require_nonnegative("power_fractions", power_fractions)

    return total_power * fractions

def component_birth_rates_from_power_fractions_s(total_beam_power_W: float, component_energies_J: ArrayLike, power_fractions: ArrayLike) -> np.ndarray:
    """Return component particle birth rates for power fraction beam components"""
    energies = _require_positive("component_energies_J", component_energies_J)
    powers = component_powers_from_fractions_W(total_beam_power_W=total_beam_power_W, power_fractions=power_fractions)
    if powers.shape != energies.shape:
        raise ValueError("component_energies_J and power_fractions must have the same shape")
    return particle_birth_rate_from_power_s(beam_power_W=powers, beam_energy_J=energies)

def source_rate_density_from_birth_rate_m3_s(birth_rate_s: ArrayLike, volume_m3: float):
    """Return uniform volume source rate density Γ_b = Ndot_b / V"""
    birth_rate = _require_nonnegative("birth_rate_s", birth_rate_s)
    volume = float(_require_positive("volume_m3", volume_m3))
    return birth_rate / volume

def beam_lambda_birth(pitch_angle_rad: ArrayLike, B_tilde_birth: ArrayLike):
    """Return birth invariant Λ_b = sin^2(θ_b) / B_tilde_birth"""
    pitch = _require_finite("pitch_angle_rad", pitch_angle_rad)
    B_tilde = _require_positive("B_tilde_birth", B_tilde_birth)
    
    return lambda_from_pitch_angle(pitch_angle_rad=pitch, B_tilde=B_tilde)

def beam_source_definition_from_birth_rate(beam_energy_J: float, birth_rate_s: float, pitch_angle_rad: float, B_tilde_birth: float, source_volume_m3: float, fast_ion_mass_kg: float) -> BeamSourceDefinition:
    """Build the Λ space source definition from a particle birth rate"""
    speed = float(beam_speed_from_energy_m_s(beam_energy_J=beam_energy_J, fast_ion_mass_kg=fast_ion_mass_kg))
    source_rate_density = float(source_rate_density_from_birth_rate_m3_s(birth_rate_s=birth_rate_s, volume_m3=source_volume_m3))
    Lambda_b = float(beam_lambda_birth(pitch_angle_rad=pitch_angle_rad, B_tilde_birth=B_tilde_birth))
    
    return BeamSourceDefinition(
        energy_J=float(_require_positive("beam_energy_J", beam_energy_J)),
        speed_m_s=speed,
        birth_rate_s=float(_require_nonnegative("birth_rate_s", birth_rate_s)),
        source_rate_density_m3_s=source_rate_density,
        pitch_angle_rad=float(_require_finite("pitch_angle_rad", pitch_angle_rad)),
        B_tilde_birth=float(_require_positive("B_tilde_birth", B_tilde_birth)),
        Lambda_birth=Lambda_b,
    )

def beam_source_definition_from_power(beam_energy_J: float, beam_power_W: float, pitch_angle_rad: float, B_tilde_birth: float, source_volume_m3: float, fast_ion_mass_kg: float) -> BeamSourceDefinition:
    """Build the Λ space source definition from beam power"""
    birth_rate = float(particle_birth_rate_from_power_s(beam_power_W=beam_power_W, beam_energy_J=beam_energy_J))
    return beam_source_definition_from_birth_rate(beam_energy_J=beam_energy_J, birth_rate_s=birth_rate, pitch_angle_rad=pitch_angle_rad, B_tilde_birth=B_tilde_birth, source_volume_m3=source_volume_m3, fast_ion_mass_kg=fast_ion_mass_kg)

def beam_source_definition_from_current(beam_energy_J: float, beam_current_A: float, charge_state: float, pitch_angle_rad: float, B_tilde_birth: float, source_volume_m3: float, fast_ion_mass_kg: float) -> BeamSourceDefinition:
    """Build the Λ space source definition from equivalent ion current"""
    birth_rate = float(particle_birth_rate_from_current_s(beam_current_A=beam_current_A, charge_state=charge_state))
    return beam_source_definition_from_birth_rate(beam_energy_J=beam_energy_J, birth_rate_s=birth_rate, pitch_angle_rad=pitch_angle_rad, B_tilde_birth=B_tilde_birth, source_volume_m3=source_volume_m3, fast_ion_mass_kg=fast_ion_mass_kg)

def multi_energy_beam_source_definition_from_power_fractions(total_beam_power_W: float, component_energies_J: ArrayLike, power_fractions: ArrayLike, pitch_angles_rad: ArrayLike, B_tilde_birth_values: ArrayLike, source_volume_m3: float, fast_ion_mass_kg: float) -> MultiEnergyBeamSourceDefinition:
    """Build source definitions for energy components of one neutral beam"""
    energies = _require_positive("component_energies_J", component_energies_J)
    pitch_angles = _require_finite("pitch_angles_rad", pitch_angles_rad)
    B_values = _require_positive("B_tilde_birth_values", B_tilde_birth_values)
    birth_rates = component_birth_rates_from_power_fractions_s(total_beam_power_W=total_beam_power_W, component_energies_J=energies, power_fractions=power_fractions)
    if not (energies.shape == birth_rates.shape == pitch_angles.shape == B_values.shape):
        raise ValueError("component_energies_J, birth rates, pitch_angles_rad, and B_tilde_birth_values must have the same shape")
    components: list[BeamSourceDefinition] = []
    for E, rate, pitch, B_birth in zip(energies, birth_rates, pitch_angles, B_values):
        components.append(beam_source_definition_from_birth_rate(beam_energy_J=float(E), birth_rate_s=float(rate), pitch_angle_rad=float(pitch), B_tilde_birth=float(B_birth), source_volume_m3=source_volume_m3, fast_ion_mass_kg=fast_ion_mass_kg))
    total_birth = float(np.sum(birth_rates))
    total_density = float(source_rate_density_from_birth_rate_m3_s(birth_rate_s=total_birth, volume_m3=source_volume_m3))
    
    return MultiEnergyBeamSourceDefinition(component_definitions=tuple(components), total_birth_rate_s=total_birth, total_source_rate_density_m3_s=total_density)