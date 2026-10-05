"""
Atomic rate interfaces for neutral beam attenuation
    File converts supplied atomic data into the effective attenuation coefficients

Cross section gives
    κ_c(s, E_b) = sum_t n_t(s) * σ_{c, t}(E_b)

If using rate coefficients K_c(E_b) then attenuation rate is
    κ_c = sm_t n_t * K_c / v_b

Channels are kept separate
    ionization -> net plasma fueling and fast ion birth
    charge exchange -> fast ion birth source sink bookkeeping
    other loss -> primary neutral removed
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.beam_source_definition import beam_speed_from_energy_m_s

@dataclass(frozen=True)
class TabulatedCrossSection:
    """
    Immutable energy dependent cross section table
    
    energies_J and cross_sections_m2 follow the same table axis
    use_log_interpolation selects linear interpolation in logarithmic energy and value coordinates
    """
    energies_J: np.ndarray
    cross_sections_m2: np.ndarray
    use_log_interpolation: bool = False

@dataclass(frozen=True)
class TabulatedRateCoefficient:
    """
    Immutable energy dependent rate coefficient table
    
    energies_J and rate_coefficients_m3_s follow the same table axis
    use_log_interpolation selects linear interpolation in logarithmic energy and value coordinates
    """
    energies_J: np.ndarray
    rate_coefficients_m3_s: np.ndarray
    use_log_interpolation: bool = False

@dataclass(frozen=True)
class AtomicChannelData:
    """
    Atomic data for one target species and one attenuation channel
    
    A channel may contain a cross section a rate coefficient or both
    When both are present their equivalent cross section contributions are added
    """
    cross_section: TabulatedCrossSection | None = None
    rate_coefficient: TabulatedRateCoefficient | None = None

@dataclass(frozen=True)
class AtomicTargetProfile:
    """
    Target density profile and channel data for one background species
    
    The density array follows the attenuation grid supplied by the caller
    Each optional channel is evaluated against that same target density
    """
    density_profile_m3: np.ndarray
    ionization: AtomicChannelData | None = None
    charge_exchange: AtomicChannelData | None = None
    other_loss: AtomicChannelData | None = None

@dataclass(frozen=True)
class BeamAtomicAttenuationRates:
    """
    Line attenuation coefficients for aggregate and resolved atomic channels
    
    All stored rates have units m⁻¹
    Single energy calculations normally use axial arrays
    Multi energy calculations use component by axial cell arrays with component index first
    """
    ionization_rate_per_m: np.ndarray
    charge_exchange_rate_per_m: np.ndarray
    other_loss_rate_per_m: np.ndarray
    electron_ionization_rate_per_m: np.ndarray | None = None
    deuterium_ionization_rate_per_m: np.ndarray | None = None
    tritium_ionization_rate_per_m: np.ndarray | None = None
    deuterium_charge_exchange_rate_per_m: np.ndarray | None = None
    tritium_charge_exchange_rate_per_m: np.ndarray | None = None

def tabulated_cross_section(energies_J: ArrayLike, cross_sections_m2: ArrayLike, use_log_interpolation: bool = False):
    """Build a cross section table"""
    return TabulatedCrossSection(energies_J=np.asarray(energies_J, dtype=float), cross_sections_m2=np.asarray(cross_sections_m2, dtype=float), use_log_interpolation=bool(use_log_interpolation))

def tabulated_rate_coefficient(energies_J: ArrayLike, rate_coefficients_m3_s: ArrayLike, use_log_interpolation: bool = False):
    """Build a rate coefficient table"""
    return TabulatedRateCoefficient(energies_J=np.asarray(energies_J, dtype=float), rate_coefficients_m3_s=np.asarray(rate_coefficients_m3_s, dtype=float), use_log_interpolation=bool(use_log_interpolation))

def atomic_channel_from_cross_section(energies_J: ArrayLike, cross_sections_m2: ArrayLike, use_log_interpolation: bool = False):
    """Atomic channel data from a tabulated cross section"""
    return AtomicChannelData(cross_section=tabulated_cross_section(energies_J=energies_J, cross_sections_m2=cross_sections_m2, use_log_interpolation=use_log_interpolation))

def atomic_channel_from_rate_coefficient(energies_J: ArrayLike, rate_coefficients_m3_s: ArrayLike, use_log_interpolation: bool = False):
    """Atomic channel data from a tabulated rate coefficient"""
    return AtomicChannelData(rate_coefficient=tabulated_rate_coefficient(energies_J=energies_J, rate_coefficients_m3_s=rate_coefficients_m3_s, use_log_interpolation=use_log_interpolation))

def interpolate_table_values(energies_J: ArrayLike, values: ArrayLike, query_energy_J: ArrayLike, use_log_interpolation: bool = False):
    """Interpolate tabulated energy data"""
    energies = np.asarray(energies_J, dtype=float)
    table_values = np.asarray(values, dtype=float)
    query = np.asarray(query_energy_J, dtype=float)
    if use_log_interpolation:
        log_values = np.interp(np.log(query), np.log(energies), np.log(table_values))
        return np.exp(log_values)

    return np.interp(query, energies, table_values)

def evaluate_cross_section_m2(table: TabulatedCrossSection, energy_J: ArrayLike):
    """Evaluate σ(E) [m^2]"""
    return interpolate_table_values(table.energies_J, table.cross_sections_m2, energy_J, use_log_interpolation=table.use_log_interpolation)

def evaluate_rate_coefficient_m3_s(table: TabulatedRateCoefficient, energy_J: ArrayLike):
    """Evaluate K(E) [m^3/s]"""
    return interpolate_table_values(table.energies_J, table.rate_coefficients_m3_s, energy_J, use_log_interpolation=table.use_log_interpolation)

def rate_coefficient_to_equivalent_cross_section_m2(rate_coefficient_m3_s: ArrayLike, beam_speed_m_s: ArrayLike):
    """Equivalent path length cross section σ_eff = K / v_b"""
    K = np.asarray(rate_coefficient_m3_s, dtype=float)
    speed = np.asarray(beam_speed_m_s, dtype=float)

    return K / speed

def evaluate_channel_cross_section_m2(channel_data: AtomicChannelData | None, energy_J: float, beam_particle_mass_kg: float):
    """Evaluate a channels effective cross section at beam energy"""
    if channel_data is None:
        return 0.0
    sigma_total = 0.0
    if channel_data.cross_section is not None:
        sigma_total += evaluate_cross_section_m2(channel_data.cross_section, energy_J)
    if channel_data.rate_coefficient is not None:
        speed = beam_speed_from_energy_m_s(energy_J, beam_particle_mass_kg)
        rate_coefficient = evaluate_rate_coefficient_m3_s(channel_data.rate_coefficient, energy_J)
        sigma_total += rate_coefficient_to_equivalent_cross_section_m2(rate_coefficient, speed)

    return sigma_total

def attenuation_rate_from_cross_section_per_m(target_density_profile_m3: ArrayLike, cross_section_m2: ArrayLike):
    """κ = n * σ"""
    density = np.asarray(target_density_profile_m3, dtype=float)
    sigma = np.asarray(cross_section_m2, dtype=float)

    return density * sigma

def target_attenuation_rates_per_m(target_profile: AtomicTargetProfile, energy_J: float, beam_particle_mass_kg: float):
    """Channel attenuation rates from one target density profile"""
    density = np.asarray(target_profile.density_profile_m3, dtype=float)
    sigma_ion = evaluate_channel_cross_section_m2(target_profile.ionization, energy_J, beam_particle_mass_kg)
    sigma_cx = evaluate_channel_cross_section_m2(target_profile.charge_exchange, energy_J, beam_particle_mass_kg)
    sigma_other = evaluate_channel_cross_section_m2(target_profile.other_loss, energy_J, beam_particle_mass_kg)

    return BeamAtomicAttenuationRates(
        ionization_rate_per_m=attenuation_rate_from_cross_section_per_m(density, sigma_ion),
        charge_exchange_rate_per_m=attenuation_rate_from_cross_section_per_m(density, sigma_cx),
        other_loss_rate_per_m=attenuation_rate_from_cross_section_per_m(density, sigma_other),
    )

def sum_beam_atomic_attenuation_rates(rate_sets: tuple[BeamAtomicAttenuationRates, ...] | list[BeamAtomicAttenuationRates]):
    """Sum channel attenuation rates from several target species"""
    ion = np.sum([rates.ionization_rate_per_m for rates in rate_sets], axis=0)
    cx = np.sum([rates.charge_exchange_rate_per_m for rates in rate_sets], axis=0)
    other = np.sum([rates.other_loss_rate_per_m for rates in rate_sets], axis=0)

    return BeamAtomicAttenuationRates(ionization_rate_per_m=ion, charge_exchange_rate_per_m=cx, other_loss_rate_per_m=other,)

def beam_atomic_attenuation_rates_per_m(target_profiles: tuple[AtomicTargetProfile, ...] | list[AtomicTargetProfile], energy_J: float, beam_particle_mass_kg: float):
    """Total channel attenuation rates from all supplied targets"""
    rate_sets = tuple(target_attenuation_rates_per_m(target_profile=target, energy_J=energy_J, beam_particle_mass_kg=beam_particle_mass_kg,) for target in target_profiles)
    return sum_beam_atomic_attenuation_rates(rate_sets)

def multi_energy_atomic_attenuation_rates_per_m(target_profiles: tuple[AtomicTargetProfile, ...] | list[AtomicTargetProfile], component_energies_J: ArrayLike, beam_particle_mass_kg: float):
    """
    Evaluate aggregate attenuation coefficients for several beam component energies
    
    The returned arrays have shape component by target profile cell
    """
    energies = np.asarray(component_energies_J, dtype=float)
    rates = tuple(beam_atomic_attenuation_rates_per_m(target_profiles=target_profiles, energy_J=float(energy), beam_particle_mass_kg=beam_particle_mass_kg,) for energy in energies)
    
    return BeamAtomicAttenuationRates(
        ionization_rate_per_m=np.asarray([rate.ionization_rate_per_m for rate in rates], dtype=float),
        charge_exchange_rate_per_m=np.asarray([rate.charge_exchange_rate_per_m for rate in rates], dtype=float),
        other_loss_rate_per_m=np.asarray([rate.other_loss_rate_per_m for rate in rates], dtype=float),
    )

def atomic_target_profile(density_profile_m3: ArrayLike, ionization: AtomicChannelData | None = None, charge_exchange: AtomicChannelData | None = None, other_loss: AtomicChannelData | None = None):
    """Build an atomic target profile"""
    return AtomicTargetProfile(density_profile_m3=np.asarray(density_profile_m3, dtype=float), ionization=ionization, charge_exchange=charge_exchange, other_loss=other_loss)