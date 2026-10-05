"""
Qualified ground state atomic model for neutral deuterium and tritium beam attenuation

Ion impact ionization uses the embedded CollisionDB D107356 ground state hydrogen data
Charge exchange and electron impact ionization use the analytic ground state fits already encoded in this file
Heavy target rates average sigma times relative speed over Maxwellian deuterium or tritium targets
Electron ionization uses a Maxwellian electron average with the neutral treated as stationary for that average
"""
from __future__ import annotations
from functools import lru_cache
from io import StringIO
import hashlib
import numpy as np
from numpy.typing import ArrayLike
from scipy.constants import atomic_mass
from scipy.integrate import quad
from source_model_revamp.beam.atomic_rates import BeamAtomicAttenuationRates
from source_model_revamp.constants import ELECTRON_MASS_KG, KEV_TO_J
from source_model_revamp.fbis.beam_source_definition import beam_speed_from_energy_m_s
from source_model_revamp.fbis.species import DEUTERON, TRITON, IonSpecies

HYDROGENIC_DT_GROUND_STATE_MODEL = "hydrogenic_dt_ground_state"
COLLISIONDB_D107356_QID = "D107356"
COLLISIONDB_D107356_DOI = "10.1140/epjd/e2019-100380-x"
COLLISIONDB_D107356_METHOD = "TC-BGM"
COLLISIONDB_D107356_ENERGY_FRAME = "target"
COLLISIONDB_D107356_ENERGY_UNITS = "eV/u"
COLLISIONDB_D107356_CROSS_SECTION_UNITS = "cm2"
HYDROGENIC_DT_GROUND_STATE_REFERENCE = ("CollisionDB_D107356_Leung_Kirchner_2019_H1s_proton_ionization, Janev_Reiter_Samm_2003_ground_state_H_H_charge_exchange_and_electron_ionization")
SUPPORTED_COMPONENT_ENERGY_MIN_KEV = 10.0
SUPPORTED_COMPONENT_ENERGY_MAX_KEV = 200.0
SUPPORTED_ION_TEMPERATURE_MAX_KEV = 10.0
COLLISIONDB_D107356_MIN_KEV_PER_U = 1.0
COLLISIONDB_D107356_MAX_KEV_PER_U = 300.0
COLLISIONDB_D107356_DATA_SHA256 = "9205bf77d754e20b7c36dc0413cf0f6a58e51cc975c22e68719bef782a7188f0"
COLLISIONDB_D107356_SOURCE_PDF_SHA256 = "ca3531c6e892fd154be62788de96296f05745dc6d23e59695f93866891fd2656"
IONIZATION_LOW_ENERGY_TRUNCATION_MAX_RELATIVE_RATE = 3.0e-7
IONIZATION_HIGH_ENERGY_TRUNCATION_MAX_RELATIVE_RATE = 3.0e-6
_COLLISIONDB_D107356_DATA = """1.00e+03 3.91e-20
1.50e+03 1.7e-19
1.75e+03 2.29e-19
2.00e+03 3.14e-19
2.50e+03 6.18e-19
3.00e+03 1.11e-18
4.00e+03 2.59e-18
5.00e+03 4.38e-18
6.00e+03 6.51e-18
8.00e+03 1.28e-17
1.00e+04 1.99e-17
1.20e+04 2.78e-17
1.50e+04 4.13e-17
2.00e+04 6.84e-17
3.00e+04 1.23e-16
4.00e+04 1.59e-16
5.00e+04 1.72e-16
6.00e+04 1.73e-16
8.00e+04 1.6e-16
1.00e+05 1.39e-16
1.50e+05 1.05e-16
2.00e+05 8.24e-17
3.00e+05 5.79e-17
"""

if hashlib.sha256(_COLLISIONDB_D107356_DATA.encode()).hexdigest() != COLLISIONDB_D107356_DATA_SHA256:
    raise RuntimeError("CollisionDB D107356 embedded data hash mismatch")

_COLLISIONDB_D107356_TABLE = np.loadtxt(StringIO(_COLLISIONDB_D107356_DATA), dtype=float)
_COLLISIONDB_D107356_ENERGIES_KEV_PER_U = 1.0e-3 * _COLLISIONDB_D107356_TABLE[:, 0]
_COLLISIONDB_D107356_CROSS_SECTIONS_M2 = 1.0e-4 * _COLLISIONDB_D107356_TABLE[:, 1]
_JANEV_CHARGE_EXCHANGE_COEFFICIENTS = np.asarray([3.2345, 2.3588e2, 2.3713, 3.8371e-2, 3.8068e-6, 1.1832e-10], dtype=float)
_JANEV_ELECTRON_IONIZATION_COEFFICIENTS = np.asarray([0.18450, -0.032226, -0.034539, 1.4003, -2.8115, 2.2986], dtype=float)

def _finite_positive_scalar(name: str, value: float) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    
    return result

def _finite_nonnegative_profile(name: str, values: ArrayLike) -> np.ndarray:
    result = np.asarray(values, dtype=float)
    if result.ndim != 1:
        raise ValueError(f"{name} must be one dimensional")
    if np.any(~np.isfinite(result)) or np.any(result < 0.0):
        raise ValueError(f"{name} must be finite and nonnegative")
    
    return result

def _supported_species(species: IonSpecies) -> IonSpecies:
    if species.species_id not in {DEUTERON.species_id, TRITON.species_id}:
        raise ValueError("neutral beam projectile species must be deuterium or tritium")
    
    return species

def _scalar_or_array(original: ArrayLike, values: np.ndarray) -> float | np.ndarray:
    if np.asarray(original).ndim == 0:
        return float(values.reshape(-1)[0])
    
    return values


def collisiondb_d107356_ground_state_proton_impact_ionization_cross_section_m2(relative_energy_keV_per_u: ArrayLike) -> float | np.ndarray:
    """
    Evaluate the embedded CollisionDB D107356 ground state ionization cross section
    
    Input is relative collision energy in keV per atomic mass unit
    The table is interpolated in logarithmic energy and cross section coordinates
    Cross section is set to zero outside the embedded 1 to 300 keV per atomic mass unit interval
    """
    energy_input = np.asarray(relative_energy_keV_per_u, dtype=float)
    energy = np.atleast_1d(energy_input)
    if np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
        raise ValueError("relative_energy_keV_per_u must be finite and nonnegative")
    result = np.zeros_like(energy, dtype=float)
    active = (energy >= _COLLISIONDB_D107356_ENERGIES_KEV_PER_U[0]) & (energy <= _COLLISIONDB_D107356_ENERGIES_KEV_PER_U[-1])
    if np.any(active):
        result[active] = np.exp(np.interp(np.log(energy[active]), np.log(_COLLISIONDB_D107356_ENERGIES_KEV_PER_U), np.log(_COLLISIONDB_D107356_CROSS_SECTIONS_M2)))

    return _scalar_or_array(relative_energy_keV_per_u, result)

def hydrogenic_ground_state_ion_impact_ionization_cross_section_m2(relative_energy_keV_per_u: ArrayLike) -> float | np.ndarray:
    """Return the selected hydrogenic ground state ion impact ionization cross section"""
    return collisiondb_d107356_ground_state_proton_impact_ionization_cross_section_m2(relative_energy_keV_per_u)

def janev_ground_state_charge_exchange_cross_section_m2(relative_energy_keV_per_u: ArrayLike) -> float | np.ndarray:
    """
    Evaluate the implemented ground state hydrogenic charge exchange fit
    
    Input is relative collision energy in keV per atomic mass unit
    Zero relative energy returns zero cross section
    """
    energy_input = np.asarray(relative_energy_keV_per_u, dtype=float)
    energy = np.atleast_1d(energy_input)
    if np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
        raise ValueError("relative_energy_keV_per_u must be finite and nonnegative")
    result = np.zeros_like(energy, dtype=float)
    active = energy > 0.0
    if np.any(active):
        work = energy[active]
        a1, a2, a3, a4, a5, a6 = _JANEV_CHARGE_EXCHANGE_COEFFICIENTS
        result[active] = 1.0e-20 * a1 * np.log(a2 / work + a3) / (1.0 + a4 * work + a5 * work**3.5 + a6 * work**5.4)

    return _scalar_or_array(relative_energy_keV_per_u, result)

def janev_ground_state_electron_impact_ionization_cross_section_m2(collision_energy_keV: ArrayLike) -> float | np.ndarray:
    """
    Evaluate the implemented ground state electron impact ionization fit
    
    Input collision energy is in keV
    The cross section is zero at and below the 13.6 eV ionization threshold
    """
    energy_input = np.asarray(collision_energy_keV, dtype=float)
    energy = np.atleast_1d(energy_input)
    if np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
        raise ValueError("collision_energy_keV must be finite and nonnegative")
    result = np.zeros_like(energy, dtype=float)
    threshold_eV = 13.6
    energy_eV = energy * 1.0e3
    active = energy_eV > threshold_eV
    if np.any(active):
        work = energy_eV[active]
        x = 1.0 - threshold_eV / work
        a1, a2, a3, a4, a5, a6 = _JANEV_ELECTRON_IONIZATION_COEFFICIENTS
        polynomial = a2 * x + a3 * x**2 + a4 * x**3 + a5 * x**4 + a6 * x**5
        sigma_cm2 = 1.0e-13 * (a1 * np.log(work / threshold_eV) + polynomial) / (threshold_eV * work)
        result[active] = np.maximum(sigma_cm2, 0.0) * 1.0e-4

    return _scalar_or_array(collision_energy_keV, result)

def _relative_speed_probability_factor(y: float, q: float) -> float:
    """
    Return the dimensionless relative speed probability factor for a fixed projectile and Maxwellian target
    
    y is relative speed divided by target thermal speed and q is beam speed divided by target thermal speed
    """
    if y <= 0.0:
        return 0.0
    return float(y * np.exp(-(y - q) ** 2) * (-np.expm1(-4.0 * y * q)) / (np.sqrt(np.pi) * q))

@lru_cache(maxsize=512)
def heavy_particle_maxwellian_rate_coefficient_m3_s(projectile_energy_J: float, projectile_mass_kg: float, target_mass_kg: float, target_temperature_J: float, channel: str) -> float:
    """
    Return the Maxwellian target average of sigma times relative speed for one heavy particle channel
    
    The target thermal speed is sqrt of 2 T divided by target mass
    The cross section is evaluated from relative energy per atomic mass unit
    Ionization quadrature includes the embedded CollisionDB energy nodes as integration break points
    Results are cached for repeated beam component and target states
    """
    energy = _finite_positive_scalar("projectile_energy_J", projectile_energy_J)
    projectile_mass = _finite_positive_scalar("projectile_mass_kg", projectile_mass_kg)
    target_mass = _finite_positive_scalar("target_mass_kg", target_mass_kg)
    temperature = _finite_positive_scalar("target_temperature_J", target_temperature_J)
    channel_name = str(channel).strip().lower()
    if channel_name not in {"ionization", "charge_exchange"}:
        raise ValueError("channel must be ionization or charge_exchange")
    beam_speed = float(beam_speed_from_energy_m_s(energy, projectile_mass))
    target_thermal_speed = float(np.sqrt(2.0 * temperature / target_mass))
    q = beam_speed / target_thermal_speed
    upper = q + 10.0

    def integrand(y: float) -> float:
        """Return the relative speed integrand for the selected heavy particle channel"""
        relative_speed = y * target_thermal_speed
        relative_energy_keV_per_u = 0.5 * atomic_mass * relative_speed**2 / KEV_TO_J
        if channel_name == "ionization":
            sigma = hydrogenic_ground_state_ion_impact_ionization_cross_section_m2(relative_energy_keV_per_u)
        else:
            sigma = janev_ground_state_charge_exchange_cross_section_m2(relative_energy_keV_per_u)
        return _relative_speed_probability_factor(y, q) * relative_speed * float(sigma)

    # Ionization table nodes are supplied as break points because the embedded cross section is piecewise interpolated
    quadrature_points = None
    if channel_name == "ionization":
        points = []
        for threshold_keV_per_u in _COLLISIONDB_D107356_ENERGIES_KEV_PER_U:
            point = np.sqrt(2.0 * threshold_keV_per_u * KEV_TO_J / atomic_mass) / target_thermal_speed
            if 0.0 < point < upper:
                points.append(float(point))
        quadrature_points = tuple(sorted(set(points)))
    value, error = quad(integrand, 0.0, upper, epsabs=1.0e-30, epsrel=1.0e-8, limit=200, points=quadrature_points)
    if not np.isfinite(value) or value < 0.0 or not np.isfinite(error):
        raise ValueError("heavy particle rate coefficient must be finite and nonnegative")
    
    return float(value)

@lru_cache(maxsize=512)
def heavy_particle_maxwellian_target_energy_moment_J(projectile_energy_J: float, projectile_mass_kg: float, target_mass_kg: float, target_temperature_J: float, channel: str) -> float:
    """
    Return the reaction weighted mean kinetic energy of a Maxwellian target ion
    
    A 48 point angular quadrature and Maxwellian speed integral weight the target energy by sigma times relative speed
    The independently integrated rate is required to reproduce the primary heavy particle rate coefficient before the moment is returned
    """
    energy = _finite_positive_scalar("projectile_energy_J", projectile_energy_J)
    projectile_mass = _finite_positive_scalar("projectile_mass_kg", projectile_mass_kg)
    target_mass = _finite_positive_scalar("target_mass_kg", target_mass_kg)
    temperature = _finite_positive_scalar("target_temperature_J", target_temperature_J)
    channel_name = str(channel).strip().lower()
    if channel_name not in {"ionization", "charge_exchange"}:
        raise ValueError("channel must be ionization or charge_exchange")
    beam_speed = float(beam_speed_from_energy_m_s(energy, projectile_mass))
    target_thermal_speed = float(np.sqrt(2.0 * temperature / target_mass))
    mu_nodes, mu_weights = np.polynomial.legendre.leggauss(48)
    normalization = 2.0 / np.sqrt(np.pi)

    def angular_rate(x: float) -> float:
        """Return the polar angle integral of sigma times relative speed at target speed x"""
        target_speed = x * target_thermal_speed
        relative_speed = np.sqrt(np.maximum(beam_speed**2 + target_speed**2 - 2.0 * beam_speed * target_speed * mu_nodes, 0.0))
        relative_energy_keV_per_u = 0.5 * atomic_mass * relative_speed**2 / KEV_TO_J
        if channel_name == "ionization":
            sigma = np.asarray(hydrogenic_ground_state_ion_impact_ionization_cross_section_m2(relative_energy_keV_per_u), dtype=float)
        else:
            sigma = np.asarray(janev_ground_state_charge_exchange_cross_section_m2(relative_energy_keV_per_u), dtype=float)
        return float(np.sum(mu_weights * relative_speed * sigma))

    def rate_integrand(x: float) -> float:
        """Return the Maxwellian speed integrand for the reaction rate coefficient"""
        return normalization * x**2 * np.exp(-x**2) * angular_rate(x)

    def energy_integrand(x: float) -> float:
        """Return the reaction rate integrand weighted by target kinetic energy"""
        return temperature * x**2 * rate_integrand(x)

    rate, rate_error = quad(rate_integrand, 0.0, 10.0, epsabs=1.0e-30, epsrel=1.0e-8, limit=200)
    energy_rate, energy_error = quad(energy_integrand, 0.0, 10.0, epsabs=1.0e-45, epsrel=1.0e-8, limit=200)
    reference_rate = heavy_particle_maxwellian_rate_coefficient_m3_s(energy, projectile_mass, target_mass, temperature, channel_name)
    rate_relative_error = abs(rate - reference_rate) / max(abs(reference_rate), np.finfo(float).tiny)
    if not np.isfinite(rate) or rate <= 0.0 or not np.isfinite(energy_rate) or energy_rate < 0.0 or not np.isfinite(rate_error) or not np.isfinite(energy_error):
        raise ValueError("heavy particle target energy moment must be finite and nonnegative")
    if rate_relative_error > 5.0e-7:
        raise ValueError("heavy particle target energy quadrature does not reproduce the attenuation rate coefficient")
    
    return float(energy_rate / rate)

@lru_cache(maxsize=64)
def electron_maxwellian_ionization_rate_coefficient_m3_s(electron_temperature_J: float) -> float:
    """
    Return the Maxwellian electron impact ionization rate coefficient for a stationary neutral
    
    The energy integral uses collision energy equal to x times electron temperature
    Results are cached by electron temperature
    """
    temperature = _finite_positive_scalar("electron_temperature_J", electron_temperature_J)
    temperature_keV = temperature / KEV_TO_J
    speed_scale = np.sqrt(8.0 * temperature / (np.pi * ELECTRON_MASS_KG))

    def integrand(x: float) -> float:
        """Return the dimensionless Maxwellian energy integrand for electron ionization"""
        sigma = janev_ground_state_electron_impact_ionization_cross_section_m2(x * temperature_keV)
        return float(sigma) * x * np.exp(-x)

    integral, error = quad(integrand, 0.0, 80.0, epsabs=1.0e-30, epsrel=1.0e-9, limit=200)
    value = speed_scale * integral
    if not np.isfinite(value) or value < 0.0 or not np.isfinite(error):
        raise ValueError("electron ionization rate coefficient must be finite and nonnegative")
    
    return float(value)

def hydrogenic_dt_atomic_rates_per_m(*, component_energies_J: ArrayLike, projectile_species: IonSpecies, electron_density_profile_m3: ArrayLike, deuterium_density_profile_m3: ArrayLike, tritium_density_profile_m3: ArrayLike, electron_temperature_J: float, ion_temperature_J: float, model_name: str = HYDROGENIC_DT_GROUND_STATE_MODEL) -> BeamAtomicAttenuationRates:
    """
    Build resolved ground state attenuation coefficients for D or T neutral beam components
    
    Beam component energies must lie within the qualified 10 to 200 keV total energy range and ion temperature must not exceed 10 keV
    Electron ionization uses the Maxwellian electron rate coefficient
    Deuterium and tritium ionization and charge exchange use Maxwellian heavy target averages
    Dividing each volumetric reaction coefficient n K by beam speed converts it to line attenuation in m⁻¹
    Returned arrays have shape component by axial cell and other loss is zero in this model
    """
    species = _supported_species(projectile_species)
    model = str(model_name).strip().lower().replace("-", "_")
    if model != HYDROGENIC_DT_GROUND_STATE_MODEL:
        raise ValueError(f"beam atomic data model must be {HYDROGENIC_DT_GROUND_STATE_MODEL!r}")
    energies = np.atleast_1d(np.asarray(component_energies_J, dtype=float))
    if energies.ndim != 1 or np.any(~np.isfinite(energies)) or np.any(energies <= 0.0):
        raise ValueError("component_energies_J must contain positive finite values")
    energies_keV = energies / KEV_TO_J
    if np.any(energies_keV < SUPPORTED_COMPONENT_ENERGY_MIN_KEV) or np.any(energies_keV > SUPPORTED_COMPONENT_ENERGY_MAX_KEV):
        raise ValueError("beam component total energies must lie from 10 to 200 keV")
    electron_density = _finite_nonnegative_profile("electron_density_profile_m3", electron_density_profile_m3)
    deuterium_density = _finite_nonnegative_profile("deuterium_density_profile_m3", deuterium_density_profile_m3)
    tritium_density = _finite_nonnegative_profile("tritium_density_profile_m3", tritium_density_profile_m3)
    if electron_density.shape != deuterium_density.shape or electron_density.shape != tritium_density.shape:
        raise ValueError("electron and thermal ion density profiles must have matching shapes")
    electron_temperature = _finite_positive_scalar("electron_temperature_J", electron_temperature_J)
    ion_temperature = _finite_positive_scalar("ion_temperature_J", ion_temperature_J)
    if ion_temperature / KEV_TO_J > SUPPORTED_ION_TEMPERATURE_MAX_KEV:
        raise ValueError(f"ion_temperature_J exceeds the qualified {SUPPORTED_ION_TEMPERATURE_MAX_KEV:g} keV range")
    electron_rate_coefficient = electron_maxwellian_ionization_rate_coefficient_m3_s(electron_temperature)
    rate_shape = (energies.size, electron_density.size)
    electron_ionization = np.zeros(rate_shape, dtype=float)
    deuterium_ionization = np.zeros(rate_shape, dtype=float)
    tritium_ionization = np.zeros(rate_shape, dtype=float)
    deuterium_charge_exchange = np.zeros(rate_shape, dtype=float)
    tritium_charge_exchange = np.zeros(rate_shape, dtype=float)
    for index, energy_J in enumerate(energies):
        beam_speed = float(beam_speed_from_energy_m_s(float(energy_J), species.mass_kg))
        # Dividing n K by beam speed converts a temporal reaction rate to attenuation per path length
        electron_ionization[index] = electron_density * electron_rate_coefficient / beam_speed
        deuterium_ionization_coefficient = heavy_particle_maxwellian_rate_coefficient_m3_s(float(energy_J), species.mass_kg, DEUTERON.mass_kg, ion_temperature, "ionization")
        tritium_ionization_coefficient = heavy_particle_maxwellian_rate_coefficient_m3_s(float(energy_J), species.mass_kg, TRITON.mass_kg, ion_temperature, "ionization")
        deuterium_charge_exchange_coefficient = heavy_particle_maxwellian_rate_coefficient_m3_s(float(energy_J), species.mass_kg, DEUTERON.mass_kg, ion_temperature, "charge_exchange")
        tritium_charge_exchange_coefficient = heavy_particle_maxwellian_rate_coefficient_m3_s(float(energy_J), species.mass_kg, TRITON.mass_kg, ion_temperature, "charge_exchange")
        deuterium_ionization[index] = deuterium_density * deuterium_ionization_coefficient / beam_speed
        tritium_ionization[index] = tritium_density * tritium_ionization_coefficient / beam_speed
        deuterium_charge_exchange[index] = deuterium_density * deuterium_charge_exchange_coefficient / beam_speed
        tritium_charge_exchange[index] = tritium_density * tritium_charge_exchange_coefficient / beam_speed
    ionization = electron_ionization + deuterium_ionization + tritium_ionization
    charge_exchange = deuterium_charge_exchange + tritium_charge_exchange
    zeros = np.zeros_like(ionization)

    return BeamAtomicAttenuationRates(
        ionization_rate_per_m=ionization,
        charge_exchange_rate_per_m=charge_exchange,
        other_loss_rate_per_m=zeros,
        electron_ionization_rate_per_m=electron_ionization,
        deuterium_ionization_rate_per_m=deuterium_ionization,
        tritium_ionization_rate_per_m=tritium_ionization,
        deuterium_charge_exchange_rate_per_m=deuterium_charge_exchange,
        tritium_charge_exchange_rate_per_m=tritium_charge_exchange,
    )
