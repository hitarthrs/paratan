"""
Non Maxwellian deuterium and tritium target rates for primary neutral beam attenuation

The heavy target rate replaces a Maxwellian ion approximation with solved local gyrotropic D and T distributions
For each beam component the code gyroaverages sigma times relative speed and integrates that kernel over the local target velocity distribution
Electron impact ionization remains a Maxwellian electron contribution

The resulting reaction rate divided by beam speed is a line attenuation coefficient in m⁻¹
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.constants import atomic_mass
from source_model_revamp.beam.atomic_rates import BeamAtomicAttenuationRates
from source_model_revamp.beam.hydrogenic_dt_atomic_data import electron_maxwellian_ionization_rate_coefficient_m3_s, hydrogenic_ground_state_ion_impact_ionization_cross_section_m2, janev_ground_state_charge_exchange_cross_section_m2
from source_model_revamp.constants import KEV_TO_J
from source_model_revamp.fbis.beam_source_definition import beam_speed_from_energy_m_s
from source_model_revamp.fbis.species import DEUTERON, TRITON, IonSpecies
from source_model_revamp.fbis.velocity_space_grid import gyrotropic_velocity_cell_volumes
from source_model_revamp.integration.pipeline_types import KineticStageResult

@dataclass(frozen=True)
class KineticTargetAtomicRateState:
    """
    Resolved attenuation state evaluated against solved local ion distributions
    
    rates contains aggregate and species resolved line coefficients
    Deuterium and tritium density profiles are velocity integrals of the supplied kinetic distributions
    Charge exchange mean energies are reaction weighted target kinetic energies for each beam component and axial cell
    """
    rates: BeamAtomicAttenuationRates
    deuterium_density_profile_m3: np.ndarray
    tritium_density_profile_m3: np.ndarray
    deuterium_charge_exchange_target_mean_energy_J: np.ndarray
    tritium_charge_exchange_target_mean_energy_J: np.ndarray
    target_model: str
    gyroangle_points: int

def kinetic_species_distribution(kinetic: KineticStageResult, species_id: str, axial_count: int) -> tuple[object, np.ndarray] | None:
    """
    Return one solved local gyrotropic species distribution when available
    
    The required array shape is axial cell by local speed cell by pitch cell
    Missing grids or distributions return None so absent target species contribute no heavy particle rate
    """
    grids = kinetic.local_speed_grid_by_species
    distributions = kinetic.local_distribution_z_v_pitch_by_species
    if grids is None or distributions is None:
        return None
    grid = grids.get(species_id)
    values = distributions.get(species_id)
    if grid is None or values is None:
        return None
    array = np.asarray(values, dtype=float)
    expected = (axial_count, grid.centers_m_s.size, kinetic.pitch_grid.centers.size)
    if array.shape != expected:
        raise ValueError(f"kinetic {species_id} target distribution must have shape {expected}")
    if np.any(~np.isfinite(array)) or np.any(array < 0.0):
        raise ValueError(f"kinetic {species_id} target distribution must be finite and nonnegative")
  
    return grid, array

@lru_cache(maxsize=128)
def _gyroaveraged_sigma_g_cached(speed_centers_m_s: tuple[float, ...], pitch_centers: tuple[float, ...], beam_speed_m_s: float, beam_axial_cosine: float, channel: str, gyroangle_points: int) -> np.ndarray:
    """
    Build an immutable gyrophase averaged sigma times relative speed kernel
    
    The target pitch coordinate is the cosine relative to the magnetic axis
    Gyrotropy allows the beam transverse direction to define the zero gyrophase without loss of generality
    The returned array has shape target speed by target pitch
    """
    target_speed = np.asarray(speed_centers_m_s, dtype=float)[:, None, None]
    target_pitch = np.asarray(pitch_centers, dtype=float)[None, :, None]
    # Midpoint gyrophase samples average the gyrotropic target around the magnetic axis
    phi = (np.arange(gyroangle_points, dtype=float) + 0.5) * (2.0 * np.pi / gyroangle_points)
    transverse_sine = float(np.sqrt(max(0.0, 1.0 - beam_axial_cosine**2)))
    direction_cosine = target_pitch * beam_axial_cosine + np.sqrt(np.maximum(0.0, 1.0 - target_pitch**2)) * transverse_sine * np.cos(phi)[None, None, :]
    relative_speed = np.sqrt(np.maximum(0.0, beam_speed_m_s**2 + target_speed**2 - 2.0 * beam_speed_m_s * target_speed * direction_cosine))
    relative_energy_keV_per_u = 0.5 * atomic_mass * relative_speed**2 / KEV_TO_J
    if channel == "ionization":
        sigma = hydrogenic_ground_state_ion_impact_ionization_cross_section_m2(relative_energy_keV_per_u)
    elif channel == "charge_exchange":
        sigma = janev_ground_state_charge_exchange_cross_section_m2(relative_energy_keV_per_u)
    else:
        raise ValueError("channel must be ionization or charge_exchange")
    kernel = np.mean(np.asarray(sigma, dtype=float) * relative_speed, axis=2)
    if np.any(~np.isfinite(kernel)) or np.any(kernel < 0.0):
        raise ValueError("kinetic target sigma g kernel must be finite and nonnegative")
    kernel.setflags(write=False)
   
    return kernel

def gyroaveraged_sigma_g(speed_grid: object, pitch_centers: np.ndarray, beam_speed_m_s: float, beam_axial_cosine: float, channel: str, gyroangle_points: int) -> np.ndarray:
    """Validate inputs and return the cached gyrophase averaged sigma times relative speed kernel"""
    speed_values = np.asarray(speed_grid.centers_m_s, dtype=float)
    pitch_values = np.asarray(pitch_centers, dtype=float)
    speed = float(beam_speed_m_s)
    axial_cosine = float(np.clip(beam_axial_cosine, -1.0, 1.0))
    count = int(gyroangle_points)
    if speed_values.ndim != 1 or np.any(~np.isfinite(speed_values)) or np.any(speed_values < 0.0):
        raise ValueError("speed grid centers must be finite and nonnegative")
    if pitch_values.ndim != 1 or np.any(~np.isfinite(pitch_values)) or np.any(np.abs(pitch_values) > 1.0):
        raise ValueError("pitch centers must be finite and within minus one to one")
    if not np.isfinite(speed) or speed <= 0.0:
        raise ValueError("beam speed must be positive and finite")
    if count < 1:
        raise ValueError("gyroangle points must be positive")
  
    return _gyroaveraged_sigma_g_cached(tuple(float(value) for value in speed_values), tuple(float(value) for value in pitch_values), speed, axial_cosine, str(channel), count)

def _species_rates(distribution_z_v_pitch: np.ndarray, speed_grid: object, pitch_grid: object, component_energies_J: np.ndarray, projectile_species: IonSpecies, target_species: IonSpecies, beam_axial_cosine: float, gyroangle_points: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Integrate ionization and charge exchange coefficients over one solved target species
    
    The gyrotropic velocity cell measure converts the local distribution to density
    For each beam energy the velocity integral of f sigma g divided by beam speed gives attenuation in m⁻¹
    The charge exchange target energy is averaged with the same reaction weight
    """
    weights = gyrotropic_velocity_cell_volumes(speed_grid, pitch_grid)
    density = np.sum(distribution_z_v_pitch * weights[None, :, :], axis=(1, 2))
    ionization = np.zeros((component_energies_J.size, distribution_z_v_pitch.shape[0]), dtype=float)
    charge_exchange = np.zeros_like(ionization)
    charge_exchange_mean_energy = np.zeros_like(ionization)
    target_energy = 0.5 * target_species.mass_kg * np.asarray(speed_grid.centers_m_s, dtype=float) ** 2
    for index, energy_J in enumerate(component_energies_J):
        beam_speed = float(beam_speed_from_energy_m_s(float(energy_J), projectile_species.mass_kg))
        ion_kernel = gyroaveraged_sigma_g(speed_grid, pitch_grid.centers, beam_speed, beam_axial_cosine, "ionization", gyroangle_points)
        cx_kernel = gyroaveraged_sigma_g(speed_grid, pitch_grid.centers, beam_speed, beam_axial_cosine, "charge_exchange", gyroangle_points)
        ion_weight = distribution_z_v_pitch * weights[None, :, :] * ion_kernel[None, :, :]
        cx_weight = distribution_z_v_pitch * weights[None, :, :] * cx_kernel[None, :, :]
        # The velocity integral gives events per volume per time and division by beam speed gives attenuation per length
        ionization[index] = np.sum(ion_weight, axis=(1, 2)) / beam_speed
        charge_exchange[index] = np.sum(cx_weight, axis=(1, 2)) / beam_speed
        cx_normalization = np.sum(cx_weight, axis=(1, 2))
        cx_energy_numerator = np.sum(cx_weight * target_energy[None, :, None], axis=(1, 2))
        charge_exchange_mean_energy[index] = np.divide(cx_energy_numerator, cx_normalization, out=np.zeros_like(cx_normalization), where=cx_normalization > 0.0)
   
    return ionization, charge_exchange, density, charge_exchange_mean_energy

def kinetic_target_atomic_rates_per_m(*, component_energies_J: np.ndarray, projectile_species: IonSpecies, electron_density_profile_m3: np.ndarray, electron_temperature_J: float, kinetic_target_state: KineticStageResult, beam_direction_unit: np.ndarray, gyroangle_points: int) -> KineticTargetAtomicRateState:
    """
    Evaluate primary beam attenuation against solved local D and T kinetic states
    
    Electron impact ionization uses the Maxwellian electron rate coefficient on the supplied electron density profile
    Available D and T distributions provide non Maxwellian ionization charge exchange density and target energy moments
    Returned heavy particle arrays have shape component by axial cell
    """
    energies = np.asarray(component_energies_J, dtype=float)
    electron_density = np.asarray(electron_density_profile_m3, dtype=float)
    direction = np.asarray(beam_direction_unit, dtype=float)
    if energies.ndim != 1 or np.any(~np.isfinite(energies)) or np.any(energies <= 0.0):
        raise ValueError("component_energies_J must contain positive finite values")
    if electron_density.ndim != 1 or np.any(~np.isfinite(electron_density)) or np.any(electron_density < 0.0):
        raise ValueError("electron_density_profile_m3 must be finite and nonnegative")
    if direction.shape != (3,) or np.any(~np.isfinite(direction)) or np.linalg.norm(direction) <= 0.0:
        raise ValueError("beam_direction_unit must be one finite nonzero vector")
    nphi = int(gyroangle_points)
    if nphi < 1:
        raise ValueError("gyroangle_points must be positive")
    axial_cosine = float(direction[2] / np.linalg.norm(direction))
    electron_coefficient = electron_maxwellian_ionization_rate_coefficient_m3_s(float(electron_temperature_J))
    electron_ionization = np.zeros((energies.size, electron_density.size), dtype=float)
    deuterium_ionization = np.zeros_like(electron_ionization)
    tritium_ionization = np.zeros_like(electron_ionization)
    deuterium_charge_exchange = np.zeros_like(electron_ionization)
    tritium_charge_exchange = np.zeros_like(electron_ionization)
    deuterium_charge_exchange_mean_energy = np.zeros_like(electron_ionization)
    tritium_charge_exchange_mean_energy = np.zeros_like(electron_ionization)
    deuterium_density = np.zeros(electron_density.shape, dtype=float)
    tritium_density = np.zeros(electron_density.shape, dtype=float)
    for index, energy_J in enumerate(energies):
        beam_speed = float(beam_speed_from_energy_m_s(float(energy_J), projectile_species.mass_kg))
        electron_ionization[index] = electron_density * electron_coefficient / beam_speed
    deuterium_state = kinetic_species_distribution(kinetic_target_state, DEUTERON.species_id, electron_density.size)
    if deuterium_state is not None:
        deuterium_grid, deuterium_distribution = deuterium_state
        deuterium_ionization, deuterium_charge_exchange, deuterium_density, deuterium_charge_exchange_mean_energy = _species_rates(deuterium_distribution, deuterium_grid, kinetic_target_state.pitch_grid, energies, projectile_species, DEUTERON, axial_cosine, nphi)
    tritium_state = kinetic_species_distribution(kinetic_target_state, TRITON.species_id, electron_density.size)
    if tritium_state is not None:
        tritium_grid, tritium_distribution = tritium_state
        tritium_ionization, tritium_charge_exchange, tritium_density, tritium_charge_exchange_mean_energy = _species_rates(tritium_distribution, tritium_grid, kinetic_target_state.pitch_grid, energies, projectile_species, TRITON, axial_cosine, nphi)
    ionization = electron_ionization + deuterium_ionization + tritium_ionization
    charge_exchange = deuterium_charge_exchange + tritium_charge_exchange
    zeros = np.zeros_like(ionization)
    
    return KineticTargetAtomicRateState(
        rates=BeamAtomicAttenuationRates(
            ionization_rate_per_m=ionization,
            charge_exchange_rate_per_m=charge_exchange,
            other_loss_rate_per_m=zeros,
            electron_ionization_rate_per_m=electron_ionization,
            deuterium_ionization_rate_per_m=deuterium_ionization,
            tritium_ionization_rate_per_m=tritium_ionization,
            deuterium_charge_exchange_rate_per_m=deuterium_charge_exchange,
            tritium_charge_exchange_rate_per_m=tritium_charge_exchange,
        ),
        deuterium_density_profile_m3=deuterium_density,
        tritium_density_profile_m3=tritium_density,
        deuterium_charge_exchange_target_mean_energy_J=deuterium_charge_exchange_mean_energy,
        tritium_charge_exchange_target_mean_energy_J=tritium_charge_exchange_mean_energy,
        target_model="solved_local_gyrotropic_D_T_sigma_g_average",
        gyroangle_points=nphi,
    )

__all__ = ["KineticTargetAtomicRateState", "gyroaveraged_sigma_g", "kinetic_species_distribution", "kinetic_target_atomic_rates_per_m"]