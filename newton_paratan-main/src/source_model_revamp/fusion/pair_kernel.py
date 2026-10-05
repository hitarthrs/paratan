"""Bosch Hale Table IV reactant pair kernel and Table VII thermal validation"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.mirror_invariant_grid import LambdaGrid,require_mirror_distribution_shape
from source_model_revamp.fbis.mirror_pitch_remap import conservative_local_pitch_distribution_from_lambda_distribution
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid, isotropic_maxwellian_distribution
from source_model_revamp.fusion.bosch_hale import bosch_hale_table_iv_fit_domain_keV
from source_model_revamp.fusion.cross_sections import bosch_hale_cross_section_m2_from_J, bosch_hale_thermal_reactivity_for_reaction_m3_s, reduced_mass_kg
from source_model_revamp.fusion.reactions import J_TO_KEV, FusionReaction
from source_model_revamp.fusion.reactivity import gyrotropic_pair_integral_kernel_and_state_moments_profiles, population_pair_counting_factor, reaction_rate_density_from_reactivity

@dataclass(frozen=True)
class FusionPairKernelResult:
    """Table IV pair rate, reactant energy moments, and fit domain partition

    Axial rate and power arrays have shape (n_z,)
    """
    reaction: FusionReaction
    rate_density_m3_s: np.ndarray
    reactant_a_kinetic_energy_removal_density_W_m3: np.ndarray
    reactant_b_kinetic_energy_removal_density_W_m3: np.ndarray
    inside_fit_domain_rate_density_m3_s: np.ndarray
    outside_fit_domain_rate_density_m3_s: np.ndarray
    fit_domain_E_cm_keV: tuple[float, float]
    identical_population: bool
    num_gyroangle_points: int
    max_pair_states_per_batch: int
    kernel_model: str = "bosch_hale_table_iv_total_cross_section_quadrature"

    @property
    def outside_fit_domain_fraction_profile(self) -> np.ndarray:
        """Fraction of the local reaction rate from E_cm outside the Table IV fit interval"""
        return np.divide(self.outside_fit_domain_rate_density_m3_s, self.rate_density_m3_s, out=np.zeros_like(self.rate_density_m3_s), where=self.rate_density_m3_s > 0.0)

@dataclass(frozen=True)
class ThermalReactivityValidationResult:
    """Table IV Maxwellian quadrature compared with the independent Table VII thermal fit"""
    reaction: FusionReaction
    table_iv_rate_density_m3_s: np.ndarray
    table_vii_rate_density_m3_s: np.ndarray
    relative_difference_profile: np.ndarray
    temperature_keV_profile: np.ndarray
    table_vii_temperature_domain_keV: tuple[float, float]
    table_vii_domain_covered: bool

def _broadcast_profile(values: ArrayLike, size: int, name: str) -> np.ndarray:
    """Broadcast one scalar or validate one axial profile of the requested length"""
    array = np.asarray(values, dtype=float)
    if array.ndim == 0:
        result = np.full(size, float(array), dtype=float)
    elif array.ndim == 1 and array.size == size:
        result = array.astype(float, copy=False)
    else:
        raise ValueError(f"{name} must be scalar or have length {size}")
    if np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    
    return result

def isotropic_maxwellian_gyrotropic_stack(speed_grid: SpeedGrid, pitch_grid: PitchGrid, density_profile_m3: ArrayLike, temperature_energy_J: ArrayLike, particle_mass_kg: float) -> np.ndarray:
    """Build isotropic Maxwellian f(z, v, ξ) with shape (n_z, n_speed, n_pitch)"""
    density = np.asarray(density_profile_m3, dtype=float)
    if density.ndim == 0:
        density = density[None]
    if density.ndim != 1 or density.size == 0:
        raise ValueError("density_profile_m3 must be a nonempty 1D profile")
    if np.any(~np.isfinite(density)) or np.any(density < 0.0):
        raise ValueError("density_profile_m3 must be finite and nonnegative")
    temperature = _broadcast_profile(temperature_energy_J, density.size, "temperature_energy_J")
    if np.any(temperature <= 0.0):
        raise ValueError("temperature_energy_J must be positive")
    
    return np.stack([np.repeat(isotropic_maxwellian_distribution(speed_m_s=speed_grid.centers_m_s, density_m3=float(density[index]), temperature_energy_J=float(temperature[index]), particle_mass_kg=float(particle_mass_kg))[:, None], pitch_grid.centers.size, axis=1) for index in range(density.size)], axis=0,)

def local_pitch_stack_from_lambda_distribution(speed_grid: SpeedGrid, lambda_grid: LambdaGrid, distribution_v_lambda: ArrayLike, pitch_grid: PitchGrid, B_tilde_profile: ArrayLike) -> np.ndarray:
    """Conservatively map one f(v, Λ) invariant distribution to local f(z, v, ξ) planes"""
    distribution = require_mirror_distribution_shape(distribution_v_lambda, speed_grid, lambda_grid, name="distribution_v_lambda")
    if np.any(distribution < 0.0):
        raise ValueError("distribution_v_lambda must be nonnegative")
    B_profile = np.asarray(B_tilde_profile, dtype=float)
    if (B_profile.ndim != 1 or B_profile.size == 0 or np.any(~np.isfinite(B_profile)) or np.any(B_profile <= 0.0)):
        raise ValueError("B_tilde_profile must be a finite positive 1D profile")
    
    return np.stack([conservative_local_pitch_distribution_from_lambda_distribution(lambda_grid=lambda_grid, distribution_v_lambda=distribution, pitch_grid=pitch_grid, B_tilde=float(B_tilde), fill_value=0.0) for B_tilde in B_profile], axis=0,)

def bosch_hale_table_iv_pair_kernel_profiles(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distributions_a_z_v_xi: ArrayLike, mass_a_kg: float, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distributions_b_z_v_xi: ArrayLike, mass_b_kg: float, reaction: FusionReaction, *, identical_population: bool = False, num_gyroangle_points: int = 32, max_pair_states_per_batch: int = 262144) -> FusionPairKernelResult:
    """Evaluate the Table IV pair rate, reactant kinetic energy moments, and fit domain coverage

    The total cross section is evaluated at E_cm = 0.5 * μ * g^2
    identical_population applies the 1/2 pair counting factor after the shared quadrature
    """
    reduced_mass = reduced_mass_kg(mass_a_kg, mass_b_kg)
    fit_min_keV, fit_max_keV = bosch_hale_table_iv_fit_domain_keV(reaction.key)

    def table_iv_sigma_g(relative_speed):
        """Evaluate the Table IV σ(E_cm) g kernel"""
        speed = np.asarray(relative_speed, dtype=float)
        energy_J = 0.5 * reduced_mass * speed**2
        return bosch_hale_cross_section_m2_from_J(reaction.key, energy_J) * speed

    def inside_fit_domain_sigma_g(relative_speed):
        """Evaluate σ(E_cm) g only inside the stated Table IV E_cm interval"""
        speed = np.asarray(relative_speed, dtype=float)
        energy_J = 0.5 * reduced_mass * speed**2
        energy_keV = energy_J * J_TO_KEV
        inside = (energy_keV >= fit_min_keV) & (energy_keV <= fit_max_keV)
        return table_iv_sigma_g(speed) * inside

    def outside_fit_domain_sigma_g(relative_speed):
        """Evaluate σ(E_cm) g only outside the stated Table IV E_cm interval"""
        speed = np.asarray(relative_speed, dtype=float)
        energy_J = 0.5 * reduced_mass * speed**2
        energy_keV = energy_J * J_TO_KEV
        outside = (energy_keV < fit_min_keV) | (energy_keV > fit_max_keV)
        return table_iv_sigma_g(speed) * outside

    energy_a_J = 0.5 * float(mass_a_kg) * np.asarray(speed_grid_a.centers_m_s, dtype=float) ** 2
    energy_b_J = 0.5 * float(mass_b_kg) * np.asarray(speed_grid_b.centers_m_s, dtype=float) ** 2
    raw, moments_a, moments_b = gyrotropic_pair_integral_kernel_and_state_moments_profiles(speed_grid_a, pitch_grid_a, distributions_a_z_v_xi, speed_grid_b, pitch_grid_b, distributions_b_z_v_xi, (table_iv_sigma_g, inside_fit_domain_sigma_g, outside_fit_domain_sigma_g), num_gyroangle_points=num_gyroangle_points, max_pair_states_per_batch=max_pair_states_per_batch, state_moment_weights_a=(energy_a_J,), state_moment_weights_b=(energy_b_J,))
    pair_factor = population_pair_counting_factor(identical_population)
    total = pair_factor * raw[0]
    reactant_a_energy = pair_factor * moments_a[0, 0]
    reactant_b_energy = pair_factor * moments_b[0, 0]
    inside = pair_factor * raw[1]
    outside = pair_factor * raw[2]
    # Shared quadrature requires inside plus outside contributions to close to the total
    closure_scale = np.maximum(total, np.finfo(float).tiny)
    closure_error = np.abs(total - inside - outside) / closure_scale
    if np.any(closure_error > 5.0e-13):
        raise RuntimeError("Table IV fit domain partition does not close")

    return FusionPairKernelResult(
        reaction=reaction,
        rate_density_m3_s=total,
        reactant_a_kinetic_energy_removal_density_W_m3=reactant_a_energy,
        reactant_b_kinetic_energy_removal_density_W_m3=reactant_b_energy,
        inside_fit_domain_rate_density_m3_s=inside,
        outside_fit_domain_rate_density_m3_s=outside,
        fit_domain_E_cm_keV=(fit_min_keV, fit_max_keV),
        identical_population=bool(identical_population),
        num_gyroangle_points=int(num_gyroangle_points),
        max_pair_states_per_batch=int(max_pair_states_per_batch),
    )

def thermal_table_vii_validation(reaction: FusionReaction, table_iv_rate_density_m3_s: ArrayLike, density_a_m3: ArrayLike, density_b_m3: ArrayLike, temperature_energy_J: ArrayLike, *, identical_population: bool) -> ThermalReactivityValidationResult:
    """Compare Table IV Maxwellian quadrature with the Table VII thermal fit

    Table VII is checked over its 0.2 to 100 keV D T and D D temperature interval
    """
    table_iv = np.asarray(table_iv_rate_density_m3_s, dtype=float)
    if table_iv.ndim != 1 or np.any(~np.isfinite(table_iv)) or np.any(table_iv < 0.0):
        raise ValueError("table_iv_rate_density_m3_s must be finite and nonnegative")
    density_a = _broadcast_profile(density_a_m3, table_iv.size, "density_a_m3")
    density_b = _broadcast_profile(density_b_m3, table_iv.size, "density_b_m3")
    temperature = _broadcast_profile(temperature_energy_J, table_iv.size, "temperature_energy_J")
    if np.any(density_a < 0.0) or np.any(density_b < 0.0):
        raise ValueError("thermal densities must be nonnegative")
    if np.any(temperature <= 0.0):
        raise ValueError("temperature_energy_J must be positive")
    table_vii_reactivity = bosch_hale_thermal_reactivity_for_reaction_m3_s(reaction.key, temperature)
    table_vii = reaction_rate_density_from_reactivity(density_a, density_b, table_vii_reactivity, identical_population=identical_population)
    denominator = np.maximum(table_vii, np.finfo(float).tiny)
    relative_difference = (table_iv - table_vii) / denominator
    temperature_keV = temperature * J_TO_KEV
    domain = (0.2, 100.0)

    return ThermalReactivityValidationResult(
        reaction=reaction,
        table_iv_rate_density_m3_s=table_iv,
        table_vii_rate_density_m3_s=table_vii,
        relative_difference_profile=relative_difference,
        temperature_keV_profile=temperature_keV,
        table_vii_temperature_domain_keV=domain,
        table_vii_domain_covered=bool(np.all((temperature_keV >= domain[0]) & (temperature_keV <= domain[1])))
    )