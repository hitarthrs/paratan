"""Gyrotropic fusion reactivity and reaction rate integrals"""
from __future__ import annotations
from collections.abc import Callable
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid, gyrotropic_velocity_cell_volumes, require_gyrotropic_distribution_shape
from source_model_revamp.fusion.cross_sections import sigma_g_from_cross_section
from source_model_revamp.fusion.reactions import FusionReaction
from source_model_revamp.velocity_geometry import relative_speed_m_s

_MAX_PAIR_STATES_PER_BATCH = 262144

def population_pair_counting_factor(identical_population: bool = False) -> float:
    """Return one half only when both reactants are drawn from one represented population"""
    return 0.5 if bool(identical_population) else 1.0

def reaction_rate_density_from_reactivity(density_a_m3: ArrayLike, density_b_m3: ArrayLike, reactivity_m3_s: ArrayLike, *, identical_population: bool = False):
    """Evaluate n_a n_b reactivity with identical population pair counting

    The supplied reactivity has units m^3/s and the result has units m^−3 s^−1
    """
    n_a = np.asarray(density_a_m3, dtype=float)
    n_b = np.asarray(density_b_m3, dtype=float)
    reactivity = np.asarray(reactivity_m3_s, dtype=float)
    if np.any(n_a < 0.0) or np.any(n_b < 0.0):
        raise ValueError("densities must be nonnegative")

    return (population_pair_counting_factor(identical_population) * n_a * n_b * reactivity)

def thermal_reaction_rate_density(density_a_m3: ArrayLike, density_b_m3: ArrayLike, temperature_energy_J: ArrayLike, thermal_reactivity_function_m3_s: Callable[[ArrayLike], ArrayLike], *, identical_population: bool = False):
    """Evaluate a Maxwellian reaction rate density with a supplied thermal reactivity"""
    reactivity = np.asarray(thermal_reactivity_function_m3_s(temperature_energy_J), dtype=float)
   
    return reaction_rate_density_from_reactivity(density_a_m3=density_a_m3, density_b_m3=density_b_m3, reactivity_m3_s=reactivity, identical_population=identical_population)

def gyrotropic_pair_integral_sigma_g(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distribution_a_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distribution_b_v_xi: ArrayLike, sigma_g_function: Callable[[ArrayLike], ArrayLike], num_gyroangle_points: int = 32) -> float:
    """Integrate f_a f_b σ(g) g over both gyrotropic velocity spaces"""
    result = gyrotropic_pair_integral_sigma_g_profiles(speed_grid_a, pitch_grid_a, np.asarray(distribution_a_v_xi, dtype=float)[None, :, :], speed_grid_b, pitch_grid_b, np.asarray(distribution_b_v_xi, dtype=float)[None, :, :], sigma_g_function, num_gyroangle_points=num_gyroangle_points)
   
    return float(result[0])

def _require_distribution_stack(values: ArrayLike, speed_grid: SpeedGrid, pitch_grid: PitchGrid, name: str) -> np.ndarray:
    """Validate f(z, v, ξ) with shape (n_z, n_speed, n_pitch)"""
    array = np.asarray(values, dtype=float)
    expected = (speed_grid.centers_m_s.size, pitch_grid.centers.size)
    if array.ndim == 2:
        array = array[None, :, :]
    if array.ndim != 3 or array.shape[1:] != expected:
        raise ValueError(f"{name} must have shape (n_z, {expected[0]}, {expected[1]})")
    if np.any(~np.isfinite(array)) or np.any(array < 0.0):
        raise ValueError(f"{name} must be finite and nonnegative")
    
    return array

def _kernel_array(kernel_function: Callable[[ArrayLike], ArrayLike], relative_speed: np.ndarray, name: str) -> np.ndarray:
    """Evaluate one nonnegative pair kernel with scalar fallback and broadcast validation"""
    try:
        values = np.asarray(kernel_function(relative_speed), dtype=float)
        if values.ndim == 0:
            values = np.full(relative_speed.shape, float(values), dtype=float)
        else:
            values = np.broadcast_to(values, relative_speed.shape).astype(float, copy=False)
    except (TypeError, ValueError):
        values = np.asarray([kernel_function(float(value)) for value in relative_speed], dtype=float)
    if (values.shape != relative_speed.shape or np.any(~np.isfinite(values)) or np.any(values < 0.0)):
        raise ValueError(f"{name} must return finite nonnegative values")
    
    return values

def _state_moment_weights(values: tuple[ArrayLike, ...], speed_grid: SpeedGrid, pitch_grid: PitchGrid, name: str) -> tuple[np.ndarray, ...]:
    """Validate nonnegative state moment weights on the (v, ξ) grid"""
    expected = (speed_grid.centers_m_s.size, pitch_grid.centers.size)
    result = []
    for index, item in enumerate(values):
        array = np.asarray(item, dtype=float)
        if array.shape == (expected[0],):
            array = np.repeat(array[:, None], expected[1], axis=1)
        if array.shape != expected or np.any(~np.isfinite(array)) or np.any(array < 0.0):
            raise ValueError(f"{name}[{index}] must be finite nonnegative with shape {expected} or ({expected[0]},)")
        result.append(array.reshape(-1))
    
    return tuple(result)

def gyrotropic_pair_integral_kernel_and_state_moments_profiles(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distributions_a_z_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distributions_b_z_v_xi: ArrayLike, kernel_functions: tuple[Callable[[ArrayLike], ArrayLike], ...], num_gyroangle_points: int = 32, *, max_pair_states_per_batch: int = _MAX_PAIR_STATES_PER_BATCH, state_moment_weights_a: tuple[ArrayLike, ...] = (), state_moment_weights_b: tuple[ArrayLike, ...] = ()) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate shared pair kernels and reactant state weighted moments

    Kernel output has shape (n_kernel, n_z)
    Reactant moment output has shape (n_moment, n_kernel, n_z)
    The relative gyrophase is integrated by midpoint averaging
    """
    if not kernel_functions:
        raise ValueError("at least one pair kernel function is required")
    stack_a = _require_distribution_stack(distributions_a_z_v_xi, speed_grid_a, pitch_grid_a, "distributions_a_z_v_xi")
    stack_b = _require_distribution_stack(distributions_b_z_v_xi, speed_grid_b, pitch_grid_b, "distributions_b_z_v_xi")
    if stack_a.shape[0] != stack_b.shape[0]:
        raise ValueError("reactant distribution stacks must have the same axial size")
    n_z = stack_a.shape[0]
    nphi = int(num_gyroangle_points)
    if nphi < 1:
        raise ValueError("num_gyroangle_points must be positive")
    max_states = int(max_pair_states_per_batch)
    if max_states < nphi:
        raise ValueError("max_pair_states_per_batch must be at least num_gyroangle_points")
    moment_weights_a = _state_moment_weights(state_moment_weights_a, speed_grid_a, pitch_grid_a, "state_moment_weights_a")
    moment_weights_b = _state_moment_weights(state_moment_weights_b, speed_grid_b, pitch_grid_b, "state_moment_weights_b")
    # Velocity cell measures already contain the gyrotropic 2π azimuthal factor
    weights_a = gyrotropic_velocity_cell_volumes(speed_grid_a, pitch_grid_a).reshape(-1)
    weights_b = gyrotropic_velocity_cell_volumes(speed_grid_b, pitch_grid_b).reshape(-1)
    weighted_a = stack_a.reshape(n_z, -1) * weights_a[None, :]
    weighted_b = stack_b.reshape(n_z, -1) * weights_b[None, :]
    active_a = np.flatnonzero(np.any(weighted_a != 0.0, axis=0))
    active_b = np.flatnonzero(np.any(weighted_b != 0.0, axis=0))
    results = np.zeros((len(kernel_functions), n_z), dtype=float)
    moments_a = np.zeros((len(moment_weights_a), len(kernel_functions), n_z), dtype=float)
    moments_b = np.zeros((len(moment_weights_b), len(kernel_functions), n_z), dtype=float)
    if active_a.size == 0 or active_b.size == 0:
        return results, moments_a, moments_b
    speed_a = np.repeat(np.asarray(speed_grid_a.centers_m_s, dtype=float), pitch_grid_a.centers.size)
    pitch_a = np.tile(np.asarray(pitch_grid_a.centers, dtype=float), speed_grid_a.centers_m_s.size)
    speed_b = np.repeat(np.asarray(speed_grid_b.centers_m_s, dtype=float), pitch_grid_b.centers.size)
    pitch_b = np.tile(np.asarray(pitch_grid_b.centers, dtype=float), speed_grid_b.centers_m_s.size)
    gyroangles = (np.arange(nphi, dtype=float) + 0.5) * (2.0 * np.pi / nphi)
    # Chunk state pairs while retaining every relative gyrophase point for each pair
    b_chunk_size = min(active_b.size, max(1, int(np.sqrt(max_states / nphi))),)
    a_chunk_size = min(active_a.size, max(1, max_states // (b_chunk_size * nphi)))
    for a_start in range(0, active_a.size, a_chunk_size):
        a_indices = active_a[a_start : a_start + a_chunk_size]
        for b_start in range(0, active_b.size, b_chunk_size):
            b_indices = active_b[b_start : b_start + b_chunk_size]
            states_per_a = b_indices.size * nphi
            a_pair_indices = np.repeat(a_indices, states_per_a)
            b_pair_indices = np.tile(np.repeat(b_indices, nphi), a_indices.size)
            phi = np.tile(np.tile(gyroangles, b_indices.size), a_indices.size)
            relative_speed = relative_speed_m_s(speed_a[a_pair_indices], pitch_a[a_pair_indices], speed_b[b_pair_indices], pitch_b[b_pair_indices], phi)
            kernel_values = np.stack([_kernel_array(function, relative_speed, f"kernel_functions[{index}]") for index, function in enumerate(kernel_functions)], axis=0,)
            gyroaveraged = np.sum((kernel_values / nphi).reshape(len(kernel_functions), a_indices.size, b_indices.size, nphi), axis=3,)
            weighted_a_chunk = weighted_a[:, a_indices]
            weighted_b_chunk = weighted_b[:, b_indices]
            results += np.einsum("za,kab,zb->kz", weighted_a_chunk, gyroaveraged, weighted_b_chunk, optimize=True)
            for index, state_weight in enumerate(moment_weights_a):
                moments_a[index] += np.einsum("za,kab,zb->kz", weighted_a_chunk * state_weight[a_indices][None, :], gyroaveraged, weighted_b_chunk, optimize=True)
            for index, state_weight in enumerate(moment_weights_b):
                moments_b[index] += np.einsum("za,kab,zb->kz", weighted_a_chunk, gyroaveraged, weighted_b_chunk * state_weight[b_indices][None, :], optimize=True)

    return results, moments_a, moments_b

def gyrotropic_pair_integral_kernel_profiles(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distributions_a_z_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distributions_b_z_v_xi: ArrayLike, kernel_functions: tuple[Callable[[ArrayLike], ArrayLike], ...], num_gyroangle_points: int = 32, *, max_pair_states_per_batch: int = _MAX_PAIR_STATES_PER_BATCH) -> np.ndarray:
    """Evaluate one or more shared gyrotropic pair kernels for each axial plane"""
    results, _, _ = gyrotropic_pair_integral_kernel_and_state_moments_profiles(speed_grid_a, pitch_grid_a, distributions_a_z_v_xi, speed_grid_b, pitch_grid_b, distributions_b_z_v_xi, kernel_functions, num_gyroangle_points=num_gyroangle_points, max_pair_states_per_batch=max_pair_states_per_batch)
   
    return results

def gyrotropic_pair_integral_sigma_g_profiles(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distributions_a_z_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distributions_b_z_v_xi: ArrayLike, sigma_g_function: Callable[[ArrayLike], ArrayLike], num_gyroangle_points: int = 32, *, max_pair_states_per_batch: int = _MAX_PAIR_STATES_PER_BATCH) -> np.ndarray:
    """Evaluate one shared σ(g) g kernel for multiple axial planes"""
  
    return gyrotropic_pair_integral_kernel_profiles(speed_grid_a, pitch_grid_a, distributions_a_z_v_xi, speed_grid_b, pitch_grid_b, distributions_b_z_v_xi, (sigma_g_function,), num_gyroangle_points=num_gyroangle_points, max_pair_states_per_batch=max_pair_states_per_batch)[0]

def gyrotropic_reaction_rate_density(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distribution_a_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distribution_b_v_xi: ArrayLike, reduced_mass_kg: float, cross_section_function_m2: Callable[[ArrayLike], ArrayLike], reaction: FusionReaction, *, identical_population: bool = False, num_gyroangle_points: int = 32) -> float:
    """Evaluate one reaction rate density from two local gyrotropic distributions"""
    def sigma_g(relative_speed):
        """Evaluate σ(E_cm) g for the supplied relative speed"""
        return sigma_g_from_cross_section(relative_speed, reduced_mass_kg, cross_section_function_m2)
    raw_integral = gyrotropic_pair_integral_sigma_g(speed_grid_a=speed_grid_a, pitch_grid_a=pitch_grid_a, distribution_a_v_xi=distribution_a_v_xi, speed_grid_b=speed_grid_b, pitch_grid_b=pitch_grid_b, distribution_b_v_xi=distribution_b_v_xi, sigma_g_function=sigma_g, num_gyroangle_points=num_gyroangle_points)
   
    return population_pair_counting_factor(identical_population) * raw_integral

def gyrotropic_reaction_rate_density_profile(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distributions_a_z_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distributions_b_z_v_xi: ArrayLike, reduced_mass_kg: float, cross_section_function_m2: Callable[[ArrayLike], ArrayLike], reaction: FusionReaction, *, identical_population: bool = False, num_gyroangle_points: int = 32, max_pair_states_per_batch: int = _MAX_PAIR_STATES_PER_BATCH) -> np.ndarray:
    """Evaluate an axial reaction rate density profile with one pair traversal"""
    def sigma_g(relative_speed):
        """Evaluate σ(E_cm) g for the supplied relative speed"""
        return sigma_g_from_cross_section(relative_speed, reduced_mass_kg, cross_section_function_m2)
    raw = gyrotropic_pair_integral_sigma_g_profiles(speed_grid_a, pitch_grid_a, distributions_a_z_v_xi, speed_grid_b, pitch_grid_b, distributions_b_z_v_xi, sigma_g, num_gyroangle_points=num_gyroangle_points, max_pair_states_per_batch=max_pair_states_per_batch)
  
    return population_pair_counting_factor(identical_population) * raw

def gyrotropic_reactivity_average(speed_grid_a: SpeedGrid, pitch_grid_a: PitchGrid, distribution_a_v_xi: ArrayLike, speed_grid_b: SpeedGrid, pitch_grid_b: PitchGrid, distribution_b_v_xi: ArrayLike, reduced_mass_kg: float, cross_section_function_m2: Callable[[ArrayLike], ArrayLike], num_gyroangle_points: int = 32) -> float:
    """Evaluate <σg> by dividing the gyrotropic pair integral by n_a n_b"""
    f_a = require_gyrotropic_distribution_shape(distribution_a_v_xi, speed_grid_a, pitch_grid_a, name="distribution_a_v_xi", allow_negative=False)
    f_b = require_gyrotropic_distribution_shape(distribution_b_v_xi, speed_grid_b, pitch_grid_b, name="distribution_b_v_xi", allow_negative=False)
    volume_a = gyrotropic_velocity_cell_volumes(speed_grid_a, pitch_grid_a)
    volume_b = gyrotropic_velocity_cell_volumes(speed_grid_b, pitch_grid_b)
    density_a = float(np.sum(f_a * volume_a))
    density_b = float(np.sum(f_b * volume_b))
    if density_a <= 0.0 or density_b <= 0.0:
        return 0.0
    def sigma_g(relative_speed):
        """Evaluate σ(E_cm) g for the supplied relative speed"""
        return sigma_g_from_cross_section(relative_speed, reduced_mass_kg, cross_section_function_m2)
    integral = gyrotropic_pair_integral_sigma_g(speed_grid_a, pitch_grid_a, f_a, speed_grid_b, pitch_grid_b, f_b, sigma_g, num_gyroangle_points=num_gyroangle_points)
    
    return integral / (density_a * density_b)