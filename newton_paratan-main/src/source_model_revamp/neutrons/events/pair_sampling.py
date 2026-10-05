"""
Reactant pair sampling on represented gyrotropic velocity cells

Velocity cells are sampled from their represented particle density measure and the fusion event importance factor is `sigma(E_cm) * g`
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid, gyrotropic_velocity_cell_volumes
from source_model_revamp.fusion.cross_sections import bosch_hale_cross_section_m2_from_J, reduced_mass_kg
from source_model_revamp.neutrons.events.types import CorrelatedNeutronEventComponentSpec

@dataclass(frozen=True)
class SampledReactantPairs:
    """
    Sampled reactant pair state and reaction importance diagnostics
    
    Velocity arrays have shape `(n_event, 3)` and `sigma_g_m3_s` contains `sigma(E_cm) * g` for each sampled pair
    """
    velocity_a_m_s: np.ndarray
    velocity_b_m_s: np.ndarray
    center_of_mass_energy_J: np.ndarray
    sigma_g_m3_s: np.ndarray
    raw_rate_density_estimate_m3_s: float
    effective_sample_size: float
    projectile_exchange_count: int

def _sample_cell_indices(distribution_v_xi: np.ndarray, speed_grid: SpeedGrid, pitch_grid: PitchGrid, event_count: int, rng: np.random.Generator) -> tuple[np.ndarray, float]:
    """
    Sample flattened speed and pitch cells with probability proportional to `f d^3v`
    
    The returned scalar density is the represented integral of the input distribution
    """
    distribution = np.asarray(distribution_v_xi, dtype=float)
    expected = (speed_grid.centers_m_s.size, pitch_grid.centers.size)
    if distribution.shape != expected:
        raise ValueError(f"distribution_v_xi must have shape {expected}, got {distribution.shape}")
    if np.any(~np.isfinite(distribution)) or np.any(distribution < 0.0):
        raise ValueError("distribution_v_xi must be finite and nonnegative")
    cell_population = distribution * gyrotropic_velocity_cell_volumes(speed_grid, pitch_grid)
    density = float(np.sum(cell_population))
    if density <= 0.0:
        raise ValueError("reactant distribution has zero represented density")
    # Sampling from normalized f d^3v leaves sigma * g as the conditional reaction importance
    probability = cell_population.reshape(-1) / density
    index = rng.choice(probability.size, size=int(event_count), p=probability)

    return np.asarray(index, dtype=np.int64), density

def _sample_velocity_vectors(speed_grid: SpeedGrid, pitch_grid: PitchGrid, flat_cell_index: np.ndarray, gyro_angle_rad: np.ndarray) -> np.ndarray:
    """Build Cartesian velocity vectors at sampled speed and pitch cell centers for supplied gyro angles"""
    pitch_count = pitch_grid.centers.size
    speed_index = flat_cell_index // pitch_count
    pitch_index = flat_cell_index % pitch_count
    speed = speed_grid.centers_m_s[speed_index]
    pitch = pitch_grid.centers[pitch_index]
    perpendicular = speed * np.sqrt(np.maximum(1.0 - pitch**2, 0.0))

    return np.column_stack((perpendicular * np.cos(gyro_angle_rad), perpendicular * np.sin(gyro_angle_rad), speed * pitch))

def sample_reactant_pairs(component: CorrelatedNeutronEventComponentSpec, axial_index: int, event_count: int, rng: np.random.Generator) -> SampledReactantPairs:
    """
    Sample reactant pairs for one component and axial stratum
    
    Absolute gyro angle is uniform and relative gyro angle is sampled from the midpoint nodes used by the corresponding fusion pair kernel
    Identical D D populations can randomly exchange projectile and target ordering to symmetrize the evaluated angular axis
    """
    axial = int(axial_index)
    count = int(event_count)
    if count < 1:
        raise ValueError("event_count must be positive")
    if axial < 0 or axial >= component.physical_rate_density_m3_s.size:
        raise IndexError("axial_index is outside the component profile")
    index_a, density_a = _sample_cell_indices(component.distributions_a_z_v_xi[axial], component.speed_grid_a, component.pitch_grid_a, count, rng)
    index_b, density_b = _sample_cell_indices(component.distributions_b_z_v_xi[axial], component.speed_grid_b, component.pitch_grid_b, count, rng)
    base_azimuth = rng.uniform(0.0, 2.0 * np.pi, count)
    gyro_index = rng.integers(0, component.num_gyroangle_points, size=count,)
    # Reuse the midpoint relative gyro angle nodes represented by the deterministic fusion pair kernel
    relative_azimuth = (gyro_index.astype(float) + 0.5) * (2.0 * np.pi / component.num_gyroangle_points)
    velocity_a = _sample_velocity_vectors(component.speed_grid_a, component.pitch_grid_a, index_a, base_azimuth)
    velocity_b = _sample_velocity_vectors(component.speed_grid_b, component.pitch_grid_b, index_b, base_azimuth + relative_azimuth)
    projectile_exchange_count = 0
    if component.projectile_assignment_model == "random_exchange_symmetrized":
        exchange = rng.random(count) < 0.5
        projectile_exchange_count = int(np.count_nonzero(exchange))
        if np.any(exchange):
            temporary = velocity_a[exchange].copy()
            velocity_a[exchange] = velocity_b[exchange]
            velocity_b[exchange] = temporary
    relative_velocity = velocity_a - velocity_b
    relative_speed = np.linalg.norm(relative_velocity, axis=1)
    reduced_mass = reduced_mass_kg(component.mass_a_kg, component.mass_b_kg)
    center_of_mass_energy = 0.5 * reduced_mass * relative_speed**2
    sigma_g = (bosch_hale_cross_section_m2_from_J(component.reaction.key, center_of_mass_energy) * relative_speed)
    if np.any(~np.isfinite(sigma_g)) or np.any(sigma_g < 0.0):
        raise ValueError("sampled Bosch Hale sigma g values are invalid")
    pair_factor = 0.5 if component.identical_population else 1.0
    raw_rate = pair_factor * density_a * density_b * float(np.mean(sigma_g))
    importance_sum = float(np.sum(sigma_g))
    importance_square_sum = float(np.sum(sigma_g**2))
    effective_sample_size = (importance_sum**2 / importance_square_sum if importance_square_sum > 0.0 else 0.0)

    return SampledReactantPairs(
        velocity_a_m_s=velocity_a,
        velocity_b_m_s=velocity_b,
        center_of_mass_energy_J=center_of_mass_energy,
        sigma_g_m3_s=sigma_g,
        raw_rate_density_estimate_m3_s=raw_rate,
        effective_sample_size=effective_sample_size,
        projectile_exchange_count=projectile_exchange_count,
    )