"""
Relativistic two body neutron kinematics shared by the deterministic spectrum and correlated event paths

Supported neutron branches are D + T → n + α and D + D → n + ³He

Reactant velocities are supplied in the lab frame
Emission directions are defined in the center of momentum frame and are Lorentz transformed with the two body final state into the lab frame
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.constants import (SPEED_OF_LIGHT_M_S, NEUTRON_MASS_KG, DEUTERON_MASS_KG, TRITON_MASS_KG, HELION_MASS_KG, ALPHA_PARTICLE_MASS_KG)
from source_model_revamp.fusion.reactions import DD_NEUTRON, DT_NEUTRON, FusionReaction, J_TO_MEV

@dataclass(frozen=True)
class TwoBodyReactionKinematics:
    """
    Mass metadata for a two body reaction a + b → n + residual
    
    Masses are in kg and the associated FusionReaction supplies the reaction identity used by the fusion and neutron modules
    """
    reaction: FusionReaction
    reactant_a_mass_kg: float
    reactant_b_mass_kg: float
    neutron_mass_kg: float
    residual_mass_kg: float
    residual_label: str

    @property
    def rest_mass_Q_J(self) -> float:
        """Return the rest mass Q value `(m_a + m_b − m_n − m_r) c²` in J"""
        return (self.reactant_a_mass_kg + self.reactant_b_mass_kg - self.neutron_mass_kg - self.residual_mass_kg) * SPEED_OF_LIGHT_M_S**2
    @property
    def rest_mass_Q_MeV(self) -> float:
        """Return the rest mass Q value in MeV"""
        return self.rest_mass_Q_J * J_TO_MEV

DT_NEUTRON_KINEMATICS = TwoBodyReactionKinematics(reaction=DT_NEUTRON, reactant_a_mass_kg=DEUTERON_MASS_KG, reactant_b_mass_kg=TRITON_MASS_KG, neutron_mass_kg=NEUTRON_MASS_KG, residual_mass_kg=ALPHA_PARTICLE_MASS_KG, residual_label="alpha")
DD_NEUTRON_KINEMATICS = TwoBodyReactionKinematics(reaction=DD_NEUTRON, reactant_a_mass_kg=DEUTERON_MASS_KG, reactant_b_mass_kg=DEUTERON_MASS_KG, neutron_mass_kg=NEUTRON_MASS_KG, residual_mass_kg=HELION_MASS_KG, residual_label="He3")
NEUTRON_KINEMATICS_BY_REACTION_KEY = {"dt_n": DT_NEUTRON_KINEMATICS, "dd_n": DD_NEUTRON_KINEMATICS}

def neutron_kinematics_for_reaction(reaction: str | FusionReaction) -> TwoBodyReactionKinematics:
    """Return the supported neutron kinematics record for a reaction key or FusionReaction"""
    key = reaction.key if isinstance(reaction, FusionReaction) else str(reaction)
    normalized = key.strip().lower()
    try:
        return NEUTRON_KINEMATICS_BY_REACTION_KEY[normalized]
    except KeyError as exc:
        raise ValueError(f"Unsupported neutron producing reaction {key!r}, supported keys are {sorted(NEUTRON_KINEMATICS_BY_REACTION_KEY)}") from exc

def as_vector3(values: ArrayLike, name: str = "vector") -> np.ndarray:
    """Return a finite Cartesian vector with shape `(3,)`"""
    vec = np.asarray(values, dtype=float)
    if vec.shape != (3,):
        vec = np.reshape(vec, (3,))
    if not np.all(np.isfinite(vec)):
        raise ValueError(f"{name} must contain only finite values")
    
    return vec

def unit_vector(vector: ArrayLike) -> np.ndarray:
    """Return the unit vector in the supplied Cartesian direction"""
    vec = as_vector3(vector)
    norm = float(np.linalg.norm(vec))
    if norm <= 0.0:
        raise ValueError("cannot normalize a zero vector")
    
    return vec / norm

def velocity_vector_from_speed_pitch_gyro(speed_m_s: float, pitch_cosine: float, gyro_angle_rad: float = 0.0) -> np.ndarray:
    """
    Build a local velocity vector for B along the z axis
    
    The pitch cosine is `ξ = v_parallel / v` and the gyro angle sets the perpendicular direction about the magnetic field
    """
    speed = float(speed_m_s)
    xi = float(pitch_cosine)
    angle = float(gyro_angle_rad)
    if not np.isfinite(speed) or speed < 0.0:
        raise ValueError("speed_m_s must be finite and nonnegative")
    if not np.isfinite(xi) or xi < -1.0 or xi > 1.0:
        raise ValueError("pitch_cosine must lie inside [-1, 1]")
    if not np.isfinite(angle):
        raise ValueError("gyro_angle_rad must be finite")
    v_perp = speed * np.sqrt(max(1.0 - xi**2, 0.0))

    return np.array([v_perp * np.cos(angle), v_perp * np.sin(angle), speed * xi], dtype=float)

def relative_velocity_vector_m_s(velocity_a_m_s: ArrayLike, velocity_b_m_s: ArrayLike) -> np.ndarray:
    """Return the relative velocity `v_a − v_b` in m s⁻¹"""
    return as_vector3(velocity_a_m_s, "velocity_a_m_s") - as_vector3(velocity_b_m_s, "velocity_b_m_s")

def relative_speed_m_s(velocity_a_m_s: ArrayLike, velocity_b_m_s: ArrayLike) -> float:
    """Return the relative speed `|v_a − v_b|` in m s⁻¹"""
    return float(np.linalg.norm(relative_velocity_vector_m_s(velocity_a_m_s, velocity_b_m_s)))

def reduced_mass_kg(mass_a_kg: float, mass_b_kg: float) -> float:
    """Return the reduced mass `μ = m_a m_b / (m_a + m_b)` in kg"""
    ma = float(mass_a_kg)
    mb = float(mass_b_kg)
    if ma <= 0.0 or mb <= 0.0 or not np.isfinite(ma) or not np.isfinite(mb):
        raise ValueError("masses must be positive and finite")
    
    return ma * mb / (ma + mb)

def center_of_mass_kinetic_energy_J(reactant_a_mass_kg: float, velocity_a_m_s: ArrayLike, reactant_b_mass_kg: float, velocity_b_m_s: ArrayLike) -> float:
    """
    Return the nonrelativistic relative kinetic energy `E_cm = 1 / 2 μ |v_a − v_b|²` in J
    
    This energy is used to evaluate the total fusion cross section
    """
    mu = reduced_mass_kg(reactant_a_mass_kg, reactant_b_mass_kg)
    g = relative_speed_m_s(velocity_a_m_s, velocity_b_m_s)

    return 0.5 * mu * g**2

def lorentz_gamma_from_velocity(velocity_m_s: ArrayLike) -> float:
    """Return `γ = 1 / sqrt(1 − |v|² / c²)` for one lab velocity"""
    speed2 = float(np.dot(as_vector3(velocity_m_s), as_vector3(velocity_m_s)))
    beta2 = speed2 / SPEED_OF_LIGHT_M_S**2
    if beta2 >= 1.0:
        raise ValueError("reactant speed must be smaller than the speed of light")
    
    return float(1.0 / np.sqrt(1.0 - beta2))

def total_relativistic_energy_J(mass_kg: float, velocity_m_s: ArrayLike) -> float:
    """Return the total relativistic energy `E = γ m c²` in J"""
    mass = float(mass_kg)
    if mass <= 0.0 or not np.isfinite(mass):
        raise ValueError("mass_kg must be positive and finite")
    
    return lorentz_gamma_from_velocity(velocity_m_s) * mass * SPEED_OF_LIGHT_M_S**2

def relativistic_momentum_kg_m_s(mass_kg: float, velocity_m_s: ArrayLike) -> np.ndarray:
    """Return the relativistic momentum `p = γ m v` in kg m s⁻¹"""
    mass = float(mass_kg)
    if mass <= 0.0 or not np.isfinite(mass):
        raise ValueError("mass_kg must be positive and finite")
    velocity = as_vector3(velocity_m_s)

    return lorentz_gamma_from_velocity(velocity) * mass * velocity

def invariant_s_J2(total_energy_J: float, total_momentum_kg_m_s: ArrayLike) -> float:
    """Return the invariant `s = E_tot² − |p_tot c|²` in J²"""
    momentum = as_vector3(total_momentum_kg_m_s, "total_momentum_kg_m_s")
    value = float(total_energy_J) ** 2 - float(np.dot(momentum, momentum)) * SPEED_OF_LIGHT_M_S**2
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("invariant s must be positive and finite")
    
    return value

def cm_product_total_energy_J(s_J2: float, product_mass_kg: float, residual_mass_kg: float) -> float:
    """Return the center of momentum total energy of one two body final product in J"""
    s = float(s_J2)
    if s <= 0.0 or not np.isfinite(s):
        raise ValueError("s_J2 must be positive and finite")
    sqrt_s = np.sqrt(s)
    m1c2 = float(product_mass_kg) * SPEED_OF_LIGHT_M_S**2
    m2c2 = float(residual_mass_kg) * SPEED_OF_LIGHT_M_S**2

    return float((s + m1c2**2 - m2c2**2) / (2.0 * sqrt_s))

def cm_product_momentum_magnitude_kg_m_s(s_J2: float, product_mass_kg: float, residual_mass_kg: float) -> float:
    """Return the common center of momentum momentum magnitude of the two final products"""
    s = float(s_J2)
    sqrt_s = np.sqrt(s)
    m1c2 = float(product_mass_kg) * SPEED_OF_LIGHT_M_S**2
    m2c2 = float(residual_mass_kg) * SPEED_OF_LIGHT_M_S**2
    factor = (s - (m1c2 + m2c2) ** 2) * (s - (m1c2 - m2c2) ** 2)
    if factor < -1.0e-20 * s**2:
        raise ValueError("reaction is below threshold for the supplied reactant velocities")
    
    return float(np.sqrt(max(factor, 0.0)) / (2.0 * sqrt_s * SPEED_OF_LIGHT_M_S))

def lorentz_boost_momentum_energy_to_lab(momentum_cm_kg_m_s: ArrayLike, energy_cm_J: float, beta_cm_vector: ArrayLike) -> tuple[np.ndarray, float]:
    """
    Boost one center of momentum four momentum into the lab frame
    
    `beta_cm_vector` is the dimensionless center of momentum velocity vector `v_cm / c`
    """
    p_star = as_vector3(momentum_cm_kg_m_s, "momentum_cm_kg_m_s")
    beta = as_vector3(beta_cm_vector, "beta_cm_vector")
    beta2 = float(np.dot(beta, beta))
    if beta2 >= 1.0:
        raise ValueError("CM beta must be smaller than one")
    if beta2 == 0.0:
        return p_star, float(energy_cm_J)
    gamma = 1.0 / np.sqrt(1.0 - beta2)
    beta_dot_p = float(np.dot(beta, p_star))
    p_lab = p_star + (((gamma - 1.0) * beta_dot_p / beta2) + gamma * float(energy_cm_J) / SPEED_OF_LIGHT_M_S) * beta
    energy_lab = gamma * (float(energy_cm_J) + beta_dot_p * SPEED_OF_LIGHT_M_S)

    return p_lab, float(energy_lab)

def neutron_lab_momentum_energy_relativistic(reaction_kinematics: TwoBodyReactionKinematics, velocity_a_m_s: ArrayLike, velocity_b_m_s: ArrayLike, cm_emission_direction: ArrayLike) -> tuple[np.ndarray, float]:
    """
    Return one neutron lab momentum and total energy from relativistic two body kinematics
    
    The supplied emission direction is interpreted in the center of momentum frame
    """
    kin = reaction_kinematics
    v_a = as_vector3(velocity_a_m_s, "velocity_a_m_s")
    v_b = as_vector3(velocity_b_m_s, "velocity_b_m_s")
    direction = unit_vector(cm_emission_direction)
    E_a = total_relativistic_energy_J(kin.reactant_a_mass_kg, v_a)
    E_b = total_relativistic_energy_J(kin.reactant_b_mass_kg, v_b)
    p_a = relativistic_momentum_kg_m_s(kin.reactant_a_mass_kg, v_a)
    p_b = relativistic_momentum_kg_m_s(kin.reactant_b_mass_kg, v_b)
    E_tot = E_a + E_b
    p_tot = p_a + p_b
    # Build the incoming four momentum invariant before solving the two body final state
    s = invariant_s_J2(E_tot, p_tot)
    E_star = cm_product_total_energy_J(s, kin.neutron_mass_kg, kin.residual_mass_kg)
    p_star_mag = cm_product_momentum_magnitude_kg_m_s(s, kin.neutron_mass_kg, kin.residual_mass_kg)
    beta_cm = p_tot * SPEED_OF_LIGHT_M_S**2 / E_tot

    return lorentz_boost_momentum_energy_to_lab(p_star_mag * direction, E_star, beta_cm / SPEED_OF_LIGHT_M_S)

def neutron_lab_kinetic_energy_relativistic_J(reaction_kinematics: TwoBodyReactionKinematics, velocity_a_m_s: ArrayLike, velocity_b_m_s: ArrayLike, cm_emission_direction: ArrayLike) -> float:
    """Return one neutron lab kinetic energy in J from relativistic two body kinematics"""
    _, total_energy = neutron_lab_momentum_energy_relativistic(reaction_kinematics, velocity_a_m_s, velocity_b_m_s, cm_emission_direction)
    kinetic = total_energy - reaction_kinematics.neutron_mass_kg * SPEED_OF_LIGHT_M_S**2
    return float(max(kinetic, 0.0))

def neutron_lab_direction(reaction_kinematics: TwoBodyReactionKinematics, velocity_a_m_s: ArrayLike, velocity_b_m_s: ArrayLike, cm_emission_direction: ArrayLike) -> np.ndarray:
    """Return one neutron lab direction unit vector from relativistic two body kinematics"""
    momentum, _ = neutron_lab_momentum_energy_relativistic(reaction_kinematics, velocity_a_m_s, velocity_b_m_s, cm_emission_direction)

    return unit_vector(momentum)

def neutron_lab_events_relativistic_batch(reaction_kinematics: TwoBodyReactionKinematics, velocity_a_m_s: ArrayLike, velocity_b_m_s: ArrayLike, cm_emission_directions: ArrayLike, *, include_lab_directions: bool = False) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Evaluate all requested center of momentum emission directions for a batch of reactant pairs
    
    Input reactant velocities have shape `(n_pair, 3)` and directions have shape `(n_direction, 3)`
    Returned kinetic energies have shape `(n_pair, n_direction)`
    Lab directions are returned with shape `(n_pair, n_direction, 3)` when requested
    """
    kin = reaction_kinematics
    v_a = np.asarray(velocity_a_m_s, dtype=float)
    v_b = np.asarray(velocity_b_m_s, dtype=float)
    directions = np.asarray(cm_emission_directions, dtype=float)
    if v_a.ndim == 1:
        v_a = np.reshape(v_a, (1, 3))
    if v_b.ndim == 1:
        v_b = np.reshape(v_b, (1, 3))
    if directions.ndim == 1:
        directions = np.reshape(directions, (1, 3))
    if v_a.ndim != 2 or v_a.shape[1] != 3:
        raise ValueError("velocity_a_m_s must have shape (n, 3)")
    if v_b.shape != v_a.shape:
        raise ValueError("velocity_b_m_s must match velocity_a_m_s shape")
    if directions.ndim != 2 or directions.shape[1] != 3 or directions.shape[0] < 1:
        raise ValueError("cm_emission_directions must have shape (m, 3)")
    if not np.all(np.isfinite(v_a)) or not np.all(np.isfinite(v_b)) or not np.all(np.isfinite(directions)):
        raise ValueError("velocities and emission directions must contain only finite values")
    direction_norms = np.linalg.norm(directions, axis=1)
    if np.any(direction_norms <= 0.0):
        raise ValueError("cm_emission_directions must not contain a zero vector")
    unit_directions = directions / direction_norms[:, None]

    c = SPEED_OF_LIGHT_M_S
    speed2_a = np.einsum("ij,ij->i", v_a, v_a)
    speed2_b = np.einsum("ij,ij->i", v_b, v_b)
    beta2_a = speed2_a / c**2
    beta2_b = speed2_b / c**2
    if np.any(beta2_a >= 1.0) or np.any(beta2_b >= 1.0):
        raise ValueError("reactant speed must be smaller than the speed of light")
    gamma_a = 1.0 / np.sqrt(1.0 - beta2_a)
    gamma_b = 1.0 / np.sqrt(1.0 - beta2_b)
    E_a = gamma_a * kin.reactant_a_mass_kg * c**2
    E_b = gamma_b * kin.reactant_b_mass_kg * c**2
    p_a = gamma_a[:, None] * kin.reactant_a_mass_kg * v_a
    p_b = gamma_b[:, None] * kin.reactant_b_mass_kg * v_b
    E_total = E_a + E_b
    p_total = p_a + p_b
    s = E_total**2 - np.einsum("ij,ij->i", p_total, p_total) * c**2
    if np.any(~np.isfinite(s)) or np.any(s <= 0.0):
        raise ValueError("invariant s must be positive and finite")

    sqrt_s = np.sqrt(s)
    neutron_rest_energy = kin.neutron_mass_kg * c**2
    residual_rest_energy = kin.residual_mass_kg * c**2
    E_star = (s + neutron_rest_energy**2 - residual_rest_energy**2) / (2.0 * sqrt_s)
    factor = (s - (neutron_rest_energy + residual_rest_energy) ** 2) * (s - (neutron_rest_energy - residual_rest_energy) ** 2)
    if np.any(factor < -1.0e-20 * s**2):
        raise ValueError("reaction is below threshold for the supplied reactant velocities")
    p_star_magnitude = np.sqrt(np.maximum(factor, 0.0)) / (2.0 * sqrt_s * c)

    beta_cm = p_total * c / E_total[:, None]
    beta2_cm = np.einsum("ij,ij->i", beta_cm, beta_cm)
    if np.any(beta2_cm >= 1.0):
        raise ValueError("CM beta must be smaller than one")
    gamma_cm = 1.0 / np.sqrt(1.0 - beta2_cm)
    beta_dot_direction = beta_cm @ unit_directions.T
    beta_dot_p = p_star_magnitude[:, None] * beta_dot_direction
    total_energy_lab = gamma_cm[:, None] * (E_star[:, None] + beta_dot_p * c)
    kinetic_energy_lab = np.maximum(total_energy_lab - neutron_rest_energy, 0.0)

    if not include_lab_directions:
        return kinetic_energy_lab, None

    p_star = p_star_magnitude[:, None, None] * unit_directions[None, :, :]
    beta2_safe = np.where(beta2_cm > 0.0, beta2_cm, 1.0)
    coefficient = ((gamma_cm[:, None] - 1.0) * beta_dot_p / beta2_safe[:, None] + gamma_cm[:, None] * E_star[:, None] / c)
    coefficient = np.where(beta2_cm[:, None] > 0.0, coefficient, 0.0)
    p_lab = p_star + coefficient[:, :, None] * beta_cm[:, None, :]
    p_norm = np.linalg.norm(p_lab, axis=2)
    if np.any(p_norm <= 0.0):
        raise ValueError("cannot normalize a zero neutron momentum")

    return kinetic_energy_lab, p_lab / p_norm[:, :, None]

def neutron_lab_kinetic_energy_endpoints_relativistic_batch(reaction_kinematics: TwoBodyReactionKinematics, velocity_a_m_s: ArrayLike, velocity_b_m_s: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    """
    Return exact lab kinetic energy endpoints for isotropic center of momentum emission
    
    For each reactant pair, the Lorentz boost gives
    `E_lab = γ_cm * (E_star + β_cm * p_star * c * cos(θ_star))`
    
    Uniform `cos(θ_star)` therefore maps to a uniform lab energy interval between the returned endpoints
    """
    kin = reaction_kinematics
    v_a = np.asarray(velocity_a_m_s, dtype=float)
    v_b = np.asarray(velocity_b_m_s, dtype=float)
    if v_a.ndim == 1:
        v_a = np.reshape(v_a, (1, 3))
    if v_b.ndim == 1:
        v_b = np.reshape(v_b, (1, 3))
    if v_a.ndim != 2 or v_a.shape[1] != 3:
        raise ValueError("velocity_a_m_s must have shape (n, 3)")
    if v_b.shape != v_a.shape:
        raise ValueError("velocity_b_m_s must match velocity_a_m_s shape")
    if not np.all(np.isfinite(v_a)) or not np.all(np.isfinite(v_b)):
        raise ValueError("velocities must contain only finite values")
    c = SPEED_OF_LIGHT_M_S
    speed2_a = np.einsum("ij,ij->i", v_a, v_a)
    speed2_b = np.einsum("ij,ij->i", v_b, v_b)
    beta2_a = speed2_a / c**2
    beta2_b = speed2_b / c**2
    if np.any(beta2_a >= 1.0) or np.any(beta2_b >= 1.0):
        raise ValueError("reactant speed must be smaller than the speed of light")
    gamma_a = 1.0 / np.sqrt(1.0 - beta2_a)
    gamma_b = 1.0 / np.sqrt(1.0 - beta2_b)
    E_a = gamma_a * kin.reactant_a_mass_kg * c**2
    E_b = gamma_b * kin.reactant_b_mass_kg * c**2
    p_a = gamma_a[:, None] * kin.reactant_a_mass_kg * v_a
    p_b = gamma_b[:, None] * kin.reactant_b_mass_kg * v_b
    E_total = E_a + E_b
    p_total = p_a + p_b
    s = E_total**2 - np.einsum("ij,ij->i", p_total, p_total) * c**2
    if np.any(~np.isfinite(s)) or np.any(s <= 0.0):
        raise ValueError("invariant s must be positive and finite")

    sqrt_s = np.sqrt(s)
    neutron_rest_energy = kin.neutron_mass_kg * c**2
    residual_rest_energy = kin.residual_mass_kg * c**2
    E_star = (s + neutron_rest_energy**2 - residual_rest_energy**2) / (2.0 * sqrt_s)
    factor = (s - (neutron_rest_energy + residual_rest_energy) ** 2) * (s - (neutron_rest_energy - residual_rest_energy) ** 2)
    if np.any(factor < -1.0e-20 * s**2):
        raise ValueError("reaction is below threshold for the supplied reactant velocities")
    p_star_magnitude = np.sqrt(np.maximum(factor, 0.0)) / (2.0 * sqrt_s * c)
    beta_cm = p_total * c / E_total[:, None]
    beta2_cm = np.einsum("ij,ij->i", beta_cm, beta_cm)
    if np.any(beta2_cm >= 1.0):
        raise ValueError("CM beta must be smaller than one")
    gamma_cm = 1.0 / np.sqrt(1.0 - beta2_cm)
    energy_half_width = gamma_cm * np.sqrt(beta2_cm) * p_star_magnitude * c
    energy_center = gamma_cm * E_star - neutron_rest_energy
    energy_min = np.maximum(energy_center - energy_half_width, 0.0)
    energy_max = np.maximum(energy_center + energy_half_width, 0.0)

    return energy_min, energy_max

def isotropic_emission_directions(num_directions: int) -> np.ndarray:
    """Return deterministic near uniform unit vectors on the sphere"""
    n = int(num_directions)
    if n < 1:
        raise ValueError("num_directions must be at least one")
    indices = np.arange(n, dtype=float)
    z = 1.0 - 2.0 * (indices + 0.5) / n
    phi = np.pi * (3.0 - np.sqrt(5.0)) * indices
    radius = np.sqrt(np.maximum(1.0 - z**2, 0.0))

    return np.column_stack((radius * np.cos(phi), radius * np.sin(phi), z))

def neutron_lab_events_relativistic_pairwise_batch(reaction_kinematics: TwoBodyReactionKinematics, velocity_a_m_s: ArrayLike, velocity_b_m_s: ArrayLike, cm_emission_directions: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    """
    Evaluate relativistic neutron kinematics with one center of momentum direction per reactant pair
    
    All three input arrays have shape `(n_pair, 3)`
    The outputs are lab kinetic energy with shape `(n_pair,)` and lab direction with shape `(n_pair, 3)`
    """
    kin = reaction_kinematics
    v_a = np.asarray(velocity_a_m_s, dtype=float)
    v_b = np.asarray(velocity_b_m_s, dtype=float)
    directions = np.asarray(cm_emission_directions, dtype=float)
    if v_a.ndim != 2 or v_a.shape[1] != 3:
        raise ValueError("velocity_a_m_s must have shape n by 3")
    if v_b.shape != v_a.shape or directions.shape != v_a.shape:
        raise ValueError("pairwise velocity and direction arrays must have equal shape")
    if np.any(~np.isfinite(v_a)) or np.any(~np.isfinite(v_b)) or np.any(~np.isfinite(directions)):
        raise ValueError("pairwise kinematic arrays must be finite")
    direction_norm = np.linalg.norm(directions, axis=1)
    if np.any(direction_norm <= 0.0):
        raise ValueError("CM emission directions must be nonzero")
    unit_directions = directions / direction_norm[:, None]
    c = SPEED_OF_LIGHT_M_S
    speed2_a = np.einsum("ij,ij->i", v_a, v_a)
    speed2_b = np.einsum("ij,ij->i", v_b, v_b)
    beta2_a = speed2_a / c**2
    beta2_b = speed2_b / c**2
    if np.any(beta2_a >= 1.0) or np.any(beta2_b >= 1.0):
        raise ValueError("reactant speed must be smaller than the speed of light")
    gamma_a = 1.0 / np.sqrt(1.0 - beta2_a)
    gamma_b = 1.0 / np.sqrt(1.0 - beta2_b)
    energy_a = gamma_a * kin.reactant_a_mass_kg * c**2
    energy_b = gamma_b * kin.reactant_b_mass_kg * c**2
    momentum_a = gamma_a[:, None] * kin.reactant_a_mass_kg * v_a
    momentum_b = gamma_b[:, None] * kin.reactant_b_mass_kg * v_b
    total_energy = energy_a + energy_b
    total_momentum = momentum_a + momentum_b
    invariant_s = total_energy**2 - np.einsum("ij,ij->i", total_momentum, total_momentum,) * c**2
    if np.any(~np.isfinite(invariant_s)) or np.any(invariant_s <= 0.0):
        raise ValueError("invariant s must be positive and finite")
    sqrt_s = np.sqrt(invariant_s)
    neutron_rest_energy = kin.neutron_mass_kg * c**2
    residual_rest_energy = kin.residual_mass_kg * c**2
    energy_star = (invariant_s + neutron_rest_energy**2 - residual_rest_energy**2) / (2.0 * sqrt_s)
    factor = (invariant_s - (neutron_rest_energy + residual_rest_energy) ** 2) * (invariant_s - (neutron_rest_energy - residual_rest_energy) ** 2)
    if np.any(factor < -1.0e-20 * invariant_s**2):
        raise ValueError("reaction is below threshold for the supplied reactant velocities")
    momentum_star_magnitude = np.sqrt(np.maximum(factor, 0.0)) / (2.0 * sqrt_s * c)
    beta_cm = total_momentum * c / total_energy[:, None]
    beta2_cm = np.einsum("ij,ij->i", beta_cm, beta_cm)
    if np.any(beta2_cm >= 1.0):
        raise ValueError("CM beta must be smaller than one")
    gamma_cm = 1.0 / np.sqrt(1.0 - beta2_cm)
    momentum_star = momentum_star_magnitude[:, None] * unit_directions
    beta_dot_p = np.einsum("ij,ij->i", beta_cm, momentum_star)
    total_energy_lab = gamma_cm * (energy_star + beta_dot_p * c)
    beta2_safe = np.where(beta2_cm > 0.0, beta2_cm, 1.0)
    coefficient = ((gamma_cm - 1.0) * beta_dot_p / beta2_safe + gamma_cm * energy_star / c)
    coefficient = np.where(beta2_cm > 0.0, coefficient, 0.0)
    momentum_lab = momentum_star + coefficient[:, None] * beta_cm
    momentum_norm = np.linalg.norm(momentum_lab, axis=1)
    if np.any(momentum_norm <= 0.0):
        raise ValueError("cannot normalize a zero neutron momentum")
    kinetic_energy = np.maximum(total_energy_lab - neutron_rest_energy, 0.0)

    return kinetic_energy, momentum_lab / momentum_norm[:, None]