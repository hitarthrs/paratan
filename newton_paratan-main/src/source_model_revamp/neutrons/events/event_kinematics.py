"""
Relativistic event kinematics adapters for correlated neutron sampling

The helpers connect sampled reactant pairs to the equivalent evaluated deuteron incident energy, center of mass emission direction, lab neutron state, and four momentum closure diagnostics
"""
from __future__ import annotations
import numpy as np
from source_model_revamp.constants import SPEED_OF_LIGHT_M_S
from source_model_revamp.nuclear_data.endf_b_viii1.incident_energy import equivalent_projectile_lab_kinetic_energy_eV_from_invariant_s_array
from source_model_revamp.neutrons.events.kinematics import directions_about_axes, pair_invariant_and_projectile_cm_direction
from source_model_revamp.neutrons.kinematics import TwoBodyReactionKinematics, neutron_lab_events_relativistic_pairwise_batch

def reactant_pair_invariant_s_J2(reaction_kinematics: TwoBodyReactionKinematics, velocity_a_m_s: np.ndarray, velocity_b_m_s: np.ndarray) -> np.ndarray:
    """Return the pair invariant `s = E_tot^2 − |p_tot c|^2` in J^2"""
    invariant, _ = pair_invariant_and_projectile_cm_direction(reaction_kinematics, velocity_a_m_s, velocity_b_m_s,)

    return invariant

def projectile_direction_in_cm(reaction_kinematics: TwoBodyReactionKinematics, velocity_a_m_s: np.ndarray, velocity_b_m_s: np.ndarray) -> np.ndarray:
    """Return the evaluated deuteron projectile direction in the reactant pair center of mass frame"""
    _, direction = pair_invariant_and_projectile_cm_direction(reaction_kinematics, velocity_a_m_s, velocity_b_m_s)
    return direction

def equivalent_deuteron_lab_energy_eV(reaction_kinematics: TwoBodyReactionKinematics, invariant_s_J2: np.ndarray) -> np.ndarray:
    """Map pair invariant `s` to the deuteron lab kinetic energy for the equivalent stationary target system"""
    return equivalent_projectile_lab_kinetic_energy_eV_from_invariant_s_array(invariant_s_J2, reaction_kinematics.reactant_a_mass_kg, reaction_kinematics.reactant_b_mass_kg)

def directions_from_axis_mu_phi(projectile_axis_cm: np.ndarray, mu_cm: np.ndarray, azimuth_rad: np.ndarray) -> np.ndarray:
    """Construct center of mass emission unit vectors from projectile axis, `mu`, and azimuth"""
    return directions_about_axes(projectile_axis_cm, mu_cm, azimuth_rad)

def correlated_neutron_lab_events(reaction_kinematics: TwoBodyReactionKinematics, velocity_a_m_s: np.ndarray, velocity_b_m_s: np.ndarray, cm_emission_direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Lorentz transform one sampled center of mass neutron direction per reactant pair into lab energy and direction"""
    return neutron_lab_events_relativistic_pairwise_batch(reaction_kinematics, velocity_a_m_s, velocity_b_m_s, cm_emission_direction)

def residual_mass_shell_relative_error(reaction_kinematics: TwoBodyReactionKinematics, velocity_a_m_s: np.ndarray, velocity_b_m_s: np.ndarray, neutron_kinetic_energy_J: np.ndarray, neutron_direction: np.ndarray) -> np.ndarray:
    """
    Return the residual product mass shell relative error after four momentum closure
    
    The diagnostic reconstructs the unobserved residual from incoming reactants and the sampled lab neutron state
    """
    kin = reaction_kinematics
    velocity_a = np.asarray(velocity_a_m_s, dtype=float)
    velocity_b = np.asarray(velocity_b_m_s, dtype=float)
    kinetic_energy = np.asarray(neutron_kinetic_energy_J, dtype=float)
    direction = np.asarray(neutron_direction, dtype=float)
    if velocity_a.ndim != 2 or velocity_a.shape[1] != 3:
        raise ValueError("velocity_a_m_s must have shape n by 3")
    if velocity_b.shape != velocity_a.shape or direction.shape != velocity_a.shape:
        raise ValueError("event kinematic arrays must have equal shape")
    if kinetic_energy.shape != (velocity_a.shape[0],):
        raise ValueError("neutron_kinetic_energy_J must have one value per event")
    c = SPEED_OF_LIGHT_M_S
    speed2_a = np.einsum("ij,ij->i", velocity_a, velocity_a)
    speed2_b = np.einsum("ij,ij->i", velocity_b, velocity_b)
    gamma_a = 1.0 / np.sqrt(1.0 - speed2_a / c**2)
    gamma_b = 1.0 / np.sqrt(1.0 - speed2_b / c**2)
    energy_a = gamma_a * kin.reactant_a_mass_kg * c**2
    energy_b = gamma_b * kin.reactant_b_mass_kg * c**2
    momentum_a = gamma_a[:, None] * kin.reactant_a_mass_kg * velocity_a
    momentum_b = gamma_b[:, None] * kin.reactant_b_mass_kg * velocity_b
    neutron_rest_energy = kin.neutron_mass_kg * c**2
    neutron_total_energy = kinetic_energy + neutron_rest_energy
    neutron_momentum_magnitude = np.sqrt(np.maximum(neutron_total_energy**2 - neutron_rest_energy**2, 0.0)) / c
    neutron_momentum = neutron_momentum_magnitude[:, None] * direction
    residual_energy = energy_a + energy_b - neutron_total_energy
    residual_momentum = momentum_a + momentum_b - neutron_momentum
    residual_invariant = residual_energy**2 - np.einsum("ij,ij->i", residual_momentum, residual_momentum,) * c**2
    target = (kin.residual_mass_kg * c**2) ** 2

    return np.abs(residual_invariant - target) / target