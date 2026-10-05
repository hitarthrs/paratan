"""Typed inputs and outputs for correlated neutron event generation"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import numpy as np
from source_model_revamp.fbis.velocity_space_grid import PitchGrid, SpeedGrid
from source_model_revamp.fusion.reactions import FusionReaction
from source_model_revamp.neutrons.kinematics import TwoBodyReactionKinematics

@dataclass(frozen=True)
class CorrelatedNeutronEventComponentSpec:
    """
    One neutron producing reactant pair on the shared axial grid
    
    Local reactant distributions use shape `(n_z, n_speed, n_pitch)` and `physical_rate_density_m3_s` supplies the fixed axial neutron rate density from the fusion calculation
    """
    label: str
    kind: str
    reaction: FusionReaction
    reaction_kinematics: TwoBodyReactionKinematics
    speed_grid_a: SpeedGrid
    pitch_grid_a: PitchGrid
    distributions_a_z_v_xi: np.ndarray
    mass_a_kg: float
    speed_grid_b: SpeedGrid
    pitch_grid_b: PitchGrid
    distributions_b_z_v_xi: np.ndarray
    mass_b_kg: float
    physical_rate_density_m3_s: np.ndarray
    identical_population: bool
    projectile_assignment_model: str
    num_gyroangle_points: int

    def __post_init__(self) -> None:
        """Validate component grids, distributions, rate profile, masses, reaction identity, and projectile assignment model"""
        distribution_a = np.asarray(self.distributions_a_z_v_xi, dtype=float)
        distribution_b = np.asarray(self.distributions_b_z_v_xi, dtype=float)
        rate = np.asarray(self.physical_rate_density_m3_s, dtype=float)
        shape_a = (rate.size, self.speed_grid_a.centers_m_s.size, self.pitch_grid_a.centers.size)
        shape_b = (rate.size, self.speed_grid_b.centers_m_s.size, self.pitch_grid_b.centers.size)
        if distribution_a.shape != shape_a:
            raise ValueError("distributions_a_z_v_xi does not match its axial and velocity grids")
        if distribution_b.shape != shape_b:
            raise ValueError("distributions_b_z_v_xi does not match its axial and velocity grids")
        if rate.ndim != 1 or rate.size == 0:
            raise ValueError("physical_rate_density_m3_s must be a nonempty 1D array")
        for name, values in (("distributions_a_z_v_xi", distribution_a), ("distributions_b_z_v_xi", distribution_b), ("physical_rate_density_m3_s", rate)):
            if np.any(~np.isfinite(values)) or np.any(values < 0.0):
                raise ValueError(f"{name} must be finite and nonnegative")
        for name, value in (("mass_a_kg", self.mass_a_kg), ("mass_b_kg", self.mass_b_kg)):
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if self.reaction_kinematics.reaction.key != self.reaction.key:
            raise ValueError("reaction_kinematics does not match the fusion reaction")
        if not np.isclose(self.reaction_kinematics.reactant_a_mass_kg, self.mass_a_kg, rtol=0.0, atol=0.0) or not np.isclose(self.reaction_kinematics.reactant_b_mass_kg, self.mass_b_kg, rtol=0.0, atol=0.0):
            raise ValueError("reactant masses do not match the reaction kinematics")
        gyroangle_count = int(self.num_gyroangle_points)
        if gyroangle_count < 1:
            raise ValueError("num_gyroangle_points must be positive")
        allowed_assignment = { "reactant_a_deuteron_projectile", "random_exchange_symmetrized"}
        if self.projectile_assignment_model not in allowed_assignment:
            raise ValueError("projectile_assignment_model must identify a supported deuteron ordering")
        object.__setattr__(self, "distributions_a_z_v_xi", distribution_a)
        object.__setattr__(self, "distributions_b_z_v_xi", distribution_b)
        object.__setattr__(self, "physical_rate_density_m3_s", rate)
        object.__setattr__(self, "num_gyroangle_points", gyroangle_count)

@dataclass(frozen=True)
class CorrelatedNeutronEventBank:
    """
    Correlated lab frame neutron source particles and audit fields
    
    All event arrays share length `n_event`
    `normalized_weights` sum to one while `physical_total_rate_s` stores the absolute neutron rate separately
    """
    positions_m: np.ndarray
    directions: np.ndarray
    energies_J: np.ndarray
    normalized_weights: np.ndarray
    physical_total_rate_s: float
    reaction_keys: np.ndarray
    component_labels: np.ndarray
    axial_cell_indices: np.ndarray
    center_of_mass_energy_J: np.ndarray
    equivalent_deuteron_lab_energy_eV: np.ndarray
    cm_emission_mu: np.ndarray
    reactant_a_velocity_m_s: np.ndarray
    reactant_b_velocity_m_s: np.ndarray
    metadata: dict[str, Any]

    def __post_init__(self) -> None:
        """Validate event array shapes, finite values, unit directions, normalized weights, and physical rate"""
        positions = np.asarray(self.positions_m, dtype=float)
        directions = np.asarray(self.directions, dtype=float)
        energies = np.asarray(self.energies_J, dtype=float)
        weights = np.asarray(self.normalized_weights, dtype=float)
        axial = np.asarray(self.axial_cell_indices, dtype=np.int64)
        center_energy = np.asarray(self.center_of_mass_energy_J, dtype=float)
        incident_energy = np.asarray(self.equivalent_deuteron_lab_energy_eV, dtype=float)
        mu = np.asarray(self.cm_emission_mu, dtype=float)
        velocity_a = np.asarray(self.reactant_a_velocity_m_s, dtype=float)
        velocity_b = np.asarray(self.reactant_b_velocity_m_s, dtype=float)
        reaction = np.asarray(self.reaction_keys)
        component = np.asarray(self.component_labels)
        if positions.ndim != 2 or positions.shape[1] != 3:
            raise ValueError("positions_m must have shape (n, 3)")
        count = positions.shape[0]
        if count < 1:
            raise ValueError("a correlated event bank must contain at least one event")
        if directions.shape != (count, 3):
            raise ValueError("directions must have shape (n, 3)")
        if velocity_a.shape != (count, 3) or velocity_b.shape != (count, 3):
            raise ValueError("reactant velocity arrays must have shape (n, 3)")
        for name, values in (("positions_m", positions), ("directions", directions), ("energies_J", energies), ("normalized_weights", weights), ("center_of_mass_energy_J", center_energy), ("equivalent_deuteron_lab_energy_eV", incident_energy), ("cm_emission_mu", mu), ("reactant_a_velocity_m_s", velocity_a), ("reactant_b_velocity_m_s", velocity_b)):
            if np.any(~np.isfinite(values)):
                raise ValueError(f"{name} must contain only finite values")
        for name, values in (("energies_J", energies), ("normalized_weights", weights), ("center_of_mass_energy_J", center_energy), ("equivalent_deuteron_lab_energy_eV", incident_energy)):
            if np.any(values < 0.0):
                raise ValueError(f"{name} must be nonnegative")
        if np.any(mu < -1.0) or np.any(mu > 1.0):
            raise ValueError("cm_emission_mu must lie inside [-1, 1]")
        if axial.shape != (count,):
            raise ValueError("axial_cell_indices must have length n")
        if np.any(axial < 0):
            raise ValueError("axial_cell_indices must be nonnegative")
        if energies.shape != (count,) or weights.shape != (count,):
            raise ValueError("energies_J and normalized_weights must have length n")
        if center_energy.shape != (count,) or incident_energy.shape != (count,):
            raise ValueError("event energy audit arrays must have length n")
        if mu.shape != (count,) or reaction.shape != (count,) or component.shape != (count,):
            raise ValueError("event labels and CM cosine must have length n")
        if not np.isfinite(self.physical_total_rate_s) or self.physical_total_rate_s <= 0.0:
            raise ValueError("physical_total_rate_s must be positive and finite")
        direction_norm = np.linalg.norm(directions, axis=1)
        if not np.allclose(direction_norm, 1.0, rtol=0.0, atol=2.0e-12):
            raise ValueError("directions must be unit vectors")
        weight_sum = float(np.sum(weights))
        if not np.isclose(weight_sum, 1.0, rtol=0.0, atol=2.0e-12):
            raise ValueError("normalized_weights must sum to one")
        object.__setattr__(self, "positions_m", positions)
        object.__setattr__(self, "directions", directions)
        object.__setattr__(self, "energies_J", energies)
        object.__setattr__(self, "normalized_weights", weights)
        object.__setattr__(self, "axial_cell_indices", axial)
        object.__setattr__(self, "center_of_mass_energy_J", center_energy)
        object.__setattr__(self, "equivalent_deuteron_lab_energy_eV", incident_energy)
        object.__setattr__(self, "cm_emission_mu", mu)
        object.__setattr__(self, "reactant_a_velocity_m_s", velocity_a)
        object.__setattr__(self, "reactant_b_velocity_m_s", velocity_b)
        object.__setattr__(self, "reaction_keys", reaction.astype("U8"))
        object.__setattr__(self, "component_labels", component.astype("U96"))

    @property
    def event_count(self) -> int:
        """Return the number of stored neutron events"""
        return int(self.energies_J.size)

    @property
    def physical_rate_weights_s(self) -> np.ndarray:
        """Return per event physical neutron rate weights in s⁻¹"""
        return self.normalized_weights * float(self.physical_total_rate_s)