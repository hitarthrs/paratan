"""
Conservative speed resolved target ion removal from primary beam charge exchange

Primary neutral charge exchange creates a fast beam ion but removes one ion from the target population
The attenuation calculation supplies the exact cell event rates
This module redistributes those events over target speed using the solved local kinetic distributions and the same charge exchange sigma times relative speed kernel

The final sink is reduced to a volume averaged speed dependent loss frequency on the modal target speed grid
"""
from __future__ import annotations
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING
import numpy as np
from source_model_revamp.beam.kinetic_target_rates import gyroaveraged_sigma_g, kinetic_species_distribution
from source_model_revamp.fbis.beam_source_definition import beam_speed_from_energy_m_s
from source_model_revamp.fbis.species import DEUTERON, TRITON, IonSpecies, ion_species
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid, gyrotropic_velocity_cell_volumes

if TYPE_CHECKING:
    from source_model_revamp.integration.pipeline_types import BeamDepositionResult, GeometryStageResult, KineticStageResult

@dataclass(frozen=True)
class ChargeExchangeTargetSinkState:
    """
    Speed resolved target ion loss state for one D or T species
    
    loss_frequency_s stores a nonnegative removal frequency in s⁻¹ on speed_grid
    reference_event_rate_s is the exact attenuation event rate that the sink must reproduce
    reference_energy_removal_W is the target kinetic energy removed per unit time
    The projectile resolved mapping must sum to the same reference event rate
    limitation records the pitch and axial information removed by the reduced sink representation
    """
    target_species: IonSpecies
    speed_grid: SpeedGrid
    loss_frequency_s: np.ndarray
    reference_event_rate_s: float
    reference_energy_removal_W: float
    reference_event_rate_by_projectile_species_s: Mapping[str, float]
    model: str
    reference_rate_identity_relative_error: float
    qualified: bool
    limitation: str | None

    def __post_init__(self) -> None:
        """
        Validate sink arrays rates power identity error and projectile resolved rate conservation
        
        The projectile rate mapping is copied into a read only mapping proxy
        """
        frequency = np.asarray(self.loss_frequency_s, dtype=float)
        expected = np.asarray(self.speed_grid.centers_m_s, dtype=float).shape
        if frequency.shape != expected or np.any(~np.isfinite(frequency)) or np.any(frequency < 0.0):
            raise ValueError("charge exchange loss frequency must be finite and nonnegative on the target speed grid")
        rate = float(self.reference_event_rate_s)
        power = float(self.reference_energy_removal_W)
        error = float(self.reference_rate_identity_relative_error)
        if not np.isfinite(rate) or rate < 0.0:
            raise ValueError("charge exchange reference event rate must be finite and nonnegative")
        if not np.isfinite(power) or power < 0.0:
            raise ValueError("charge exchange reference energy removal must be finite and nonnegative")
        if not np.isfinite(error) or error < 0.0:
            raise ValueError("charge exchange reference identity error must be finite and nonnegative")
        rates = {str(key): float(value) for key, value in dict(self.reference_event_rate_by_projectile_species_s).items()}
        if any(not np.isfinite(value) or value < 0.0 for value in rates.values()):
            raise ValueError("charge exchange projectile resolved rates must be finite and nonnegative")
        scale = max(rate, float(sum(rates.values())), 1.0)
        if abs(rate - float(sum(rates.values()))) > 1.0e-12 * scale:
            raise ValueError("charge exchange projectile resolved rates must sum to the reference event rate")
        if not self.model:
            raise ValueError("charge exchange sink model must be nonempty")
        object.__setattr__(self, "loss_frequency_s", frequency)
        object.__setattr__(self, "reference_event_rate_by_projectile_species_s", MappingProxyType(rates))

def _target_cell_event_rates(component: object, target_species_id: str) -> np.ndarray:
    """Return target species resolved charge exchange event rates for one attenuated beam component"""
    if target_species_id == DEUTERON.species_id:
        return np.asarray(component.cell_deuterium_charge_exchange_rates_s, dtype=float)
    if target_species_id == TRITON.species_id:
        return np.asarray(component.cell_tritium_charge_exchange_rates_s, dtype=float)
    raise ValueError("charge exchange target species must be deuterium or tritium")

def _modal_speed_population_particles(kinetic: KineticStageResult, target_species: IonSpecies, volume_m3: float) -> tuple[SpeedGrid, np.ndarray]:
    """
    Return modal target particle inventory resolved by speed cell
    
    The invariant distribution is integrated over eta then multiplied by speed shell volume and total plasma volume
    """
    system = kinetic.fast_ion_system_state
    if system is None:
        raise ValueError("charge exchange redistribution requires a solved fast ion system state")
    state = system.species_states.get(target_species.species_id)
    if state is None or state.modal_result is None:
        raise ValueError(f"charge exchange redistribution requires a solved {target_species.species_id} target state")
    modal = state.modal_result
    eta = np.asarray(modal.basis.physical_basis.eta_grid, dtype=float)
    distribution = np.maximum(np.asarray(modal.distribution_v_eta, dtype=float), 0.0)
    eta_integral = np.trapezoid(distribution, eta, axis=1)
    particles = float(volume_m3) * np.asarray(modal.speed_grid.shell_volumes_m3_s3, dtype=float) * eta_integral
    if np.any(~np.isfinite(particles)) or np.any(particles < 0.0):
        raise ValueError("charge exchange target speed population must be finite and nonnegative")
  
    return modal.speed_grid, particles

def _interpolated_loss_frequency(local_speed_grid: SpeedGrid, local_event_rate_by_speed_s: np.ndarray, local_inventory_by_speed: np.ndarray, modal_speed_grid: SpeedGrid, modal_population_particles: np.ndarray, reference_event_rate_s: float) -> tuple[np.ndarray, float]:
    """
    Map a local speed resolved event rate to the modal speed grid as a conservative loss frequency
    
    Local event rate divided by local particle inventory defines the initial frequency
    The frequency is interpolated onto the modal speed centers and renormalized so its modal inventory integral reproduces reference_event_rate_s
    """
    local_rate = np.asarray(local_event_rate_by_speed_s, dtype=float)
    local_inventory = np.asarray(local_inventory_by_speed, dtype=float)
    # Event rate divided by represented inventory defines the local speed dependent removal frequency
    local_frequency = np.divide(local_rate, local_inventory, out=np.zeros_like(local_rate), where=local_inventory > 0.0)
    frequency = np.interp(np.asarray(modal_speed_grid.centers_m_s, dtype=float), np.asarray(local_speed_grid.centers_m_s, dtype=float), local_frequency, left=0.0, right=0.0)
    reconstructed = float(np.sum(frequency * modal_population_particles))
    reference = float(reference_event_rate_s)
    if reference > 0.0:
        if reconstructed <= 0.0:
            raise ValueError("charge exchange target event rate has no represented kinetic support")
        # Renormalization preserves the exact attenuation event rate after interpolation to the modal grid
        frequency *= reference / reconstructed
    elif reconstructed > 0.0:
        frequency.fill(0.0)
    final_rate = float(np.sum(frequency * modal_population_particles))
    error = abs(final_rate - reference) / max(abs(reference), 1.0)
 
    return frequency, float(error)

def build_charge_exchange_target_sink_states(*, depositions: Sequence[BeamDepositionResult], kinetic_target_state: KineticStageResult, geometry: GeometryStageResult, gyroangle_points: int) -> dict[str, ChargeExchangeTargetSinkState]:
    """
    Build conservative target ion charge exchange sinks for available solved D and T states
    
    Each attenuation cell event rate is distributed across local target speed according to the charge exchange reaction weight
    Events from all beam components and projectile species are accumulated before interpolation to the modal target grid
    The returned frequency is qualified when its reconstructed total event rate matches the attenuation reference within the implemented 1e−12 identity threshold
    """
    nphi = int(gyroangle_points)
    if nphi < 1:
        raise ValueError("gyroangle_points must be positive")
    volumes = np.asarray(geometry.cell_volumes_m3, dtype=float)
    if np.any(~np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("charge exchange redistribution requires positive finite cell volumes")
    states: dict[str, ChargeExchangeTargetSinkState] = {}
    for target_species in (DEUTERON, TRITON):
        resolved = kinetic_species_distribution(kinetic_target_state, target_species.species_id, volumes.size)
        if resolved is None:
            continue
        local_speed_grid, local_distribution = resolved
        velocity_weights = gyrotropic_velocity_cell_volumes(local_speed_grid, kinetic_target_state.pitch_grid)
        local_inventory_by_speed = np.sum(local_distribution * velocity_weights[None, :, :] * volumes[:, None, None], axis=(0, 2))
        local_event_rate_by_speed = np.zeros(local_speed_grid.centers_m_s.shape, dtype=float)
        event_rate_by_projectile: dict[str, float] = {}
        energy_removed_W = 0.0
        target_energy_J = 0.5 * target_species.mass_kg * np.asarray(local_speed_grid.centers_m_s, dtype=float) ** 2
        for deposition in depositions:
            projectile = ion_species(deposition.projectile_species)
            axial_cosine = float(deposition.path_geometry.direction_unit[2] / np.linalg.norm(deposition.path_geometry.direction_unit))
            projectile_rate = 0.0
            for component_index, attenuation_component in enumerate(deposition.attenuation.component_results):
                cell_events = _target_cell_event_rates(attenuation_component, target_species.species_id)
                component_rate = float(np.sum(cell_events))
                if component_rate <= 0.0:
                    continue
                beam_speed = float(beam_speed_from_energy_m_s(float(deposition.component_energies_J[component_index]), projectile.mass_kg))
                kernel = gyroaveraged_sigma_g(local_speed_grid, kinetic_target_state.pitch_grid.centers, beam_speed, axial_cosine, "charge_exchange", nphi)
                reaction_weight = local_distribution * velocity_weights[None, :, :] * kernel[None, :, :]
                normalization = np.sum(reaction_weight, axis=(1, 2))
                speed_weight = np.sum(reaction_weight, axis=2)
                # The local reaction kernel partitions each cell event rate over target speed without changing its total
                speed_probability = np.divide(speed_weight, normalization[:, None], out=np.zeros_like(speed_weight), where=normalization[:, None] > 0.0)
                represented_cells = cell_events > 0.0
                if np.any(represented_cells & (normalization <= 0.0)):
                    raise ValueError("charge exchange attenuation event has no represented target velocity support")
                event_by_speed = np.sum(cell_events[:, None] * speed_probability, axis=0)
                local_event_rate_by_speed += event_by_speed
                projectile_rate += component_rate
                energy_removed_W += float(np.sum(event_by_speed * target_energy_J))
            event_rate_by_projectile[projectile.species_id] = event_rate_by_projectile.get(projectile.species_id, 0.0) + projectile_rate
        reference_rate = float(np.sum(local_event_rate_by_speed))
        if reference_rate <= 0.0:
            continue
        modal_speed_grid, modal_population = _modal_speed_population_particles(kinetic_target_state, target_species, float(geometry.volume_m3))
        frequency, identity_error = _interpolated_loss_frequency(local_speed_grid, local_event_rate_by_speed, local_inventory_by_speed, modal_speed_grid, modal_population, reference_rate)
        states[target_species.species_id] = ChargeExchangeTargetSinkState(
            target_species=target_species,
            speed_grid=modal_speed_grid,
            loss_frequency_s=frequency,
            reference_event_rate_s=reference_rate,
            reference_energy_removal_W=energy_removed_W,
            reference_event_rate_by_projectile_species_s=event_rate_by_projectile,
            model="reaction_weighted_speed_resolved_pitch_averaged_primary_beam_charge_exchange_loss",
            reference_rate_identity_relative_error=identity_error,
            qualified=bool(identity_error <= 1.0e-12),
            limitation="pitch dependence and axial localization are reduced to a speed resolved volume averaged loss frequency",
        )
  
    return states

__all__ = ["ChargeExchangeTargetSinkState", "build_charge_exchange_target_sink_states"]