"""Stationary NBI supported D and T particle source and loss accounting"""
from __future__ import annotations
from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping
import numpy as np
from source_model_revamp.integration.pipeline_types import BeamEnsembleResult, KineticStageResult

@dataclass(frozen=True)
class StationarySpeciesParticleBalance:
    """
    Particle balance for one kinetic ion species
    
    The identity is `ionization + charge exchange gain − charge exchange target sink = terminal loss`
    All rates are in s⁻¹
    """
    species_id: str
    ionization_source_rate_s: float
    charge_exchange_gain_rate_s: float
    charge_exchange_target_sink_rate_s: float
    net_nbi_source_rate_s: float
    terminal_loss_rate_s: float
    residual_rate_s: float
    relative_residual: float
    converged: bool

    def as_metadata(self) -> dict[str, object]:
        """Return the species particle balance as a serializable metadata record"""
        return {
            "species_id": self.species_id,
            "ionization_source_rate_s": self.ionization_source_rate_s,
            "charge_exchange_gain_rate_s": self.charge_exchange_gain_rate_s,
            "charge_exchange_target_sink_rate_s": self.charge_exchange_target_sink_rate_s,
            "net_nbi_source_rate_s": self.net_nbi_source_rate_s,
            "terminal_loss_rate_s": self.terminal_loss_rate_s,
            "residual_rate_s": self.residual_rate_s,
            "relative_residual": self.relative_residual,
            "converged": self.converged,
        }

@dataclass(frozen=True)
class StationaryParticleBalanceState:
    """Immutable collection of separate D and T stationary particle balances and their shared relative tolerance"""
    by_species: Mapping[str, StationarySpeciesParticleBalance]
    relative_tolerance: float
    converged: bool

    def __post_init__(self) -> None:
        """Require at least one species, validate the tolerance, and freeze the species mapping"""
        records = dict(self.by_species)
        if not records:
            raise ValueError("stationary particle balance requires at least one species")
        if not np.isfinite(self.relative_tolerance) or self.relative_tolerance < 0.0:
            raise ValueError("stationary particle balance tolerance must be finite and nonnegative")
        object.__setattr__(self, "by_species", MappingProxyType(records))

    def as_metadata(self) -> dict[str, object]:
        """Return system particle balance metadata including the roles of same isotope and cross isotope charge exchange"""
        return {
            "nbi_supported_particle_balance_model": "ionization_plus_charge_exchange_gain_minus_reaction_weighted_target_sink_equals_terminal_loss",
            "nbi_supported_particle_balance_relative_tolerance": self.relative_tolerance,
            "nbi_supported_particle_balance_converged": self.converged,
            "nbi_supported_particle_balance_by_species": {key: value.as_metadata() for key, value in self.by_species.items()},
            "nbi_supported_net_fueling_authority": "ionization_plus_cross_isotope_charge_exchange_transfer",
            "same_isotope_charge_exchange_inventory_role": "velocity_space_redistribution_only",
            "cross_isotope_charge_exchange_inventory_role": "target_isotope_to_projectile_isotope_transfer",
        }

def _beam_species_rates(beam: BeamEnsembleResult) -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    """Accumulate beam ionization births, charge exchange fast births, and reaction weighted target pumpout rates by ion species"""
    ionization: dict[str, float] = {}
    charge_exchange_gain: dict[str, float] = {}
    charge_exchange_sink: dict[str, float] = {"deuterium": 0.0, "tritium": 0.0}
    for result in beam.beam_results_by_id.values():
        species_id = str(result.projectile_species)
        metadata = result.metadata
        ionization[species_id] = ionization.get(species_id, 0.0) + float(metadata.get("beam_ionization_birth_rate_s", 0.0))
        charge_exchange_gain[species_id] = charge_exchange_gain.get(species_id, 0.0) + float(metadata.get("beam_charge_exchange_fast_birth_rate_s", 0.0))
        charge_exchange_sink["deuterium"] += float(metadata.get("beam_deuterium_target_pumpout_rate_s", 0.0))
        charge_exchange_sink["tritium"] += float(metadata.get("beam_tritium_target_pumpout_rate_s", 0.0))
 
    return ionization, charge_exchange_gain, charge_exchange_sink

def build_stationary_particle_balance(*, beam: BeamEnsembleResult, kinetic: KineticStageResult, relative_tolerance: float) -> StationaryParticleBalanceState:
    """
    Build stationary D and T particle identities from the beam source and kinetic terminal losses
    
    The normalized residual uses the larger of total source side, total sink side, or one particle per second as its scale
    """
    tolerance = float(relative_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("stationary particle balance tolerance must be finite and nonnegative")
    system = kinetic.fast_ion_system_state
    if system is None:
        raise ValueError("stationary particle balance requires a fast ion system state")
    ionization, charge_exchange_gain, charge_exchange_sink = _beam_species_rates(beam)
    species_ids = tuple(sorted(set(ionization) | set(charge_exchange_gain) | set(charge_exchange_sink) | set(system.species_states)))
    balances: dict[str, StationarySpeciesParticleBalance] = {}
    for species_id in species_ids:
        state = system.species_states.get(species_id)
        if state is None:
            terminal_loss = 0.0
        elif state.modal_result is not None:
            terminal_loss = float(state.modal_result.ion_particle_loss_rate_s)
        else:
            terminal_loss = float(state.prompt_only_particle_loss_rate_s)
        source_ionization = float(ionization.get(species_id, 0.0))
        source_charge_exchange = float(charge_exchange_gain.get(species_id, 0.0))
        target_sink = float(charge_exchange_sink.get(species_id, 0.0))
        # Cross isotope charge exchange transfers inventory while same isotope exchange only redistributes velocity
        net_source = source_ionization + source_charge_exchange - target_sink
        residual = net_source - terminal_loss
        scale = max(abs(source_ionization + source_charge_exchange), abs(target_sink + terminal_loss), 1.0)
        relative = abs(residual) / scale
        balances[species_id] = StationarySpeciesParticleBalance(
            species_id=species_id,
            ionization_source_rate_s=source_ionization,
            charge_exchange_gain_rate_s=source_charge_exchange,
            charge_exchange_target_sink_rate_s=target_sink,
            net_nbi_source_rate_s=net_source,
            terminal_loss_rate_s=terminal_loss,
            residual_rate_s=residual,
            relative_residual=relative,
            converged=bool(relative <= tolerance),
        )
   
    return StationaryParticleBalanceState(by_species=balances, relative_tolerance=tolerance, converged=all(value.converged for value in balances.values()))

__all__ = [
    "StationarySpeciesParticleBalance",
    "StationaryParticleBalanceState",
    "build_stationary_particle_balance",
]