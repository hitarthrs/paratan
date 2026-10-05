"""Define species requests and shared system states for the modal FBIS solver"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING
import numpy as np
from source_model_revamp.constants import ELECTRON_CHARGE_C
from source_model_revamp.beam.source_from_attenuation import AttenuatedMultiEnergyBeamSource
from source_model_revamp.fbis.species import IonSpecies
from source_model_revamp.fbis.collision_parameters import FBISCollisionParameterState
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid
from source_model_revamp.electrostatic.current_balance import AmbipolarCurrentBalance
from source_model_revamp.fbis.modal.types import ModalElectrostaticProfile, ModalFBISBasis, ModalFBISResult

if TYPE_CHECKING:
    from source_model_revamp.beam.charge_exchange_redistribution import ChargeExchangeTargetSinkState

@dataclass(frozen=True)
class FastIonSpeciesRequest:
    """Input source, species, speed grid, and optional charge exchange target sink for one fast ion species"""
    species: IonSpecies
    speed_grid: SpeedGrid
    attenuated_source: AttenuatedMultiEnergyBeamSource
    charge_exchange_sink_state: ChargeExchangeTargetSinkState | None = None

    def __post_init__(self) -> None:
        """Validate source activity, source array shape, species identity, and shared speed grid requirements"""
        rate = float(self.attenuated_source.total_birth_rate_s)
        if not np.isfinite(rate) or rate < 0.0:
            raise ValueError("fast ion species source rate must be nonnegative and finite")
        source = self.attenuated_source.total_source_v_lambda_per_s
        if rate > 0.0 and source is None:
            raise ValueError("an active fast ion species source requires a spatial velocity source array")
        if source is not None:
            values = np.asarray(source, dtype=float)
            if values.ndim != 3 or values.shape[1] != self.speed_grid.centers_m_s.size:
                raise ValueError("fast ion species source must have shape (n_z, n_speed, n_lambda) on the species speed grid")
            if np.any(~np.isfinite(values)) or np.any(values < 0.0):
                raise ValueError("fast ion species source must be finite and nonnegative")
            if rate <= 0.0 and np.any(values > 0.0):
                raise ValueError("a zero rate fast ion species source cannot contain nonzero source density")
        sink = self.charge_exchange_sink_state
        if sink is not None:
            if sink.target_species != self.species:
                raise ValueError("charge exchange target sink species must match the fast ion species request")
            if sink.speed_grid is not self.speed_grid:
                raise ValueError("charge exchange target sink must use the fast ion species speed grid object")
            if float(sink.reference_event_rate_s) > 0.0 and rate <= 0.0:
                raise ValueError("an active charge exchange target sink requires an active kinetic species source")

@dataclass(frozen=True)
class FastIonSpeciesState:
    """Solved, inactive, or prompt only state for one fast ion species"""
    species: IonSpecies
    speed_grid: SpeedGrid
    active: bool
    source_particle_rate_s: float
    modal_result: ModalFBISResult | None
    collision_state: FBISCollisionParameterState | None = None
    full_device_local_distribution_z_v_pitch: np.ndarray | None = None
    prompt_only_particle_loss_rate_s: float = 0.0
    prompt_only_midplane_kinetic_power_loss_W: float = 0.0
    prompt_only_ion_wall_power_loss_W: float | None = None
    status: str = "active"

    def __post_init__(self) -> None:
        """Validate consistency between activity, modal result, collision state, and prompt only losses"""
        rate = float(self.source_particle_rate_s)
        if not np.isfinite(rate) or rate < 0.0:
            raise ValueError("fast ion species state source rate must be nonnegative and finite")
        if self.active and self.modal_result is None:
            raise ValueError("an active fast ion species state requires a modal result")
        if not self.active and self.modal_result is not None:
            raise ValueError("an inactive fast ion species state cannot contain a modal result")
        if self.active and self.collision_state is None:
            raise ValueError("an active fast ion species state requires a collision state")
        if not self.active and self.collision_state is not None:
            raise ValueError("an inactive fast ion species state cannot contain a collision state")
        if self.collision_state is not None and self.collision_state.fast_ion_species != self.species:
            raise ValueError("fast ion species state and collision state species must match")
        if self.collision_state is not None and self.collision_state.pairwise_collision_state.test_species != self.species:
            raise ValueError("fast ion species state and pairwise collision test species must match")
        if self.modal_result is not None and self.modal_result.species != self.species:
            raise ValueError("fast ion species state and modal result species must match")
        if self.modal_result is not None and self.modal_result.speed_grid is not self.speed_grid:
            raise ValueError("fast ion species state and modal result must use the same speed grid object")
        prompt_rate = float(self.prompt_only_particle_loss_rate_s)
        prompt_power = float(self.prompt_only_midplane_kinetic_power_loss_W)
        prompt_wall_power = self.prompt_only_ion_wall_power_loss_W
        if not np.isfinite(prompt_rate) or prompt_rate < 0.0:
            raise ValueError("prompt only particle loss rate must be finite and nonnegative")
        if not np.isfinite(prompt_power) or prompt_power < 0.0:
            raise ValueError("prompt only midplane kinetic power must be finite and nonnegative")
        if prompt_wall_power is not None and (not np.isfinite(float(prompt_wall_power)) or float(prompt_wall_power) < 0.0):
            raise ValueError("prompt only ion wall power must be finite and nonnegative")
        if self.modal_result is not None and (prompt_rate > 0.0 or prompt_power > 0.0 or prompt_wall_power is not None):
            raise ValueError("prompt only loss fields cannot duplicate an active modal result")
        if prompt_rate > rate + 128.0 * np.finfo(float).eps * max(rate, 1.0):
            raise ValueError("prompt only particle loss rate cannot exceed the species source rate")

@dataclass(frozen=True)
class FastIonSystemState:
    """Species indexed fast ion states sharing one basis and electrostatic closure where active"""
    species_states: Mapping[str, FastIonSpeciesState]
    shared_basis: ModalFBISBasis | None = None
    shared_current_balance: AmbipolarCurrentBalance | None = None
    shared_electrostatic_profile: ModalElectrostaticProfile | None = None
    electron_wall_power_loss_W: float | None = None
    status: str = "active"
    metadata: Mapping[str, object] | None = None

    def __post_init__(self) -> None:
        """Validate shared object identity across active species and canonicalize mappings"""
        states = dict(self.species_states)
        active_states = tuple(state for state in states.values() if state.active)
        if active_states and self.shared_basis is None:
            raise ValueError("an active fast ion system state requires a shared basis")
        for key, state in states.items():
            if str(key) != state.species.species_id:
                raise ValueError("fast ion system state keys must match species identifiers")
            if state.modal_result is not None and state.modal_result.basis is not self.shared_basis:
                raise ValueError("active fast ion species states must reuse the shared basis object")
            if self.shared_electrostatic_profile is not None and state.modal_result is not None and state.modal_result.electrostatic_profile is not self.shared_electrostatic_profile:
                raise ValueError("active fast ion species states must reuse the shared electrostatic profile object")
            if self.shared_current_balance is not None and state.modal_result is not None and state.modal_result.current_balance is not self.shared_current_balance:
                raise ValueError("active fast ion species states must reuse the shared current balance object")
        electron_power = self.electron_wall_power_loss_W
        if electron_power is not None and (not np.isfinite(float(electron_power)) or float(electron_power) < 0.0):
            raise ValueError("fast ion system electron wall power must be finite and nonnegative")
        object.__setattr__(self, "species_states", states)
        object.__setattr__(self, "metadata", {} if self.metadata is None else dict(self.metadata))

    @property
    def active_species_ids(self) -> tuple[str, ...]:
        """Return sorted identifiers for species with an active confined modal state"""
        return tuple(sorted(key for key, state in self.species_states.items() if state.active))

    @property
    def source_species_ids(self) -> tuple[str, ...]:
        """Return sorted identifiers for species with a nonzero physical source"""
        return tuple(sorted(key for key, state in self.species_states.items() if float(state.source_particle_rate_s) > 0.0))

    @property
    def total_fast_ion_inventory_particles(self) -> float:
        """Return the summed active fast ion particle inventory"""
        return float(sum(state.modal_result.inventory_particles for state in self.species_states.values() if state.modal_result is not None))

    @property
    def total_fast_ion_charge_inventory_C(self) -> float:
        """Return the summed active fast ion charge inventory in C"""
        return float(ELECTRON_CHARGE_C * sum(state.species.charge_number * state.modal_result.inventory_particles for state in self.species_states.values() if state.modal_result is not None))

    @property
    def total_ion_particle_loss_rate_s(self) -> float:
        """Return the summed ion particle loss rate in s⁻¹ including prompt only states"""
        return float(sum(float(state.modal_result.ion_particle_loss_rate_s) if state.modal_result is not None else float(state.prompt_only_particle_loss_rate_s) for state in self.species_states.values()))

    @property
    def total_ion_loss_current_A(self) -> float:
        """Return the summed ion loss current in A including prompt only states"""
        return float(ELECTRON_CHARGE_C * sum(state.species.charge_number * (float(state.modal_result.ion_particle_loss_rate_s) if state.modal_result is not None else float(state.prompt_only_particle_loss_rate_s)) for state in self.species_states.values()))

    @property
    def total_ion_wall_power_W(self) -> float | None:
        """Return summed ion wall power in W when every active contribution is available"""
        powers: list[float] = []
        for state in self.species_states.values():
            if state.modal_result is not None:
                if state.modal_result.ion_wall_power_loss_W is None:
                    return None
                powers.append(float(state.modal_result.ion_wall_power_loss_W))
            elif float(state.prompt_only_particle_loss_rate_s) > 0.0:
                if state.prompt_only_ion_wall_power_loss_W is None:
                    return None
                powers.append(float(state.prompt_only_ion_wall_power_loss_W))
        return float(sum(powers))

    def state_for(self, species_id: str) -> FastIonSpeciesState:
        """Return the state for a canonical species identifier"""
        try:
            return self.species_states[str(species_id).strip().lower()]
        except KeyError as exc:
            raise ValueError(f"fast ion species state {species_id!r} is unavailable") from exc
