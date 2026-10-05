"""
Represented electron energy balance assembly for one closed kinetic state

The stage combines Eq 59 ion to electron energy transfer, prescribed external electron heating, electron wall kinetic loss, and ambipolar ion acceleration power
"""
from __future__ import annotations
from dataclasses import replace
from source_model_revamp.constants import ELECTRON_CHARGE_C
from source_model_revamp.coupling.electron_energy_balance import ElectronEnergyBalanceState, build_electron_energy_balance_state, electron_energy_balance_metadata
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.pipeline_types import ExpanderStageResult, KineticStageResult

def _unavailable_metadata(reason: str) -> dict[str, object]:
    """Return an explicit unavailable electron energy balance metadata state and failure reason"""
    return {
        "physical_electron_energy_equation_active": True,
        "electron_energy_balance_model": "classical_represented_confined_plasma_electron_energy_balance",
        "electron_energy_balance_scope": "classical_represented_confined_plasma",
        "represented_electron_energy_terms_available": False,
        "represented_electron_energy_terms_evaluable": False,
        "represented_electron_energy_terms_qualified": False,
        "electron_energy_balance_failure_reason": str(reason),
        "electron_energy_balance_converged": False,
        "global_plasma_power_balance_claimed": False,
        "global_device_power_balance_claimed": False,
    }

def _mark_fast_system(kinetic: KineticStageResult, metadata: dict[str, object]) -> KineticStageResult:
    """Attach electron energy balance metadata to the kinetic result, shared fast ion system, and each available modal species result"""
    system = kinetic.fast_ion_system_state
    if system is None:
        return replace(kinetic, metadata={**kinetic.metadata, **metadata})
    states = {}
    for species_id, state in system.species_states.items():
        modal = state.modal_result
        if modal is None:
            states[species_id] = state
            continue
        states[species_id] = replace(state, modal_result=replace(modal, metadata={**modal.metadata, **metadata}))
    updated_system = replace(system, species_states=states, metadata={**system.metadata, **metadata})
    
    return replace(kinetic, fast_ion_system_state=updated_system, metadata={**kinetic.metadata, **metadata})

def build_electron_energy_balance_stage(config: SourceModelRunConfig, kinetic: KineticStageResult, *, electron_temperature_J: float, expander: ExpanderStageResult | None=None) -> tuple[KineticStageResult, ElectronEnergyBalanceState | None]:
    """
    Assemble the represented electron energy balance at one electron temperature
    
    When an expander state is available, terminal ion current and terminal electrostatic power replace the throat level acceleration terms
    Unresolved prompt loss uses the throat barrier as a diagnostic fallback and is marked unqualified
    """
    fast_system = kinetic.fast_ion_system_state
    if fast_system is None:
        metadata = _unavailable_metadata("fast ion system state is unavailable")
        return _mark_fast_system(kinetic, metadata), None
    terminal_current = None
    terminal_power = None
    terminal_prompt_fallback = False
    expander_system = None if expander is None else expander.system_state
    if expander_system is not None:
        current_balance = fast_system.shared_current_balance
        if current_balance is None:
            metadata = _unavailable_metadata("terminal electron energy balance requires the shared current balance")
            return _mark_fast_system(kinetic, metadata), None
        prompt_rate = float(expander_system.prompt_unresolved_particle_rate_s)
        barrier = float(current_balance.electron_balance.barrier_energy_J)
        terminal_current = float(expander_system.terminal_current_A + ELECTRON_CHARGE_C * prompt_rate)
        terminal_power = float(expander_system.terminal_electrostatic_power_W + prompt_rate * barrier)
        terminal_prompt_fallback = prompt_rate > 0.0
    try:
        state = build_electron_energy_balance_state(electron_temperature_J=float(electron_temperature_J), fast_ion_system_state=fast_system, external_electron_heating_W=float(config.power_balance.external_electron_heating_W), ambipolar_power_relative_tolerance=float(config.kinetic_electrostatic.total_current_balance_relative_tolerance), current_balance_relative_tolerance=float(config.kinetic_electrostatic.total_current_balance_relative_tolerance), terminal_ion_current_A=terminal_current, terminal_independent_ion_acceleration_power_W=terminal_power)
    except (ValueError, RuntimeError, FloatingPointError) as error:
        metadata = _unavailable_metadata(f"electron energy balance assembly failed: {error}")
        return _mark_fast_system(kinetic, metadata), None
    metadata = electron_energy_balance_metadata(state)
    metadata.update({'electron_energy_balance_converged': False, 'electron_energy_balance_external_source_model': 'prescribed_external_electron_heating', 'electron_energy_balance_fast_source_model': 'direct_pair_resolved_Eq59_electron_collision_energy_moment', 'electron_energy_balance_loss_model': 'exact_Egedal_electron_wall_kinetic_power_plus_terminal_ambipolar_ion_acceleration' if expander_system is not None else 'exact_Egedal_electron_wall_kinetic_power_plus_ambipolar_ion_acceleration', 'electron_energy_balance_uses_terminal_ion_current': expander_system is not None, 'electron_energy_balance_uses_terminal_ion_power': expander_system is not None, 'electron_energy_balance_prompt_terminal_power_uses_throat_fallback': terminal_prompt_fallback, 'electron_energy_balance_prompt_terminal_power_fallback_is_qualified': False})
    return _mark_fast_system(kinetic, metadata), state

__all__ = ["build_electron_energy_balance_stage"]
