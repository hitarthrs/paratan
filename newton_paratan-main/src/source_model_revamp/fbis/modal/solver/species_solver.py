"""Single species modal FBIS entry point routed through the shared system solver"""
from __future__ import annotations
from typing import Any
from source_model_revamp.fbis.modal.species_state import FastIonSpeciesRequest, FastIonSpeciesState
from source_model_revamp.fbis.modal.solver.system_solver import solve_fast_ion_system

def solve_fast_ion_species(*, request: FastIonSpeciesRequest, **solver_arguments: Any) -> FastIonSpeciesState:
    """Solve one requested fast ion species through the shared system solver path"""
    collision_state = solver_arguments.pop("collision_state", None)
    warm_start_state = solver_arguments.pop("_eq59_warm_start_state", None)
    compatibility_policy = solver_arguments.pop("_eq59_warm_start_compatibility_policy", None)
    if compatibility_policy is not None and compatibility_policy not in {"exact_physics", "closure_continuation"}:
        raise ValueError("_eq59_warm_start_compatibility_policy must be exact_physics or closure_continuation")
    system = solve_fast_ion_system(
        requests={request.species.species_id: request},
        collision_states_by_species={} if collision_state is None else {request.species.species_id: collision_state},
        eq59_warm_start_states_by_species={} if warm_start_state is None else {request.species.species_id: warm_start_state},
        eq59_warm_start_compatibility_policy_by_species={} if compatibility_policy is None else {request.species.species_id: compatibility_policy},
        **solver_arguments,
    )

    return system.state_for(request.species.species_id)
