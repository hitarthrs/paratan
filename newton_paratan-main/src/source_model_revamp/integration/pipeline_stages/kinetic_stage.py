"""Thin integration wrapper for the modal FBIS kinetic stage"""
from __future__ import annotations
from typing import Any
from source_model_revamp.integration.pipeline_types import BeamEnsembleResult, GeometryStageResult, KineticStageResult
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.integration.modal_stage import build_modal_kinetic_stage
from source_model_revamp.integration.evaluation import OperatingPointEvaluationMode

def build_kinetic_stage(config: SourceModelRunConfig, geometry: GeometryStageResult, beam: BeamEnsembleResult, *, electron_temperature_J: float | None = None, collision_state_recompute_policy: str | None = None, collision_state_is_self_consistent_electron_temperature: bool = False, modal_basis: Any | None = None, eq59_warm_start_state: object | None = None, eq59_warm_start_state_by_species: dict[str, object] | None = None, eq70_initial_potential_energy_J: Any | None = None, initial_representative_fast_density_m3: float | None = None, initial_fast_D_midpoint_density_m3: float | None = None, initial_fast_D_confined_average_density_m3: float | None = None, initial_fast_T_midpoint_density_m3: float | None = None, initial_fast_T_confined_average_density_m3: float | None = None, initial_electron_midpoint_density_m3: float | None = None, initial_electron_confined_average_density_m3: float | None = None, initial_electron_parent_n0_m3: float | None = None, initial_electron_collision_density_m3: float | None = None, initial_eq42_density_profile: object | None = None, scalar_collision_operator_warm_start_by_species: dict[str, tuple[object, object, object, object]] | None = None, scalar_external_fast_field_warm_start_by_test_species: dict[str, tuple[object, ...]] | None = None, scalar_wall_barrier_energy_warm_start_J: float | None = None, evaluation_mode: OperatingPointEvaluationMode = OperatingPointEvaluationMode.FIXED_TEMPERATURE_FINAL) -> KineticStageResult:
    """
    Run the modal FBIS kinetic stage with the supplied operating point continuation state
    
    The wrapper forwards scalar densities, Eq 59 warm states, Eq 70 potential, Eq 42 profile, collision operator fields, wall barrier state, and evaluation mode unchanged
    """
    return build_modal_kinetic_stage(
        config,
        geometry,
        beam,
        electron_temperature_J=electron_temperature_J,
        collision_state_recompute_policy=collision_state_recompute_policy,
        collision_state_is_self_consistent_electron_temperature=(collision_state_is_self_consistent_electron_temperature),
        modal_basis=modal_basis,
        eq59_warm_start_state=eq59_warm_start_state,
        eq59_warm_start_state_by_species=eq59_warm_start_state_by_species,
        eq70_initial_potential_energy_J=eq70_initial_potential_energy_J,
        initial_representative_fast_density_m3=initial_representative_fast_density_m3,
        initial_fast_D_midpoint_density_m3=initial_fast_D_midpoint_density_m3,
        initial_fast_D_confined_average_density_m3=initial_fast_D_confined_average_density_m3,
        initial_fast_T_midpoint_density_m3=initial_fast_T_midpoint_density_m3,
        initial_fast_T_confined_average_density_m3=initial_fast_T_confined_average_density_m3,
        initial_electron_midpoint_density_m3=initial_electron_midpoint_density_m3,
        initial_electron_confined_average_density_m3=initial_electron_confined_average_density_m3,
        initial_electron_parent_n0_m3=initial_electron_parent_n0_m3,
        initial_electron_collision_density_m3=initial_electron_collision_density_m3,
        initial_eq42_density_profile=initial_eq42_density_profile,
        scalar_collision_operator_warm_start_by_species=scalar_collision_operator_warm_start_by_species,
        scalar_external_fast_field_warm_start_by_test_species=scalar_external_fast_field_warm_start_by_test_species,
        scalar_wall_barrier_energy_warm_start_J=scalar_wall_barrier_energy_warm_start_J,
        evaluation_mode=evaluation_mode,
    )
