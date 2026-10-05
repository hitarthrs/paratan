"""Compact run assessment built from stage metadata after direct state validation"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any
import numpy as np

@dataclass(frozen=True)
class AssessmentRecord:
    """
    One convergence or applicability assessment
    
    `passed` is True for a passed check, False for a failed check, and None when the check is not applicable or cannot be evaluated from the available evidence
    """
    name: str
    passed: bool | None
    value: object = None
    tolerance: object = None
    reason: str | None = None

    def to_metadata(self) -> dict[str, object]:
        """Return the compact serializable assessment record"""
        return {"name": self.name, "passed": self.passed, "value": self.value, "tolerance": self.tolerance, "reason": self.reason}

@dataclass(frozen=True)
class RunAssessment:
    """
    Canonical convergence and applicability records for one completed source model run
    
    A None assessment state does not fail the corresponding aggregate
    """
    convergence: tuple[AssessmentRecord, ...]
    applicability: tuple[AssessmentRecord, ...]

    @staticmethod
    def _aggregate(records: tuple[AssessmentRecord, ...]) -> bool:
        """Return False only when at least one applicable assessment explicitly failed"""
        return not any(record.passed is False for record in records)

    @property
    def calculation_valid(self) -> bool:
        """Return True because invalid calculated states are rejected before RunAssessment is constructed"""
        return True

    @property
    def numerically_converged(self) -> bool:
        """Return the aggregate numerical convergence state"""
        return self._aggregate(self.convergence)

    @property
    def model_applicable(self) -> bool:
        """Return the aggregate model applicability state"""
        return self._aggregate(self.applicability)

    @property
    def failed_convergence_checks(self) -> tuple[str, ...]:
        """Return names of convergence assessments that explicitly failed"""
        return tuple(record.name for record in self.convergence if record.passed is False)

    @property
    def applicability_limitations(self) -> tuple[str, ...]:
        """Return names of applicability assessments that explicitly failed"""
        return tuple(record.name for record in self.applicability if record.passed is False)

    @property
    def result_status(self) -> str:
        """Return the combined convergence and applicability status string"""
        if self.numerically_converged and self.model_applicable:
            return "qualified"
        if not self.numerically_converged and self.model_applicable:
            return "unconverged"
        if self.numerically_converged:
            return "applicability_limited"
        return "unconverged_and_applicability_limited"

CONVERGENCE_CHECKS = ('electron_temperature_convergence', 'source_projection_convergence', 'eigenbasis_convergence', 'eq42_density_convergence', 'collision_composition_convergence', 'beam_density_fixed_point_convergence', 'eq59_convergence', 'eq70_convergence', 'local_reconstruction_convergence', 'expander_potential_convergence', 'terminal_current_convergence')
APPLICABILITY_CHECKS = ("equivalent_end_assumption", "collision_model_applicability", "expander_fast_boundary_applicability", "expander_velocity_mapping", "fast_ion_fusion_burnup_applicability", "nbi_supported_stationary_closure", "bosch_hale_domain_coverage")

def _literal_boolean(value: object) -> bool | None:
    """Return a Python boolean only for literal boolean evidence and otherwise return None"""
    return bool(value) if isinstance(value, (bool, np.bool_)) else None

def _failure_reason(items: list[tuple[str, str | None]]) -> str | None:
    """Join named assessment failures into one compact reason string"""
    return "; ".join(f"{name}: {reason or 'evidence unavailable'}" for name, reason in items) or None

_ITERATIVE_DENSITY_CLOSURE_MODELS = frozenset(('nbi_supported_stationary',))

def _density_closure_model(evidence: Mapping[str, Any]) -> str:
    """Return the normalized modal density closure model from stage evidence"""
    return str(evidence.get("modal_density_closure_model", "")).strip().lower()

def _iterative_density_closure_active(evidence: Mapping[str, Any]) -> bool:
    """Return whether the active density closure requires the coupled iterative operating point path"""
    return _density_closure_model(evidence) in _ITERATIVE_DENSITY_CLOSURE_MODELS

def _electron_temperature(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """
    Assess self consistent electron temperature convergence
    
    The check combines the absolute or relative represented power residual with preflight, scalar search, final qualification, final residual consistency, and diagnostic root evidence
    """
    if str(evidence.get("electron_temperature_mode", "")).strip().lower() != "self_consistent_electron_energy":
        return AssessmentRecord("electron_temperature_convergence", None, reason="not applicable")
    settings = evidence.get("electron_temperature_solve_numerics_settings")
    controls = settings if isinstance(settings, Mapping) else {}
    residual_W = evidence.get("electron_temperature_power_residual_W")
    relative_residual = evidence.get("electron_temperature_relative_power_residual")
    absolute_tolerance = controls.get("residual_tolerance_W")
    relative_tolerance = controls.get("relative_residual_tolerance")
    residual_passed = False
    try:
        residual_passed = bool(
            (np.isfinite(float(residual_W)) and np.isfinite(float(absolute_tolerance)) and abs(float(residual_W)) <= float(absolute_tolerance))
            or (np.isfinite(float(relative_residual)) and np.isfinite(float(relative_tolerance)) and abs(float(relative_residual)) <= float(relative_tolerance))
        )
    except (TypeError, ValueError):
        residual_passed = False
    skipped = evidence.get("electron_temperature_search_skipped") is True
    requirements = [
        ("electron_temperature_power_balance_closure", evidence.get("electron_temperature_solve_converged") is True and evidence.get("electron_temperature_power_balance_converged") is True and residual_passed),
        ("electron_temperature_applicability_preflight", evidence.get("electron_temperature_preflight_passed") is True),
        ("electron_temperature_final_qualification", evidence.get("electron_temperature_final_qualification_passed") is True),
        ("electron_temperature_final_residual_consistency", evidence.get("electron_temperature_final_residual_consistent_with_trial") is True),
        ("electron_temperature_diagnostic_root_convergence", evidence.get("electron_temperature_diagnostic_root_converged") is True),
    ]
    if not skipped:
        requirements[2:2] = [
            ("electron_temperature_search_bracket_found", evidence.get("electron_temperature_search_bracket_found") is True),
            ("electron_temperature_search_convergence", evidence.get("electron_temperature_search_converged") is True),
        ]
    failed = [(name, f"{name}_not_passed") for name, passed in requirements if not passed]

    return AssessmentRecord(
        "electron_temperature_convergence",
        not failed,
        value={
            "solved_electron_temperature_keV": evidence.get("electron_temperature_solved_keV"),
            "power_residual_W": residual_W,
            "relative_power_residual": relative_residual,
            "search_bracket_found": evidence.get("electron_temperature_search_bracket_found"),
            "search_converged": evidence.get("electron_temperature_search_converged"),
            "final_qualification_passed": evidence.get("electron_temperature_final_qualification_passed"),
            "final_residual_consistent_with_trial": evidence.get("electron_temperature_final_residual_consistent_with_trial"),
        },
        tolerance={"power_residual_W": absolute_tolerance, "relative_power_residual": relative_tolerance},
        reason=_failure_reason(failed),
    )

def _source_projection(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """Assess conservation of the modal source projection rate"""
    passed = evidence.get("modal_source_projection_rate_converged") is True
   
    return AssessmentRecord("source_projection_convergence", passed, value=evidence.get("modal_source_projection_relative_error"), tolerance=evidence.get("modal_source_projection_relative_tolerance"), reason=None if passed else "modal_source_projection_rate: modal_source_projection_rate_not_converged")

def _eigenbasis(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """
    Assess physical eigenbasis grid convergence and retained physical mode convergence when that second check applies
    
    The retained mode test is omitted when fewer than two physical modes make it inapplicable
    """
    basis_assessed = evidence.get("kinetic_basis_convergence_assessed") is True
    basis_passed = bool(basis_assessed and evidence.get("kinetic_basis_converged") is True)
    retained_applicable = bool(evidence.get("kinetic_retained_mode_count_applicable", int(evidence.get("modal_n_physical_modes", 0) or 0) >= 2))
    names = ["physical_eigenbasis_grid_convergence"]
    values: list[object] = [evidence.get("modal_basis_max_eigenvalue_relative_change")]
    tolerances: list[object] = [evidence.get("modal_basis_eigenvalue_relative_tolerance")]
    states = [basis_passed]
    failures = [] if basis_passed else [(names[0], "modal_basis_grid_convergence_failed" if basis_assessed else "modal_basis_grid_convergence_not_assessed")]
    if retained_applicable:
        retained_assessed = evidence.get("modal_retained_mode_count_assessed") is True
        retained_passed = bool(retained_assessed and evidence.get("modal_retained_mode_count_converged") is True)
        names.append("retained_physical_mode_count_convergence")
        values.append({"selected_mode_count": evidence.get("modal_retained_mode_selected_count"), "final_distribution_relative_change": evidence.get("modal_retained_mode_final_distribution_relative_change"), "final_inventory_relative_change": evidence.get("modal_retained_mode_final_inventory_relative_change")})
        tolerances.append(evidence.get("modal_retained_mode_relative_tolerance"))
        states.append(retained_passed)
        if not retained_passed:
            failures.append((names[-1], "modal_retained_mode_count_convergence_failed" if retained_assessed else "modal_retained_mode_count_convergence_not_assessed"))
    value = values[0] if len(values) == 1 else dict(zip(names, values))
    tolerance = tolerances[0] if len(tolerances) == 1 else dict(zip(names, tolerances))
  
    return AssessmentRecord("eigenbasis_convergence", all(states), value=value, tolerance=tolerance, reason=_failure_reason(failures))

def _eq42(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """
    Assess Eq 42 density weighting closure
    
    The check requires the density ratio to be coupled to the modal basis, the density profile to satisfy mirror symmetry, and the density basis outer iteration to converge when active
    """
    if evidence.get("modal_eq42_nonuniform_density_weighting_required") is not True:
        return AssessmentRecord("eq42_density_convergence", None, reason="not applicable")
 
    coupled = evidence.get("modal_eq42_density_weighting_coupled") is True
    outer_applicable = evidence.get("modal_eq42_outer_iteration_applicable") is True
    outer = evidence.get("modal_eq42_outer_iteration_converged") is True
    symmetric = evidence.get("modal_eq42_density_profile_symmetric") is True
    closure = bool(coupled and symmetric and (not outer_applicable or outer))
    values = {
        "eq42_density_weighted_bounce_operator": {"model": evidence.get("modal_eq42_density_weighting_model"), "density_weighting_coupled": coupled, "nonuniform_density_weighting_active": evidence.get("modal_eq42_nonuniform_density_weighting_active"), "eta_mapping_density_weighted": evidence.get("modal_eq42_eta_mapping_density_weighted")},
        "eq42_density_profile_symmetry": evidence.get("modal_eq42_density_profile_symmetry_relative_error"),
        "eq42_density_closure": {"density_weighting_coupled": coupled, "outer_iteration_applicable": outer_applicable, "outer_iteration_converged": outer, "profile_symmetric": symmetric, "final_density_shape_relative_change": evidence.get("modal_eq42_density_shape_final_relative_change")},
    }
    tolerances = {
        "eq42_density_weighted_bounce_operator": True,
        "eq42_density_profile_symmetry": evidence.get("modal_eq42_density_profile_symmetry_tolerance"),
        "eq42_density_closure": {"density_shape_relative_tolerance": evidence.get("modal_eq42_density_shape_relative_tolerance"), "basis_eigenvalue_relative_tolerance": evidence.get("modal_eq42_basis_eigenvalue_relative_tolerance"), "basis_eigenfunction_minimum_overlap": evidence.get("modal_eq42_basis_eigenfunction_overlap_tolerance")},
    }
    failures = []
    if not coupled:
        failures.append(("eq42_density_weighted_bounce_operator", "eq42_density_ratio_not_coupled_to_modal_basis"))
    if outer_applicable:
        values = {"eq42_density_weighted_bounce_operator": values["eq42_density_weighted_bounce_operator"], "eq42_density_basis_outer_iteration": evidence.get("modal_eq42_density_shape_final_relative_change"), **{key: value for key, value in values.items() if key != "eq42_density_weighted_bounce_operator"}}
        tolerances = {"eq42_density_weighted_bounce_operator": True, "eq42_density_basis_outer_iteration": tolerances["eq42_density_closure"], **{key: value for key, value in tolerances.items() if key not in {"eq42_density_weighted_bounce_operator", "eq42_density_closure"}}, "eq42_density_closure": tolerances["eq42_density_closure"]}
        if not outer:
            failures.append(("eq42_density_basis_outer_iteration", str(evidence.get("modal_eq42_outer_iteration_failure_reason") or "eq42_density_basis_outer_iteration_not_converged")))
    if not symmetric:
        failures.append(("eq42_density_profile_symmetry", "eq42_density_profile_exceeds_symmetric_modal_well_tolerance"))
    if not closure:
        failures.append(("eq42_density_closure", "eq42_density_closure_not_converged_or_applicable"))
   
    return AssessmentRecord("eq42_density_convergence", not failures, value=values, tolerance=tolerances, reason=_failure_reason(failures))

def _collision_convergence(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """
    Assess collision composition and scalar density closure convergence
    
    For the NBI supported stationary closure the check also requires reduced fast D and fast T cross species collision fields when more than one fast species is active
    """
    model = _density_closure_model(evidence)
    source_state = str(evidence.get("kinetic_source_state", "")).strip().lower()
  
    applicable = bool(_iterative_density_closure_active(evidence) and evidence.get("modal_density_closure_applicable", source_state != "no_confined_fbis_source") is True)
    if not applicable:
        return AssessmentRecord("collision_composition_convergence", None, reason="not applicable")
  
    closure = evidence.get("modal_density_closure_converged") is True
    residual = evidence.get("modal_density_closure_relative_error")
    tolerance = evidence.get("modal_density_closure_relative_tolerance")
    values: dict[str, object] = {"collision_density_closure": residual}
    tolerances: dict[str, object] = {"collision_density_closure": tolerance}
    failures = []
    if not closure:
        failures.append(("collision_density_closure", "modal_collision_density_closure_not_converged"))
    if model == 'nbi_supported_stationary':
        active_species = tuple((str(value) for value in evidence.get('active_fast_species', ()) if str(value)))
        cross_species_required = len(set(active_species)) > 1
        cross_species_available = evidence.get('fast_D_fast_T_cross_collisions_available') is True
        values['reduced_fast_cross_species_collisions'] = {'required': cross_species_required, 'available': cross_species_available}
        tolerances['reduced_fast_cross_species_collisions'] = True
        if cross_species_required and (not cross_species_available):
            failures.append(('reduced_fast_cross_species_collisions', 'fast_D_fast_T_cross_collision_state_unavailable'))
  
    return AssessmentRecord("collision_composition_convergence", not failures, value=values, tolerance=tolerances, reason=_failure_reason(failures))

def _beam_density(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """
    Assess the beam and density fixed point for iterative plasma closure
    
    The record combines density closure, total birth rate, deposited power, axial birth shape, component deposition, candidate confirmation, and final consistency evidence
    """
    active = _iterative_density_closure_active(evidence)
    source_state = str(evidence.get("kinetic_source_state", "")).strip().lower()
  
    if not active or evidence.get("modal_density_closure_applicable", source_state != "no_confined_fbis_source") is not True:
        return AssessmentRecord("beam_density_fixed_point_convergence", None, reason="not applicable")
  
    controls = evidence.get("beam_density_coupling_numerics_settings")
    controls = controls if isinstance(controls, Mapping) else {}
    states = {
        "beam_birth_profile_convergence": evidence.get("beam_axial_birth_profile_converged") is True,
        "beam_component_deposition_convergence": evidence.get("beam_energy_component_deposition_converged") is True,
        "final_beam_kinetic_consistency_recompute": evidence.get("final_consistency_recompute_passed") is True,
        "beam_density_fixed_point_convergence": evidence.get("beam_density_fixed_point_converged") is True,
    }
    values = {
        "beam_birth_profile_convergence": {"pointwise_relative_change": evidence.get("beam_birth_profile_pointwise_relative_change"), "volume_L2_relative_change": evidence.get("beam_birth_profile_volume_L2_relative_change"), "absolute_reference_change": evidence.get("beam_birth_profile_absolute_reference_change")},
        "beam_component_deposition_convergence": {"pointwise_relative_change": evidence.get("beam_component_pointwise_relative_change"), "L2_relative_change": evidence.get("beam_component_L2_relative_change"), "absolute_reference_change": evidence.get("beam_component_absolute_reference_change")},
        "final_beam_kinetic_consistency_recompute": {"candidate_fixed_point_passed": evidence.get("candidate_fixed_point_passed"), "final_consistency_recompute_passed": evidence.get("final_consistency_recompute_passed")},
        "beam_density_fixed_point_convergence": {"density_pointwise_relative_change": evidence.get("density_pointwise_relative_change"), "density_volume_L2_relative_change": evidence.get("density_volume_L2_relative_change"), "density_absolute_reference_change": evidence.get("density_absolute_reference_change"), "birth_rate_relative_change": evidence.get("beam_birth_rate_relative_change"), "deposited_power_relative_change": evidence.get("beam_deposited_power_relative_change"), "history": evidence.get("beam_density_coupling_history")},
    }
    tolerances = {
        "beam_birth_profile_convergence": controls.get("axial_birth_profile_relative_tolerance"),
        "beam_component_deposition_convergence": controls.get("energy_component_deposition_relative_tolerance"),
        "final_beam_kinetic_consistency_recompute": True,
        "beam_density_fixed_point_convergence": {"pointwise": controls.get("electron_profile_pointwise_relative_tolerance"), "volume_L2": controls.get("electron_profile_volume_L2_relative_tolerance"), "absolute_reference": controls.get("electron_profile_absolute_reference_tolerance"), "birth_rate": controls.get("total_birth_rate_relative_tolerance"), "deposited_power": controls.get("deposited_power_relative_tolerance")},
    }
    messages = {
        "beam_birth_profile_convergence": "beam_birth_profile_not_converged",
        "beam_component_deposition_convergence": "beam_component_deposition_not_converged",
        "final_beam_kinetic_consistency_recompute": "final_beam_kinetic_consistency_recompute_not_converged",
        "beam_density_fixed_point_convergence": "beam_density_fixed_point_not_converged",
    }
    failures = [(name, messages[name]) for name, passed in states.items() if not passed]
  
    return AssessmentRecord("beam_density_fixed_point_convergence", not failures, value=values, tolerance=tolerances, reason=_failure_reason(failures))

def _eq59(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """
    Assess Eq 59 nonlinear convergence, speed resolution, and speed domain convergence
    
    The record retains nonlinear residual, high speed tail fractions, selected speed grid, and the largest nested speed refinement change
    """
    if evidence.get("kinetic_eq59_convergence_applicable") is not True:
        return AssessmentRecord("eq59_convergence", None, reason="not applicable")
  
    speed_changes = []
    for key in ("modal_rosenbluth_speed_common_domain_distribution_relative_change_history", "modal_rosenbluth_speed_inventory_relative_change_history", "modal_rosenbluth_speed_effective_energy_relative_change_history", "modal_rosenbluth_speed_particle_loss_relative_change_history", "modal_rosenbluth_speed_power_loss_relative_change_history", "modal_rosenbluth_speed_source_to_loss_ratio_relative_change_history"):
        for value in evidence.get(key, ()) or ():
            try:
                numeric = float(value)
            except (TypeError, ValueError):
                continue
            if np.isfinite(numeric):
                speed_changes.append(abs(numeric))
    passed = bool(evidence.get("kinetic_eq59_convergence_assessed") is True and evidence.get("kinetic_eq59_converged") is True)
    value = {"nonlinear_relative_residual": evidence.get("modal_rosenbluth_final_relative_residual"), "high_speed_tail_population_fraction": evidence.get("modal_rosenbluth_high_speed_tail_population_fraction"), "high_speed_tail_energy_fraction": evidence.get("modal_rosenbluth_high_speed_tail_energy_fraction"), "speed_resolution_assessed": evidence.get("modal_rosenbluth_speed_resolution_assessed"), "speed_resolution_converged": evidence.get("modal_rosenbluth_speed_resolution_converged"), "speed_domain_assessed": evidence.get("modal_rosenbluth_speed_domain_assessed"), "speed_domain_converged": evidence.get("modal_rosenbluth_speed_domain_converged"), "maximum_nested_speed_relative_change": max(speed_changes) if speed_changes else None, "selected_speed_cell_count": evidence.get("modal_rosenbluth_selected_speed_cell_count"), "selected_speed_domain_upper_m_s": evidence.get("modal_rosenbluth_selected_speed_domain_upper_m_s"), "selected_grid_matched": evidence.get("modal_rosenbluth_speed_selected_grid_matched")}
    tolerance = {"relative_tolerance": evidence.get("modal_rosenbluth_relative_tolerance"), "high_speed_tail_population_tolerance": evidence.get("modal_rosenbluth_high_speed_tail_population_tolerance"), "high_speed_tail_energy_tolerance": evidence.get("modal_rosenbluth_high_speed_tail_energy_tolerance"), "speed_convergence_relative_tolerance": evidence.get("modal_rosenbluth_speed_convergence_relative_tolerance")}
  
    return AssessmentRecord("eq59_convergence", passed, value=value, tolerance=tolerance, reason=None if passed else "eq59_nonlinear_and_speed_domain_convergence: modal_eq59_nonlinear_or_speed_domain_convergence_failed")

def _eq70(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """
    Assess Eq 70 quasineutrality and confined electron ion inventory consistency
    
    The interior profile check is omitted for a zero potential reference or a source state without a confined FBIS population
    """
    source_state = str(evidence.get("kinetic_source_state", "")).strip().lower()
    zero_reference = str(evidence.get("electrostatic_feedback_model", "")).strip().lower() == "zero_potential_magnetic_reference"
    eq70_applicable = evidence.get("kinetic_eq70_profile_applicable", source_state != "no_confined_fbis_source" and not zero_reference) is True
    active_density = _iterative_density_closure_active(evidence)
    density_applicable = bool(active_density and evidence.get("modal_density_closure_applicable", source_state != "no_confined_fbis_source") is True)
    if not eq70_applicable and not density_applicable:
        return AssessmentRecord("eq70_convergence", None, reason="not applicable")
  
    values = {}
    tolerances = {}
    failures = []
    if density_applicable:
        inventory_passed = evidence.get("electron_ion_confined_inventory_check_passed") is True
        values["electron_ion_confined_inventory_consistency"] = evidence.get("electron_ion_confined_inventory_relative_error")
        tolerances["electron_ion_confined_inventory_consistency"] = evidence.get("density_inventory_relative_tolerance", evidence.get("modal_phi_z_relative_tolerance"))
        if not inventory_passed:
            failures.append(("electron_ion_confined_inventory_consistency", "electron_ion_confined_inventory_identity_failed"))
    if eq70_applicable:
        passed = evidence.get("modal_phi_z_interior_converged", evidence.get("modal_phi_z_converged", False)) is True
        values["eq70_interior_quasineutrality_convergence"] = {"maximum_supported_cell_relative_residual": evidence.get("modal_phi_z_max_relative_quasineutrality_error"), "maximum_absolute_residual_normalized_to_midplane": evidence.get("modal_phi_z_maximum_absolute_density_residual_normalized_to_reference"), "volume_integrated_absolute_mismatch_fraction": evidence.get("modal_phi_z_volume_integrated_absolute_particle_mismatch_fraction"), "exact_midplane_residual": evidence.get("modal_phi_z_exact_midplane_quasineutrality_residual"), "interior_converged": passed}
        tolerances["eq70_interior_quasineutrality_convergence"] = evidence.get("modal_phi_z_relative_tolerance")
        if not passed:
            failures.append(("eq70_interior_quasineutrality_convergence", "modal_eq70_interior_profile_invalid"))
    value = next(iter(values.values())) if len(values) == 1 else values
    tolerance = next(iter(tolerances.values())) if len(tolerances) == 1 else tolerances
   
    return AssessmentRecord("eq70_convergence", not failures, value=value, tolerance=tolerance, reason=_failure_reason(failures))

def _local_reconstruction(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """
    Assess local physical speed reconstruction convergence
    
    Species resolved evidence is used when available, otherwise the legacy aggregate density, energy, inventory, and overflow diagnostics are assessed
    """
    species_records = evidence.get("species_local_reconstruction_convergence")
    if isinstance(species_records, Mapping) and species_records:
        if evidence.get("kinetic_local_speed_refinement_applicable") is not True:
            return AssessmentRecord("local_reconstruction_convergence", None, reason="not applicable")
  
        active_species = tuple(str(value) for value in evidence.get("active_fast_species", ()) if str(value))
        required_species = active_species or tuple(str(value) for value in species_records)
        values: dict[str, object] = {}
        tolerances: dict[str, object] = {}
        failures = []
        for species_id in required_species:
            item = species_records.get(species_id)
            if not isinstance(item, Mapping):
                values[species_id] = None
                tolerances[species_id] = None
                failures.append((species_id, "local_reconstruction_convergence_evidence_unavailable"))
                continue
            assessed = item.get("assessed") is True
            converged = item.get("converged") is True
            values[species_id] = dict(item)
            tolerances[species_id] = item.get("relative_tolerance")
            if not assessed:
                failures.append((species_id, "local_physical_speed_grid_convergence_not_assessed"))
            elif not converged:
                failures.append((species_id, str(item.get("failure_reason") or "local_physical_speed_grid_convergence_failed")))
        return AssessmentRecord("local_reconstruction_convergence", not failures, value=values, tolerance=tolerances, reason=_failure_reason(failures))
    applicable = bool(evidence.get("kinetic_local_speed_refinement_applicable", False) and (evidence.get("modal_local_distribution_shape_z_v_pitch") is not None or evidence.get("modal_local_distribution_z_v_pitch") is not None))
 
    if not applicable:
        return AssessmentRecord("local_reconstruction_convergence", None, reason="not applicable")
  
    assessed = evidence.get("modal_local_speed_refinement_assessed") is True
    passed = bool(assessed and evidence.get("modal_local_speed_refinement_converged") is True)
    value = {"density_profile_relative_change": evidence.get("modal_local_speed_final_density_relative_change"), "energy_density_profile_relative_change": evidence.get("modal_local_speed_final_energy_density_relative_change"), "particle_inventory_relative_change": evidence.get("modal_local_speed_final_particle_inventory_relative_change"), "energy_inventory_relative_change": evidence.get("modal_local_speed_final_energy_inventory_relative_change"), "domain_sufficient": evidence.get("modal_local_speed_domain_sufficient"), "overflow_particle_fraction": evidence.get("modal_local_speed_overflow_particle_fraction"), "overflow_energy_fraction": evidence.get("modal_local_speed_overflow_energy_fraction"), "tail_clipped": evidence.get("modal_local_speed_tail_clipped")}
    reason = None if passed else str(evidence.get("modal_local_speed_refinement_failure_reason") or "local_physical_speed_grid_convergence_not_assessed")
  
    return AssessmentRecord("local_reconstruction_convergence", passed, value=value, tolerance=evidence.get("modal_local_speed_refinement_relative_tolerance"), reason=None if passed else f"local_physical_speed_grid_convergence: {reason}")

def _terminal_current(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """Assess terminal ion electron current closure, component identities, end symmetry, and unresolved prompt fraction"""
    measured = _literal_boolean(evidence.get("terminal_current_balance_qualified"))
    passed = measured
    value = {"target_ion_current_A": evidence.get("terminal_current_balance_target_ion_current_A"), "electron_current_A": evidence.get("terminal_current_balance_electron_current_A"), "current_relative_residual": evidence.get("terminal_current_balance_current_relative_residual"), "final_component_relative_change": evidence.get("terminal_current_balance_final_consistency_relative_change"), "component_identity_relative_error": evidence.get("terminal_current_balance_component_current_identity_relative_error"), "side_sum_identity_relative_error": evidence.get("terminal_current_balance_side_sum_identity_relative_error"), "end_asymmetry_relative_error": evidence.get("terminal_current_end_asymmetry_relative_error"), "prompt_unresolved_fraction": evidence.get("terminal_current_prompt_unresolved_fraction")}
    tolerance = {"relative_tolerance": evidence.get("terminal_current_balance_relative_tolerance"), "end_asymmetry_relative_tolerance": evidence.get("terminal_current_end_asymmetry_relative_tolerance"), "prompt_unresolved_fraction_tolerance": evidence.get("terminal_current_prompt_unresolved_fraction_tolerance")}
    reason = None if passed is True else str(evidence.get("terminal_current_balance_failure_reason") or ("terminal_current_balance_not_evaluated" if passed is None else "terminal_current_balance_not_qualified"))
 
    return AssessmentRecord("terminal_current_convergence", passed, value=value, tolerance=tolerance, reason=None if passed is True else f"terminal_ion_electron_current_balance: {reason}")

def _equivalent_end(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """
    Assess the symmetric end and equivalent end assumptions
    
    When total current balance applies the record also includes the side resolved total ion current asymmetry check
    """
    names = ["symmetric_end_model_applicability"]
    states: list[bool | None] = [evidence.get("kinetic_symmetric_end_model_applicable") is True]
    values = [evidence.get("kinetic_left_right_mirror_ratio_relative_difference")]
    tolerances = [evidence.get("modal_loss_convention_relative_tolerance")]
    reasons = [] if states[0] else [(names[0], "asymmetric_branch_kinetics_not_implemented")]
    total_applicable = evidence.get("total_current_balance_applicable", False) is True
    if total_applicable:
        equivalent = evidence.get("equivalent_end_model_applicable", evidence.get("total_current_balance_end_model_applicable", False)) is True
        names.append("equivalent_end_model_applicability")
        states.append(equivalent)
        values.append(evidence.get("equivalent_end_model_applicable"))
        tolerances.append(None)
        if not equivalent:
            reasons.append((names[-1], "equivalent_end_model_not_applicable"))
        asymmetry_applicable = evidence.get("total_current_balance_end_asymmetry_applicable") is True
        asymmetry = evidence.get("total_current_balance_end_asymmetry_check_passed") is True if asymmetry_applicable else None
        names.append("total_ion_current_end_asymmetry")
        states.append(asymmetry)
        values.append(evidence.get("total_current_balance_end_asymmetry_relative_error") if asymmetry_applicable else None)
        tolerances.append(evidence.get("total_current_balance_end_asymmetry_relative_tolerance") if asymmetry_applicable else None)
        if asymmetry is False:
            reasons.append((names[-1], "side_resolved_total_ion_current_is_asymmetric"))
        elif asymmetry is None and not reasons:
            reasons.append((names[-1], "evidence unavailable"))
    passed = False if any(state is False for state in states) else None if any(state is None for state in states) else True
    value = values[0] if len(values) == 1 else dict(zip(names, values))
    tolerance = tolerances[0] if len(tolerances) == 1 else dict(zip(names, tolerances))
  
    return AssessmentRecord("equivalent_end_assumption", passed, value=value, tolerance=tolerance, reason=_failure_reason(reasons))

def _collision_applicability(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """Assess the implemented multispecies collision composition and Coulomb logarithm applicability evidence"""
    definitions = (
        ("multispecies_beta_m_applicability", True, "mixed_D_T_beta_m_mass_rule_unavailable"),
        ("fast_ion_collision_decomposition_applicability", float(evidence.get("collision_fast_D_density_m3", 0.0) or 0.0) > 0.0, "thermal_background_plus_non_Maxwellian_fast_D_collision_decomposition_unavailable"),
        ("thermal_screening_composition_applicability", True, "thermal_screening_composition_unavailable"),
        ("ion_fast_coulomb_log_applicability", True, "mixed_field_ion_fast_Coulomb_log_unavailable"),
        ("collision_physics_applicability", True, "collision_physics_reference_closure_unavailable"),
    )
    selected = [(name, _literal_boolean(evidence.get(name)), message) for name, applicable, message in definitions if applicable]
 
    if not selected:
        return AssessmentRecord("collision_model_applicability", None, reason="not applicable")
 
    limitations = tuple(f"{name}: {message}" for name, passed, message in selected if passed is not True)
    passed = False if any(value is False for _, value, _ in selected) else None if any(value is None for _, value, _ in selected) else True
 
    return AssessmentRecord("collision_model_applicability", passed, value={"limitations": limitations, "collision_scope": evidence.get("collision_scope")}, tolerance=True, reason="; ".join(limitations) or None)

def fast_ion_fusion_burnup_assessment(fusion: Mapping[str, Any]) -> AssessmentRecord:
    """Assess whether fast ion fusion burnup is represented by an Eq 59 sink or is negligible against the configured source and terminal loss scale"""
    burnup_D = float(fusion.get("fusion_fast_deuterium_burnup_rate_s", 0.0))
    burnup_T = float(fusion.get("fusion_fast_tritium_burnup_rate_s", 0.0))
    applicable = max(0.0, burnup_D) + max(0.0, burnup_T) > 0.0
  
    if not applicable:
        return AssessmentRecord("fast_ion_fusion_burnup_applicability", None, reason="not applicable")
 
    sink = fusion.get('fusion_fast_population_burnup_in_Eq59_sink') is True
    fusion_authority = fusion.get("fusion_fast_population_burnup_negligibility_tolerance_available") is True
    tolerance_available = fusion_authority
    assessed = fusion.get('fusion_fast_population_burnup_applicability_assessed') is True
    negligible = fusion.get('fusion_fast_population_burnup_applicability_passed') is True
    evaluated = bool(sink or tolerance_available and assessed)
    passed = bool(sink or tolerance_available and assessed and negligible) if evaluated else None
    if passed is True:
        reason = None
    elif evaluated:
        reason = "fast_ion_fusion_burnup_applicability: fast_ion_fusion_burnup_not_negligible"
    elif tolerance_available:
        reason = "fast_ion_fusion_burnup_applicability: fast_ion_fusion_burnup_not_assessable_against_source_and_terminal_loss"
    else:
        reason = "fast_ion_fusion_burnup_applicability: fast_ion_fusion_sink_not_in_Eq59_and_negligibility_tolerance_unavailable"
    tolerance = fusion.get('fusion_fast_population_burnup_negligibility_tolerance')
    maximum_fraction = fusion.get('fusion_fast_population_burnup_maximum_relative_fraction')
    value = {"fast_D_burnup_rate_s": burnup_D, "fast_T_burnup_rate_s": burnup_T, "Eq59_fusion_sink_active": sink, "negligibility_tolerance_available": tolerance_available, "negligibility_assessed": assessed, "negligibility_tolerance": tolerance, "maximum_relative_fraction": maximum_fraction, "negligibility_passed": negligible}
  
    return AssessmentRecord("fast_ion_fusion_burnup_applicability", passed, value=value, tolerance=tolerance, reason=reason)

def _nbi_supported_stationary_closure(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """
    Assess the final population authority for the NBI supported stationary closure
    
    The check requires kinetic only final populations, zero startup seed authority, stationary particle balance, final consistency, self consistent electron temperature, active charge exchange target sinks, and required fast cross species collision fields
    """
    if str(evidence.get("fusion_plasma_closure_model", "")).strip().lower() != "nbi_supported_stationary":
        return AssessmentRecord("nbi_supported_stationary_closure", None, reason="not applicable")
 
    active_species = tuple(evidence.get("active_fast_species", ()) or ())
    collision_required = "deuterium" in active_species and "tritium" in active_species
    seed_fraction = evidence.get("startup_seed_final_authority_fraction")
    try:
        seed_authority_zero = bool(np.isfinite(float(seed_fraction)) and float(seed_fraction) == 0.0)
    except (TypeError, ValueError):
        seed_authority_zero = False
    balances = evidence.get("nbi_supported_particle_balance_by_species", {})
    cx_activity = evidence.get("fusion_charge_exchange_target_sink_active_by_species", {})
    balances = balances if isinstance(balances, Mapping) else {}
    cx_activity = cx_activity if isinstance(cx_activity, Mapping) else {}
    cx_sink_required = tuple(species_id for species_id, record in balances.items() if isinstance(record, Mapping) and float(record.get("charge_exchange_target_sink_rate_s", 0.0) or 0.0) > 0.0)
    cx_sink_active = all(cx_activity.get(species_id) is True for species_id in cx_sink_required)
    checks = {
        "kinetic_populations_only": evidence.get("fusion_nbi_supported_kinetic_population_only") is True,
        "startup_seed_authority_zero": seed_authority_zero,
        "stationary_particle_balance": evidence.get("nbi_supported_particle_balance_converged") is True,
        "final_consistency_recompute": evidence.get("final_consistency_recompute_passed") is True,
        "self_consistent_electron_temperature": str(evidence.get("electron_temperature_mode", "")).strip().lower() == "self_consistent_electron_energy",
        "charge_exchange_target_sink": cx_sink_active,
        "reduced_fast_cross_collisions": not collision_required or evidence.get("fast_D_fast_T_cross_collisions_available") is True,
    }
    passed = all(checks.values())
    failed = tuple(name for name, value in checks.items() if not value)
  
    return AssessmentRecord("nbi_supported_stationary_closure", passed, value=checks, tolerance=True, reason=None if passed else "; ".join(failed))

def _bosch_hale(evidence: Mapping[str, Any]) -> AssessmentRecord:
    """Assess active fusion rate coverage of the Bosch and Hale Table IV center of mass energy fit intervals"""
    if evidence.get("fusion_bosch_hale_fit_domain_coverage_required", True) is not True:
        return AssessmentRecord("bosch_hale_domain_coverage", None, reason="not applicable")
  
    assessed = evidence.get("fusion_bosch_hale_fit_domain_coverage_assessed") is True
    passed = evidence.get("fusion_bosch_hale_fit_domain_coverage_passed") is True if assessed else None
    value = {"cross_section_fit_domains_E_cm_keV": evidence.get("fusion_bosch_hale_cross_section_fit_domains_E_cm_keV"), "active_channels": evidence.get("active_fusion_channels"), "coverage_assessment_method": evidence.get("fusion_bosch_hale_fit_domain_coverage_assessment_method"), "rate_weighted_out_of_domain_fraction": evidence.get("fusion_bosch_hale_rate_weighted_out_of_fit_domain_fraction"), "rate_weighted_out_of_domain_fraction_tolerance": evidence.get("fusion_bosch_hale_rate_weighted_out_of_fit_domain_fraction_tolerance"), "authoritative_tolerance_available": evidence.get("fusion_bosch_hale_fit_domain_authoritative_tolerance_available")}
    tolerance = {"rate_weighted_out_of_domain_fraction": evidence.get("fusion_bosch_hale_rate_weighted_out_of_fit_domain_fraction_tolerance")}
    reason = None if passed is True else "bosch_hale_fit_domain_coverage: bosch_hale_fit_domain_coverage_not_assessed" if not assessed else "bosch_hale_fit_domain_coverage: bosch_hale_fit_domain_coverage_failed"
  
    return AssessmentRecord("bosch_hale_domain_coverage", passed, value=value, tolerance=tolerance, reason=reason)

def kinetic_numerical_convergence(metadata: Mapping[str, Any]) -> tuple[bool, tuple[str, ...]]:
    """Return the aggregate inner kinetic convergence state and names of the applicable kinetic checks"""
    records = (_source_projection(metadata), _eigenbasis(metadata), _eq42(metadata), _collision_convergence(metadata), _eq59(metadata), _eq70(metadata), _local_reconstruction(metadata))
    applicable = tuple(record for record in records if record.passed is not None)
 
    return all(record.passed is True for record in applicable), tuple(record.name for record in applicable)

def build_run_assessment(*, operating_point: object, kinetic: object, expander: object, fusion: object) -> RunAssessment:
    """
    Build the final RunAssessment from operating point, kinetic, expander, and fusion metadata
    
    Stage evidence is merged in order and converted into the canonical convergence and applicability records
    """
    stages = (operating_point, kinetic, expander, fusion)
    # Later stage metadata can refine evidence assembled by earlier stages
    evidence: dict[str, Any] = {}
    for stage in stages:
        metadata = getattr(stage, "metadata", None)
        if isinstance(metadata, Mapping):
            evidence.update(metadata)
    expander_value = _literal_boolean(evidence.get("expander_potential_converged"))
    convergence = (_electron_temperature(evidence), _source_projection(evidence), _eigenbasis(evidence), _eq42(evidence), _collision_convergence(evidence), _beam_density(evidence), _eq59(evidence), _eq70(evidence), _local_reconstruction(evidence), AssessmentRecord('expander_potential_convergence', expander_value, value=evidence.get('expander_potential_converged'), reason=None if expander_value is True else 'expander_potential_converged is false' if expander_value is False else 'expander_potential_converged is unavailable'), _terminal_current(evidence))
    expander_measured = _literal_boolean(evidence.get("expander_population_classification_passed"))
    expander_applicable = expander_measured is not None
    fast_boundary = _literal_boolean(evidence.get("expander_fixed_fast_boundary_applicable")) if expander_applicable else None
    local_closed = _literal_boolean(evidence.get("expander_velocity_mapping_complete")) if expander_applicable else None
    applicability = (
        _equivalent_end(evidence),
        _collision_applicability(evidence),
        AssessmentRecord("expander_fast_boundary_applicability", fast_boundary, value={"fast_nonterminal_fraction": evidence.get("expander_fast_nonterminal_fraction"), "fixed_fast_boundary_applicable": evidence.get("expander_fixed_fast_boundary_applicable")}, tolerance=evidence.get("expander_fast_nonterminal_fraction_tolerance"), reason=None if fast_boundary is True else "fixed fast boundary is not applicable" if fast_boundary is False else "not applicable"),
        AssessmentRecord("expander_velocity_mapping", local_closed, value={"velocity_domain_overflow_fraction": evidence.get("expander_velocity_domain_overflow_inventory_fraction_max"), "velocity_mapping_complete": evidence.get("expander_velocity_mapping_complete")}, tolerance=evidence.get("expander_source_disconnected_population_fraction_tolerance"), reason=None if local_closed is True else "expander velocity domain overflow is unresolved" if local_closed is False else "not applicable"),
        fast_ion_fusion_burnup_assessment(getattr(fusion, "metadata", {})),
        _nbi_supported_stationary_closure(evidence),
        _bosch_hale(evidence),
    )
  
    return RunAssessment(convergence, applicability)

__all__ = [
    "APPLICABILITY_CHECKS",
    "AssessmentRecord",
    "CONVERGENCE_CHECKS",
    "RunAssessment",
    "build_run_assessment",
    "fast_ion_fusion_burnup_assessment",
    "kinetic_numerical_convergence",
]