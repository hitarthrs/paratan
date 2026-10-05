"""Shared convergence and status metadata helpers for modal kinetic stage results"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable

KINETIC_STATUS_METADATA_KEYS = (
    "exact_midplane_quasineutrality_check_passed",
    "modal_source_projection_rate_converged",
    "kinetic_min_distribution_value",
    "kinetic_geometry_boundary_check_passed",
    "kinetic_symmetric_end_model_applicable",
    "density_roles_separated",
    "electron_ion_confined_inventory_check_passed",
    "kinetic_convergence_failure_reason",
)

def ordered_active_failure_reasons(states: Iterable[tuple[str, bool]]) -> tuple[str, ...]:
    """Return active convergence failure reasons once each in their supplied priority order"""
    reasons: list[str] = []
    for reason, active in states:
        value = str(reason).strip()
        if active and value and value not in reasons:
            reasons.append(value)
 
    return tuple(reasons)

def convergence_failure_summary(reasons: Iterable[str], *, converged: bool) -> tuple[str | None, tuple[str, ...], int]:
    """Return the primary convergence failure reason, complete ordered reason tuple, and additional reason count"""
    ordered = tuple(dict.fromkeys(str(reason).strip() for reason in reasons if str(reason).strip()))
    if converged:
        return None, ordered, 0
    if not ordered:
        ordered = ("kinetic_stage_not_converged_without_component_reason",)
    additional = len(ordered) - 1
    summary = ordered[0] if additional == 0 else f"{ordered[0]} (+{additional} additional reasons)"
   
    return summary, ordered, additional

@dataclass(frozen=True)
class KineticMetadataContract:
    """Common convergence and applicability fields emitted by the single species and coupled species modal stage paths"""
    exact_midplane_quasineutrality_check_passed: bool
    modal_source_projection_rate_converged: bool
    modal_source_projection_relative_error: float
    modal_source_projection_relative_tolerance: float
    kinetic_min_distribution_value: float
    kinetic_geometry_boundary_check_passed: bool
    kinetic_symmetric_end_model_applicable: bool
    density_roles_separated: bool
    electron_ion_confined_inventory_check_passed: bool
    kinetic_convergence_status: str
    kinetic_convergence_failure_reasons: tuple[str, ...] = ()

    def as_metadata(self) -> dict[str, object]:
        """Return the stable metadata representation and summarized convergence failure state"""
        converged = self.kinetic_convergence_status == "converged"
        summary, reasons, additional = convergence_failure_summary(self.kinetic_convergence_failure_reasons, converged=converged)
      
        return {
            "exact_midplane_quasineutrality_check_passed": bool(self.exact_midplane_quasineutrality_check_passed),
            "modal_source_projection_rate_converged": bool(self.modal_source_projection_rate_converged),
            "modal_source_projection_relative_error": float(self.modal_source_projection_relative_error),
            "modal_source_projection_relative_tolerance": float(self.modal_source_projection_relative_tolerance),
            "kinetic_min_distribution_value": float(self.kinetic_min_distribution_value),
            "kinetic_geometry_boundary_check_passed": bool(self.kinetic_geometry_boundary_check_passed),
            "kinetic_symmetric_end_model_applicable": bool(self.kinetic_symmetric_end_model_applicable),
            "density_roles_separated": bool(self.density_roles_separated),
            "electron_ion_confined_inventory_check_passed": bool(self.electron_ion_confined_inventory_check_passed),
            "kinetic_convergence_status": self.kinetic_convergence_status,
            "kinetic_convergence_failure_reason": summary,
            "kinetic_convergence_failure_reasons": reasons,
            "kinetic_convergence_additional_failure_count": additional,
        }
