"""Operating point evaluation modes and optional progress event types"""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from typing import Callable

class OperatingPointEvaluationMode(str, Enum):
    """
    Numerical qualification policy for one operating point evaluation
    
    Temperature trials can skip final only refinement studies while fixed temperature final and final qualification evaluations run the configured qualification checks
    """
    TEMPERATURE_TRIAL = "temperature_trial"
    FINAL_QUALIFICATION = "final_qualification"
    FIXED_TEMPERATURE_FINAL = "fixed_temperature_final"

    @classmethod
    def parse(cls, value: "OperatingPointEvaluationMode | str") -> "OperatingPointEvaluationMode":
        """Return a validated evaluation mode from an enum value or string"""
        if isinstance(value, cls):
            return value
        try:
            return cls(str(value))
        except ValueError as exc:
            choices = ", ".join(item.value for item in cls)
            raise ValueError(f"operating point evaluation mode must be one of {choices}") from exc

    @property
    def runs_final_qualification(self) -> bool:
        """Return whether final only numerical qualification studies must run"""
        return self is not OperatingPointEvaluationMode.TEMPERATURE_TRIAL

FINAL_ONLY_QUALIFICATION_CHECK_NAMES = ("physical_eigenbasis_grid_convergence", "retained_physical_mode_count_convergence", "eq59_nonlinear_and_speed_domain_convergence", "local_physical_speed_grid_convergence")

@dataclass(frozen=True)
class OperatingPointProgressEvent:
    """
    One optional progress update from the temperature and operating point solvers
    
    All fields are diagnostic and do not alter solver state
    """
    phase: str
    evaluation_mode: OperatingPointEvaluationMode
    temperature_trial_index: int | None = None
    temperature_keV: float | None = None
    warm_start_source_temperature_keV: float | None = None
    beam_iteration: int | None = None
    beam_iteration_limit: int | None = None
    eq42_iteration: int | None = None
    scalar_density_iteration: int | None = None
    elapsed_stage_s: float | None = None
    current_residual_W: float | None = None
    current_relative_residual: float | None = None
    residual_tolerance_W: float | None = None
    relative_residual_tolerance: float | None = None
    power_balance_converged: bool | None = None
    beam_profile_change: float | None = None
    message: str | None = None

OperatingPointProgressCallback = Callable[[OperatingPointProgressEvent], None]

def emit_progress(callback: OperatingPointProgressCallback | None, event: OperatingPointProgressEvent) -> None:
    """Send one progress event when a callback is configured and otherwise remain silent"""
    if callback is not None:
        callback(event)

__all__ = [
    "FINAL_ONLY_QUALIFICATION_CHECK_NAMES",
    "OperatingPointEvaluationMode",
    "OperatingPointProgressCallback",
    "OperatingPointProgressEvent",
    "emit_progress",
]