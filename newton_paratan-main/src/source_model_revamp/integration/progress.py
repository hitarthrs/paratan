"""Readable console progress for long operating point searches"""
from __future__ import annotations
from source_model_revamp.integration.evaluation import OperatingPointProgressEvent

def _short_message(message: str | None, limit: int = 180) -> str:
    if not message:
        return ""
    compact = " ".join(str(message).split())
    compact = compact.split("; eq61_diagnostic_json=", 1)[0]
    if len(compact) <= limit:
        return compact
    return compact[: limit - 3].rstrip() + "..."

class OperatingPointConsoleProgressReporter:
    """Print bounded stage level progress events"""

    def __init__(self, *, verbose: bool = False) -> None:
        self._last_trial: int | None = None
        self._verbose = bool(verbose)

    def __call__(self, event: OperatingPointProgressEvent) -> None:
        if event.phase == "temperature_trial_started":
            self._last_trial = event.temperature_trial_index
            print(f"Te trial {event.temperature_trial_index}: " f"Te = {event.temperature_keV:.6g} keV")
            if self._verbose:
                print(f"evaluation mode = {event.evaluation_mode.value}")
        elif event.phase == "temperature_trial_cache_hit":
            print(f"Te cache hit: Te = {event.temperature_keV:.6g} keV")
        elif event.phase == "temperature_trial_invalid":
            reason = _short_message(event.message)
            suffix = "" if not reason else f", reason = {reason}"
            print(f"Te trial invalid: Te = {event.temperature_keV:.6g} keV, " f"runtime = {event.elapsed_stage_s:.3f} s{suffix}")
            if self._verbose and event.message and reason != event.message:
                print(f"full reason = {event.message}")
        elif event.phase == "temperature_trial_completed":
            relative = ("unavailable" if event.current_relative_residual is None else f"{event.current_relative_residual:.6g}")
            passed = ("unavailable" if event.power_balance_converged is None else str(bool(event.power_balance_converged)))
            print(f"Te trial complete: relative residual = {relative}, " f"passed = {passed}, runtime = {event.elapsed_stage_s:.3f} s")
            if self._verbose:
                absolute_tolerance = ("unavailable" if event.residual_tolerance_W is None else f"{event.residual_tolerance_W:.6g} W")
                relative_tolerance = ("unavailable" if event.relative_residual_tolerance is None else f"{event.relative_residual_tolerance:.6g}")
                print(f"residual = {event.current_residual_W:.6g} W, " f"tolerance = |residual| <= {absolute_tolerance} or " f"|relative| <= {relative_tolerance}")
        elif event.phase == "beam_density_iteration_started":
            if self._verbose:
                print(f"beam iteration {event.beam_iteration} of " f"{event.beam_iteration_limit}")
        elif event.phase == "beam_density_iteration_completed":
            residual = ("unavailable" if event.beam_profile_change is None else f"{event.beam_profile_change:.6g}")
            print(f"beam iteration {event.beam_iteration}/{event.beam_iteration_limit}: " f"profile change = {residual}, runtime = {event.elapsed_stage_s:.3f} s")
            if self._verbose and event.message:
                print(f"{event.message}")
        elif event.phase == "kinetic_stage_completed":
            if self._verbose:
                print(f"kinetic stage runtime = {event.elapsed_stage_s:.3f} s, " f"{event.message}")
        elif event.phase == "final_beam_kinetic_consistency_completed":
            print("final beam and kinetic consistency recompute complete, " f"runtime = {event.elapsed_stage_s:.3f} s")
            if self._verbose and event.message:
                print(f"{event.message}")
        elif event.phase == "final_qualification_started":
            print(f"Final qualification: Te = {event.temperature_keV:.6g} keV")
        elif event.phase == "final_qualification_completed":
            relative = ("unavailable" if event.current_relative_residual is None else f"{event.current_relative_residual:.6g}")
            message = _short_message(event.message)
            status = "complete" if not message else message
            print(f"Final qualification {status}: relative residual = {relative}, " f"runtime = {event.elapsed_stage_s:.3f} s")
            if self._verbose:
                absolute_tolerance = ("unavailable" if event.residual_tolerance_W is None else f"{event.residual_tolerance_W:.6g} W")
                relative_tolerance = ("unavailable" if event.relative_residual_tolerance is None else f"{event.relative_residual_tolerance:.6g}")
                print(f"residual = {event.current_residual_W}, " f"tolerance = |residual| <= {absolute_tolerance} or " f"|relative| <= {relative_tolerance}")

__all__ = ["OperatingPointConsoleProgressReporter"]