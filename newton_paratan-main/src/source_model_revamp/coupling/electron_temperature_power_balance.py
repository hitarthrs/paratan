"""Solve the self consistent electron temperature power balance

Each temperature evaluation rebuilds an operating point with temperature dependent collisions electrostatics and end loss powers"""
from __future__ import annotations
from dataclasses import dataclass
from math import isfinite
from time import perf_counter
from typing import Any, Callable
import numpy as np
from scipy.optimize import brentq
from source_model_revamp.fbis.modal.types import ModalPhysicalDistributionError
from source_model_revamp.integration.evaluation import OperatingPointEvaluationMode, OperatingPointProgressCallback, OperatingPointProgressEvent, emit_progress

class ElectronTemperatureStaticApplicabilityError(ValueError):
    """Signal that temperature independent physics prevents residual admission"""
    def __init__(self, message: str, *, blocking_check_names: tuple[str, ...]) -> None:
        """Store the blocking applicability check identifiers"""
        super().__init__(message)
        self.blocking_check_names = tuple(blocking_check_names)

class ElectronTemperatureTrialRejectedError(RuntimeError):
    """Signal that one temperature trial cannot enter the scalar root search"""
    def __init__(self, message: str, *, rejection_kind: str = "trial_rejected", blocking_check_names: tuple[str, ...] = ()) -> None:
        """Store the trial rejection class and blocking check identifiers"""
        super().__init__(message)
        self.rejection_kind = str(rejection_kind)
        self.blocking_check_names = tuple(str(item) for item in blocking_check_names)

class _InvalidRootTrial(RuntimeError):
    """Carry one rejected temperature encountered inside a valid Brent bracket"""
    def __init__(self, electron_temperature_keV: float, error: str) -> None:
        """Store the rejected temperature and evaluator error"""
        super().__init__(error)
        self.electron_temperature_keV = float(electron_temperature_keV)
        self.error = str(error)

class _ConvergedRootTrial(RuntimeError):
    """Stop Brent immediately when a trial already satisfies the power residual tolerance"""
    def __init__(self, evaluation: "ElectronTemperatureTrialEvaluation") -> None:
        """Store the converged trial evaluation"""
        super().__init__("electron temperature trial satisfied the configured power balance tolerance")
        self.evaluation = evaluation

@dataclass(frozen=True)
class ElectronTemperatureSolveConfig:
    """Numerical controls for the geometric bracket scan and scalar Te root"""
    bracket_min_keV: float = 0.5
    bracket_max_keV: float = 100.0
    bracket_scan_points: int = 9
    max_iterations: int = 30
    residual_tolerance_W: float = 1.0e-6
    relative_residual_tolerance: float = 1.0e-3
    root_temperature_relative_tolerance: float = 1.0e-4
    initial_trial_temperature_keV: float | None = None

@dataclass(frozen=True)
class ElectronTemperatureTrialEvaluation:
    """One admitted Te trial with balance terms and optional operating point state"""
    electron_temperature_keV: float
    terms: Any
    state: Any = None
    evaluation_mode: OperatingPointEvaluationMode = (OperatingPointEvaluationMode.TEMPERATURE_TRIAL)
    runtime_s: float | None = None

    @property
    def residual_W(self) -> float:
        """Return the signed heating minus loss residual in W"""
        return float(self.terms.residual_W)

    @property
    def relative_residual(self) -> float:
        """Return the residual normalized by the represented power scale"""
        return float(self.terms.relative_residual)

@dataclass(frozen=True)
class ElectronTemperatureSolveResult:
    """Result of the admitted trial scan and optional Brent root solve"""
    converged: bool
    electron_temperature_keV: float | None
    residual_W: float | None
    relative_residual: float | None
    evaluations: tuple[ElectronTemperatureTrialEvaluation, ...]
    final_evaluation: ElectronTemperatureTrialEvaluation | None
    iterations: int
    failure_reason: str | None = None
    failed_trials: tuple[tuple[float, str], ...] = ()
    search_method: str = "scipy_brentq"
    valid_bracket_keV: tuple[float, float] | None = None

@dataclass(frozen=True)
class ElectronTemperatureFinalQualificationComparison:
    """Compare the selected trial residual with the rebuilt final operating point"""
    electron_temperature_keV: float
    trial_residual_W: float
    final_residual_W: float | None
    absolute_difference_W: float | None
    relative_difference: float | None
    final_residual_within_tolerance: bool
    residual_consistent_with_trial: bool
    required_final_evidence_passed: bool
    operating_point_numerically_converged: bool
    relative_difference_normalization_power_W: float | None = None

    @property
    def passed(self) -> bool:
        """Return whether every final temperature qualification condition passed"""
        return bool(self.final_residual_within_tolerance and self.residual_consistent_with_trial and self.required_final_evidence_passed and self.operating_point_numerically_converged)

    def as_metadata(self) -> dict[str, Any]:
        """Return stable final temperature qualification metadata"""
        return {
            "electron_temperature_trial_candidate_residual_W": self.trial_residual_W,
            "electron_temperature_final_qualified_residual_W": self.final_residual_W,
            "electron_temperature_final_residual_difference_W": self.absolute_difference_W,
            "electron_temperature_final_residual_difference_relative": self.relative_difference,
            "electron_temperature_final_residual_difference_normalization_power_W": self.relative_difference_normalization_power_W,
            "electron_temperature_final_residual_difference_normalization_model": "maximum_absolute_trial_and_final_electron_heating_or_loss_power",
            "electron_temperature_final_residual_within_tolerance": self.final_residual_within_tolerance,
            "electron_temperature_final_residual_consistent_with_trial": self.residual_consistent_with_trial,
            "electron_temperature_final_operating_point_numerically_converged": self.operating_point_numerically_converged,
            "electron_temperature_final_qualification_passed": self.passed,
        }

def _require_positive_finite(value: float, name: str) -> float:
    """Return one validated positive finite scalar"""
    scalar = float(value)
    if not np.isfinite(scalar) or scalar <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return scalar

def _require_nonnegative_finite(value: float, name: str) -> float:
    """Return one validated nonnegative finite scalar"""
    scalar = float(value)
    if not np.isfinite(scalar) or scalar < 0.0:
        raise ValueError(f"{name} must be nonnegative and finite")
    return scalar

def _residual_converged(evaluation: ElectronTemperatureTrialEvaluation, config: ElectronTemperatureSolveConfig) -> bool:
    """Test the absolute or relative represented power residual tolerance"""
    return bool(abs(evaluation.residual_W) <= float(config.residual_tolerance_W) or abs(evaluation.relative_residual) <= float(config.relative_residual_tolerance))

def _evaluate_checked(evaluator: Callable[[float], ElectronTemperatureTrialEvaluation], electron_temperature_keV: float) -> ElectronTemperatureTrialEvaluation:
    """Evaluate one Te trial and validate its returned temperature and residuals"""
    evaluation = evaluator(float(electron_temperature_keV))
    if not isinstance(evaluation, ElectronTemperatureTrialEvaluation):
        raise TypeError("electron temperature evaluator must return ElectronTemperatureTrialEvaluation")
    if not isfinite(evaluation.residual_W) or not isfinite(evaluation.relative_residual):
        raise ValueError("electron temperature evaluator returned nonfinite residuals")
    if not np.isclose(evaluation.electron_temperature_keV, electron_temperature_keV, rtol=1.0e-12, atol=0.0):
        raise ValueError("electron temperature evaluator returned the wrong temperature")
    return evaluation

def _best_evaluation(evaluations: tuple[ElectronTemperatureTrialEvaluation, ...]) -> ElectronTemperatureTrialEvaluation | None:
    """Select the valid trial with the smallest relative then absolute residual"""
    if not evaluations:
        return None
    return min(evaluations, key=lambda item: (abs(item.relative_residual), abs(item.residual_W)))

def solve_electron_temperature_power_balance(evaluator: Callable[[float], ElectronTemperatureTrialEvaluation], config: ElectronTemperatureSolveConfig, *, progress_callback: OperatingPointProgressCallback | None = None, clock: Callable[[], float] = perf_counter) -> ElectronTemperatureSolveResult:
    """Solve R(Te) = 0 using admitted operating points geometric scanning and Brent bracketing"""
    lo = _require_positive_finite(config.bracket_min_keV, "bracket_min_keV")
    hi = _require_positive_finite(config.bracket_max_keV, "bracket_max_keV")
    if hi <= lo:
        raise ValueError("bracket_max_keV must be greater than bracket_min_keV")
    scan_points = int(config.bracket_scan_points)
    if scan_points < 2:
        raise ValueError("bracket_scan_points must be at least two")
    _require_nonnegative_finite(config.residual_tolerance_W, "residual_tolerance_W")
    _require_nonnegative_finite(config.relative_residual_tolerance, "relative_residual_tolerance")
    root_tolerance = _require_positive_finite(config.root_temperature_relative_tolerance, "root_temperature_relative_tolerance")
    initial_trial_temperature_keV = config.initial_trial_temperature_keV
    if initial_trial_temperature_keV is not None:
        initial_trial_temperature_keV = _require_positive_finite(
            initial_trial_temperature_keV,
            "initial_trial_temperature_keV",
        )
        if initial_trial_temperature_keV < lo or initial_trial_temperature_keV > hi:
            raise ValueError(
                "initial_trial_temperature_keV must lie within the configured bracket"
            )
    mode = OperatingPointEvaluationMode.TEMPERATURE_TRIAL
    # Exact binary temperature keys prevent repeated operating point rebuilds
    cache: dict[str, ElectronTemperatureTrialEvaluation | str] = {}
    ordered_valid: list[ElectronTemperatureTrialEvaluation] = []
    failed_trials: list[tuple[float, str]] = []
    trial_request_index = 0
    terminal_check_names: tuple[str, ...] = ()

    def result(*, converged: bool, final: ElectronTemperatureTrialEvaluation | None, failure_reason: str | None, bracket: tuple[float, float] | None = None, iterations: int = 0) -> ElectronTemperatureSolveResult:
        """Build one immutable scalar temperature solve result"""
        return ElectronTemperatureSolveResult(
            converged=converged,
            electron_temperature_keV=(None if final is None else float(final.electron_temperature_keV)),
            residual_W=None if final is None else float(final.residual_W),
            relative_residual=(None if final is None else float(final.relative_residual)),
            evaluations=tuple(ordered_valid),
            final_evaluation=final,
            iterations=int(iterations),
            failure_reason=failure_reason,
            failed_trials=tuple(failed_trials),
            valid_bracket_keV=bracket,
        )

    def evaluate(temperature_keV: float) -> ElectronTemperatureTrialEvaluation | None:
        """Evaluate or reuse one exact Te trial and record admission status"""
        nonlocal trial_request_index, terminal_check_names
        temperature = _require_positive_finite(temperature_keV, "electron_temperature_keV")
        trial_request_index += 1
        key = float(temperature).hex()
        cached = cache.get(key)
        if cached is not None:
            emit_progress(progress_callback, OperatingPointProgressEvent(phase="temperature_trial_cache_hit", evaluation_mode=mode, temperature_trial_index=trial_request_index, temperature_keV=temperature),)
            return cached if isinstance(cached, ElectronTemperatureTrialEvaluation) else None
        emit_progress(progress_callback, OperatingPointProgressEvent(phase="temperature_trial_started", evaluation_mode=mode, temperature_trial_index=trial_request_index, temperature_keV=temperature),)
        trial_started = clock()
        try:
            candidate = _evaluate_checked(evaluator, temperature)
        except ElectronTemperatureStaticApplicabilityError as exc:
            terminal_check_names = tuple(exc.blocking_check_names)
            error = f"{type(exc).__name__}: {exc}"
            cache[key] = error
            failed_trials.append((temperature, error))
            emit_progress(progress_callback, OperatingPointProgressEvent(phase="temperature_trial_invalid", evaluation_mode=mode, temperature_trial_index=trial_request_index, temperature_keV=temperature, elapsed_stage_s=max(0.0, float(clock() - trial_started)), message=error),)
            return None
        except (ModalPhysicalDistributionError, ElectronTemperatureTrialRejectedError) as exc:
            error = f"{type(exc).__name__}: {exc}"
            cache[key] = error
            failed_trials.append((temperature, error))
            emit_progress(progress_callback, OperatingPointProgressEvent(phase="temperature_trial_invalid", evaluation_mode=mode, temperature_trial_index=trial_request_index, temperature_keV=temperature, elapsed_stage_s=max(0.0, float(clock() - trial_started)), message=error),)
            return None
        runtime = max(0.0, float(clock() - trial_started))
        if candidate.runtime_s is None:
            candidate = ElectronTemperatureTrialEvaluation(electron_temperature_keV=candidate.electron_temperature_keV, terms=candidate.terms, state=candidate.state, evaluation_mode=candidate.evaluation_mode, runtime_s=runtime)
        cache[key] = candidate
        ordered_valid.append(candidate)
        emit_progress(
            progress_callback,
            OperatingPointProgressEvent(
                phase="temperature_trial_completed",
                evaluation_mode=mode,
                temperature_trial_index=trial_request_index,
                temperature_keV=temperature,
                elapsed_stage_s=runtime,
                current_residual_W=float(candidate.residual_W),
                current_relative_residual=float(candidate.relative_residual),
                residual_tolerance_W=float(config.residual_tolerance_W),
                relative_residual_tolerance=float(config.relative_residual_tolerance),
                power_balance_converged=_residual_converged(candidate, config),
            ),
        )
        
        return candidate

    if initial_trial_temperature_keV is not None:
        initial_evaluation = evaluate(initial_trial_temperature_keV)
        if terminal_check_names:
            return result(
                converged=False,
                final=_best_evaluation(tuple(ordered_valid)),
                failure_reason="temperature_model_statically_unavailable_after_one_operating_point",
            )
        if initial_evaluation is not None and _residual_converged(initial_evaluation, config):
            return result(
                converged=True,
                final=initial_evaluation,
                failure_reason=None,
                bracket=(initial_trial_temperature_keV, initial_trial_temperature_keV),
            )

    scan_temperatures = tuple(float(value) for value in np.geomspace(lo, hi, scan_points))
    bracket: tuple[ElectronTemperatureTrialEvaluation, ElectronTemperatureTrialEvaluation] | None = None
    previous_scan_temperature: float | None = None
    previous_scan_evaluation: ElectronTemperatureTrialEvaluation | None = None
    admitted_scan_evaluation_count = 0

    # Never infer a sign change across a rejected lower temperature interval
    def recover_lower_invalid_gap(
        lower_invalid_temperature_keV: float,
        upper_evaluation: ElectronTemperatureTrialEvaluation,
    ) -> tuple[ElectronTemperatureTrialEvaluation | None, tuple[ElectronTemperatureTrialEvaluation, ElectronTemperatureTrialEvaluation] | None]:
        """Refine an invalid lower scan gap until a valid root or sign bracket is found"""
        lower_temperature = float(lower_invalid_temperature_keV)
        upper = upper_evaluation
        refinement_limit = max(int(config.max_iterations), 1)
        for _ in range(refinement_limit):
            upper_temperature = float(upper.electron_temperature_keV)
            relative_span = (upper_temperature - lower_temperature) / max(upper_temperature, np.finfo(float).tiny)
            if relative_span <= root_tolerance:
                break
            probe_temperature = float(np.sqrt(lower_temperature * upper_temperature))
            if not lower_temperature < probe_temperature < upper_temperature:
                break
            probe = evaluate(probe_temperature)
            if terminal_check_names:
                return None, None
            if probe is None:
                lower_temperature = probe_temperature
                continue
            if _residual_converged(probe, config):
                return probe, None
            if probe.residual_W * upper.residual_W < 0.0:
                return None, (probe, upper)
            upper = probe
        return None, None

    for temperature in scan_temperatures:
        evaluation = evaluate(temperature)
        if terminal_check_names:
            return result(converged=False, final=_best_evaluation(tuple(ordered_valid)), failure_reason="temperature_model_statically_unavailable_after_one_operating_point")
        if evaluation is not None and _residual_converged(evaluation, config):
            return result(converged=True, final=evaluation, failure_reason=None, bracket=(temperature, temperature))
        if evaluation is not None and previous_scan_evaluation is not None and previous_scan_evaluation.residual_W * evaluation.residual_W < 0.0:
            bracket = (previous_scan_evaluation, evaluation)
            break
        if (
            evaluation is not None
            and evaluation.residual_W < 0.0
            and admitted_scan_evaluation_count == 0
            and previous_scan_temperature is not None
            and previous_scan_evaluation is None
        ):
            recovered, recovered_bracket = recover_lower_invalid_gap(previous_scan_temperature, evaluation)
            if terminal_check_names:
                return result(converged=False, final=_best_evaluation(tuple(ordered_valid)), failure_reason="temperature_model_statically_unavailable_after_one_operating_point")
            if recovered is not None:
                return result(converged=True, final=recovered, failure_reason=None, bracket=(recovered.electron_temperature_keV, recovered.electron_temperature_keV))
            if recovered_bracket is not None:
                bracket = recovered_bracket
                break
        if evaluation is not None:
            admitted_scan_evaluation_count += 1
        previous_scan_temperature = float(temperature)
        previous_scan_evaluation = evaluation
    if bracket is None:
        best = _best_evaluation(tuple(ordered_valid))
        failure = "no_valid_electron_temperature_trials" if best is None else "no_power_balance_sign_change_between_adjacent_valid_scan_temperatures"
        return result(converged=False, final=best, failure_reason=failure)

    left, right = bracket
    bracket_keV = (float(left.electron_temperature_keV), float(right.electron_temperature_keV))

    root_start_evaluation_count = len(ordered_valid)

    def residual(temperature_keV: float) -> float:
        """Return the admitted power residual used by Brent"""
        evaluation = evaluate(float(temperature_keV))
        if terminal_check_names:
            raise ElectronTemperatureStaticApplicabilityError("temperature model became statically unavailable during Brent solve", blocking_check_names=terminal_check_names)
        if evaluation is None:
            key = float(temperature_keV).hex()
            error = cache.get(key)
            raise _InvalidRootTrial(float(temperature_keV), str(error or "invalid_temperature_trial"))
        if _residual_converged(evaluation, config):
            raise _ConvergedRootTrial(evaluation)
        return float(evaluation.residual_W)

    try:
        root_temperature_keV, root_info = brentq(residual, bracket_keV[0], bracket_keV[1], rtol=root_tolerance, maxiter=int(config.max_iterations), full_output=True, disp=False)
    except _ConvergedRootTrial as exc:
        root_evaluations = max(len(ordered_valid) - root_start_evaluation_count, 1)
        return result(converged=True, final=exc.evaluation, failure_reason=None, bracket=bracket_keV, iterations=root_evaluations)
    except ElectronTemperatureStaticApplicabilityError:
        return result(converged=False, final=_best_evaluation(tuple(ordered_valid)), failure_reason="temperature_model_statically_unavailable_after_one_operating_point", bracket=bracket_keV)
    except _InvalidRootTrial as exc:
        return result(converged=False, final=_best_evaluation(tuple(ordered_valid)), failure_reason=("invalid_temperature_trial_inside_power_balance_bracket:" f"Te={exc.electron_temperature_keV:.16g}keV:{exc.error}"), bracket=bracket_keV)
    except ValueError as exc:
        return result(converged=False, final=_best_evaluation(tuple(ordered_valid)), failure_reason=f"scipy_brentq_failed:{exc}", bracket=bracket_keV)

    final = evaluate(float(root_temperature_keV))
    if final is None:
        return result(converged=False, final=_best_evaluation(tuple(ordered_valid)), failure_reason="scipy_brentq_returned_invalid_temperature_trial", bracket=bracket_keV, iterations=int(root_info.iterations))
    if not root_info.converged:
        return result(converged=False, final=final, failure_reason="scipy_brentq_max_iterations_reached", bracket=bracket_keV, iterations=int(root_info.iterations))
    if not _residual_converged(final, config):
        return result(converged=False, final=final, failure_reason="scipy_brentq_temperature_converged_but_power_residual_not_within_tolerance", bracket=bracket_keV, iterations=int(root_info.iterations))
    return result(converged=True, final=final, failure_reason=None, bracket=bracket_keV, iterations=int(root_info.iterations))

__all__ = [
    "ElectronTemperatureSolveConfig",
    "ElectronTemperatureFinalQualificationComparison",
    "ElectronTemperatureSolveResult",
    "ElectronTemperatureStaticApplicabilityError",
    "ElectronTemperatureTrialRejectedError",
    "ElectronTemperatureTrialEvaluation",
    "solve_electron_temperature_power_balance",
]
