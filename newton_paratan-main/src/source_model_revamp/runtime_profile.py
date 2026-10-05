"""Lightweight wall clock runtime profiling helpers"""
from __future__ import annotations
from collections.abc import Mapping, MutableMapping
from typing import Any

def new_runtime_profile() -> dict[str, Any]:
    """Return one mutable runtime accumulator"""
    return {"top_level_s": {}, "top_level_calls": {}, "nested_s": {}, "nested_calls": {}}

def record_runtime(profile: MutableMapping[str, Any] | None, key: str, elapsed_s: float, *, nested: bool = False) -> None:
    """Accumulate one measured wall clock interval"""
    if profile is None:
        return
    seconds = max(0.0, float(elapsed_s))
    values_key = "nested_s" if nested else "top_level_s"
    calls_key = "nested_calls" if nested else "top_level_calls"
    values = profile.setdefault(values_key, {})
    calls = profile.setdefault(calls_key, {})
    values[str(key)] = float(values.get(str(key), 0.0)) + seconds
    calls[str(key)] = int(calls.get(str(key), 0)) + 1

def record_eq70_confined_interpolation_audit_call(profile: MutableMapping[str, Any] | None, record: Mapping[str, Any]) -> None:
    """Append one Eq 70 confined interpolation audit record"""
    if profile is None:
        return
    calls = profile.setdefault("eq70_confined_interpolation_audit_calls", [])
    if not isinstance(calls, list):
        raise TypeError("Eq 70 confined interpolation audit call storage must be a list")
    calls.append(dict(record))

def finalize_runtime_profile(profile: Mapping[str, Any] | None, *, total_s: float) -> dict[str, Any]:
    """Return one JSON ready runtime profile"""
    source = {} if profile is None else dict(profile)
    top_level = {str(key): float(value) for key, value in dict(source.get("top_level_s", {})).items()}
    top_calls = {str(key): int(value) for key, value in dict(source.get("top_level_calls", {})).items()}
    nested = {str(key): float(value) for key, value in dict(source.get("nested_s", {})).items()}
    nested_calls = {str(key): int(value) for key, value in dict(source.get("nested_calls", {})).items()}
    eq70_confined_interpolation_audit_calls = tuple(dict(item) for item in source.get("eq70_confined_interpolation_audit_calls", ()) if isinstance(item, Mapping))
    measured_total = max(0.0, float(total_s))
    accounted = float(sum(top_level.values()))

    return {
        "schema": "source_model_wall_clock_profile_v1",
        "total_s": measured_total,
        "top_level_s": top_level,
        "top_level_calls": top_calls,
        "top_level_accounted_s": accounted,
        "top_level_unaccounted_s": max(0.0, measured_total - accounted),
        "nested_s": nested,
        "nested_calls": nested_calls,
        "nested_component_semantics": "nested intervals are contained inside one or more top level intervals and must not be summed with top level timings",
        "eq70_confined_interpolation_audit_schema": "eq70_confined_interpolation_audit_v1" if eq70_confined_interpolation_audit_calls else None,
        "eq70_confined_interpolation_audit_call_count": len(eq70_confined_interpolation_audit_calls),
        "eq70_confined_interpolation_audit_calls": eq70_confined_interpolation_audit_calls,
    }

def compact_runtime_summary(profile: Mapping[str, Any] | None, *, limit: int = 5) -> str | None:
    """Return a compact ranked nested component summary"""
    if not isinstance(profile, Mapping):
        return None
    nested = profile.get("nested_s")
    if not isinstance(nested, Mapping) or not nested:
        return None
    ranked = sorted(((str(key), float(value)) for key, value in nested.items()), key=lambda item: item[1], reverse=True)[: max(int(limit), 1)]

    return ", ".join(f"{key}={value:.3f}s" for key, value in ranked)

__all__ = [
    "compact_runtime_summary",
    "finalize_runtime_profile",
    "new_runtime_profile",
    "record_eq70_confined_interpolation_audit_call",
    "record_runtime",
]
