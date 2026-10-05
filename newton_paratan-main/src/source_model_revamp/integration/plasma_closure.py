"""Helpers that interpret the role of the configured plasma closure"""
from __future__ import annotations
from typing import Any
from source_model_revamp.integration.config.constants import NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL

def plasma_closure_model(value: Any) -> str:
    """Return the normalized closure model string from a configuration object or raw value"""
    model = getattr(value, "model", value)
    return str(model).strip().lower()

def uses_startup_seed_only(value: Any) -> bool:
    """Return whether configured D and T profiles are initialization data rather than final density authority"""
    return plasma_closure_model(value) == NBI_SUPPORTED_STATIONARY_PLASMA_CLOSURE_MODEL