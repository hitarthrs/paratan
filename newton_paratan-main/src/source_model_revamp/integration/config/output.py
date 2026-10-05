"""Optional source model reporting configuration"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any
from source_model_revamp.integration.config.common import _bool, _check_unknown, _mapping

@dataclass(frozen=True)
class OutputConfig:
    """Controls optional run reporting outputs that do not alter the solved physics state"""
    write_iteration_histories: bool = False

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any] | None, *, strict: bool) -> "OutputConfig":
        """Parse optional reporting controls from the source model output block"""
        values = _mapping(data, "source_model.output")
        _check_unknown(values, {"write_iteration_histories"}, "source_model.output", strict)
        return cls(write_iteration_histories=_bool(values.get("write_iteration_histories", False)))
