"""Electron temperature closure and represented electron energy balance configuration"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass, fields
from typing import Any
from source_model_revamp.integration.config.common import canonical_model_name, _check_unknown, _number, _integer
from source_model_revamp.integration.config.constants import *

@dataclass(frozen=True)
class PowerBalanceConfig:
    """
    Electron temperature closure and scalar energy balance numerical controls
    
    The bracket temperatures are in keV and power residual tolerances are in W
    """
    electron_temperature_mode: str = "fixed_closure"
    external_electron_heating_W: float = 0.0
    bracket_min_keV: float = 0.5
    bracket_max_keV: float = 200.0
    bracket_scan_points: int = 11
    max_iterations: int = 40
    residual_tolerance_W: float = 1.0e-3
    relative_residual_tolerance: float = 1.0e-4
    root_temperature_relative_tolerance: float = 1.0e-4
    temperature_search_initial_trial_keV: float | None = None

    def __post_init__(self) -> None:
        """Require an increasing temperature bracket and keep any initial trial inside that bracket"""
        if self.bracket_max_keV <= self.bracket_min_keV:
            raise ValueError("power_balance.bracket_max_keV must be greater than bracket_min_keV")
        if self.temperature_search_initial_trial_keV is not None:
            if not self.bracket_min_keV <= self.temperature_search_initial_trial_keV <= self.bracket_max_keV:
                raise ValueError("power_balance.temperature_search_initial_trial_keV must lie within the configured bracket")

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, strict: bool = False) -> "PowerBalanceConfig":
        """Parse electron temperature mode, external electron heating, root bracket, and convergence controls"""
        allowed = {item.name for item in fields(cls)}
        _check_unknown(data, allowed, "source_model.power_balance", strict)
        optional = {
            "electron_temperature_mode",
            "external_electron_heating_W",
            "temperature_search_initial_trial_keV",
        }
        missing = sorted((allowed - optional) - set(data))
        if missing:
            raise ValueError("source_model.power_balance is missing required fields: " + ", ".join(missing))
        mode = canonical_model_name(data.get("electron_temperature_mode"), POWER_BALANCE_MODE_ALIASES, "fixed_closure")

        return cls(
            electron_temperature_mode=mode,
            external_electron_heating_W=_number(data.get("external_electron_heating_W"), "power_balance.external_electron_heating_W", nonnegative=True, default=0.0),
            bracket_min_keV=_number(data.get("bracket_min_keV"), "power_balance.bracket_min_keV", positive=True),
            bracket_max_keV=_number(data.get("bracket_max_keV"), "power_balance.bracket_max_keV", positive=True),
            bracket_scan_points=_integer(data.get("bracket_scan_points"), "power_balance.bracket_scan_points", minimum=2),
            max_iterations=_integer(data.get("max_iterations"), "power_balance.max_iterations", minimum=1),
            residual_tolerance_W=_number(data.get("residual_tolerance_W"), "power_balance.residual_tolerance_W", nonnegative=True),
            relative_residual_tolerance=_number(data.get("relative_residual_tolerance"), "power_balance.relative_residual_tolerance", nonnegative=True),
            root_temperature_relative_tolerance=_number(data.get("root_temperature_relative_tolerance"), "power_balance.root_temperature_relative_tolerance", positive=True),
            temperature_search_initial_trial_keV=(None if data.get("temperature_search_initial_trial_keV") is None else _number(data.get("temperature_search_initial_trial_keV"), "power_balance.temperature_search_initial_trial_keV", positive=True)),
        )

    def solve_numerics_metadata(self) -> dict[str, Any]:
        """Return the resolved scalar electron temperature solve controls"""
        metadata = {
            "bracket_min_keV": self.bracket_min_keV,
            "bracket_max_keV": self.bracket_max_keV,
            "bracket_scan_points": self.bracket_scan_points,
            "max_iterations": self.max_iterations,
            "residual_tolerance_W": self.residual_tolerance_W,
            "relative_residual_tolerance": self.relative_residual_tolerance,
            "root_temperature_relative_tolerance": self.root_temperature_relative_tolerance,
            "temperature_search_method": "scipy_brentq",
        }
        if self.temperature_search_initial_trial_keV is not None:
            metadata["temperature_search_initial_trial_keV"] = self.temperature_search_initial_trial_keV
        return metadata
