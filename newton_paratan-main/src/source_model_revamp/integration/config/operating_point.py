"""Beam attenuation and solved density fixed point numerical configuration"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass, fields
from typing import Any
import numpy as np
from source_model_revamp.integration.config.common import _check_unknown, _integer, _number

@dataclass(frozen=True)
class BeamDensityCouplingConfig:
    """
    Numerical controls for the beam and density fixed point
    
    The convergence tests include electron density profile changes, total beam birth rate, deposited power, axial birth shape, and energy component deposition
    """
    max_iterations: int = 12
    electron_profile_relative_tolerance: float = 1.0e-3
    total_birth_rate_relative_tolerance: float = 1.0e-3
    deposited_power_relative_tolerance: float = 1.0e-3
    axial_birth_profile_relative_tolerance: float = 1.0e-3
    energy_component_deposition_relative_tolerance: float = 1.0e-3
    target_density_relaxation: float = 0.5
    electron_profile_volume_L2_relative_tolerance: float | None = None
    electron_profile_absolute_reference_tolerance: float | None = None
    density_floor_absolute_m3: float = 0.0
    density_floor_reference_fraction: float = 1.0e-12
    birth_profile_floor_absolute_m3_s: float = 0.0
    birth_profile_floor_reference_fraction: float = 1.0e-12
    component_power_floor_absolute_W: float = 0.0
    component_power_floor_reference_fraction: float = 1.0e-12

    def __post_init__(self) -> None:
        """Validate fixed point tolerances, relaxation, iteration count, and absolute or reference scaled floors"""
        if self.max_iterations < 2:
            raise ValueError("beam density coupling requires at least two iterations")
        for name, value in (
            ("electron_profile_relative_tolerance", self.electron_profile_relative_tolerance),
            ("total_birth_rate_relative_tolerance", self.total_birth_rate_relative_tolerance),
            ("deposited_power_relative_tolerance", self.deposited_power_relative_tolerance),
            ("axial_birth_profile_relative_tolerance", self.axial_birth_profile_relative_tolerance),
            ("energy_component_deposition_relative_tolerance", self.energy_component_deposition_relative_tolerance),
            ("electron_profile_volume_L2_relative_tolerance", self.profile_volume_L2_tolerance),
            ("electron_profile_absolute_reference_tolerance", self.profile_absolute_reference_tolerance),
        ):
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if not np.isfinite(self.target_density_relaxation) or not 0.0 < self.target_density_relaxation <= 1.0:
            raise ValueError("target_density_relaxation must lie inside (0, 1]")
        if not np.isfinite(self.density_floor_absolute_m3) or self.density_floor_absolute_m3 < 0.0:
            raise ValueError("density_floor_absolute_m3 must be finite and nonnegative")
        if not np.isfinite(self.density_floor_reference_fraction) or self.density_floor_reference_fraction < 0.0:
            raise ValueError("density_floor_reference_fraction must be finite and nonnegative")
        for name, value in (
            ("birth_profile_floor_absolute_m3_s", self.birth_profile_floor_absolute_m3_s),
            ("birth_profile_floor_reference_fraction", self.birth_profile_floor_reference_fraction),
            ("component_power_floor_absolute_W", self.component_power_floor_absolute_W),
            ("component_power_floor_reference_fraction", self.component_power_floor_reference_fraction),
        ):
            if not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")

    @property
    def profile_volume_L2_tolerance(self) -> float:
        """Return the dedicated volume L2 tolerance or fall back to the pointwise electron profile tolerance"""
        if self.electron_profile_volume_L2_relative_tolerance is None:
            return float(self.electron_profile_relative_tolerance)
        return float(self.electron_profile_volume_L2_relative_tolerance)

    @property
    def profile_absolute_reference_tolerance(self) -> float:
        """Return the dedicated absolute reference tolerance or fall back to the pointwise electron profile tolerance"""
        if self.electron_profile_absolute_reference_tolerance is None:
            return float(self.electron_profile_relative_tolerance)
        return float(self.electron_profile_absolute_reference_tolerance)

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, strict: bool = False) -> "BeamDensityCouplingConfig":
        """Parse the beam density fixed point controls and optional profile metric tolerances"""
        allowed = {item.name for item in fields(cls)}
        _check_unknown(data, allowed, "source_model.beam_density_coupling", strict)
        optional = {"electron_profile_volume_L2_relative_tolerance", "electron_profile_absolute_reference_tolerance"}
        missing = sorted((allowed - optional) - set(data))
        if missing:
            raise ValueError("source_model.beam_density_coupling is missing required fields: " + ", ".join(missing))
        volume_tolerance = data.get("electron_profile_volume_L2_relative_tolerance")
        absolute_tolerance = data.get("electron_profile_absolute_reference_tolerance")

        return cls(
            max_iterations=_integer(data.get("max_iterations"), "beam_density_coupling.max_iterations", minimum=2),
            electron_profile_relative_tolerance=_number(data.get("electron_profile_relative_tolerance"), "beam_density_coupling.electron_profile_relative_tolerance", nonnegative=True),
            total_birth_rate_relative_tolerance=_number(data.get("total_birth_rate_relative_tolerance"), "beam_density_coupling.total_birth_rate_relative_tolerance", nonnegative=True),
            deposited_power_relative_tolerance=_number(data.get("deposited_power_relative_tolerance"), "beam_density_coupling.deposited_power_relative_tolerance", nonnegative=True),
            axial_birth_profile_relative_tolerance=_number(data.get("axial_birth_profile_relative_tolerance"), "beam_density_coupling.axial_birth_profile_relative_tolerance", nonnegative=True),
            energy_component_deposition_relative_tolerance=_number(data.get("energy_component_deposition_relative_tolerance"), "beam_density_coupling.energy_component_deposition_relative_tolerance", nonnegative=True),
            target_density_relaxation=_number(data.get("target_density_relaxation"), "beam_density_coupling.target_density_relaxation", positive=True),
            electron_profile_volume_L2_relative_tolerance=(None if volume_tolerance is None else _number(volume_tolerance, "beam_density_coupling.electron_profile_volume_L2_relative_tolerance", nonnegative=True)),
            electron_profile_absolute_reference_tolerance=(None if absolute_tolerance is None else _number(absolute_tolerance, "beam_density_coupling.electron_profile_absolute_reference_tolerance", nonnegative=True)),
            density_floor_absolute_m3=_number(data.get("density_floor_absolute_m3"), "beam_density_coupling.density_floor_absolute_m3", nonnegative=True),
            density_floor_reference_fraction=_number(data.get("density_floor_reference_fraction"), "beam_density_coupling.density_floor_reference_fraction", nonnegative=True),
            birth_profile_floor_absolute_m3_s=_number(data.get("birth_profile_floor_absolute_m3_s"), "beam_density_coupling.birth_profile_floor_absolute_m3_s", nonnegative=True),
            birth_profile_floor_reference_fraction=_number(data.get("birth_profile_floor_reference_fraction"), "beam_density_coupling.birth_profile_floor_reference_fraction", nonnegative=True),
            component_power_floor_absolute_W=_number(data.get("component_power_floor_absolute_W"), "beam_density_coupling.component_power_floor_absolute_W", nonnegative=True),
            component_power_floor_reference_fraction=_number(data.get("component_power_floor_reference_fraction"), "beam_density_coupling.component_power_floor_reference_fraction", nonnegative=True),
        )

    def as_metadata(self) -> dict[str, object]:
        """Return the fixed point controls and resolved fallback tolerances used by the operating point stage"""
        return {
            "max_iterations": self.max_iterations,
            "electron_profile_relative_tolerance": self.electron_profile_relative_tolerance,
            "electron_profile_pointwise_relative_tolerance": self.electron_profile_relative_tolerance,
            "electron_profile_volume_L2_relative_tolerance": self.profile_volume_L2_tolerance,
            "electron_profile_absolute_reference_tolerance": self.profile_absolute_reference_tolerance,
            "density_floor_absolute_m3": self.density_floor_absolute_m3,
            "density_floor_reference_fraction": self.density_floor_reference_fraction,
            "density_floor_rule": "maximum_of_named_absolute_floor_and_reference_scaled_floor",
            "birth_profile_floor_absolute_m3_s": self.birth_profile_floor_absolute_m3_s,
            "birth_profile_floor_reference_fraction": self.birth_profile_floor_reference_fraction,
            "birth_profile_floor_rule": "maximum_of_named_absolute_floor_and_reference_scaled_floor",
            "component_power_floor_absolute_W": self.component_power_floor_absolute_W,
            "component_power_floor_reference_fraction": self.component_power_floor_reference_fraction,
            "component_power_floor_rule": "maximum_of_named_absolute_floor_and_total_power_reference_scaled_floor",
            "total_birth_rate_relative_tolerance": self.total_birth_rate_relative_tolerance,
            "deposited_power_relative_tolerance": self.deposited_power_relative_tolerance,
            "axial_birth_profile_relative_tolerance": self.axial_birth_profile_relative_tolerance,
            "energy_component_deposition_relative_tolerance": self.energy_component_deposition_relative_tolerance,
            "target_density_relaxation": self.target_density_relaxation,
        }
