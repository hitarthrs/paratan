"""Check configuration support for the self consistent electron temperature closure"""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from typing import Any
import numpy as np
from source_model_revamp.fbis.modal.lost_ion_distribution import BALDWIN_1972_THROAT_DENSITY_CLOSURE, EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY

class ElectronTemperaturePreflightDisposition(str, Enum):
    """Configuration level disposition for the self consistent temperature solve"""
    APPLICABLE = "statically_applicable"
    UNAVAILABLE = "statically_unavailable"
    REQUIRES_OPERATING_POINT = "requires_one_operating_point_evaluation"

@dataclass(frozen=True)
class ElectronTemperaturePreflightResult:
    """Configuration checks that determine whether Te trials may be attempted"""
    disposition: ElectronTemperaturePreflightDisposition
    blocking_gate_ids: tuple[str, ...] = ()
    required_operating_point_gate_ids: tuple[str, ...] = ()
    evaluated: bool = True

    @property
    def passed(self) -> bool:
        """Return whether no configuration level check blocks the temperature solve"""
        return self.disposition is not ElectronTemperaturePreflightDisposition.UNAVAILABLE

    def as_metadata(self) -> dict[str, object]:
        """Return stable preflight metadata"""
        return {
            "electron_temperature_preflight_evaluated": bool(self.evaluated),
            "electron_temperature_preflight_passed": self.passed,
            "electron_temperature_preflight_disposition": self.disposition.value,
            "electron_temperature_preflight_blocking_gate_ids": list(self.blocking_gate_ids),
            "electron_temperature_preflight_required_operating_point_gate_ids": list(self.required_operating_point_gate_ids),
        }

def evaluate_self_consistent_temperature_preflight(config: Any) -> ElectronTemperaturePreflightResult:
    """Separate unavailable configuration choices from checks requiring one solved operating point"""
    unavailable: list[str] = []
    requires: list[str] = []
    kinetic = getattr(config, "kinetic_electrostatic", None)
    closure = getattr(config, "plasma_closure", None)
    power = getattr(config, "power_balance", None)
    if kinetic is None or closure is None or power is None:
        return ElectronTemperaturePreflightResult(
            disposition=ElectronTemperaturePreflightDisposition.REQUIRES_OPERATING_POINT,
            required_operating_point_gate_ids=("operating_point_density_model",),
        )
    density_model = str(getattr(kinetic, "modal_density_closure_model", "")).strip().lower()
    nbi_supported_density_model = "nbi_supported_stationary"
    supported_density_models = {nbi_supported_density_model}
    if density_model not in supported_density_models:
        unavailable.append("operating_point_density_model")
    loss_model = str(getattr(kinetic, "ion_loss_closure_model", "")).strip().lower()
    lost_temperature_model = str(getattr(kinetic, "lost_ion_parallel_temperature_model", "")).strip().lower()
    if loss_model != EGEDAL_HOT_FIXED_MAGNETIC_BOUNDARY or lost_temperature_model != BALDWIN_1972_THROAT_DENSITY_CLOSURE:
        unavailable.append("electron_and_ion_end_loss_power_terms")
    if str(getattr(power, "electron_temperature_mode", "")).strip().lower() != "self_consistent_electron_energy":
        unavailable.append("electron_temperature_power_balance_closure")
    deuterium = float(getattr(closure, "background_deuterium_midplane_density_m3", 0.0))
    tritium = float(getattr(closure, "background_tritium_midplane_density_m3", 0.0))
    if not np.isfinite(deuterium) or not np.isfinite(tritium):
        unavailable.append("operating_point_density_model")
    # These checks require collision and species states from one solved operating point
    if density_model == nbi_supported_density_model:
        requires.extend(('multispecies_beta_m_applicability', 'ion_fast_coulomb_log_applicability', 'fast_ion_collision_decomposition_applicability'))
    unavailable = list(dict.fromkeys(unavailable))
    requires = [check_name for check_name in dict.fromkeys(requires) if check_name not in unavailable]
    if unavailable:
        disposition = ElectronTemperaturePreflightDisposition.UNAVAILABLE
    elif requires:
        disposition = ElectronTemperaturePreflightDisposition.REQUIRES_OPERATING_POINT
    else:
        disposition = ElectronTemperaturePreflightDisposition.APPLICABLE

    return ElectronTemperaturePreflightResult(
        disposition=disposition,
        blocking_gate_ids=tuple(unavailable),
        required_operating_point_gate_ids=tuple(requires),
    )

__all__ = ["ElectronTemperaturePreflightDisposition", "ElectronTemperaturePreflightResult", "evaluate_self_consistent_temperature_preflight"]
