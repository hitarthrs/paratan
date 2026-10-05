"""Validation of end loss power terms used by the electron temperature energy balance"""
from __future__ import annotations
import numpy as np

def end_loss_power_availability_metadata(metadata: dict[str, object], *, fixed_boundary_selected: bool, current_balance_applicable: bool, current_balance_valid: bool, loss_convention_valid: bool, full_device_lost_population_applicable: bool, full_device_lost_population_valid: bool) -> dict[str, object]:
    """
    Validate whether modal end loss powers are complete and internally consistent
    
    The checks require the fixed magnetic ion loss path, valid Eq 68 current balance, valid total device loss convention, valid full device lost population, and the identities
    `P_ion_wall = P_ion_midplane + P_ion_barrier`
    `P_total_wall = P_ion_wall + P_electron_wall`
    """
    failures: list[str] = []
    model_applicable = bool(fixed_boundary_selected and current_balance_applicable and full_device_lost_population_applicable)
    if not fixed_boundary_selected:
        failures.append("selected_ion_loss_model_is_not_egedal_hot_fixed_magnetic_boundary")
    if metadata.get("fixed_boundary_lost_reconstruction_complete") is not True:
        failures.append("fixed_boundary_lost_reconstruction_not_complete")
    if not current_balance_applicable:
        failures.append("electron_current_balance_not_applicable")
    elif not current_balance_valid:
        failures.append("electron_current_balance_not_valid")
    if not loss_convention_valid:
        failures.append("total_device_loss_convention_not_valid")
    if not full_device_lost_population_applicable:
        failures.append("full_device_lost_population_not_applicable")
    elif not full_device_lost_population_valid:
        failures.append("full_device_lost_population_not_valid")
    power_keys = ("modal_ion_midplane_kinetic_power_loss_W", "modal_ion_barrier_power_loss_W", "modal_ion_wall_power_loss_W", "modal_electron_wall_power_loss_W", "modal_total_wall_power_loss_W")
    powers: dict[str, float] = {}
    for key in power_keys:
        value = metadata.get(key)
        try:
            parsed = float(value)
        except (TypeError, ValueError):
            failures.append(f"{key}_missing_or_not_numeric")
            continue
        if not np.isfinite(parsed):
            failures.append(f"{key}_not_finite")
        elif parsed < 0.0:
            failures.append(f"{key}_negative")
        powers[key] = parsed
    identity_tolerance = 64.0 * np.finfo(float).eps
    ion_identity_relative_error: float | None = None
    total_identity_relative_error: float | None = None
    ion_keys = ("modal_ion_midplane_kinetic_power_loss_W", "modal_ion_barrier_power_loss_W", "modal_ion_wall_power_loss_W",)
    if all(key in powers for key in ion_keys):
        expected_ion = powers[ion_keys[0]] + powers[ion_keys[1]]
        ion_scale = max(abs(expected_ion), abs(powers[ion_keys[2]]), np.finfo(float).tiny,)
        ion_identity_relative_error = (powers[ion_keys[2]] - expected_ion) / ion_scale
        if abs(ion_identity_relative_error) > identity_tolerance:
            failures.append("ion_wall_power_identity_not_satisfied")
    total_keys = ("modal_ion_wall_power_loss_W", "modal_electron_wall_power_loss_W", "modal_total_wall_power_loss_W",)
    if all(key in powers for key in total_keys):
        expected_total = powers[total_keys[0]] + powers[total_keys[1]]
        total_scale = max(abs(expected_total), abs(powers[total_keys[2]]), np.finfo(float).tiny,)
        total_identity_relative_error = (powers[total_keys[2]] - expected_total) / total_scale
        if abs(total_identity_relative_error) > identity_tolerance:
            failures.append("total_wall_power_identity_not_satisfied")
    failures = list(dict.fromkeys(failures))
    available = not failures
   
    return {
        "modal_end_loss_power_terms_model_applicable": model_applicable,
        "modal_end_loss_power_terms_available": available,
        "modal_end_loss_power_unavailability_reason": (None if available else ";".join(failures)),
        "modal_end_loss_power_validation_failures": failures,
        "modal_end_loss_power_validation_model": ("Egedal_fixed_boundary_ion_wall_plus_Eq68_Eq69_electron_wall"),
        "modal_end_loss_power_validation_dependencies": [ "fixed_boundary_lost_reconstruction", "eq68_current_balance", "total_device_loss_convention", "full_device_lost_population"],
        "modal_end_loss_power_identity_relative_tolerance": identity_tolerance,
        "modal_ion_wall_power_identity_relative_error": (ion_identity_relative_error),
        "modal_total_wall_power_identity_relative_error": (total_identity_relative_error),
    }