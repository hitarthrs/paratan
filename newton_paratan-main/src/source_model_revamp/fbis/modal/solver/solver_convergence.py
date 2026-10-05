"""Convergence metrics shared by the modal FBIS solver orchestration"""
from __future__ import annotations
import numpy as np
from source_model_revamp.fbis.velocity_space_grid import SpeedGrid
from source_model_revamp.fbis.modal.types import ModalBeamComponentResult
from source_model_revamp.fbis.modal.types import ModalFBISBasis
from source_model_revamp.fbis.modal.solver.solver_loss_rates import _modal_mode_eta_integrals

def _safe_ratio(numerator: float | None, denominator: float | None) -> float | None:
    """Return numerator / denominator when both values are finite and the denominator is positive"""
    if numerator is None or denominator is None:
        return None
    n = float(numerator)
    d = float(denominator)
    if not np.isfinite(n) or not np.isfinite(d) or d <= 0.0:
        return None
    
    return n / d

def _upper_speed_tail_population_fraction(speed_grid: SpeedGrid, cell_population_weights: np.ndarray, *, tail_start_fraction: float = 0.95) -> float:
    """Return the particle fraction above a selected fraction of the finite speed domain"""
    weights = np.asarray(cell_population_weights, dtype=float)
    lower = np.asarray(speed_grid.faces_m_s[:-1], dtype=float)
    upper = np.asarray(speed_grid.faces_m_s[1:], dtype=float)
    if weights.shape != lower.shape or np.any(~np.isfinite(weights)) or np.any(weights < 0.0):
        raise ValueError("cell_population_weights must be finite, nonnegative, and match the speed grid")
    fraction = float(tail_start_fraction)
    if not np.isfinite(fraction) or not 0.0 < fraction < 1.0:
        raise ValueError("tail_start_fraction must lie inside (0, 1)")
    total = float(np.sum(weights))
    if total <= 0.0:
        return 0.0
    threshold = fraction * float(speed_grid.faces_m_s[-1])
    clipped_lower = np.maximum(lower, threshold)
    # The v³ shell measure retains the partial cell cut by the tail threshold
    full_measure = upper**3 - lower**3
    tail_measure = np.maximum(upper**3 - clipped_lower**3, 0.0)
    partial = np.divide(tail_measure, full_measure, out=np.zeros_like(tail_measure), where=full_measure > 0.0)

    return float(np.sum(weights * partial) / total)

def _upper_speed_tail_energy_fraction(speed_grid: SpeedGrid, eta_integral_by_speed: np.ndarray, *, particle_mass_kg: float, tail_start_fraction: float = 0.95) -> float:
    """Return the kinetic energy fraction above a selected fraction of the finite speed domain"""
    eta_integral = np.asarray(eta_integral_by_speed, dtype=float)
    lower = np.asarray(speed_grid.faces_m_s[:-1], dtype=float)
    upper = np.asarray(speed_grid.faces_m_s[1:], dtype=float)
    if (eta_integral.shape != lower.shape or np.any(~np.isfinite(eta_integral)) or np.any(eta_integral < 0.0)):
        raise ValueError("eta_integral_by_speed must be finite, nonnegative, and match the speed grid")
    fraction = float(tail_start_fraction)
    if not np.isfinite(fraction) or not 0.0 < fraction < 1.0:
        raise ValueError("tail_start_fraction must lie inside (0, 1)")
    mass = float(particle_mass_kg)
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("particle_mass_kg must be positive and finite")
    full_weights = (2.0 * np.pi * mass / 5.0 * eta_integral * (upper**5 - lower**5))
    total = float(np.sum(full_weights))
    if total <= 0.0:
        return 0.0
    threshold = fraction * float(speed_grid.faces_m_s[-1])
    clipped_lower = np.maximum(lower, threshold)
    # The v⁵ energy measure retains the partial cell cut by the tail threshold
    tail_weights = (2.0 * np.pi * mass / 5.0 * eta_integral * np.maximum(upper**5 - clipped_lower**5, 0.0))

    return float(np.sum(tail_weights) / total)

def _projected_modal_source_rate_from_coefficients(*, basis: ModalFBISBasis, component_results: list[ModalBeamComponentResult], volume_m3: float) -> tuple[float, list[float]]:
    """Integrate modal source coefficients over η and return component and total particle rates"""
    mode_eta_integrals = _modal_mode_eta_integrals(basis)
    component_rates: list[float] = []
    for component in component_results:
        coeff = np.asarray(component.modal_source_coefficients, dtype=float)
        component_rate_density = 4.0 * np.pi * float(np.sum(coeff * mode_eta_integrals))
        component_rates.append(component_rate_density * float(volume_m3))
        
    return float(np.sum(component_rates)), component_rates