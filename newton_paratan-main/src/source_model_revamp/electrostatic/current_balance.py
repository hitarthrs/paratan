"""Ambipolar electron and ion current balance for the Egedal FBIS closure

Ion particle losses are converted to positive current magnitudes with
I_i = e Σ_s Z_s Γ_s

Electron losses use the exact tail refilling kernel E1(y) with
y = E_barrier / T_e
The barrier is inverted so the electron loss current matches the positive ion loss current

Two electron density normalizations are supported
The actual midplane density path couples the electron tail integral to the Eq 70 midplane confined fraction
The parent Maxwellian density path holds n0 fixed after the electrostatic profile supplies that normalization

Currents are stored as positive loss magnitudes in A
Energies are in J densities are in m⁻³ and particle loss rates are in s⁻¹
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from scipy.constants import elementary_charge
from scipy.optimize import brentq
from scipy.special import gammainc
from source_model_revamp.losses.electron_loss_boundary import electron_wall_potential_from_barrier_V
from source_model_revamp.losses.electron_loss_rate import EGEDAL_TAIL_REFILLING_MODEL, egedal_eq68_asymptotic_relative_correction, egedal_eq68_asymptotic_tail_refilling_kernel, egedal_eq69_asymptotic_relative_correction, egedal_exact_mean_wall_kinetic_energy_J, egedal_exact_tail_refilling_kernel, egedal_tail_refilling_loss_prefactor_s

@dataclass(frozen=True)
class IonCurrentLoss:
    """Positive ion loss current state

    ion_particle_loss_rates_s and ion_charge_numbers carry species on axis zero
    species_currents_A = e Z_s Γ_s
    total_ion_current_A = Σ_s e Z_s Γ_s
    total_ion_particle_loss_rate_s = Σ_s Γ_s
    """
    ion_particle_loss_rates_s: float | np.ndarray
    ion_charge_numbers: float | np.ndarray
    species_currents_A: float | np.ndarray
    total_ion_current_A: float | np.ndarray
    total_ion_particle_loss_rate_s: float | np.ndarray

@dataclass(frozen=True)
class ElectronCurrentBalance:
    """Electron barrier state that matches a specified positive loss current

    target_current_A and electron_current_A are positive current magnitudes
    current_residual_A = I_e − I_target
    normalized_barrier = y = E_barrier / T_e
    wall_potential_relative_to_midplane_V follows the electron confining potential sign
    electron_parent_maxwellian_n0_m3 is the Eq 70 parent Maxwellian density n0
    electron_collision_density_m3 records the density associated with the supplied collision state
    mean_wall_kinetic_energy_J comes from the exact tail energy integral
    Eq 68 and Eq 69 correction fields compare the exact integrals with their high barrier limits
    """
    target_current_A: float | np.ndarray
    electron_current_A: float | np.ndarray
    current_residual_A: float | np.ndarray
    barrier_energy_J: float | np.ndarray
    normalized_barrier: float | np.ndarray
    wall_potential_relative_to_midplane_V: float | np.ndarray
    electron_particle_loss_rate_s: float | np.ndarray
    electron_wall_power_loss_W: float | np.ndarray
    electron_model: str
    number_of_ends: float
    electron_collision_frequency_s: float | np.ndarray | None = None
    electron_loss_prefactor_s: float | np.ndarray | None = None
    loss_geometry_factor: float | np.ndarray | None = None
    loss_cone_slope_dI_dlambda: float | np.ndarray | None = None
    electron_midplane_density_m3: float | np.ndarray | None = None
    electron_parent_maxwellian_n0_m3: float | np.ndarray | None = None
    electron_collision_density_m3: float | np.ndarray | None = None
    mean_wall_kinetic_energy_J: float | np.ndarray | None = None
    exact_tail_kernel: float | np.ndarray | None = None
    eq68_asymptotic_tail_kernel: float | np.ndarray | None = None
    eq68_asymptotic_relative_correction: float | np.ndarray | None = None
    eq69_asymptotic_relative_correction: float | np.ndarray | None = None

@dataclass(frozen=True)
class AmbipolarCurrentBalance:
    """Combined ion and electron current balance state

    current_residual_A = I_e − I_i
    electron_wall_potential_relative_to_midplane_V stores Φ_wall − Φ_midplane
    """
    ion_current_loss: IonCurrentLoss
    electron_balance: ElectronCurrentBalance
    current_residual_A: float | np.ndarray
    electron_wall_potential_relative_to_midplane_V: float | np.ndarray

def _ion_species_loss_currents_A(ion_particle_loss_rates_s: ArrayLike, ion_charge_numbers: ArrayLike) -> np.ndarray:
    """Return I_s = e Z_s Γ_s as positive species loss current magnitudes

    Inputs broadcast through NumPy and species is expected on axis zero
    Negative rates and charge numbers are clipped to zero
    """
    rates = np.maximum(np.asarray(ion_particle_loss_rates_s, dtype=float), 0.0)
    charges = np.maximum(np.asarray(ion_charge_numbers, dtype=float), 0.0)
    return elementary_charge * charges * rates

def build_ion_current_loss(ion_particle_loss_rates_s: ArrayLike, ion_charge_numbers: ArrayLike) -> IonCurrentLoss:
    """Build species and total ion loss current bookkeeping

    I_s = e Z_s Γ_s
    Species is summed over axis zero and any remaining broadcast axes are retained
    """
    rates = np.maximum(np.asarray(ion_particle_loss_rates_s, dtype=float), 0.0)
    charges = np.maximum(np.asarray(ion_charge_numbers, dtype=float), 0.0)
    currents = _ion_species_loss_currents_A(rates, charges)

    return IonCurrentLoss(ion_particle_loss_rates_s=rates, ion_charge_numbers=charges, species_currents_A=currents, total_ion_current_A=np.sum(currents, axis=0), total_ion_particle_loss_rate_s=np.sum(rates, axis=0))

def egedal_tail_refilling_current_balance_for_actual_midplane_density(*, target_current_A: float, electron_midplane_density_m3: float, electron_collision_density_m3: float, electron_temperature_J: float, midplane_area_m2: float, half_length_m: float, mirror_ratio: float, loss_geometry_factor: float, loss_cone_slope_dI_dlambda: float, electron_collision_frequency_s: float, number_of_ends: float = 2.0) -> ElectronCurrentBalance:
    """Solve I_e(y) = I_target with the Eq 70 midplane density normalization

    y = E_barrier / T_e
    n_e(0) = n0 P(3/2, y)
    n0 = n_e(0) / P(3/2, y)

    P is the regularized lower incomplete gamma function represented by gammainc
    For each trial y the parent Maxwellian density n0 is updated before the exact E1(y) tail rate is evaluated
    The root therefore couples the electron tail integral to the Eq 70 midplane confined fraction

    target_current_A is a positive loss current magnitude
    number_of_ends scales the summed electron loss rate explicitly
    """
    target_current = float(target_current_A)
    midplane_density = float(electron_midplane_density_m3)
    collision_density = float(electron_collision_density_m3)
    temperature = float(electron_temperature_J)
    collision_frequency = float(electron_collision_frequency_s)
    ends = float(number_of_ends)

    if not np.isfinite(target_current) or target_current < 0.0:
        raise ValueError("target_current_A must be finite and nonnegative")
    if not np.isfinite(midplane_density) or midplane_density <= 0.0:
        raise ValueError("electron_midplane_density_m3 must be positive and finite")
    if not np.isfinite(collision_density) or collision_density <= 0.0:
        raise ValueError("electron_collision_density_m3 must be positive and finite")
    if not np.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("electron_temperature_J must be positive and finite")
    if not np.isfinite(collision_frequency) or collision_frequency <= 0.0:
        raise ValueError("electron_collision_frequency_s must be positive and finite")
    if not np.isfinite(ends) or ends <= 0.0:
        raise ValueError("number_of_ends must be positive and finite")

    # Separate the trial parent density from the supplied collision frequency
    unit_density_prefactor = float(egedal_tail_refilling_loss_prefactor_s(electron_density_midplane_m3=1.0, electron_temperature_J=temperature, midplane_area_m2=midplane_area_m2, half_length_m=half_length_m, mirror_ratio=mirror_ratio, loss_geometry_factor=loss_geometry_factor, loss_cone_slope_dI_dlambda=loss_cone_slope_dI_dlambda, electron_collision_frequency_s=collision_frequency, number_of_ends=ends))
    target_rate = target_current / elementary_charge
    def state_at_normalized_barrier(value: float) -> tuple[float, float]:
        """Return n0 and the electron particle loss rate for one trial y"""
        normalized_barrier = float(value)
        confined_fraction = float(gammainc(1.5, normalized_barrier))
        parent_density = midplane_density / max(confined_fraction, np.finfo(float).tiny)
        particle_rate = unit_density_prefactor * parent_density * float(egedal_exact_tail_refilling_kernel(normalized_barrier))
        return parent_density, particle_rate
    # Zero target current is the infinite barrier limit
    if target_rate == 0.0:
        normalized_barrier = float("inf")
        parent_density = midplane_density
        achieved_rate = 0.0
        barrier_energy = float("inf")
    else:
        lower = 1.0e-12
        upper = 1.0
        def residual(value: float) -> float:
            """Return Γ_e(y) − Γ_target for root finding"""
            return state_at_normalized_barrier(value)[1] - target_rate
        # Increase y until the electron loss falls below the target
        while residual(upper) > 0.0 and upper < 2048.0:
            upper *= 2.0
        if residual(upper) > 0.0:
            raise ValueError("unable to bracket the coupled Eq 68 and Eq 70 barrier")
        normalized_barrier = float(brentq(residual, lower, upper, rtol=8.0 * np.finfo(float).eps,))
        parent_density, achieved_rate = state_at_normalized_barrier(normalized_barrier)
        barrier_energy = normalized_barrier * temperature
    achieved_current = elementary_charge * achieved_rate
    loss_prefactor = unit_density_prefactor * parent_density
    mean_wall_energy = float(egedal_exact_mean_wall_kinetic_energy_J(normalized_barrier=normalized_barrier, electron_temperature_J=temperature))
    exact_tail_kernel = float(egedal_exact_tail_refilling_kernel(normalized_barrier))
    eq68_tail_kernel = float(egedal_eq68_asymptotic_tail_refilling_kernel(normalized_barrier))

    return ElectronCurrentBalance(
        target_current_A=target_current,
        electron_current_A=achieved_current,
        current_residual_A=achieved_current - target_current,
        barrier_energy_J=barrier_energy,
        normalized_barrier=normalized_barrier,
        wall_potential_relative_to_midplane_V=(electron_wall_potential_from_barrier_V(barrier_energy)),
        electron_particle_loss_rate_s=achieved_rate,
        electron_wall_power_loss_W=achieved_rate * mean_wall_energy,
        electron_model=EGEDAL_TAIL_REFILLING_MODEL,
        number_of_ends=ends,
        electron_collision_frequency_s=collision_frequency,
        electron_loss_prefactor_s=loss_prefactor,
        loss_geometry_factor=float(loss_geometry_factor),
        loss_cone_slope_dI_dlambda=float(loss_cone_slope_dI_dlambda),
        electron_midplane_density_m3=midplane_density,
        electron_parent_maxwellian_n0_m3=parent_density,
        electron_collision_density_m3=collision_density,
        mean_wall_kinetic_energy_J=mean_wall_energy,
        exact_tail_kernel=exact_tail_kernel,
        eq68_asymptotic_tail_kernel=eq68_tail_kernel,
        eq68_asymptotic_relative_correction=float(egedal_eq68_asymptotic_relative_correction(normalized_barrier)),
        eq69_asymptotic_relative_correction=float(egedal_eq69_asymptotic_relative_correction(normalized_barrier)),
    )

def egedal_tail_refilling_current_balance_for_parent_maxwellian_density(*, target_current_A: float, electron_parent_maxwellian_n0_m3: float, electron_midplane_density_m3: float, electron_collision_density_m3: float, electron_temperature_J: float, midplane_area_m2: float, half_length_m: float, mirror_ratio: float, loss_geometry_factor: float, loss_cone_slope_dI_dlambda: float, electron_collision_frequency_s: float, exact_midplane_potential_drop_J: float = 0.0, number_of_ends: float = 2.0) -> ElectronCurrentBalance:
    """Solve I_e(y) = I_target with a fixed Eq 70 parent Maxwellian density

    y = E_barrier / T_e and n0 is held fixed while the exact E1(y) tail rate is inverted
    This path is used after the electrostatic profile supplies the Eq 70 parent Maxwellian normalization

    exact_midplane_potential_drop_J shifts the absolute wall barrier to the wall potential relative to the selected midplane reference
    Φ_wall − Φ_midplane = −(E_barrier − ΔE_midplane) / e

    target_current_A is a positive loss current magnitude
    number_of_ends scales the summed electron loss rate explicitly
    """
    target_current = float(target_current_A)
    parent_density = float(electron_parent_maxwellian_n0_m3)
    midplane_density = float(electron_midplane_density_m3)
    collision_density = float(electron_collision_density_m3)
    temperature = float(electron_temperature_J)
    collision_frequency = float(electron_collision_frequency_s)
    midplane_drop = float(exact_midplane_potential_drop_J)
    ends = float(number_of_ends)
    if not np.isfinite(target_current) or target_current < 0.0:
        raise ValueError("target_current_A must be finite and nonnegative")
    if not np.isfinite(parent_density) or parent_density <= 0.0:
        raise ValueError("electron_parent_maxwellian_n0_m3 must be positive and finite")
    if not np.isfinite(midplane_density) or midplane_density <= 0.0:
        raise ValueError("electron_midplane_density_m3 must be positive and finite")
    if not np.isfinite(collision_density) or collision_density <= 0.0:
        raise ValueError("electron_collision_density_m3 must be positive and finite")
    if not np.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("electron_temperature_J must be positive and finite")
    if not np.isfinite(collision_frequency) or collision_frequency <= 0.0:
        raise ValueError("electron_collision_frequency_s must be positive and finite")
    if not np.isfinite(midplane_drop) or midplane_drop < 0.0:
        raise ValueError("exact_midplane_potential_drop_J must be finite and nonnegative")
    if not np.isfinite(ends) or ends <= 0.0:
        raise ValueError("number_of_ends must be positive and finite")
    # Hold the Eq 70 parent Maxwellian density fixed while y is varied
    loss_prefactor = float(egedal_tail_refilling_loss_prefactor_s(electron_density_midplane_m3=parent_density, electron_temperature_J=temperature, midplane_area_m2=midplane_area_m2, half_length_m=half_length_m, mirror_ratio=mirror_ratio, loss_geometry_factor=loss_geometry_factor, loss_cone_slope_dI_dlambda=loss_cone_slope_dI_dlambda, electron_collision_frequency_s=collision_frequency, number_of_ends=ends))
    target_rate = target_current / elementary_charge
    # Zero target current is the infinite barrier limit
    if target_rate == 0.0:
        normalized_barrier = float("inf")
        barrier_energy = float("inf")
        achieved_rate = 0.0
    else:
        lower = 1.0e-12
        upper = 1.0
        def residual(value: float) -> float:
            """Return Γ_e(y) − Γ_target for root finding"""
            return loss_prefactor * float(egedal_exact_tail_refilling_kernel(float(value))) - target_rate
        # Increase y until the electron loss falls below the target
        while residual(upper) > 0.0 and upper < 2048.0:
            upper *= 2.0
        if residual(upper) > 0.0:
            raise ValueError("unable to bracket the parent density Eq 68 wall barrier")
        normalized_barrier = float(brentq(residual, lower, upper, rtol=8.0 * np.finfo(float).eps))
        barrier_energy = normalized_barrier * temperature
        achieved_rate = loss_prefactor * float(egedal_exact_tail_refilling_kernel(normalized_barrier))
    achieved_current = elementary_charge * achieved_rate
    mean_wall_energy = float(egedal_exact_mean_wall_kinetic_energy_J(normalized_barrier=normalized_barrier, electron_temperature_J=temperature))
    exact_tail_kernel = float(egedal_exact_tail_refilling_kernel(normalized_barrier))
    eq68_tail_kernel = float(egedal_eq68_asymptotic_tail_refilling_kernel(normalized_barrier))
    wall_relative_to_midplane = -(barrier_energy - midplane_drop) / elementary_charge
   
    return ElectronCurrentBalance(
        target_current_A=target_current,
        electron_current_A=achieved_current,
        current_residual_A=achieved_current - target_current,
        barrier_energy_J=barrier_energy,
        normalized_barrier=normalized_barrier,
        wall_potential_relative_to_midplane_V=wall_relative_to_midplane,
        electron_particle_loss_rate_s=achieved_rate,
        electron_wall_power_loss_W=achieved_rate * mean_wall_energy,
        electron_model=EGEDAL_TAIL_REFILLING_MODEL,
        number_of_ends=ends,
        electron_collision_frequency_s=collision_frequency,
        electron_loss_prefactor_s=loss_prefactor,
        loss_geometry_factor=float(loss_geometry_factor),
        loss_cone_slope_dI_dlambda=float(loss_cone_slope_dI_dlambda),
        electron_midplane_density_m3=midplane_density,
        electron_parent_maxwellian_n0_m3=parent_density,
        electron_collision_density_m3=collision_density,
        mean_wall_kinetic_energy_J=mean_wall_energy,
        exact_tail_kernel=exact_tail_kernel,
        eq68_asymptotic_tail_kernel=eq68_tail_kernel,
        eq68_asymptotic_relative_correction=float(egedal_eq68_asymptotic_relative_correction(normalized_barrier)),
        eq69_asymptotic_relative_correction=float(egedal_eq69_asymptotic_relative_correction(normalized_barrier)),
    )