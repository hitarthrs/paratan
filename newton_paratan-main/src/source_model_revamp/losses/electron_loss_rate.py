"""Electron end loss rate models for mirror confinement

maxwellian_surface_flux_reference
    Reference surface streaming model over an electrostatic barrier

egedal_tail_refilling
    Collisional tail refilling model from Egedal 2022
    Eq 66 supplies the electron collision scale
    The exact particle and wall energy integrals preceding Eq 68 and Eq 69 are evaluated directly
    Eq 68 and Eq 69 are retained as high barrier diagnostics

Electron temperatures and barrier energies are in joules
Particle rates are particles per second unless a volume density suffix is present
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from scipy.constants import elementary_charge, electron_mass, epsilon_0
from scipy.special import exp1
from source_model_revamp.losses.electron_loss_boundary import maxwellian_flux_tail_factor, mean_total_energy_at_wall_for_escaping_electrons_J, mean_total_energy_midplane_for_escaping_electrons_J, normalized_electron_barrier

ELECTRON_CHARGE_MAGNITUDE_C = elementary_charge
ELECTRON_MASS_KG = electron_mass
VACUUM_PERMITTIVITY_F_PER_M = epsilon_0
EGEDAL_TAIL_REFILLING_MODEL = "egedal_tail_refilling"
MAXWELLIAN_SURFACE_FLUX_REFERENCE_MODEL = "maxwellian_surface_flux_reference"
EGEDAL_EQ68_DEFAULT_NUMBER_OF_ENDS = 1.0

@dataclass(frozen=True)
class ElectronSurfaceLossFlux:
    """Reference Maxwellian electron loss through explicit end surfaces
    
    Flux densities use particles m^-2 s^-1
    Particle rates use particles s^-1
    Current is stored as a positive loss magnitude in A
    Kinetic powers are in W and mean energies are in J
    """
    one_end_particle_flux_m2_s: float | np.ndarray
    total_particle_flux_m2_s: float | np.ndarray
    total_particle_loss_rate_s: float | np.ndarray
    total_current_loss_A: float | np.ndarray
    midplane_kinetic_power_loss_W: float | np.ndarray
    wall_kinetic_power_loss_W: float | np.ndarray
    mean_midplane_kinetic_energy_J: float | np.ndarray
    mean_wall_kinetic_energy_J: float | np.ndarray

@dataclass(frozen=True)
class ElectronVolumeLossRate:
    """Reference Maxwellian electron loss normalized by plasma volume
    
    Rate densities use particles m^-3 s^-1
    Current density uses A m^-3 and power densities use W m^-3
    Mean escaping electron energies remain in J
    """
    particle_loss_rate_m3_s: float | np.ndarray
    current_loss_A_m3: float | np.ndarray
    midplane_kinetic_power_loss_W_m3: float | np.ndarray
    wall_kinetic_power_loss_W_m3: float | np.ndarray
    mean_midplane_kinetic_energy_J: float | np.ndarray
    mean_wall_kinetic_energy_J: float | np.ndarray

@dataclass(frozen=True)
class ElectronTailRefillingLossRate:
    """Generic high barrier collisional electron tail refilling rate density
    
    The state stores particle current and wall power rate densities together with y = E_barrier / T_e and exp(-y)
    """
    particle_loss_rate_m3_s: float | np.ndarray
    current_loss_A_m3: float | np.ndarray
    wall_kinetic_power_loss_W_m3: float | np.ndarray
    normalized_barrier: float | np.ndarray
    tail_factor: float | np.ndarray

def _as_float_array(value: ArrayLike):
    """Convert an array like input to a NumPy float array without changing its shape"""
    return np.asarray(value, dtype=float)

def _scalarize(value: np.ndarray):
    """Return a Python scalar for a zero dimensional array and preserve higher dimensional arrays"""
    arr = np.asarray(value)

    return arr.item() if arr.ndim == 0 else arr

def _require_nonnegative(name: str, value: ArrayLike):
    """Validate a finite nonnegative scalar or array and return it as float values"""
    arr = _as_float_array(value)
    if np.any(~np.isfinite(arr)) or np.any(arr < 0.0):
        raise ValueError(f"{name} must be finite and nonnegative")
    
    return arr

def _require_positive(name: str, value: ArrayLike):
    """Validate a finite strictly positive scalar or array and return it as float values"""
    arr = _as_float_array(value)
    if np.any(~np.isfinite(arr)) or np.any(arr <= 0.0):
        raise ValueError(f"{name} must be finite and positive")
    
    return arr

def _require_positive_allow_infinity(name: str, value: ArrayLike):
    """Validate positive values while allowing positive infinity for the zero loss limit"""
    arr = _as_float_array(value)
    if np.any(np.isnan(arr)) or np.any(arr <= 0.0):
        raise ValueError(f"{name} must be positive and may be infinite")
    
    return arr

def _require_number_of_ends(number_of_ends: float) -> float:
    """Validate and return the explicit positive finite number of loss ends"""
    ends = float(number_of_ends)
    if not np.isfinite(ends) or ends <= 0.0:
        raise ValueError("number_of_ends must be finite and positive")
    
    return ends

def _require_coulomb_log_electron(value: float | None) -> float:
    """Require the positive finite electron Coulomb logarithm when no collision frequency is supplied"""
    if value is None:
        raise ValueError("coulomb_log_electron must be supplied by the backend collision parameter state when electron_collision_frequency_s is not supplied")
    scalar = float(value)
    if not np.isfinite(scalar) or scalar <= 0.0:
        raise ValueError("coulomb_log_electron must be positive and finite") 
    
    return scalar

def _egedal_scaled_exponential_integral_e1_asymptotic(normalized_barrier: np.ndarray) -> np.ndarray:
    """Evaluate the large y asymptotic series for exp(y) E1(y)
    
    The series is terminated when the next term is at floating point precision or begins to grow
    """
    y = np.asarray(normalized_barrier, dtype=float)
    inverse = 1.0 / y
    term = np.ones_like(y)
    total = np.ones_like(y)
    for order in range(1, 64):
        term *= -float(order) * inverse
        candidate = total + term
        if np.all(np.abs(term) <= np.finfo(float).eps * np.maximum(np.abs(candidate), 1.0)):
            total = candidate
            break
        if order > 1 and np.any(np.abs(term) > np.abs(total)):
            break
        total = candidate

    return inverse * total

def egedal_scaled_exponential_integral_e1(normalized_barrier: ArrayLike):
    """Return exp(y) E1(y) without overflow at large y
    
    Direct evaluation is used for finite y < 64 and the asymptotic series is used for larger finite y
    The positive infinity limit is zero
    """
    y = _require_positive_allow_infinity("normalized_barrier", normalized_barrier)
    result = np.zeros_like(y, dtype=float)
    finite = np.isfinite(y)
    # Direct evaluation is safe before the exponential becomes numerically large
    direct = finite & (y < 64.0)
    # Larger finite barriers use the asymptotic form of exp(y) E1(y)
    asymptotic = finite & ~direct
    result[direct] = np.exp(y[direct]) * exp1(y[direct])
    if np.any(asymptotic):
        result[asymptotic] = _egedal_scaled_exponential_integral_e1_asymptotic(y[asymptotic])

    return _scalarize(result)

def egedal_exact_tail_refilling_kernel(normalized_barrier: ArrayLike):
    """Return the exact Egedal particle tail kernel E1(y)
    
    This kernel multiplies the Eq 68 geometric and collision prefactor before the high barrier approximation is taken
    """
    y = _require_positive_allow_infinity("normalized_barrier", normalized_barrier)
    result = np.zeros_like(y, dtype=float)
    finite = np.isfinite(y)
    result[finite] = exp1(y[finite])

    return _scalarize(result)

def egedal_eq68_asymptotic_tail_refilling_kernel(normalized_barrier: ArrayLike):
    """Return the high barrier particle kernel used in Egedal Eq 68
    
    kernel = exp(-y) / y
    """
    y = _require_positive_allow_infinity("normalized_barrier", normalized_barrier)
    result = np.zeros_like(y, dtype=float)
    finite = np.isfinite(y)
    result[finite] = np.exp(-y[finite]) / y[finite]

    return _scalarize(result)

def egedal_exact_mean_wall_kinetic_energy_J(normalized_barrier: ArrayLike, electron_temperature_J: ArrayLike):
    """Return the exact mean electron kinetic energy delivered to the wall
    
    <K_wall> = T_e * [1 / (exp(y) E1(y)) - y]
    The high barrier limit approaches T_e as used by Eq 69
    """
    y = _require_positive_allow_infinity("normalized_barrier", normalized_barrier)
    T_e = _require_positive("electron_temperature_J", electron_temperature_J)
    y, T_e = np.broadcast_arrays(y, T_e)
    ratio = np.ones_like(y, dtype=float)
    finite = np.isfinite(y)
    if np.any(finite):
        scaled = np.asarray(egedal_scaled_exponential_integral_e1(y[finite]), dtype=float)
        ratio[finite] = 1.0 / scaled - y[finite]

    return _scalarize(T_e * ratio)

def egedal_eq68_asymptotic_relative_correction(normalized_barrier: ArrayLike):
    """Return the relative difference between the exact particle kernel and Eq 68
    
    correction = abs(y * exp(y) * E1(y) - 1)
    """
    y = _require_positive_allow_infinity("normalized_barrier", normalized_barrier)
    result = np.zeros_like(y, dtype=float)
    finite = np.isfinite(y)
    if np.any(finite):
        scaled = np.asarray(egedal_scaled_exponential_integral_e1(y[finite]), dtype=float)
        result[finite] = np.abs(y[finite] * scaled - 1.0)

    return _scalarize(result)

def egedal_eq69_asymptotic_relative_correction(normalized_barrier: ArrayLike):
    """Return the relative difference between exact wall energy and Eq 69
    
    correction = abs(1 / [exp(y) E1(y)] - y - 1)
    """
    y = _require_positive_allow_infinity("normalized_barrier", normalized_barrier)
    result = np.zeros_like(y, dtype=float)
    finite = np.isfinite(y)
    if np.any(finite):
        scaled = np.asarray(egedal_scaled_exponential_integral_e1(y[finite]), dtype=float)
        exact_ratio = 1.0 / scaled - y[finite]
        result[finite] = np.abs(exact_ratio - 1.0)

    return _scalarize(result)

def one_end_maxwellian_electron_flux_density_m2_s(electron_density_m3: ArrayLike, electron_temperature_J: ArrayLike, barrier_energy_J: ArrayLike, electron_mass_kg: float = ELECTRON_MASS_KG):
    """Return the one sided Maxwellian particle flux over a parallel energy barrier
    
    Gamma_e = n_e * sqrt(T_e / (2π m_e)) * exp(-E_barrier / T_e)
    The result has units particles m^-2 s^-1
    """
    n_e = _require_nonnegative("electron_density_m3", electron_density_m3)
    T_e = _require_positive("electron_temperature_J", electron_temperature_J)
    barrier = _require_nonnegative("barrier_energy_J", barrier_energy_J)
    if not np.isfinite(electron_mass_kg) or electron_mass_kg <= 0.0:
        raise ValueError("electron_mass_kg must be positive")
    tail = maxwellian_flux_tail_factor(barrier_energy_J=barrier, electron_temperature_J=T_e)

    return n_e * np.sqrt(T_e / (2.0 * np.pi * electron_mass_kg)) * tail

def total_maxwellian_electron_flux_density_m2_s(electron_density_m3: ArrayLike, electron_temperature_J: ArrayLike, barrier_energy_J: ArrayLike, number_of_ends: float = 2.0, electron_mass_kg: float = ELECTRON_MASS_KG):
    """Sum the reference Maxwellian flux density over an explicit number of identical ends"""
    if float(number_of_ends) <= 0.0:
        raise ValueError("number_of_ends must be positive")
    
    return float(number_of_ends) * one_end_maxwellian_electron_flux_density_m2_s(electron_density_m3=electron_density_m3, electron_temperature_J=electron_temperature_J, barrier_energy_J=barrier_energy_J, electron_mass_kg=electron_mass_kg)

def electron_current_from_particle_loss_rate_A(particle_loss_rate_s: ArrayLike):
    """Return the positive electron current loss magnitude
    
    I = e * dN/dt
    """
    return ELECTRON_CHARGE_MAGNITUDE_C * _as_float_array(particle_loss_rate_s)

def _safe_rate_energy_product(rate: ArrayLike, energy: ArrayLike):
    """Return rate*energy while treating zero rate times infinite energy as zero power
    
    A nonzero rate with infinite energy retains an infinite power result with the rate sign
    """
    rate_arr, energy_arr = np.broadcast_arrays(_as_float_array(rate), _as_float_array(energy))
    out = np.zeros(rate_arr.shape, dtype=float)
    nonzero_rate = rate_arr != 0.0
    finite_energy = np.isfinite(energy_arr)
    valid = nonzero_rate & finite_energy
    np.multiply(rate_arr, energy_arr, out=out, where=valid)
    infinite_power = nonzero_rate & ~finite_energy
    if np.any(infinite_power):
        out[infinite_power] = np.where(rate_arr[infinite_power] > 0.0, np.inf, -np.inf)

    return _scalarize(out)

def maxwellian_electron_surface_loss_flux(electron_density_m3: ArrayLike, electron_temperature_J: ArrayLike, barrier_energy_J: ArrayLike, end_area_m2: ArrayLike, number_of_ends: float = 2.0, include_perpendicular_wall_energy: bool = False):
    """Build the reference Maxwellian surface particle current and power losses
    
    end_area_m2 converts flux density to particle rate
    number_of_ends scales the one end flux explicitly
    The wall energy may optionally include the perpendicular Maxwellian contribution
    """
    area = _require_positive("end_area_m2", end_area_m2)
    one_end_flux = one_end_maxwellian_electron_flux_density_m2_s(electron_density_m3=electron_density_m3, electron_temperature_J=electron_temperature_J, barrier_energy_J=barrier_energy_J)
    total_flux = float(number_of_ends) * one_end_flux
    particle_rate = total_flux * area
    current = electron_current_from_particle_loss_rate_A(particle_rate)
    mean_mid = mean_total_energy_midplane_for_escaping_electrons_J(barrier_energy_J=barrier_energy_J, electron_temperature_J=electron_temperature_J)
    mean_wall = mean_total_energy_at_wall_for_escaping_electrons_J(electron_temperature_J=electron_temperature_J, include_perpendicular_energy=include_perpendicular_wall_energy)

    return ElectronSurfaceLossFlux(
        one_end_particle_flux_m2_s=one_end_flux,
        total_particle_flux_m2_s=total_flux,
        total_particle_loss_rate_s=particle_rate,
        total_current_loss_A=current,
        midplane_kinetic_power_loss_W=_safe_rate_energy_product(particle_rate, mean_mid),
        wall_kinetic_power_loss_W=_safe_rate_energy_product(particle_rate, mean_wall),
        mean_midplane_kinetic_energy_J=mean_mid,
        mean_wall_kinetic_energy_J=mean_wall,
    )

def maxwellian_electron_volume_loss_rate(electron_density_m3: ArrayLike, electron_temperature_J: ArrayLike, barrier_energy_J: ArrayLike, end_area_m2: ArrayLike, volume_m3: ArrayLike, number_of_ends: float = 2.0, include_perpendicular_wall_energy: bool = False):
    """Convert the reference Maxwellian surface loss state into rate and power densities
    
    All extensive rates and powers are divided by volume_m3 while mean energies are unchanged
    """
    surface = maxwellian_electron_surface_loss_flux(electron_density_m3=electron_density_m3, electron_temperature_J=electron_temperature_J, barrier_energy_J=barrier_energy_J, end_area_m2=end_area_m2, number_of_ends=number_of_ends, include_perpendicular_wall_energy=include_perpendicular_wall_energy)
    volume = _require_positive("volume_m3", volume_m3)

    return ElectronVolumeLossRate(
        particle_loss_rate_m3_s=surface.total_particle_loss_rate_s / volume,
        current_loss_A_m3=surface.total_current_loss_A / volume,
        midplane_kinetic_power_loss_W_m3=surface.midplane_kinetic_power_loss_W / volume,
        wall_kinetic_power_loss_W_m3=surface.wall_kinetic_power_loss_W / volume,
        mean_midplane_kinetic_energy_J=surface.mean_midplane_kinetic_energy_J,
        mean_wall_kinetic_energy_J=surface.mean_wall_kinetic_energy_J,
    )

def high_barrier_tail_refilling_kernel(electron_temperature_J: ArrayLike, barrier_energy_J: ArrayLike):
    """Return the generic high barrier kernel
    
    (T_e / E_b) * exp(-E_b / T_e)
    """
    T_e = _require_positive("electron_temperature_J", electron_temperature_J)
    barrier = _require_positive("barrier_energy_J", barrier_energy_J)
    y = normalized_electron_barrier(barrier_energy_J=barrier, electron_temperature_J=T_e)

    return (T_e / barrier) * np.exp(-y)

def collisional_tail_refilling_loss_rate_density(electron_density_m3: ArrayLike, electron_temperature_J: ArrayLike, barrier_energy_J: ArrayLike, refill_rate_s: ArrayLike, geometry_factor: ArrayLike = 1.0):
    """Build a generic collisional tail refilling loss rate density
    
    particle rate density = n_e * refill_rate_s * geometry_factor * high barrier kernel
    The wall energy uses the residual parallel value T_e
    """
    n_e = _require_nonnegative("electron_density_m3", electron_density_m3)
    nu = _require_nonnegative("refill_rate_s", refill_rate_s)
    geom = _require_nonnegative("geometry_factor", geometry_factor)
    kernel = high_barrier_tail_refilling_kernel(electron_temperature_J=electron_temperature_J, barrier_energy_J=barrier_energy_J)
    particle_rate = n_e * nu * geom * kernel
    current = ELECTRON_CHARGE_MAGNITUDE_C * particle_rate
    mean_wall = mean_total_energy_at_wall_for_escaping_electrons_J(electron_temperature_J=electron_temperature_J, include_perpendicular_energy=False)
    y = normalized_electron_barrier(barrier_energy_J=barrier_energy_J, electron_temperature_J=electron_temperature_J)

    return ElectronTailRefillingLossRate(particle_loss_rate_m3_s=particle_rate, current_loss_A_m3=current, wall_kinetic_power_loss_W_m3=particle_rate * mean_wall, normalized_barrier=y, tail_factor=np.exp(-y))

def electron_wall_power_from_current_W(electron_current_A: ArrayLike, electron_temperature_J: ArrayLike):
    """Return the Egedal Eq 69 asymptotic wall power relation
    
    P_loss,e = I_e T_e/e
    """
    current = _as_float_array(electron_current_A)
    T_e = _require_positive("electron_temperature_J", electron_temperature_J)
    
    return current * T_e / ELECTRON_CHARGE_MAGNITUDE_C

def egedal_electron_electron_collision_frequency_s(electron_density_m3: ArrayLike, electron_temperature_J: ArrayLike, coulomb_log_electron: float | None = None):
    """Electron electron collision scale from Egedal Eq 66
        nu_0^ee = n * e^4 * lnΛ / (4π * ε0^2 * m_e^2 * v_t^3),
    where
        v_t = sqrt(2T_e / m_e)
    
    The supplied electron temperature is an energy in joules
    """
    n_e = _require_nonnegative("electron_density_m3", electron_density_m3)
    T_e = _require_positive("electron_temperature_J", electron_temperature_J)
    lnL = _require_coulomb_log_electron(coulomb_log_electron)
    v_t = np.sqrt(2.0 * T_e / ELECTRON_MASS_KG)
    nu = (n_e * ELECTRON_CHARGE_MAGNITUDE_C**4 * lnL / (4.0 * np.pi * VACUUM_PERMITTIVITY_F_PER_M**2 * ELECTRON_MASS_KG**2 * v_t**3))

    return _scalarize(nu)

def egedal_electron_electron_collision_frequency_from_tau_s(spitzer_slowing_time_s: ArrayLike, fast_ion_mass_kg: float, fast_ion_charge_number: float = 1.0):
    """Return the equivalent Eq 66 collision scale using the fast ion Spitzer slowing time τ_s
    
    nu_0^ee = (1 / τ_s) * (3 sqrt(π) / (4 z_f^2)) * (m_f / m_e)
    """
    tau_s = _require_positive("spitzer_slowing_time_s", spitzer_slowing_time_s)
    m_f = float(fast_ion_mass_kg)
    z_f = float(fast_ion_charge_number)
    if not np.isfinite(m_f) or m_f <= 0.0:
        raise ValueError("fast_ion_mass_kg must be positive")
    if not np.isfinite(z_f) or z_f == 0.0:
        raise ValueError("fast_ion_charge_number must be nonzero")
    nu = (1.0 / tau_s) * (3.0 * np.sqrt(np.pi) / (4.0 * z_f**2)) * (m_f / ELECTRON_MASS_KG)

    return _scalarize(nu)

def egedal_tail_refilling_loss_prefactor_s(electron_density_midplane_m3: ArrayLike, electron_temperature_J: ArrayLike, midplane_area_m2: ArrayLike, half_length_m: ArrayLike, mirror_ratio: ArrayLike, loss_geometry_factor: ArrayLike, loss_cone_slope_dI_dlambda: ArrayLike, electron_collision_frequency_s: ArrayLike | None = None, coulomb_log_electron: float | None = None, number_of_ends: float = EGEDAL_EQ68_DEFAULT_NUMBER_OF_ENDS):
    """Return the Eq 68 particle loss factor C_e summed over explicit ends
    
    Egedal Eq 68 is one end, so the default is one and a symmetric total device balance must explicitly request 2
    
    C_e = number_of_ends * (4 * n0 * A0 / sqrt(π)) * nu_0^ee * (dI1 / dLambda)|LambdaM * l * (1/R_M) * G(l)
    
    The exact particle rate is C_e * E1(y)
    If electron_collision_frequency_s is omitted the Eq 66 collision scale is evaluated from T_e and the supplied Coulomb logarithm
    """
    n0 = _require_nonnegative("electron_density_midplane_m3", electron_density_midplane_m3)
    T_e = _require_positive("electron_temperature_J", electron_temperature_J)
    A0 = _require_positive("midplane_area_m2", midplane_area_m2)
    half_length = _require_positive("half_length_m", half_length_m)
    R_M = _require_positive("mirror_ratio", mirror_ratio)
    if np.any(R_M < 1.0):
        raise ValueError("mirror_ratio must be >= 1")
    G_l = _require_nonnegative("loss_geometry_factor", loss_geometry_factor)
    slope = _require_positive("loss_cone_slope_dI_dlambda", loss_cone_slope_dI_dlambda)
    ends = _require_number_of_ends(number_of_ends)
    if electron_collision_frequency_s is None:
        nu = egedal_electron_electron_collision_frequency_s(electron_density_m3=n0, electron_temperature_J=T_e, coulomb_log_electron=coulomb_log_electron)
    else:
        nu = _require_nonnegative("electron_collision_frequency_s", electron_collision_frequency_s)
    n0, T_e, A0, half_length, R_M, G_l, slope, nu = np.broadcast_arrays(_as_float_array(n0), _as_float_array(T_e), _as_float_array(A0), _as_float_array(half_length), _as_float_array(R_M), _as_float_array(G_l), _as_float_array(slope), _as_float_array(nu))
    one_end_prefactor = (4.0 * n0 * A0 / np.sqrt(np.pi)) * nu * slope * half_length * (1.0 / R_M) * G_l
    prefactor = ends * one_end_prefactor

    return _scalarize(prefactor)
