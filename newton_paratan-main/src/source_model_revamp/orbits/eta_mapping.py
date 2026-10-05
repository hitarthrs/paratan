"""
η-Λ coordinate mapping, η coordinate construction from normalized bounce time τ_tilde_b(Λ)

Definitions/variables
    Λ: Egedal pitch invariant Λ = μ * B0 / E
    τ_tilde_b: normalized bounce time integral
    <1>_Λ: phase space normalization integral
        <1>_Λ = integral_0^1(τ_tilde_b(Λ) dΛ)
    η(Λ): orbit pitch coord
        η(Λ) = 1 - [integral)0_Λ(τ_tilde_b(Λ') dΛ')] / <1>_Λ
Boundary mapping
    Λ = 0        -> η = 1
    Λ = 1        -> η = 0
    Λ = 1 / R_M  -> η_TP
Trapped region
    0 <= η <= η_TP
Loss cone region
    η_TP < η <= 1
"""
from __future__ import annotations
from collections.abc import Callable
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.interpolate import PchipInterpolator
from source_model_revamp.orbits.bounce_averages import normalized_bounce_time_grid
from source_model_revamp.orbits.invariants import magnetic_loss_boundary_lambda

@dataclass(frozen=True)
class EtaLambdaMap:
    """
    Container for η-Λ mapping arrays and interpolators
        lambda_grid is strictly increasing on [0, 1]
        eta_grid decreases from 1 to 0  
        the derivative d_lambda_d_eta may be singular at η = 0 if the final τ_tilde_b value is zero
    """
    lambda_grid: np.ndarray
    tau_tilde_grid: np.ndarray
    eta_grid: np.ndarray
    lambda_normalization: float
    eta_of_lambda: Callable[[ArrayLike], np.ndarray]
    lambda_of_eta: Callable[[ArrayLike], np.ndarray]
    d_eta_d_lambda: Callable[[ArrayLike], np.ndarray]
    d_lambda_d_eta: Callable[[ArrayLike], np.ndarray]

def _as_lambda_tau_arrays(lambda_grid: ArrayLike, tau_tilde_grid: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    """Convert mapping inputs to float arrays without applying domain checks"""
    Lambda = np.asarray(lambda_grid, dtype=float)
    tau = np.asarray(tau_tilde_grid, dtype=float)
    return Lambda, tau

def _validate_lambda_tau_grid(lambda_grid: ArrayLike, tau_tilde_grid: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    """Validate the full domain Λ grid and nonnegative bounce time measure"""
    Lambda, tau = _as_lambda_tau_arrays(lambda_grid, tau_tilde_grid)
    if Lambda.ndim != 1 or tau.ndim != 1:
        raise ValueError("lambda_grid and tau_tilde_grid must be 1D arrays")
    if Lambda.shape != tau.shape:
        raise ValueError("lambda_grid and tau_tilde_grid must have the same shape")
    if Lambda.size < 2:
        raise ValueError("lambda_grid must contain at least two points")
    if not np.all(np.isfinite(Lambda)):
        raise ValueError("lambda_grid must contain only finite values")
    if not np.all(np.isfinite(tau)):
        raise ValueError("tau_tilde_grid must contain only finite values")
    if not np.isclose(Lambda[0], 0.0):
        raise ValueError("lambda_grid must start at 0 for the η-Λ mapping")
    if not np.isclose(Lambda[-1], 1.0):
        raise ValueError("lambda_grid must end at 1 for the η-Λ mapping")
    if not np.all(np.diff(Lambda) > 0.0):
        raise ValueError("lambda_grid must be strictly increasing")
    if np.any(tau < 0.0):
        raise ValueError("tau_tilde_grid must be nonnegative")
    
    return Lambda, tau

def lambda_phase_space_normalization(lambda_grid: ArrayLike, tau_tilde_grid: ArrayLike) -> float:
    """Return <1>_Λ = ∫_0^1 (τ_tilde_b(Λ) dΛ)"""
    Lambda, tau = _validate_lambda_tau_grid(lambda_grid, tau_tilde_grid)
    return float(trapezoid(tau, Lambda))

def cumulative_bounce_integral(lambda_grid: ArrayLike, tau_tilde_grid: ArrayLike) -> np.ndarray:
    """Return C(Λ) = ∫_0^Λ (τ_tilde_b(Λ') dΛ')"""
    Lambda, tau = _validate_lambda_tau_grid(lambda_grid, tau_tilde_grid)
    return cumulative_trapezoid(tau, Lambda, initial=0.0)

def eta_grid_from_tau_tilde(lambda_grid: ArrayLike, tau_tilde_grid: ArrayLike) -> np.ndarray:
    """Return η(Λ)= 1 - C(Λ) / <1>_Λ on a full domain Λ grid"""
    Lambda, tau = _validate_lambda_tau_grid(lambda_grid, tau_tilde_grid)
    C = cumulative_bounce_integral(lambda_grid=Lambda, tau_tilde_grid=tau)
    norm = lambda_phase_space_normalization(lambda_grid=Lambda, tau_tilde_grid=tau)
    eta = 1.0 - C / norm
    eta[0] = 1.0
    eta[-1] = 0.0

    return eta

def d_eta_d_lambda_grid(lambda_grid: ArrayLike, tau_tilde_grid: ArrayLike) -> np.ndarray:
    """Return dη/dΛ = -τ_tilde_b(Λ)/<1>_Λ on the Λ grid"""
    Lambda, tau = _validate_lambda_tau_grid(lambda_grid, tau_tilde_grid)
    norm = lambda_phase_space_normalization(lambda_grid=Lambda, tau_tilde_grid=tau)

    return -tau / norm

def build_eta_lambda_map(lambda_grid: ArrayLike, tau_tilde_grid: ArrayLike) -> EtaLambdaMap:
    """Build η(Λ), Λ(η), and derivative interpolators from τ_tilde_b(Λ)"""
    Lambda, tau = _validate_lambda_tau_grid(lambda_grid, tau_tilde_grid)
    norm = lambda_phase_space_normalization(lambda_grid=Lambda, tau_tilde_grid=tau)
    eta = eta_grid_from_tau_tilde(lambda_grid=Lambda, tau_tilde_grid=tau)
    d_eta_grid = d_eta_d_lambda_grid(lambda_grid=Lambda, tau_tilde_grid=tau)
    eta_of_lambda_interp = PchipInterpolator(Lambda, eta, extrapolate=False)
    lambda_of_eta_interp = PchipInterpolator(eta[::-1], Lambda[::-1], extrapolate=False)
    d_eta_d_lambda_interp = PchipInterpolator(Lambda, d_eta_grid, extrapolate=False)
    tau_of_lambda_interp = PchipInterpolator(Lambda, tau, extrapolate=False)

    def eta_of_lambda(values: ArrayLike) -> np.ndarray:
        return eta_of_lambda_interp(values)
    def lambda_of_eta(values: ArrayLike) -> np.ndarray:
        return lambda_of_eta_interp(values)
    def d_eta_d_lambda(values: ArrayLike) -> np.ndarray:
        return d_eta_d_lambda_interp(values)
    def d_lambda_d_eta(values: ArrayLike) -> np.ndarray:
        Lambda_values = lambda_of_eta_interp(values)
        tau_values = tau_of_lambda_interp(Lambda_values)
        
        return -norm / tau_values

    return EtaLambdaMap(lambda_grid=Lambda, tau_tilde_grid=tau, eta_grid=eta, lambda_normalization=norm, eta_of_lambda=eta_of_lambda, lambda_of_eta=lambda_of_eta, d_eta_d_lambda=d_eta_d_lambda, d_lambda_d_eta=d_lambda_d_eta)

def build_eta_lambda_map_from_bounce(lambda_grid: ArrayLike, B_tilde_function: Callable[[float], float], mirror_ratio: float, zeta_throat: float = 1.0) -> EtaLambdaMap:
    """Build η-Λ map from a magnetic geometry through τ_tilde_b(Λ)"""
    tau = normalized_bounce_time_grid(lambda_values=lambda_grid, B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, zeta_throat=zeta_throat)
    return build_eta_lambda_map(lambda_grid=lambda_grid, tau_tilde_grid=tau)

def eta_trapped_passing_boundary(eta_lambda_map: EtaLambdaMap, mirror_ratio: float) -> float:
    """Return η_TP = η(1/R_M)"""
    lambda_boundary = magnetic_loss_boundary_lambda(mirror_ratio)
    return float(eta_lambda_map.eta_of_lambda(lambda_boundary))

def lambda_from_eta(eta_values: ArrayLike, eta_lambda_map: EtaLambdaMap) -> np.ndarray:
    """Evaluate Λ(η)"""
    return eta_lambda_map.lambda_of_eta(eta_values)

def eta_from_lambda(lambda_values: ArrayLike, eta_lambda_map: EtaLambdaMap) -> np.ndarray:
    """Evaluate η(Λ)"""
    return eta_lambda_map.eta_of_lambda(lambda_values)

def d_eta_d_lambda(lambda_values: ArrayLike, eta_lambda_map: EtaLambdaMap) -> np.ndarray:
    """Evaluate dη/dΛ"""
    return eta_lambda_map.d_eta_d_lambda(lambda_values)

def d_lambda_d_eta(eta_values: ArrayLike, eta_lambda_map: EtaLambdaMap) -> np.ndarray:
    """Evaluate dΛ/dη"""
    return eta_lambda_map.d_lambda_d_eta(eta_values)