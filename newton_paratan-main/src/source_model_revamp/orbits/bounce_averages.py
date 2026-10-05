"""
Magnetic only bounce motion and orbit averages
    Files contains the orbit limit and bounce average helpers used by the orbit averaged Lorentz operator
    This file is magnetic only, no ambipolar potential or electrostatic loss boundary, those are added down the line in the model's calculations as corrections

Physics relations
-----------------
    v_parallel^2 / v^2 = 1 - Λ * B_tilde(ζ)

    trapped-particle magnetic turning point:
        B_tilde(ζ_b) = 1 / Λ

    normalized bounce-time integral:
        τ_tilde_b(Λ) = ∫_0^ζ_b dζ / sqrt(1 - Λ * B_tilde(ζ))

    full bounce period:
        τ_b = 4 * L * τ_tilde_b / v

    Egedal Eq. 42 spatial average:
        <Q>_ζ = (1/τ_tilde_b) * ∫_0^ζ_b Q(ζ) * [n(ζ)/<n>] dζ / sqrt(1 - Λ * B_tilde(ζ))
"""
from __future__ import annotations
from collections.abc import Callable
import numpy as np
from numpy.typing import ArrayLike
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from source_model_revamp.orbits.invariants import magnetic_loss_boundary_lambda, speed_from_energy_m_s, turning_point_B_tilde

_GAUSS_LEGENDRE_ORDER = 8
_GAUSS_LEGENDRE_PANELS = 256
_GAUSS_NODES, _GAUSS_WEIGHTS = leggauss(_GAUSS_LEGENDRE_ORDER)
_PANEL_LEFT = np.arange(_GAUSS_LEGENDRE_PANELS, dtype=float) / _GAUSS_LEGENDRE_PANELS
_PANEL_RIGHT = np.arange(1, _GAUSS_LEGENDRE_PANELS + 1, dtype=float) / _GAUSS_LEGENDRE_PANELS
_QUADRATURE_U = (0.5 * (_PANEL_RIGHT - _PANEL_LEFT)[:, None] * _GAUSS_NODES[None, :] + 0.5 * (_PANEL_LEFT + _PANEL_RIGHT)[:, None]).ravel()
_QUADRATURE_WEIGHTS = (0.5 * (_PANEL_RIGHT - _PANEL_LEFT)[:, None] * np.broadcast_to(_GAUSS_WEIGHTS[None, :], (_GAUSS_LEGENDRE_PANELS, _GAUSS_LEGENDRE_ORDER))).ravel()
_SEPARATRIX_GAUSS_NODES, _SEPARATRIX_GAUSS_WEIGHTS = leggauss(16)
TURNING_POINT_QUADRATURE_MODEL = "quadratic_endpoint_transform_composite_gauss_legendre"
SEPARATRIX_BOUNCE_TIME_MODEL = "cell_integrated_logarithmic_separatrix"

def _require_lambda_domain(Lambda: float) -> float:
    """Return Λ as a float after enforcing the magnetic only domain 0 <= Λ <= 1"""
    Lambda = float(Lambda)
    if Lambda < 0.0 or Lambda > 1.0:
        raise ValueError("Lambda must satisfy 0 <= Lambda <= 1 for magnetic only orbit averages")
    
    return Lambda

def positive_orbit_limit_zeta(Lambda: float, B_tilde_function: Callable[[float], float], mirror_ratio: float, zeta_throat: float = 1.0) -> float:
    """
    Positive integration limit for a magnetic orbit
        For trapped particles, this is the positive magnetic turning point satisfyingB_tilde(ζ_b) = 1/Λ  
        For passing/loss cone particles, no central cell magnetic turning point exists, so the integration limit is the throat
    """
    Lambda = _require_lambda_domain(Lambda)
    lambda_boundary = magnetic_loss_boundary_lambda(mirror_ratio)
    if Lambda <= lambda_boundary:
        return float(zeta_throat)
    if Lambda == 1.0:
        return 0.0
    target_B_tilde = turning_point_B_tilde(Lambda)

    return float(brentq(lambda zeta: B_tilde_function(zeta) - target_B_tilde, 0.0, zeta_throat,))

def parallel_speed_fraction(zeta: float, Lambda: float, B_tilde_function: Callable[[float], float]) -> float:
    """Return |v_parallel| / v = sqrt(1 - Λ * B_tilde(ζ))"""
    return float(np.sqrt(1.0 - Lambda * B_tilde_function(zeta)))

def inverse_parallel_speed_fraction(zeta: float, Lambda: float, B_tilde_function: Callable[[float], float]) -> float:
    """Return v / |v_parallel| = 1 / sqrt(1 - Λ * B_tilde(ζ))"""
    return 1.0 / parallel_speed_fraction(zeta=zeta, Lambda=Lambda, B_tilde_function=B_tilde_function)

def quadratic_endpoint_coordinate(u: float, zeta_b: float) -> tuple[float, float]:
    """
    Map 0 <= u <= 1 to 0 <= ζ <= ζ_b with quadratic endpoint clustering

    ζ = ζ_b * (2u - u^2) gives
    ζ_b - ζ= ζ_b * (1 - u)^2 and
    dζ / du = 2 * ζ_b * (1 - u)

    For a simple magnetic turning point 
        1 - Λ * B_tilde is proportional to ζ_b - ζ
    The Jacobian therefore cancels its inverse square root exactly and leaves a bounded transformed integrand
    """
    coordinate = float(u)
    upper = float(zeta_b)
    if not 0.0 <= coordinate <= 1.0:
        raise ValueError("u must satisfy 0 <= u <= 1")
    if upper < 0.0:
        raise ValueError("zeta_b must be nonnegative")
    distance_fraction = 1.0 - coordinate
    zeta = upper * (1.0 - distance_fraction**2)
    jacobian = 2.0 * upper * distance_fraction

    return float(zeta), float(jacobian)

def throat_turning_point_order(B_tilde_function: Callable[[float], float], zeta_throat: float = 1.0) -> float:
    """Estimate p in B_throat - B(ζ_throat - δ) proportional to δ**p"""
    throat = float(zeta_throat)
    if not np.isfinite(throat) or throat <= 0.0:
        raise ValueError("zeta_throat must be finite and positive")
    B_throat = float(B_tilde_function(throat))
    orders: list[float] = []
    for relative_distance in (1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2):
        outer_distance = max(relative_distance * throat, 1.0e-7)
        inner_distance = 0.5 * outer_distance
        inner_drop = B_throat - float(B_tilde_function(throat - inner_distance))
        outer_drop = B_throat - float(B_tilde_function(throat - outer_distance))
        if np.isfinite(inner_drop) and np.isfinite(outer_drop) and inner_drop > 0.0 and outer_drop > inner_drop:
            orders.append(float(np.log(outer_drop / inner_drop) / np.log(outer_distance / inner_distance)))
    if not orders:
        return np.nan
    
    return float(np.max(orders))

def has_logarithmic_smooth_throat_separatrix(B_tilde_function: Callable[[float], float], zeta_throat: float = 1.0) -> bool:
    """Return True when local field behavior is consistent with a quadratic throat maximum"""
    order = throat_turning_point_order(B_tilde_function=B_tilde_function, zeta_throat=zeta_throat)
    return bool(np.isfinite(order) and order >= 1.98)

def _quadratic_endpoint_weighted_integral(*, zeta_b: float, Lambda: float, B_tilde_function: Callable[[float], float], numerator_function: Callable[[float], float]) -> float:
    """Integrate numerator/sqrt(1 - Λ * B_tilde) after the quadratic transform"""
    upper = float(zeta_b)
    if upper == 0.0:
        return 0.0
    B_upper = float(B_tilde_function(upper))
    endpoint_argument = 1.0 - Lambda * B_upper
    roundoff_scale = max(1.0, abs(Lambda * B_upper))
    argument_tolerance = 128.0 * np.finfo(float).eps * roundoff_scale
    endpoint_is_turning_point = abs(endpoint_argument) <= argument_tolerance
    distance_fraction = 1.0 - _QUADRATURE_U
    zeta = upper * (1.0 - distance_fraction**2)
    jacobian = 2.0 * upper * distance_fraction
    B_values = _evaluate_function_grid(B_tilde_function, zeta, name="B_tilde_function")
    if endpoint_is_turning_point:
        argument = endpoint_argument + Lambda * (B_upper - B_values)
    else:
        argument = 1.0 - Lambda * B_values
    if float(np.min(argument)) < -argument_tolerance:
        raise ValueError("B_tilde_function produced a negative parallel energy argument inside the orbit")
    positive = argument > 0.0
    transformed_values = np.zeros_like(argument)
    if np.any(positive):
        numerator = _evaluate_function_grid(numerator_function, zeta[positive], name="bounce integral numerator")
        transformed_values[positive] = numerator * jacobian[positive] / np.sqrt(argument[positive])

    return float(np.dot(_QUADRATURE_WEIGHTS, transformed_values))

def normalized_bounce_time(Lambda: float, B_tilde_function: Callable[[float], float], mirror_ratio: float, zeta_throat: float = 1.0) -> float:
    """
    Evaluate the normalized bounce time integral τ_tilde_b(Λ)
        For trapped particles the upper limit is the turning point  
        For passingparticles it is the throat
    """
    Lambda = _require_lambda_domain(Lambda)
    lambda_boundary = magnetic_loss_boundary_lambda(mirror_ratio)
    boundary_tolerance = 64.0 * np.finfo(float).eps * max(1.0, abs(lambda_boundary))
    if abs(Lambda - lambda_boundary) <= boundary_tolerance and has_logarithmic_smooth_throat_separatrix(B_tilde_function, zeta_throat):
        return np.inf
    zeta_b = positive_orbit_limit_zeta(Lambda=Lambda, B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, zeta_throat=zeta_throat)
    
    return _quadratic_endpoint_weighted_integral(zeta_b=zeta_b, Lambda=Lambda, B_tilde_function=B_tilde_function, numerator_function=lambda zeta: 1.0)

def normalized_bounce_time_grid(lambda_values: ArrayLike, B_tilde_function: Callable[[float], float], mirror_ratio: float, zeta_throat: float = 1.0) -> np.ndarray:
    """Evaluate τ_tilde_b on a Λ grid using a cell average at 1/R_M"""
    lambdas = np.asarray(lambda_values, dtype=float)
    if lambdas.ndim != 1:
        raise ValueError("lambda_values must be a 1D array")
    boundary = magnetic_loss_boundary_lambda(mirror_ratio)
    boundary_matches = np.flatnonzero(lambdas == boundary)
    tau = np.empty_like(lambdas)
    boundary_index = int(boundary_matches[0]) if boundary_matches.size == 1 else None
    for index, Lambda in enumerate(lambdas):
        if index == boundary_index:
            continue
        tau[index] = normalized_bounce_time(Lambda=float(Lambda), B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, zeta_throat=zeta_throat)
    if boundary_index is not None:
        if boundary_index == 0 or boundary_index == lambdas.size - 1:
            raise ValueError("the exact magnetic separatrix must be an Lambda grid point")
        left = float(lambdas[boundary_index - 1])
        right = float(lambdas[boundary_index + 1])
        cell_integral = 0.0
        for interval_left, interval_right in ((left, boundary), (boundary, right)):
            mapped_lambda = (0.5 * (interval_right - interval_left) * _SEPARATRIX_GAUSS_NODES + 0.5 * (interval_left + interval_right))
            mapped_tau = np.asarray([normalized_bounce_time(Lambda=float(Lambda), B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, zeta_throat=zeta_throat,) for Lambda in mapped_lambda], dtype=float,)
            cell_integral += 0.5 * (interval_right - interval_left) * float(np.dot(_SEPARATRIX_GAUSS_WEIGHTS, mapped_tau))
        left_width = boundary - left
        right_width = right - boundary
        tau[boundary_index] = (2.0 * cell_integral - left_width * tau[boundary_index - 1] - right_width * tau[boundary_index + 1]) / (left_width + right_width)
        if not np.isfinite(tau[boundary_index]) or tau[boundary_index] <= 0.0:
            raise ValueError("failed to construct a finite positive cell averaged separatrix bounce time")
    elif boundary_matches.size > 1:
        raise ValueError("lambda_values may contain the exact magnetic separatrix at most once")
    
    return tau

def full_bounce_period_s(Lambda: float, energy_J: float, mass_kg: float, half_length_m: float, B_tilde_function: Callable[[float], float], mirror_ratio: float, zeta_throat: float = 1.0) -> float:
    """Return full bounce period τ_b = 4 * L * τ_tilde_b / v"""
    speed_m_s = speed_from_energy_m_s(energy_J=energy_J, mass_kg=mass_kg)
    tau_tilde_b = normalized_bounce_time(Lambda=Lambda, B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, zeta_throat=zeta_throat)

    return float(4.0 * half_length_m * tau_tilde_b / speed_m_s)

def spatial_average(Lambda: float, quantity_function: Callable[[float], float], B_tilde_function: Callable[[float], float], mirror_ratio: float, density_ratio_function: Callable[[float], float] | None = None, zeta_throat: float = 1.0) -> float:
    """Evaluate the Egedal Eq. 42 orbit/spatial average"""
    if density_ratio_function is None:
        density_ratio_function = lambda zeta: 1.0
    Lambda = _require_lambda_domain(Lambda)
    lambda_boundary = magnetic_loss_boundary_lambda(mirror_ratio)
    boundary_tolerance = 64.0 * np.finfo(float).eps * max(1.0, abs(lambda_boundary))
    if abs(Lambda - lambda_boundary) <= boundary_tolerance and has_logarithmic_smooth_throat_separatrix(B_tilde_function, zeta_throat):
        return float(quantity_function(zeta_throat) * density_ratio_function(zeta_throat))
    zeta_b = positive_orbit_limit_zeta(Lambda=Lambda, B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, zeta_throat=zeta_throat)
    tau_tilde_b = normalized_bounce_time(Lambda=Lambda, B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, zeta_throat=zeta_throat)
    if zeta_b == 0.0 or tau_tilde_b == 0.0:
        return float(quantity_function(0.0) * density_ratio_function(0.0))
    value = _quadratic_endpoint_weighted_integral(zeta_b=zeta_b, Lambda=Lambda, B_tilde_function=B_tilde_function, numerator_function=lambda zeta: quantity_function(zeta) * density_ratio_function(zeta))
    
    return float(value / tau_tilde_b)

def spatial_average_grid(lambda_values: ArrayLike, quantity_function: Callable[[float], float], B_tilde_function: Callable[[float], float], mirror_ratio: float, density_ratio_function: Callable[[float], float] | None = None, zeta_throat: float = 1.0) -> np.ndarray:
    """Evaluate Eq 42 spatial average on a Λ grid"""
    lambdas = np.asarray(lambda_values, dtype=float)
    return np.array([spatial_average(Lambda=float(Lambda), quantity_function=quantity_function, B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, density_ratio_function=density_ratio_function, zeta_throat=zeta_throat,) for Lambda in lambdas], dtype=float)

def _evaluate_function_grid(function: Callable[[float], float], coordinates: np.ndarray, *, name: str) -> np.ndarray:
    """Evaluate a scalar or vector callable on one finite coordinate grid"""
    values = np.asarray(coordinates, dtype=float)
    try:
        evaluated = np.asarray(function(values), dtype=float)
        if evaluated.shape == values.shape:
            result = evaluated
        elif evaluated.ndim == 0:
            result = np.full(values.shape, float(evaluated), dtype=float)
        else:
            raise ValueError
    except (TypeError, ValueError):
        result = np.fromiter((float(function(float(coordinate))) for coordinate in values), dtype=float, count=values.size)
    if result.shape != values.shape or np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must return one finite value per coordinate")
    
    return result

def _unit_density_ratio(coordinates: ArrayLike) -> float | np.ndarray:
    """Return the uniform Eq 42 density ratio on scalar or array coordinates"""
    values = np.asarray(coordinates, dtype=float)
    result = np.ones_like(values)
    return float(result) if np.ndim(coordinates) == 0 else result

def spatial_average_lorentz_weight_grid(lambda_values: ArrayLike, B_tilde_function: Callable[[float], float], mirror_ratio: float, density_ratio_function: Callable[[float], float], zeta_throat: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the two Eq 42 averages required by the full local Eq 40 operator"""
    lambdas = np.asarray(lambda_values, dtype=float)
    if lambdas.ndim != 1:
        raise ValueError("lambda_values must be a 1D array")
    inverse_B_density_average = np.empty_like(lambdas)
    density_weight_average = np.empty_like(lambdas)
    boundary = magnetic_loss_boundary_lambda(mirror_ratio)
    boundary_tolerance = 64.0 * np.finfo(float).eps * max(1.0, abs(boundary))
    smooth_separatrix = has_logarithmic_smooth_throat_separatrix(B_tilde_function, zeta_throat)

    for index, Lambda_value in enumerate(lambdas):
        Lambda = _require_lambda_domain(float(Lambda_value))
        if abs(Lambda - boundary) <= boundary_tolerance and smooth_separatrix:
            throat_density_ratio = float(density_ratio_function(float(zeta_throat)))
            throat_B_tilde = float(B_tilde_function(float(zeta_throat)))
            inverse_B_density_average[index] = throat_density_ratio / throat_B_tilde
            density_weight_average[index] = throat_density_ratio
            continue

        zeta_b = positive_orbit_limit_zeta(Lambda=Lambda, B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, zeta_throat=zeta_throat)
        if zeta_b == 0.0:
            center_density_ratio = float(density_ratio_function(0.0))
            center_B_tilde = float(B_tilde_function(0.0))
            inverse_B_density_average[index] = center_density_ratio / center_B_tilde
            density_weight_average[index] = center_density_ratio
            continue

        distance_fraction = 1.0 - _QUADRATURE_U
        zeta = zeta_b * (1.0 - distance_fraction**2)
        jacobian = 2.0 * zeta_b * distance_fraction
        B_values = _evaluate_function_grid(B_tilde_function, zeta, name="B_tilde_function",)
        if np.any(B_values <= 0.0):
            raise ValueError("B_tilde_function must remain positive along every orbit")
        B_upper = float(B_tilde_function(float(zeta_b)))
        endpoint_argument = 1.0 - Lambda * B_upper
        roundoff_scale = max(1.0, abs(Lambda * B_upper))
        argument_tolerance = 128.0 * np.finfo(float).eps * roundoff_scale

        if abs(endpoint_argument) <= argument_tolerance:
            argument = endpoint_argument + Lambda * (B_upper - B_values)
        else:
            argument = 1.0 - Lambda * B_values

        if float(np.min(argument)) < -argument_tolerance:
            raise ValueError("B_tilde_function produced a negative parallel energy argument inside the orbit")
        argument = np.maximum(argument, np.finfo(float).tiny)
        base_integrand = jacobian / np.sqrt(argument)
        density_ratio = _evaluate_function_grid(density_ratio_function, zeta, name="density_ratio_function")
        density_scale = max(float(np.max(np.abs(density_ratio))), 1.0)
        density_tolerance = 1.0e-12 * density_scale

        if float(np.min(density_ratio)) < -density_tolerance:
            raise ValueError("density_ratio_function must remain nonnegative")
        density_ratio = np.where(density_ratio < 0.0, 0.0, density_ratio)
        bounce_integral = float(np.dot(_QUADRATURE_WEIGHTS, base_integrand))
        
        if not np.isfinite(bounce_integral) or bounce_integral <= 0.0:
            raise ValueError("Eq 42 bounce integral must be finite and positive")
        density_integral = float(np.dot(_QUADRATURE_WEIGHTS, base_integrand * density_ratio))
        inverse_B_density_integral = float(np.dot(_QUADRATURE_WEIGHTS, base_integrand * density_ratio / B_values))
        density_weight_average[index] = density_integral / bounce_integral
        inverse_B_density_average[index] = (inverse_B_density_integral / bounce_integral)

    return inverse_B_density_average, density_weight_average

def spatial_average_inverse_B_tilde(Lambda: float, B_tilde_function: Callable[[float], float], mirror_ratio: float, density_ratio_function: Callable[[float], float] | None = None, zeta_throat: float = 1.0) -> float:
    """Evaluate <1/B_tilde>_z for the orbit averaged Lorentz operator"""
    return spatial_average(Lambda=Lambda, quantity_function=lambda zeta: 1.0 / B_tilde_function(zeta), B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, density_ratio_function=density_ratio_function, zeta_throat=zeta_throat)

def spatial_average_inverse_B_tilde_grid(lambda_values: ArrayLike, B_tilde_function: Callable[[float], float], mirror_ratio: float, density_ratio_function: Callable[[float], float] | None = None, zeta_throat: float = 1.0) -> np.ndarray:
    """Evaluate <1/B_tilde>_z on a Λ grid"""
    return spatial_average_grid(lambda_values=lambda_values, quantity_function=lambda zeta: 1.0 / B_tilde_function(zeta), B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, density_ratio_function=density_ratio_function, zeta_throat=zeta_throat)
