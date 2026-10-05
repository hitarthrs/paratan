"""
Lorentz pitch angle scattering coefficients in the magnetic invariant Lambda coordinate

This module represents the Egedal 2022 Lorentz operator in Lambda space
Eq 40 provides local coefficients at a specified normalized magnetic field
Eq 41 and Eq 42 provide orbit averaged coefficients with optional density weighting

The operator is stored in the form
    L[f] = A d f/dΛ + D d²f/dΛ²

Lambda and B_tilde are dimensionless
The coefficient arrays are therefore dimensionless and follow the supplied Lambda grid
"""
from __future__ import annotations
from collections.abc import Callable
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike
from source_model_revamp.orbits.bounce_averages import _unit_density_ratio, spatial_average_lorentz_weight_grid

LORENTZ_COEFFICIENT_KIND_LOCAL = "local"
LORENTZ_COEFFICIENT_KIND_ORBIT_AVERAGED = "orbit_averaged"
_ALLOWED_LORENTZ_COEFFICIENT_KINDS = {LORENTZ_COEFFICIENT_KIND_LOCAL, LORENTZ_COEFFICIENT_KIND_ORBIT_AVERAGED}

@dataclass(frozen=True)
class LorentzLambdaCoefficients:
    """
    Coefficients for a Λ space Lorentz operator, operator has the form
        L[f] = A * df/dΛ + D * d^2f/dΛ^2
    
    lambda_values is a one dimensional dimensionless Lambda grid
    first_derivative_coefficient stores A on that grid
    second_derivative_coefficient stores D on that grid
    coefficient_kind distinguishes local Eq 40 data from orbit averaged Eq 41 data
    """
    lambda_values: np.ndarray
    first_derivative_coefficient: np.ndarray
    second_derivative_coefficient: np.ndarray
    coefficient_kind: str

    def __post_init__(self) -> None:
        Lambda = np.asarray(self.lambda_values, dtype=float)
        first = np.asarray(self.first_derivative_coefficient, dtype=float)
        second = np.asarray(self.second_derivative_coefficient, dtype=float)
        kind = str(self.coefficient_kind)
        if kind not in _ALLOWED_LORENTZ_COEFFICIENT_KINDS:
            raise ValueError("coefficient_kind must be one of " f"{sorted(_ALLOWED_LORENTZ_COEFFICIENT_KINDS)}, got {kind!r}")
        for name, values in (("lambda_values", Lambda), ("first_derivative_coefficient", first), ("second_derivative_coefficient", second)):
            if values.ndim != 1:
                raise ValueError(f"{name} must be a 1D array")
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{name} must contain only finite values")
        if first.shape != Lambda.shape or second.shape != Lambda.shape:
            raise ValueError("lambda_values, first_derivative_coefficient, and second_derivative_coefficient must have matching shapes")
        
        object.__setattr__(self, "lambda_values", Lambda)
        object.__setattr__(self, "first_derivative_coefficient", first)
        object.__setattr__(self, "second_derivative_coefficient", second)
        object.__setattr__(self, "coefficient_kind", kind)

    @property
    def is_local(self) -> bool:
        """Return True for local Eq (40) coefficients"""
        return self.coefficient_kind == LORENTZ_COEFFICIENT_KIND_LOCAL
    
    @property
    def is_orbit_averaged(self) -> bool:
        """Return True for orbit averaged Eq (41) coefficients"""
        return self.coefficient_kind == LORENTZ_COEFFICIENT_KIND_ORBIT_AVERAGED

def require_orbit_averaged_lorentz_coefficients(coefficients: LorentzLambdaCoefficients, *, caller: str = "this operation") -> LorentzLambdaCoefficients:
    """
    Require orbit averaged Lorentz coefficients for physical basis operations

    Local Eq 40 coefficients and the orbit averaged Eq 41 coefficients share the smae container shape but
    they are not interchangeable in the physical mirror eigenbasis construction
    in the physical mirror eigenbasis construction
    """
    if coefficients.coefficient_kind != LORENTZ_COEFFICIENT_KIND_ORBIT_AVERAGED:
        raise ValueError(f"{caller} requires orbit averaged Lorentz coefficients " f"(coefficient_kind={LORENTZ_COEFFICIENT_KIND_ORBIT_AVERAGED!r}) " f"got coefficient_kind={coefficients.coefficient_kind!r}")
    
    return coefficients

def local_lorentz_first_derivative_coefficient(Lambda: ArrayLike, B_tilde: ArrayLike):
    """
    Egedal Eq (40) first derivative coefficient
        A_local = 4 / B_tilde - 6Λ
    """
    Lambda = np.asarray(Lambda, dtype=float)
    B_tilde = np.asarray(B_tilde, dtype=float)

    return 4.0 / B_tilde - 6.0 * Lambda

def local_lorentz_second_derivative_coefficient(Lambda: ArrayLike, B_tilde: ArrayLike):
    """
    Egedal Eq (40) second derivative coefficient
        D_local = (4 / B_tilde) * Λ * (1 - Λ * B_tilde) = 4Λ * (1 / B_tilde - Λ)
    """
    Lambda = np.asarray(Lambda, dtype=float)
    B_tilde = np.asarray(B_tilde, dtype=float)

    return 4.0 * Lambda * ((1.0 / B_tilde) - Lambda)

def local_lorentz_coefficients(Lambda: ArrayLike, B_tilde: ArrayLike) -> LorentzLambdaCoefficients:
    """Local Λ space Lorentz operator coefficients from Egedal Eq (40)"""
    Lambda_values = np.atleast_1d(np.asarray(Lambda, dtype=float))

    return LorentzLambdaCoefficients(lambda_values=Lambda_values, first_derivative_coefficient=local_lorentz_first_derivative_coefficient(Lambda=Lambda_values, B_tilde=B_tilde), second_derivative_coefficient=local_lorentz_second_derivative_coefficient(Lambda=Lambda_values, B_tilde=B_tilde), coefficient_kind=LORENTZ_COEFFICIENT_KIND_LOCAL)

def orbit_averaged_lorentz_first_derivative_coefficient(Lambda: ArrayLike, inverse_B_tilde_average: ArrayLike):
    """
    Egedal Eq (41) first derivative coefficient
        A_avg = 4 * <1/B_tilde>_z - 6Λ
    """
    Lambda = np.asarray(Lambda, dtype=float)
    inverse_B_avg = np.asarray(inverse_B_tilde_average, dtype=float)

    return 4.0 * inverse_B_avg - 6.0 * Lambda

def orbit_averaged_lorentz_second_derivative_coefficient(Lambda: ArrayLike, inverse_B_tilde_average: ArrayLike):
    """
    Egedal Eq (41) second derivative coefficient
        D_avg = 4Λ * (<1/B_tilde>_z - Λ)
    """
    Lambda = np.asarray(Lambda, dtype=float)
    inverse_B_avg = np.asarray(inverse_B_tilde_average, dtype=float)

    return 4.0 * Lambda * (inverse_B_avg - Lambda)

def orbit_averaged_lorentz_coefficients(Lambda: ArrayLike, inverse_B_tilde_average: ArrayLike) -> LorentzLambdaCoefficients:
    """Orbit averaged Λ space Lorentz coefficients from Egedal Eq (41)"""
    Lambda_values = np.atleast_1d(np.asarray(Lambda, dtype=float))

    return LorentzLambdaCoefficients(lambda_values=Lambda_values, first_derivative_coefficient=orbit_averaged_lorentz_first_derivative_coefficient(Lambda=Lambda_values, inverse_B_tilde_average=inverse_B_tilde_average), second_derivative_coefficient=orbit_averaged_lorentz_second_derivative_coefficient(Lambda=Lambda_values, inverse_B_tilde_average=inverse_B_tilde_average),  coefficient_kind=LORENTZ_COEFFICIENT_KIND_ORBIT_AVERAGED)

def density_weighted_orbit_averaged_lorentz_coefficients(Lambda: ArrayLike, inverse_B_tilde_density_average: ArrayLike, density_weight_average: ArrayLike) -> LorentzLambdaCoefficients:
    """Average every local Eq 40 coefficient using Egedal Eq 42"""
    Lambda_values = np.atleast_1d(np.asarray(Lambda, dtype=float))
    inverse_B_average = np.asarray(inverse_B_tilde_density_average, dtype=float)
    density_average = np.asarray(density_weight_average, dtype=float)
    if inverse_B_average.shape != Lambda_values.shape or density_average.shape != Lambda_values.shape:
        raise ValueError("density weighted Lorentz averages must match the Lambda grid")
    if np.any(~np.isfinite(inverse_B_average)) or np.any(~np.isfinite(density_average)):
        raise ValueError("density weighted Lorentz averages must be finite")
    if np.any(density_average < 0.0):
        raise ValueError("density weight average must be nonnegative")
    first = 4.0 * inverse_B_average - 6.0 * Lambda_values * density_average
    second = 4.0 * Lambda_values * (inverse_B_average - Lambda_values * density_average)
    scale = max(float(np.max(np.abs(second))) if second.size else 0.0, 1.0)
    tolerance = 512.0 * np.finfo(float).eps * scale
    if float(np.min(second)) < -tolerance:
        raise ValueError("density weighted orbit averaged Lorentz diffusion coefficient became negative")
    second = np.where(second < 0.0, 0.0, second)

    return LorentzLambdaCoefficients(lambda_values=Lambda_values, first_derivative_coefficient=first, second_derivative_coefficient=second, coefficient_kind=LORENTZ_COEFFICIENT_KIND_ORBIT_AVERAGED)

def orbit_averaged_lorentz_coefficients_from_bounce(lambda_values: ArrayLike, B_tilde_function: Callable[[float], float], mirror_ratio: float, density_ratio_function: Callable[[float], float] | None = None, zeta_throat: float = 1.0) -> LorentzLambdaCoefficients:
    """Build the Eq 41 operator by applying Eq 42 to the local Eq 40 coefficients"""
    Lambda = np.atleast_1d(np.asarray(lambda_values, dtype=float))
    active_density_ratio_function = _unit_density_ratio if density_ratio_function is None else density_ratio_function
    inverse_B_avg, density_average = spatial_average_lorentz_weight_grid(lambda_values=Lambda, B_tilde_function=B_tilde_function, mirror_ratio=mirror_ratio, density_ratio_function=active_density_ratio_function, zeta_throat=zeta_throat)
    return density_weighted_orbit_averaged_lorentz_coefficients(Lambda=Lambda, inverse_B_tilde_density_average=inverse_B_avg, density_weight_average=density_average)

def apply_lorentz_operator_on_lambda_grid(function_values: ArrayLike, lambda_grid: ArrayLike, coefficients: LorentzLambdaCoefficients, edge_order: int = 2) -> np.ndarray:
    """
    Apply a Λ space Lorentz operator to values on a Λ grid
        L[f] = A * df/dΛ + D * d^2f/dΛ^2
    """
    values = np.asarray(function_values, dtype=float)
    Lambda = np.asarray(lambda_grid, dtype=float)
    if values.shape != Lambda.shape:
        raise ValueError("function_values and lambda_grid must have matching shapes")
    if coefficients.lambda_values.shape != Lambda.shape:
        raise ValueError("coefficient arrays must match lambda_grid shape")
    df_dlambda = np.gradient(values, Lambda, edge_order=edge_order)
    d2f_dlambda2 = np.gradient(df_dlambda, Lambda, edge_order=edge_order)

    return (coefficients.first_derivative_coefficient * df_dlambda + coefficients.second_derivative_coefficient * d2f_dlambda2)
