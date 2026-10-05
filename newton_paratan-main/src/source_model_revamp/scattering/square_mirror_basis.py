"""
Square mirror Lorentz eigenbasis
    builds the square mirror reference eigenfunctions M_j(ξ) which are used later as a basis for the physical eigenbasis

Variables
    ξ: square mirror pitch coordinate
            ξ = v_parallel / v
        trapped interval is
            0 <= ξ <= ξ_TP
        with 
            ξ_TP = sqrt(1 - 1/R_M)
    L: pitch angle Lorentz operator
            L = d/dξ [(1 - ξ^2) d/dξ]
        eigenproblem is
            L * M_j = -λ_j * M_j
        which is the Legendre differential eq
            (1 - ξ^2) M'' - 2ξ*M' + λ * M = 0
        with
            λ = l(l + 1)
Boundary conditions
    Square mirror eigenfunctions satisfy
        M_j(0) = 1
        dM_j/dξ |_(ξ=0) = 0
        M_j(ξ_TP) = 0
Even solution
        M_l(ξ) = 2F1(-l/2, (l + 1)/2; 1/2; ξ^2)
    gives
        M_l(0) = 1
        dM_l/dξ |_{ξ=0} = 0
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import brentq
from scipy.special import hyp2f1
from numpy.polynomial.legendre import leggauss
from source_model_revamp.orbits.invariants import xi_trapped_passing_boundary

@dataclass(frozen=True)
class SquareMirrorBasis:
    """
    Dataclass for the square mirror eigenbasis

    Variables
        mirror_ratio: R_M
        xi_trapped_passing: trapped/passing boundary ξ_TP
        legendre_orders: Legendre orders ell_j
        eigenvalues: λ_j = ell_j * (ell_j +1)
        normalizations: α_j = integral_0^ξ_TP(M_j(ξ)^2 dξ)
    """
    mirror_ratio: float
    xi_trapped_passing: float
    legendre_orders: np.ndarray
    eigenvalues: np.ndarray
    normalizations: np.ndarray

def square_mirror_xi_trapped_passing(mirror_ratio: float) -> float:
    """ξ_TP = sqrt(1 - 1/R_M)"""
    return float(xi_trapped_passing_boundary(mirror_ratio))

def square_mirror_eigenvalue(legendre_order: ArrayLike):
    """λ = ell * (ell + 1)"""
    ell = np.asarray(legendre_order, dtype=float)
    
    return ell * (ell + 1.0)

def square_mirror_mode(xi: ArrayLike, legendre_order: float):
    """
    Even square mirror eigenfunction M_ell(ξ)
        M_ell(ξ) = 2F1(-ell/2, (ell + 1)/2; 1/2; ξ^2)
    even solution of the Legendre equation with
        M(0) = 1
        M'(0) = 0
    """
    xi = np.asarray(xi, dtype=float)
    ell = np.asarray(legendre_order, dtype=float)
    a = -0.5 * ell
    b = 0.5 * (ell + 1.0)
    c = 0.5

    return hyp2f1(a, b, c, xi**2)

def square_mirror_mode_derivative(xi: ArrayLike, legendre_order: float):
    """
    First derivative dM_ell/dξ of M = 2F1(a, b; c; xi^2)
    with
        a = -ell/2
        b = (ell + 1)/2
        c = 1/2
    so
        d/dx 2F1(a,b;c;x^2) = 2x * (ab/c) * 2F1(a+1,b+1;c+1;x^2)
    """
    xi = np.asarray(xi, dtype=float)
    ell = np.asarray(legendre_order, dtype=float)
    a = -0.5 * ell
    b = 0.5 * (ell + 1.0)
    c = 0.5

    return 2.0 * xi * (a * b / c) * hyp2f1(a + 1.0, b + 1.0, c + 1.0, xi**2)

def square_mirror_mode_second_derivative(xi: ArrayLike, legendre_order: float):
    """
    Second derivative d^2M_ell/dξ^2 from the Legendre equation
        (1 - ξ) M'' - 2ξ * M' + λ * M = 0
    so
        M'' = [2ξ * M' - λ * M] / (1 - ξ^2)
    """
    xi = np.asarray(xi, dtype=float)
    M = square_mirror_mode(xi=xi, legendre_order=legendre_order)
    dM = square_mirror_mode_derivative(xi=xi, legendre_order=legendre_order)
    lam = square_mirror_eigenvalue(legendre_order)

    return (2.0 * xi * dM - lam * M) / (1.0 - xi**2)

def square_mirror_boundary_residual(legendre_order: float, xi_trapped_passing: float) -> float:
    """
    Boundary residual
        M_ell(ξ_TP)
    eigenorders satisfy
        M_ell(ξ_TP) = 0
    """
    return float(square_mirror_mode(xi=xi_trapped_passing, legendre_order=legendre_order))

@lru_cache(maxsize=16)
def _normalization_quadrature(order: int) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic Gauss Legendre nodes for square mode normalization"""
    nodes, weights = leggauss(int(order))
    return np.asarray(nodes, dtype=float), np.asarray(weights, dtype=float)

def square_mirror_mode_normalization(legendre_order: float, xi_trapped_passing: float, quadrature_order: int = 384) -> float:
    """α_j = integral_0^ξ_TP M_j(ξ)^2 dξ with deterministic quadrature"""
    xi_tp = float(xi_trapped_passing)
    if not np.isfinite(xi_tp) or not 0.0 < xi_tp < 1.0:
        raise ValueError("xi_trapped_passing must lie inside (0, 1)")
    order = int(quadrature_order)
    if order < 32:
        raise ValueError("quadrature_order must be at least 32")
    nodes, weights = _normalization_quadrature(order)
    xi = 0.5 * xi_tp * (nodes + 1.0)
    values = np.asarray(square_mirror_mode(xi=xi, legendre_order=float(legendre_order)), dtype=float)
    integral = 0.5 * xi_tp * float(np.sum(weights * values**2))
    if not np.isfinite(integral) or integral <= 0.0:
        raise RuntimeError("square mirror normalization quadrature failed")
    return integral

def find_square_mirror_legendre_orders(mirror_ratio: float, n_modes: int, scan_order_max: float | None = None, scan_points: int = 20000) -> np.ndarray:
    """
    Find first n_modes for non-integer Legendre orders ell_j
    eigenorders are roots of
        M_ell(ξ_TP) = 0
    where even solution M_ell satisfies
        M(0) = 1
        M'(0) = 0
    """
    xi_tp = square_mirror_xi_trapped_passing(mirror_ratio)
    if scan_order_max is None:
        scan_order_max = 4.0 * n_modes + 10.0
    ell_scan = np.linspace(0.0, scan_order_max, scan_points)
    residuals = np.asarray(square_mirror_mode(xi=xi_tp, legendre_order=ell_scan), dtype=float)
    roots: list[float] = []
    for left, right, f_left, f_right in zip(ell_scan[:-1], ell_scan[1:], residuals[:-1], residuals[1:]):
        if len(roots) >= n_modes:
            break
        if f_left == 0.0:
            roots.append(float(left))
        if f_left * f_right < 0.0:
            root = brentq(lambda ell: square_mirror_boundary_residual(legendre_order=ell, xi_trapped_passing=xi_tp), float(left), float(right),)
            roots.append(float(root))
            
    return np.asarray(roots[:n_modes], dtype=float)

def build_square_mirror_basis(mirror_ratio: float, n_modes: int, scan_order_max: float | None = None, scan_points: int = 20000) -> SquareMirrorBasis:
    """Build square mirror basis"""
    xi_tp = square_mirror_xi_trapped_passing(mirror_ratio)
    orders = find_square_mirror_legendre_orders(mirror_ratio=mirror_ratio, n_modes=n_modes, scan_order_max=scan_order_max, scan_points=scan_points)
    eigenvalues = square_mirror_eigenvalue(orders)
    normalizations = np.array([square_mirror_mode_normalization(legendre_order=float(order), xi_trapped_passing=xi_tp) for order in orders], dtype=float)

    return SquareMirrorBasis(mirror_ratio=float(mirror_ratio), xi_trapped_passing=xi_tp, legendre_orders=orders, eigenvalues=eigenvalues, normalizations=normalizations)

def square_mirror_basis_values(xi: ArrayLike, basis: SquareMirrorBasis) -> np.ndarray:
    """Evaluate all square mirror modes on a ξ grid"""
    xi = np.asarray(xi, dtype=float)
    return np.array([square_mirror_mode(xi=xi, legendre_order=float(order)) for order in basis.legendre_orders], dtype=float)

def square_mirror_basis_derivatives(xi: ArrayLike, basis: SquareMirrorBasis) -> np.ndarray:
    """Evaluate dM_j/dξ for all square mirror modes"""
    xi = np.asarray(xi, dtype=float)
    return np.array([square_mirror_mode_derivative(xi=xi, legendre_order=float(order)) for order in basis.legendre_orders], dtype=float)

def square_mirror_basis_second_derivatives(xi: ArrayLike, basis: SquareMirrorBasis) -> np.ndarray:
    """Evaluate d^2M_j/dξ^2 for all square mirror modes"""
    xi = np.asarray(xi, dtype=float)
    return np.array([square_mirror_mode_second_derivative(xi=xi, legendre_order=float(order)) for order in basis.legendre_orders], dtype=float)