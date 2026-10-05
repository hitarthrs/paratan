"""
Physical mirror eigenbasis
    builds physical pitch/orbit eigenfunctions I_j(η) for orbit averaged Lorentz operator

Physics
    Orbit averaged Lorentz operator
        <L>_z = A(Λ) d/dΛ + D(Λ) d^2/dΛ^2
    with
        A(Λ) = 4 * <1/B_tilde>_z,n - 6Λ * <1>_z,n
        D(Λ) = 4Λ * [<1/B_tilde>_z,n - Λ * <1>_z,n]
    where Eq 42 supplies the normalized density weighting
    physical eigenproblem is
        <L>_z * I_j = -Λ_j * I_j
    physical eigenfunctions are expanded in a mapped square mirror reference basis
        I_j(η) = sum_m * V_jm * M_m(η)

Projection
    have
        I = sum_m * V_m * M_m
    project
        <L> * M_m
    onto M_l over η
        K_lm = integral_0^η_TP M_l(η) * [<L> * M_m](η) dη

        G_lm = integral_0^η_TP M_l(η) * M_m(η) dη
    The production Eq 42 form uses the conservative weak projection
        K_lm = -integral D(Λ(η)) dM_l/dΛ dM_m/dΛ dη

    solve the generalized eigenproblem

        K * V = -λ * G * V
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eig
from source_model_revamp.orbits.eta_mapping import EtaLambdaMap
from source_model_revamp.scattering.lorentz_operator import LorentzLambdaCoefficients, require_orbit_averaged_lorentz_coefficients
from source_model_revamp.scattering.square_mirror_basis import SquareMirrorBasis, square_mirror_basis_derivatives, square_mirror_basis_second_derivatives, square_mirror_basis_values

LORENTZ_GALERKIN_FORM_STRONG = "strong_eq41"
LORENTZ_GALERKIN_FORM_CONSERVATIVE_WEAK = "conservative_weak_eq39_eq42"

@dataclass(frozen=True)
class PhysicalEigenbasis:
    """
    Physical mirror eigenbasis

    Variables
        eta_grid
            Grid on 0 <= η <= η_TP
        lambda_grid
            Λ(eta_grid)
        eigenvalues
            Positive λ_j values satisfying <L> * I_j = -λ_j * I_j
        expansion_coefficients
            V_jm coefficients in I_j = sum_m V_jm * M_m
        eigenfunctions
            I_j(eta_grid)
        eigenfunction_derivatives_eta
            dI_j/dη on eta_grid
        eta_trapped_passing
            η_TP
        xi_trapped_passing
            ξ_TP from the square-mirror reference basis
        galerkin_matrix
            K_lm = integral M_l * L[M_m] dη
        mass_matrix
            G_lm = integral M_l * M_m dη
        eta_monotonic
            True when both the full eta-Lambda map and trapped eta grid have the required orientation
        eigenpair_relative_residuals
            Algebraic residuals of K v_j = -lambda_j G v_j
        galerkin_symmetry_relative_error
            Frobenius relative asymmetry of the projected operator
        gram_symmetry_relative_error and gram_condition_number
            Symmetry and conditioning diagnostics for the mass matrix
        eigenfunction_eta_norms and eigenfunction_overlap_matrix
            L2 normalization and orthogonality diagnostics on the trapped eta interval
    """
    eta_grid: np.ndarray
    lambda_grid: np.ndarray
    eigenvalues: np.ndarray
    expansion_coefficients: np.ndarray
    eigenfunctions: np.ndarray
    eigenfunction_derivatives_eta: np.ndarray
    eta_trapped_passing: float
    xi_trapped_passing: float
    galerkin_matrix: np.ndarray
    mass_matrix: np.ndarray
    eta_monotonic: bool = False
    eigenpair_relative_residuals: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    max_eigenpair_relative_residual: float = np.nan
    galerkin_symmetry_relative_error: float = np.nan
    gram_symmetry_relative_error: float = np.nan
    gram_condition_number: float = np.nan
    eigenfunction_eta_norms: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    eigenfunction_overlap_matrix: np.ndarray = field(default_factory=lambda: np.empty((0, 0), dtype=float))
    eigenfunction_boundary_residuals: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    eigenfunction_normalization_residuals: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    eigenfunction_orthogonality_residuals: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    max_eigenfunction_offdiagonal_overlap: float = np.nan
    normalization_max_abs_error: float = np.nan
    boundary_max_abs_error: float = np.nan
    galerkin_form: str = LORENTZ_GALERKIN_FORM_STRONG

def xi_basis_from_eta(eta: ArrayLike, eta_trapped_passing: float, xi_trapped_passing: float) -> np.ndarray:
    """Map physical η to the square mirror reference coordinate ξ"""
    eta = np.asarray(eta, dtype=float)
    return xi_trapped_passing * eta / eta_trapped_passing

def mapped_square_basis_values_on_eta(eta: ArrayLike, square_basis: SquareMirrorBasis, eta_trapped_passing: float) -> np.ndarray:
    """
    Evaluate mapped square basis values
        M_m(η) = M_m,ξ(ξ_TP η / η_TP)
    """
    eta = np.asarray(eta, dtype=float)
    xi = xi_basis_from_eta(eta=eta, eta_trapped_passing=eta_trapped_passing, xi_trapped_passing=square_basis.xi_trapped_passing)

    return square_mirror_basis_values(xi=xi, basis=square_basis)

def mapped_square_basis_first_eta_derivatives(eta: ArrayLike, square_basis: SquareMirrorBasis, eta_trapped_passing: float) -> np.ndarray:
    """
    Evaluate dM_m/dη for mapped square basis functions
        ξ = ξ_TP η / η_TP
        dM/dη = dM/dξ * ξ_TP / η_TP
    """
    eta = np.asarray(eta, dtype=float)
    xi_scale = square_basis.xi_trapped_passing / eta_trapped_passing
    xi = xi_basis_from_eta(eta=eta, eta_trapped_passing=eta_trapped_passing, xi_trapped_passing=square_basis.xi_trapped_passing)

    return xi_scale * square_mirror_basis_derivatives(xi=xi, basis=square_basis)

def mapped_square_basis_second_eta_derivatives(eta: ArrayLike, square_basis: SquareMirrorBasis, eta_trapped_passing: float) -> np.ndarray:
    """Evaluate d^2M_m/dη^2 for mapped square basis functions"""
    eta = np.asarray(eta, dtype=float)
    xi_scale = square_basis.xi_trapped_passing / eta_trapped_passing
    xi = xi_basis_from_eta(eta=eta, eta_trapped_passing=eta_trapped_passing, xi_trapped_passing=square_basis.xi_trapped_passing)

    return xi_scale**2 * square_mirror_basis_second_derivatives(xi=xi, basis=square_basis)

def d2_eta_dlambda2_grid(eta_lambda_map: EtaLambdaMap) -> np.ndarray:
    """
    Compute d^2η/dΛ^2 on the η-map Λ grid

    Since dη/dΛ = -τ_tilde_b / <1>_Λ
        d^2η/dΛ^2 = -(dτ_tilde_b/dΛ) / <1>_Λ
    """
    d_tau_dlambda = np.gradient(eta_lambda_map.tau_tilde_grid, eta_lambda_map.lambda_grid, edge_order=2)
    return -d_tau_dlambda / eta_lambda_map.lambda_normalization

def basis_lambda_derivatives_from_eta_derivatives(dM_deta: np.ndarray, d2M_deta2: np.ndarray, d_eta_dlambda_values: np.ndarray, d2_eta_dlambda2_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert basis derivatives from η to Λ
        dM/dΛ = dM/dη * dη/dΛ,
        d^2M/dΛ^2 = d^2M/dη^2 (dη/dΛ)^2 + dM/dη d^2η/dΛ^2
    """
    d_eta = np.asarray(d_eta_dlambda_values, dtype=float)
    d2_eta = np.asarray(d2_eta_dlambda2_values, dtype=float)
    dM_dlambda = dM_deta * d_eta[None, :]
    d2M_dlambda2 = d2M_deta2 * d_eta[None, :] ** 2 + dM_deta * d2_eta[None, :]

    return dM_dlambda, d2M_dlambda2

def apply_orbit_averaged_operator_to_basis(basis_dlambda: np.ndarray, basis_d2lambda2: np.ndarray, lorentz_coefficients: LorentzLambdaCoefficients) -> np.ndarray:
    """
    Apply the orbit averaged Λ space Lorentz operator to basis modes
        <L>_z[M_m] = A_avg * dM_m/dΛ + D_avg * d^2M_m/dΛ^2
    """
    require_orbit_averaged_lorentz_coefficients(lorentz_coefficients, caller="apply_orbit_averaged_operator_to_basis")
    d1 = np.asarray(basis_dlambda, dtype=float)
    d2 = np.asarray(basis_d2lambda2, dtype=float)
    if d1.shape != d2.shape:
        raise ValueError("basis_dlambda and basis_d2lambda2 must have matching shapes")
    if d1.ndim != 2:
        raise ValueError("basis derivative arrays must have shape (n_basis, n_eta)")
    if lorentz_coefficients.lambda_values.size != d1.shape[1]:
        raise ValueError("Lorentz coefficient grid length must match basis derivative grid length")
    A = lorentz_coefficients.first_derivative_coefficient
    D = lorentz_coefficients.second_derivative_coefficient

    return A[None, :] * d1 + D[None, :] * d2

def galerkin_matrix(eta_grid: ArrayLike, basis_values: np.ndarray, operator_basis_values: np.ndarray) -> np.ndarray:
    """Build K_lm = integral M_l(η) * L[M_m](η) dη"""
    eta = np.asarray(eta_grid, dtype=float)
    n_basis = basis_values.shape[0]
    K = np.empty((n_basis, n_basis), dtype=float)
    for l in range(n_basis):
        for m in range(n_basis):
            K[l, m] = trapezoid(basis_values[l, :] * operator_basis_values[m, :], eta)

    return K

def conservative_weak_galerkin_matrix(eta_grid: ArrayLike, basis_values: np.ndarray, basis_dlambda: np.ndarray, lorentz_coefficients: LorentzLambdaCoefficients) -> np.ndarray:
    """Build the conservative Eq 39 and Eq 42 weak Galerkin matrix"""
    require_orbit_averaged_lorentz_coefficients(lorentz_coefficients, caller="conservative_weak_galerkin_matrix")
    eta = np.asarray(eta_grid, dtype=float)
    values = np.asarray(basis_values, dtype=float)
    derivatives = np.asarray(basis_dlambda, dtype=float)
    if values.shape != derivatives.shape:
        raise ValueError("basis_values and basis_dlambda must have matching shapes")
    if values.ndim != 2 or values.shape[1] != eta.size:
        raise ValueError("basis arrays must have shape (n_basis, n_eta)")
    diffusion = np.asarray(lorentz_coefficients.second_derivative_coefficient, dtype=float)
    if diffusion.shape != eta.shape:
        raise ValueError("Lorentz diffusion coefficient must match eta_grid")
    scale = max(float(np.max(np.abs(diffusion))), 1.0)
    tolerance = 512.0 * np.finfo(float).eps * scale
    if float(np.min(diffusion)) < -tolerance:
        raise ValueError("orbit averaged Lorentz diffusion coefficient became negative")
    diffusion = np.where(diffusion < 0.0, 0.0, diffusion)
    n_basis = values.shape[0]
    matrix = np.empty((n_basis, n_basis), dtype=float)
    for left in range(n_basis):
        for right in range(left, n_basis):
            value = -trapezoid(diffusion * derivatives[left] * derivatives[right], eta)
            matrix[left, right] = value
            matrix[right, left] = value

    return matrix

def mass_matrix(eta_grid: ArrayLike, basis_values: np.ndarray) -> np.ndarray:
    """Build G_lm = integral M_l(η) * M_m(η) dη"""
    eta = np.asarray(eta_grid, dtype=float)
    n_basis = basis_values.shape[0]
    G = np.empty((n_basis, n_basis), dtype=float)
    for l in range(n_basis):
        for m in range(n_basis):
            G[l, m] = trapezoid(basis_values[l, :] * basis_values[m, :], eta)

    return G

def relative_matrix_symmetry_error(matrix: ArrayLike) -> float:
    """Return ||A - A.T||_F / ||A||_F with a finite zero matrix convention"""
    values = np.asarray(matrix, dtype=float)
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise ValueError("matrix must be square")
    scale = float(np.linalg.norm(values, ord="fro"))
    asymmetry = float(np.linalg.norm(values - values.T, ord="fro"))
    if scale == 0.0:
        return 0.0 if asymmetry == 0.0 else np.inf
    
    return float(asymmetry / scale)

def projected_eigenpair_relative_residuals(K: ArrayLike, G: ArrayLike, eigenvalues: ArrayLike, expansion_coefficients: ArrayLike) -> np.ndarray:
    """Return scale independent residuals for K * v_j = -λ_j * G * v_j"""
    galerkin = np.asarray(K, dtype=float)
    gram = np.asarray(G, dtype=float)
    values = np.asarray(eigenvalues, dtype=float)
    vectors = np.asarray(expansion_coefficients, dtype=float)
    if galerkin.ndim != 2 or gram.shape != galerkin.shape or galerkin.shape[0] != galerkin.shape[1]:
        raise ValueError("K and G must be square matrices with the same shape")
    if values.ndim != 1 or vectors.shape != (values.size, galerkin.shape[0]):
        raise ValueError("expansion_coefficients must have shape (n_eigenvalues, n_basis)")
    residuals = np.empty(values.size, dtype=float)
    for mode_index, (eigenvalue, vector) in enumerate(zip(values, vectors, strict=True)):
        K_vector = galerkin @ vector
        G_vector = gram @ vector
        residual = K_vector + eigenvalue * G_vector
        scale = np.linalg.norm(K_vector) + abs(eigenvalue) * np.linalg.norm(G_vector)
        residuals[mode_index] = np.linalg.norm(residual) / max(float(scale), np.finfo(float).tiny)

    return residuals

def eigenfunction_inner_product_matrix(eta_grid: ArrayLike, eigenfunctions: ArrayLike) -> np.ndarray:
    """Return H_ij = integral I_i(η) I_j(η) dη"""
    eta = np.asarray(eta_grid, dtype=float)
    functions = np.asarray(eigenfunctions, dtype=float)
    if eta.ndim != 1 or functions.ndim != 2 or functions.shape[1] != eta.size:
        raise ValueError("eigenfunctions must have shape (n_modes, n_eta)")
    products = functions[:, None, :] * functions[None, :, :]

    return np.asarray(trapezoid(products, eta, axis=2), dtype=float)

def normalized_eigenfunction_overlap_matrix(eta_grid: ArrayLike, eigenfunctions: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    """Return L2 norms and the normalized eigenfunction overlap matrix"""
    inner_products = eigenfunction_inner_product_matrix(eta_grid=eta_grid, eigenfunctions=eigenfunctions)
    norm_squared = np.diag(inner_products)
    if np.any(~np.isfinite(norm_squared)) or np.any(norm_squared <= 0.0):
        raise ValueError("eigenfunctions must have finite positive eta space norms")
    norms = np.sqrt(norm_squared)
    overlaps = inner_products / (norms[:, None] * norms[None, :])

    return np.asarray(norms, dtype=float), np.asarray(overlaps, dtype=float)

def solve_projected_eigenproblem(K: np.ndarray, G: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Solve K * V = -λ * G * V and return positive λ values"""
    K = np.asarray(K, dtype=float)
    G = np.asarray(G, dtype=float)
    if K.ndim != 2 or G.ndim != 2 or K.shape[0] != K.shape[1] or G.shape != K.shape:
        raise ValueError("K and G must be square matrices with the same shape")
    if np.any(~np.isfinite(K)) or np.any(~np.isfinite(G)):
        raise ValueError("K and G must contain only finite values")
    symmetric_G = 0.5 * (G + G.T)
    gram_eigenvalues = np.linalg.eigvalsh(symmetric_G)
    if np.any(gram_eigenvalues <= 0.0):
        raise ValueError("G must be positive definite")
    raw_values, raw_vectors = eig(K, G)
    raw_values = np.real_if_close(raw_values)
    raw_vectors = np.real_if_close(raw_vectors)
    if np.iscomplexobj(raw_values) or np.iscomplexobj(raw_vectors):
        raise ValueError("Projected eigenproblem produced non negligible complex values")
    eigenvalues = -np.asarray(raw_values, dtype=float)
    vectors = np.asarray(raw_vectors, dtype=float)
    order = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[order]
    vectors = vectors[:, order].T
    positive = eigenvalues > 0.0
    eigenvalues = eigenvalues[positive]
    vectors = vectors[positive, :]

    return eigenvalues, vectors

def normalize_expansion_vectors_at_eta_zero(vectors: np.ndarray, basis_values: np.ndarray) -> np.ndarray:
    """Normalize expansion vectors so the reconstructed I_j(η=0) equals 1"""
    basis_at_zero = np.asarray(basis_values, dtype=float)[:, 0]
    normalized = np.asarray(vectors, dtype=float).copy()
    for j in range(normalized.shape[0]):
        amplitude = float(np.dot(normalized[j, :], basis_at_zero))
        if not np.isfinite(amplitude) or abs(amplitude) <= np.finfo(float).tiny:
            raise ValueError("projected eigenvector cannot be normalized at eta=0")
        normalized[j, :] = normalized[j, :] / amplitude
        
    return normalized

def reconstruct_physical_eigenfunctions(expansion_coefficients: np.ndarray, basis_values: np.ndarray) -> np.ndarray:
    """Reconstruct I_j(η) = sum_m * V_jm * M_m(η)"""
    return expansion_coefficients @ basis_values

def reconstruct_physical_eigenfunction_derivatives_eta(expansion_coefficients: np.ndarray, basis_first_eta_derivatives: np.ndarray) -> np.ndarray:
    """Reconstruct dI_j/dη = sum_m * V_jm * dM_m/dη"""
    return expansion_coefficients @ basis_first_eta_derivatives

def build_physical_eigenbasis(eta_grid: ArrayLike, eta_trapped_passing: float, eta_lambda_map: EtaLambdaMap, lorentz_coefficients: LorentzLambdaCoefficients, square_basis: SquareMirrorBasis, n_modes: int, galerkin_form: str = LORENTZ_GALERKIN_FORM_STRONG) -> PhysicalEigenbasis:
    """
    Build the physical mirror Lorentz eigenbasis
        lorentz_coefficients must be the orbit averaged Eq 41 coefficients evaluated on Λ(η_grid)
    """
    if int(n_modes) <= 0:
        raise ValueError("n_modes must be positive")
    require_orbit_averaged_lorentz_coefficients(lorentz_coefficients, caller="build_physical_eigenbasis")
    eta = np.asarray(eta_grid, dtype=float)
    if eta.ndim != 1 or eta.size < 3 or np.any(~np.isfinite(eta)):
        raise ValueError("eta_grid must be a finite 1D array containing at least three points")
    if not np.isclose(eta[0], 0.0, rtol=0.0, atol=1.0e-14):
        raise ValueError("eta_grid must start at zero")
    if not np.isclose(eta[-1], eta_trapped_passing, rtol=1.0e-12, atol=1.0e-14):
        raise ValueError("eta_grid must end at eta_trapped_passing")
    physical_eta_monotonic = bool(np.all(np.diff(eta) > 0.0))
    mapping_eta = np.asarray(eta_lambda_map.eta_grid, dtype=float)
    mapping_lambda = np.asarray(eta_lambda_map.lambda_grid, dtype=float)
    mapping_eta_monotonic = bool(mapping_eta.ndim == 1 and mapping_lambda.shape == mapping_eta.shape and np.all(np.isfinite(mapping_eta)) and np.all(np.isfinite(mapping_lambda)) and np.all(np.diff(mapping_lambda) > 0.0) and np.all(np.diff(mapping_eta) < 0.0))
    eta_monotonic = bool(physical_eta_monotonic and mapping_eta_monotonic)
    if not physical_eta_monotonic:
        raise ValueError("eta_grid must be strictly increasing")
    if not mapping_eta_monotonic:
        raise ValueError("eta_lambda_map must be strictly decreasing in eta as Lambda increases")
    Lambda = eta_lambda_map.lambda_of_eta(eta)
    if np.any(~np.isfinite(Lambda)) or not np.all(np.diff(Lambda) < 0.0):
        raise ValueError("Lambda(eta_grid) must be finite and strictly decreasing")
    if lorentz_coefficients.lambda_values.shape != Lambda.shape:
        raise ValueError("Lorentz coefficient grid must match Λ(eta_grid) shape")
    if not np.allclose(lorentz_coefficients.lambda_values, Lambda, rtol=1.0e-10, atol=1.0e-12):
        raise ValueError("Lorentz coefficients must be evaluated on Λ(eta_grid)")
    d_eta = eta_lambda_map.d_eta_d_lambda(Lambda)
    d2_eta_grid_lambda = d2_eta_dlambda2_grid(eta_lambda_map=eta_lambda_map)
    d2_eta_interp = PchipInterpolator(eta_lambda_map.lambda_grid, d2_eta_grid_lambda, extrapolate=False)
    d2_eta = d2_eta_interp(Lambda)
    basis_values = mapped_square_basis_values_on_eta(eta=eta, square_basis=square_basis, eta_trapped_passing=eta_trapped_passing)
    basis_values[:, -1] = 0.0
    basis_deta = mapped_square_basis_first_eta_derivatives(eta=eta, square_basis=square_basis, eta_trapped_passing=eta_trapped_passing)
    basis_d2eta2 = mapped_square_basis_second_eta_derivatives(eta=eta, square_basis=square_basis, eta_trapped_passing=eta_trapped_passing)
    basis_dlambda, basis_d2lambda2 = basis_lambda_derivatives_from_eta_derivatives(dM_deta=basis_deta, d2M_deta2=basis_d2eta2, d_eta_dlambda_values=d_eta, d2_eta_dlambda2_values=d2_eta)
    selected_galerkin_form = str(galerkin_form)
    if selected_galerkin_form == LORENTZ_GALERKIN_FORM_CONSERVATIVE_WEAK:
        K = conservative_weak_galerkin_matrix(eta_grid=eta, basis_values=basis_values, basis_dlambda=basis_dlambda, lorentz_coefficients=lorentz_coefficients)
    elif selected_galerkin_form == LORENTZ_GALERKIN_FORM_STRONG:
        operator_basis_values = apply_orbit_averaged_operator_to_basis(basis_dlambda=basis_dlambda, basis_d2lambda2=basis_d2lambda2, lorentz_coefficients=lorentz_coefficients)
        K = galerkin_matrix(eta_grid=eta, basis_values=basis_values, operator_basis_values=operator_basis_values)
    else:
        raise ValueError("galerkin_form must be strong_eq41 or conservative_weak_eq39_eq42")
    G = mass_matrix(eta_grid=eta, basis_values=basis_values)
    eigenvalues, vectors = solve_projected_eigenproblem(K=K, G=G)
    if eigenvalues.size < n_modes:
        raise ValueError(f"Requested {n_modes} physical eigenmodes, but only " f"{eigenvalues.size} positive projected eigenvalues were found")
    vectors = vectors[:n_modes, :]
    eigenvalues = eigenvalues[:n_modes]
    vectors = normalize_expansion_vectors_at_eta_zero(vectors=vectors, basis_values=basis_values)
    eigenfunctions = reconstruct_physical_eigenfunctions(expansion_coefficients=vectors, basis_values=basis_values)
    eigenfunction_derivatives = reconstruct_physical_eigenfunction_derivatives_eta(expansion_coefficients=vectors, basis_first_eta_derivatives=basis_deta)
    eigenpair_residuals = projected_eigenpair_relative_residuals(K=K, G=G, eigenvalues=eigenvalues, expansion_coefficients=vectors)
    eigenfunction_norms, eigenfunction_overlaps = normalized_eigenfunction_overlap_matrix(eta_grid=eta, eigenfunctions=eigenfunctions)
    offdiagonal_overlaps = eigenfunction_overlaps - np.eye(eigenfunction_overlaps.shape[0])
    boundary_residuals = np.abs(eigenfunctions[:, -1])
    normalization_residuals = np.abs(eigenfunctions[:, 0] - 1.0)
    orthogonality_residuals = np.max(np.abs(offdiagonal_overlaps), axis=1)

    return PhysicalEigenbasis(
        eta_grid=eta,
        lambda_grid=Lambda,
        eigenvalues=eigenvalues,
        expansion_coefficients=vectors,
        eigenfunctions=eigenfunctions,
        eigenfunction_derivatives_eta=eigenfunction_derivatives,
        eta_trapped_passing=float(eta_trapped_passing),
        xi_trapped_passing=float(square_basis.xi_trapped_passing),
        galerkin_matrix=K,
        mass_matrix=G,
        eta_monotonic=eta_monotonic,
        eigenpair_relative_residuals=eigenpair_residuals,
        max_eigenpair_relative_residual=float(np.max(eigenpair_residuals)),
        galerkin_symmetry_relative_error=relative_matrix_symmetry_error(K),
        gram_symmetry_relative_error=relative_matrix_symmetry_error(G),
        gram_condition_number=float(np.linalg.cond(G)),
        eigenfunction_eta_norms=eigenfunction_norms,
        eigenfunction_overlap_matrix=eigenfunction_overlaps,
        eigenfunction_boundary_residuals=boundary_residuals,
        eigenfunction_normalization_residuals=normalization_residuals,
        eigenfunction_orthogonality_residuals=orthogonality_residuals,
        max_eigenfunction_offdiagonal_overlap=float(np.max(np.abs(offdiagonal_overlaps))),
        normalization_max_abs_error=float(np.max(np.abs(eigenfunctions[:, 0] - 1.0))),
        boundary_max_abs_error=float(np.max(np.abs(eigenfunctions[:, -1]))),
        galerkin_form=selected_galerkin_form,
    )