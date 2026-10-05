#!/usr/bin/env python3
"""
fbis_orbit_averaged.py
======================

Scaffolding for the **orbit-averaged** pitch eigenproblem in Egedal et al.
2022, **§4** (not §3 — §3 is the square-well reduced reactor model).

Paper map
---------
- Eq. (56)  B(z) model                          → `Bfunc` (from fbis_core)
- Eqs. (37),(47)  bounce time τ̃_b, η(Λ)         → `bounce_averages`
- Eq. (40)  local L in Λ; bounce-averaged L_z   → `orbit_averaged_L_coeffs`
- Eqs. (50)–(54)  Galerkin H → I_j(η), λ_j      → `OrbitAveragedMirror`
- Eq. (57)  n(z)/⟨n⟩ from I₁(Λ)                 → `density_from_I1`
- Fig. 10(f)  n(z)-weighted scattering average  → `density_weight=True`

Relation to the neutron-source pipeline
---------------------------------------
`fbis_core.build_midplane_f` uses **square-well / midplane** M_j(ξ).
This module builds **I_j(η)** on a real B(z).  Fig. 10 is the validation
target for this scaffolding; wiring I_j into the Sn mapper is a later step.

Eigenvalue convention
---------------------
The discrete Galerkin matrix H approximates L_z on the {M_m(η)} basis with
M(η_TP)=0.  Eigenvalues of H are ≈ −λ_j · Δη (see get_H in fbis_core).
We report λ_j = −Re(μ_j) / Δη so that the square-well limit is comparable
to λ₁^□ = ℓ₁(ℓ₁+1) from the midplane Legendre solve.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.interpolate import interp1d, make_interp_spline

from fbis_core import (
    Bfunc,
    find_eigenvalues_from_data,
    get_H,
    get_LM,
    get_M_basis,
    get_zTP,
    scan_M_table,
)

logger = logging.getLogger(__name__)


# =============================================================================
# Square-well reference λ₁ (Fig. 10 denominator)
# =============================================================================


def lambda1_square(Rm: float, *, lmax: float = 25.3, nxi: int = 80) -> float:
    """
    Square-mirror λ₁ = ℓ₁(ℓ₁+1) with M(ξ*)=0, ξ*=√(1-1/R_m).

    This is λ_{1,square} in Fig. 10(d,f).  Asymptotically ~ 2/ln(R_m).
    """
    x_star = float(np.sqrt(1.0 - 1.0 / Rm))
    lax, xax, data = scan_M_table(lmax=lmax, nxi=nxi)
    eigs = find_eigenvalues_from_data(x_star, xax, lax, data)
    if len(eigs) == 0:
        raise RuntimeError(f"No square-well eigenvalues for Rm={Rm}")
    l1 = float(eigs[0])
    return l1 * (l1 + 1.0)


def lambda1_square_approx(Rm: float) -> float:
    """Paper fit near Fig. 2(a): 2/ln(R_m) + 0.37/[ln(R_m)]^{1.3}."""
    lnR = float(np.log(Rm))
    return 2.0 / lnR + 0.37 * (1.0 / lnR) ** 1.3


# =============================================================================
# Bounce averages on a given B(z)
# =============================================================================


@dataclass
class BounceAverages:
    """τ̃_b(Λ), ⟨1/B⟩(Λ), η(Λ) for one (R_m, w, k) field."""

    Rm: float
    w: float
    k: float
    lax: np.ndarray
    z_TP: np.ndarray
    tau: np.ndarray
    avg_invB: np.ndarray
    eta: np.ndarray
    eta_TP: float
    Itau: float


def _integrate_bounce(
    Lambda: np.ndarray,
    B: Callable[[float], float],
    z_TP: np.ndarray,
    n_of_z: Callable[[np.ndarray], np.ndarray] | None = None,
    n_z: int = 120,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Bounce integrals for each Λ.

    Unweighted (n_of_z is None)
        τ(Λ) = ∫_0^{z_TP} dz / √(1-Λ B)
        ⟨1/B⟩ = ∫ (1/B)/√ / τ

    Density-weighted (Fig. 10f — scattering ∝ n(z))
        τ_n = ∫ n(z)/√ dz
        ⟨1/B⟩_n = ∫ (n/B)/√ / τ_n
    """
    tau = np.zeros_like(Lambda, dtype=float)
    avg = np.zeros_like(Lambda, dtype=float)
    for i, (L, Z) in enumerate(zip(Lambda, z_TP)):
        if Z <= 0.0:
            tau[i] = 0.0
            avg[i] = 1.0
            continue
        # Open endpoint: stay off the integrable turning-point singularity.
        z = np.linspace(0.0, float(Z), n_z, endpoint=False)
        if len(z) < 2:
            tau[i] = 0.0
            avg[i] = 1.0
            continue
        dz = float(z[1] - z[0])
        Bz = np.asarray(B(z), dtype=float)
        rad = np.sqrt(np.clip(1.0 - L * Bz, 1e-15, None))
        wgt = 1.0 / rad
        if n_of_z is not None:
            wgt = wgt * np.asarray(n_of_z(z), dtype=float)
        tau[i] = float(np.sum(wgt) * dz)
        avg[i] = float(np.sum((1.0 / Bz) * wgt) * dz / max(tau[i], 1e-30))
    return tau, avg


def bounce_averages(
    Rm: float = 16.0,
    w: float = 0.1,
    k: float = 1.3,
    N_Lambda: int = 200,
    n_z: int = 200,
    n_of_z: Callable[[np.ndarray], np.ndarray] | None = None,
) -> BounceAverages:
    """Build τ, ⟨1/B⟩, η(Λ)=1−∫τ/∫τ for Eq. (47)."""
    lax = np.linspace(1e-4, 1.0 - 1e-6, N_Lambda)
    Bf: Callable[[float], float] = lambda z: float(Bfunc(z, w=w, k=k, Rm=Rm))
    B_arr = lambda z: np.asarray(Bfunc(z, w=w, k=k, Rm=Rm), dtype=float)

    z_TP = get_zTP(lax, Bf)
    tau, avg = _integrate_bounce(lax, B_arr, z_TP, n_of_z=n_of_z, n_z=n_z)
    Itau = float(np.trapezoid(tau, lax))
    from scipy.integrate import cumulative_trapezoid

    integ = np.concatenate([[0.0], cumulative_trapezoid(tau, lax)])
    eta = 1.0 - integ / max(Itau, 1e-30)
    eta_TP = float(np.interp(1.0 / Rm, lax, eta))
    return BounceAverages(
        Rm=Rm, w=w, k=k, lax=lax, z_TP=z_TP, tau=tau,
        avg_invB=avg, eta=eta, eta_TP=eta_TP, Itau=Itau,
    )


# =============================================================================
# Orbit-averaged eigenproblem → I₁, λ₁, n(z)
# =============================================================================


@dataclass
class OrbitAveragedEigen:
    """First orbit-averaged eigenmode and density for one field."""

    Rm: float
    w: float
    k: float
    density_weighted: bool
    bounce: BounceAverages
    zax: np.ndarray
    B: np.ndarray
    eta_ax: np.ndarray
    I1: np.ndarray          # I₁(η) on eta_ax, I₁(0)=1
    l_ax: np.ndarray        # Λ(η) on eta_ax
    lambda1: float
    lambda1_sq: float
    ratio: float            # λ₁ / λ₁,square
    ne: np.ndarray          # n(z)/⟨n⟩ Eq. (57)


def density_from_I1(
    bounce: BounceAverages,
    l_ax: np.ndarray,
    I1: np.ndarray,
    zax: np.ndarray,
    B: np.ndarray,
) -> np.ndarray:
    """Eq. (57): n(z)/⟨n⟩ from bounce-weighted I₁(Λ)."""
    f_tau = make_interp_spline(bounce.lax, bounce.tau, k=3)
    f_I1 = make_interp_spline(l_ax[::-1], I1[::-1], k=3)
    L = np.linspace(1.0 / bounce.Rm, 1.0, 100)
    dL = float(np.diff(L).mean())
    Inorm = float(np.sum(f_tau(L) * f_I1(L)) * dL)
    norm = bounce.Itau / max(Inorm, 1e-30)

    ne = np.zeros_like(zax, dtype=float)
    for i, Bt in enumerate(B):
        L_max = max(1.0 / bounce.Rm + 1e-6, 1.0 / Bt - 1e-3)
        if L_max <= 1.0 / bounce.Rm:
            continue
        Lloc = np.linspace(1.0 / bounce.Rm, L_max, 80)
        dLl = float(np.diff(Lloc).mean())
        rad = np.sqrt(np.clip(1.0 - Bt * Lloc, 1e-12, None))
        ne[i] = float(np.sum(Bt * f_I1(Lloc) / rad) * dLl * norm)
    return ne


def solve_orbit_averaged_I1(
    Rm: float = 16.0,
    w: float = 0.1,
    k: float = 1.3,
    Nz: int = 401,
    N_Lambda: int = 160,
    N_eta: int = 300,
    density_weight: bool = False,
    lambda1_sq: float | None = None,
) -> OrbitAveragedEigen:
    """
    Galerkin solve for I₁(η) and λ₁ on analytic B(z).

    If density_weight=True (Fig. 10f path):
      1) unweighted solve → n(z)
      2) rebuild ⟨1/B⟩ with n(z) weight
      3) re-solve H → λ₁, I₁, n(z)
    """
    zax = np.linspace(0.0, 1.0, Nz)
    B = np.asarray(Bfunc(zax, w=w, k=k, Rm=Rm), dtype=float)

    # --- pass 1: unweighted bounce averages ---
    b0 = bounce_averages(Rm=Rm, w=w, k=k, N_Lambda=N_Lambda)
    if lambda1_sq is None:
        lambda1_sq = lambda1_square(Rm)

    def _galerkin(bounce: BounceAverages) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        from scipy.linalg import eig as geig

        spl_B = make_interp_spline(bounce.lax, bounce.avg_invB, k=3)
        M_basis, eta_ax = get_M_basis(bounce.eta_TP, N_eta=N_eta)
        # Λ(η): bounce.eta decreases with Λ; build a strictly increasing spline
        order = np.argsort(bounce.eta)
        eta_u = bounce.eta[order]
        lax_u = bounce.lax[order]
        uniq = np.concatenate([[True], np.diff(eta_u) > 1e-12])
        l_of_eta = make_interp_spline(eta_u[uniq], lax_u[uniq], k=3)
        l_ax = np.clip(np.asarray(l_of_eta(eta_ax), dtype=float), 1e-6, 1.0 - 1e-6)
        # Ensure strictly increasing Λ sample for d/dΛ splines in get_LM
        # (η↑ ⇒ Λ↓, so reverse sort for the spline inside get_LM).
        for i in range(1, len(l_ax)):
            if l_ax[i] >= l_ax[i - 1]:
                l_ax[i] = l_ax[i - 1] - 1e-12

        LM = get_LM(M_basis, spl_B, l_ax)
        # Proper generalized eigenproblem: A v = μ B v with μ = −λ
        # (the old H = ⟨M_j, L M_k⟩/||M_j||² is row-scaled and biases |λ|).
        deta = float(np.diff(eta_ax).mean())
        A = (M_basis @ LM.T) * deta
        Bmat = (M_basis @ M_basis.T) * deta
        mu, eigvec = geig(A, Bmat)
        mu = np.real(mu)
        phys = mu < -1e-8
        if not np.any(phys):
            raise RuntimeError("No negative eigenvalues for orbit-averaged L")
        ik = int(np.where(phys)[0][np.argmax(mu[phys])])  # least-negative ⇒ λ₁
        I1 = np.real(eigvec[:, ik] @ M_basis)
        if abs(I1[0]) > 0:
            I1 = I1 / I1[0]
        lam1 = float(-mu[ik])
        return eta_ax, I1, l_ax, lam1

    eta_ax, I1, l_ax, lam1 = _galerkin(b0)
    ne = density_from_I1(b0, l_ax, I1, zax, B)
    bounce_used = b0
    weighted = False

    if density_weight:
        # Interpolate n(z) for weighted bounce integrals
        n_spl = interp1d(zax, np.clip(ne, 1e-12, None), kind="linear",
                         fill_value=(ne[0], ne[-1]), bounds_error=False)
        bn = bounce_averages(
            Rm=Rm, w=w, k=k, N_Lambda=N_Lambda,
            n_of_z=lambda z: np.asarray(n_spl(z), dtype=float),
        )
        # Keep η from the unweighted geometric definition; only ⟨1/B⟩ is
        # replaced for the scattering-frequency weighting (Fig. 10f).
        bn_mix = BounceAverages(
            Rm=bn.Rm, w=bn.w, k=bn.k, lax=bn.lax, z_TP=bn.z_TP,
            tau=b0.tau, avg_invB=bn.avg_invB, eta=b0.eta,
            eta_TP=b0.eta_TP, Itau=b0.Itau,
        )
        eta_ax, I1, l_ax, lam1 = _galerkin(bn_mix)
        ne = density_from_I1(bn_mix, l_ax, I1, zax, B)
        bounce_used = bn_mix
        weighted = True

    return OrbitAveragedEigen(
        Rm=Rm, w=w, k=k, density_weighted=weighted,
        bounce=bounce_used, zax=zax, B=B,
        eta_ax=eta_ax, I1=I1, l_ax=l_ax,
        lambda1=lam1, lambda1_sq=float(lambda1_sq),
        ratio=lam1 / float(lambda1_sq),
        ne=ne,
    )


def scan_lambda_ratio_vs_Rm(
    w_values: list[float],
    Rm_values: np.ndarray,
    *,
    density_weight: bool = False,
    k: float = 1.3,
    Nz: int = 201,
    N_Lambda: int = 100,
    N_eta: int = 200,
) -> dict[float, tuple[np.ndarray, np.ndarray]]:
    """
    Fig. 10(d) or (f): λ₁/λ₁,square vs log10(R_m) for each w.

    Returns {w: (log10_Rm, ratio_array)}.
    """
    out: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    # Cache square-well λ₁(Rm)
    lam_sq_cache: dict[float, float] = {}

    for w in w_values:
        ratios = []
        logs = []
        for Rm in Rm_values:
            Rm = float(Rm)
            if Rm not in lam_sq_cache:
                lam_sq_cache[Rm] = lambda1_square(Rm)
            logger.info(
                "  w=%.3f Rm=%.3g density_weight=%s ...",
                w, Rm, density_weight,
            )
            eig = solve_orbit_averaged_I1(
                Rm=Rm, w=w, k=k, Nz=Nz, N_Lambda=N_Lambda, N_eta=N_eta,
                density_weight=density_weight,
                lambda1_sq=lam_sq_cache[Rm],
            )
            logs.append(np.log10(Rm))
            ratios.append(eig.ratio)
            logger.info("    λ1/λ1□ = %.4f", eig.ratio)
        out[w] = (np.asarray(logs), np.asarray(ratios))
    return out
