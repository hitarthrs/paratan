#!/usr/bin/env python3
"""
fbis_core.py
============

Shared physics for the FBIS (Fast Beam Ion Solver) analytic / semi-analytic
model, extracted from FBISpy.ipynb / fbis_walkthrough.py.

Paper
-----
J. Egedal, D. Endrizzi, C.B. Forest, T.K. Fowler,
"Fusion by beam ions in a low collisionality, high mirror ratio magnetic mirror",
Nucl. Fusion 62 126053 (2022).

Notes
-----
Tony Qian, latex/fbis-notes.tex — roadmap for rederiving the paper in Python.
IMPORTANT note from those notes: Figure 3 of the paper cannot be reproduced from
the j=1 eigenmode alone; the NBI source S(ξ) must be projected onto the full
{M_j} basis (the Matlab / notebook does this with a Gaussian S).

What this module is for
-----------------------
Build a **normalized axial neutron source** S_n(z) by:

  1. Solving the §2 midplane Fokker–Planck eigen-expansion for f(v, ξ)
     including multi-mode NBI pitch (injection angle).
  2. Mapping that f along an analytic B(z) using the constant of motion
     Λ = μ B_0 / E  (paper §4), with exact per-cell Λ integration of the
     axial Jacobian.
  3. Deterministic beam–beam DT reactivity sum at each z (species-aware
     D/T speed scaling; exact gyrophase / bounce-sign averaging).

This is intentionally **not** full FBIS §5 (no ambipolar Φ, no Rosenbluth
energy diffusion iteration, no separate D/T solvers). Absolute Amperes →
n [m^-3] scaling is also deferred; outputs are **shapes**.

Coordinate conventions (read this before changing units)
--------------------------------------------------------
- ξ = v_parallel / v          (pitch cosine), midplane.
- Λ = 1 - ξ^2 = μ B_0 / E     at midplane where B=B_0.
- At general z, B̃ = B(z)/B_0, an orbit with Λ reaches z iff Λ < 1/B̃,
  and locally  ξ_local^2 = 1 - Λ B̃.
- The walkthrough / notebook label the speed axis with the same numbers as
  the beam energy in keV (v0=60 ↔ E_beam=60 keV). We keep that convention:
      E_i / E_beam = (v_i / v0)^2
  so numerical `v` is a **normalized speed**, not keV itself.

Commenting convention
---------------------
Inline notes call out why a formula exists, what would break if changed,
and which paper equation is being tracked.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import mpmath
import numpy as np
from scipy.interpolate import interp1d, make_interp_spline
from scipy.optimize import brentq

logger = logging.getLogger(__name__)

OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"


# =============================================================================
# Part A — Legendre pitch eigenfunctions M_l(ξ)   [paper §2]
# =============================================================================
#
# The Lorentz pitch-angle operator L (paper Eq. 4) is the ξ-part of a
# spherical Laplacian. Eigenfunctions of L on a mirror are NOT ordinary
# integer Legendre P_l, because the boundary condition is M(ξ*)=0 at the
# loss-cone edge ξ* = sqrt(1 - 1/R_m), not M(±1)=finite alone.
#
# Solution: take a linear combination of Legendre P and Q of (generally
# non-integer) degree l, fix midplane BCs M(0)=1, M'(0)=0, then scan l until
# M(ξ*)=0. Eigenvalue of L is λ = l(l+1).


def make_M(l: float) -> tuple[float, float, Callable[[float], float]]:
    """
    Build M_l(x) = a P_l(x) + b Q_l(x) with M(0)=1 and M'(0)=0.

    Parameters
    ----------
    l : float
        Legendre degree (need not be integer).

    Returns
    -------
    a, b : float
        Coefficients of P and Q.
    M : callable
        M(x) for x in (-1, 1).
    """
    # mpmath is used because scipy.special.lqn is awkward for
    # non-integer degree + we need derivatives at 0 for the BC solve.
    P0 = float(mpmath.legenp(l, 0, 0).real)
    Q0 = float(mpmath.legenq(l, 0, 0).real)
    dP0 = float(mpmath.diff(lambda z: mpmath.legenp(l, 0, z), 0).real)
    dQ0 = float(mpmath.diff(lambda z: mpmath.legenq(l, 0, z), 0).real)

    det = P0 * dQ0 - dP0 * Q0
    a = dQ0 / det
    b = -dP0 / det

    def M(x: float) -> float:
        Px = float(mpmath.legenp(l, 0, x).real)
        Qx = float(mpmath.legenq(l, 0, x).real)
        return a * Px + b * Qx

    return a, b, M


def find_eigenvalues_from_data(
    x_star: float,
    xax: np.ndarray,
    lax: np.ndarray,
    data: np.ndarray,
    tol: float = 1e-10,
) -> np.ndarray:
    """
    Find degrees l where M_l(x_star)=0 (loss-cone eigencondition).

    `data` is optional context for callers that already scanned a table; we
    re-evaluate M_l(x_star) on `lax` for robust bracketing (nearest-column
    sign changes can disagree with the exact x_star when the ξ-grid is coarse).
    """
    # Always bracket using the exact boundary ξ*, not data[:, i_nearest].
    # Coarse Nξ grids otherwise throw "f(a) and f(b) must have different signs".
    col = np.array([make_M(l)[2](x_star) for l in lax])

    eigenvalues: list[float] = []
    for i in range(len(lax) - 1):
        if not (np.isfinite(col[i]) and np.isfinite(col[i + 1])):
            continue
        if col[i] * col[i + 1] >= 0:
            continue
        try:
            l_eig = brentq(
                lambda l: make_M(l)[2](x_star),
                lax[i],
                lax[i + 1],
                xtol=tol,
            )
            eigenvalues.append(l_eig)
        except ValueError:
            # Rare: table sign change but refined endpoints agree — skip.
            continue
    return np.asarray(eigenvalues)


def scan_M_table(
    lmin: float = 0.11,
    lmax: float = 25.3,
    nl: int = 159,
    nxi: int = 100,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precompute M_l(ξ) on a (l, ξ) grid for contouring / root bracketing."""
    lax = np.linspace(lmin, lmax, nl)
    xax = np.linspace(0.0, 1.0, nxi, endpoint=False)
    Ms = [make_M(l)[2] for l in lax]
    data = np.array([[m(xi) for xi in xax] for m in Ms])
    return lax, xax, data


# =============================================================================
# Part B — Midplane distribution f(v, ξ)   [paper §2, Eqs. 9–15]
# =============================================================================
#
# Steady Fokker–Planck with drag + pitch scattering + monoenergetic
# anisotropic source separates in the {M_j} basis:
#
#   f(v, ξ) = Σ_j  S_j * u(v)^{λ_j} / (v^3 + v_c^3) * M_j(ξ)     (Eq. 14)
#
# u(v) is the slowing-down factor. S_j are projections of the NBI pitch
# source S0(ξ). Injection angle enters ONLY through S0 → S_j.


def ufunc(v: np.ndarray, v0: float = 60.0, vc: float = 30.0, beta: float = 0.5) -> np.ndarray:
    """
    Slowing-down amplitude u(v) (paper Eq. 15).

    Parameters
    ----------
    v, v0 :
        Speed and injection speed in the notebook's normalized units
        (v0=60 ↔ E_beam=60 keV under E/E_beam=(v/v0)^2).
    vc :
        Critical speed (set by T_e); same units as v.
        Physical anchor: E_c ≈ 20 T_e (paper Eq. 16) ⇒ vc/v0 = sqrt(20 T_e / E_beam).
    beta :
        β_m = z_eff · m_i / (2 m_f)  (paper, below Eq. 3).  For a pure DT
        plasma with m_i = m_f this is 0.5.
        The combined pitch-depletion exponent in Eq. 14 is λ_j·β_m/3.
        A legacy default of 4.0 predates the u**l → u**λ correction and
        over-depletes the slowed-down spectrum by u^{λ·(4-0.5)/3}; do not
        restore it without a reason.
    """
    v = np.asarray(v, dtype=float)
    return (((v0**3 + vc**3) / (v**3 + vc**3)) * (v / v0) ** 3) ** (beta / 3.0)


def S0(xi: np.ndarray, R_b: float = 2.0, sigma: float = 0.05) -> np.ndarray:
    """
    Gaussian NBI source in pitch cosine.

    Peaked at ξ_b = sqrt(1 - 1/R_b).

    Injection angle (do NOT invert this):
      ξ = v∥/v = cosθ, so θ=90° (perp to B) ⇒ ξ=0,  θ=0° (parallel) ⇒ ξ=1.
      - Near-perpendicular NBI: R_b → 1⁺  (e.g. R_b=1.01 ⇒ ξ_b≈0.1)
      - ~45° pitch:           R_b = 2     (ξ_b = 1/√2 ≈ 0.707)
      - More parallel / sloshing toward loss cone: larger R_b → ξ_b → 1
    """
    xi = np.asarray(xi, dtype=float)
    if R_b <= 1.0:
        raise ValueError(f"R_b must be > 1 (got {R_b}); R_b→1⁺ is near-perpendicular")
    xi_b = np.sqrt(max(1.0 - 1.0 / R_b, 0.0))
    return np.exp(-((xi - xi_b) ** 2) / sigma**2)


@dataclass
class MidplaneDistribution:
    """Container for the §2 midplane solution."""

    vax: np.ndarray          # shape (Nv,)  normalized speed
    xax: np.ndarray          # shape (Nξ,)  pitch cosine ξ ∈ [0, 1)
    f: np.ndarray            # shape (Nv, Nξ)  f(v, ξ)  (relative units)
    eigs: np.ndarray         # Legendre degrees l_j  (λ_j = l_j(l_j+1))
    M_arr: np.ndarray        # shape (n_modes, Nξ)
    Sj: np.ndarray           # source projection coefficients
    x_star: float            # loss-cone edge ξ*
    Rm: float                # mirror ratio used for the loss cone
    R_b: float               # beam / injection pitch parameter
    v0: float
    vc: float
    beta: float
    E_beam_keV: float        # beam energy used for σ(E) conversion


def build_midplane_f(
    Rm: float = 10.26,
    R_b: float = 2.0,
    v0: float = 60.0,
    vc: float = 22.0,
    beta: float = 0.5,
    E_beam_keV: float = 60.0,
    n_modes: int = 10,
    Nv: int = 48,
    Nxi: int = 80,
    lmax: float = 25.3,
    source_sigma: float = 0.05,
) -> MidplaneDistribution:
    """
    Build multi-mode midplane f(v, ξ) for a given mirror ratio and NBI pitch.

    Steps
    -----
    1. Scan M_l(ξ), find eigenvalues at ξ* = sqrt(1 - 1/R_m).
    2. Evaluate the first `n_modes` eigenfunctions on the ξ grid.
    3. Project S0(ξ; R_b) → S_j.
    4. Sum Eq. 14 over j.
    """
    # ξ* is the trapped/passing boundary at the midplane for mirror ratio Rm.
    x_star = float(np.sqrt(1.0 - 1.0 / Rm))
    logger.info("Building midplane f: Rm=%.3f → ξ*=%.4f, R_b=%.3f (ξ_b=%.4f)",
                Rm, x_star, R_b, np.sqrt(max(1.0 - 1.0 / R_b, 0.0)))

    lax, xax_scan, data = scan_M_table(lmax=lmax, nxi=Nxi)
    # Use the scan's ξ grid for the distribution (consistent bracketing grid).
    xax = xax_scan
    eigs = find_eigenvalues_from_data(x_star, xax, lax, data)
    if len(eigs) == 0:
        raise RuntimeError(f"No pitch eigenvalues found for Rm={Rm} (ξ*={x_star})")

    n_use = min(n_modes, len(eigs))
    eigs = eigs[:n_use]
    logger.info("Using %d pitch eigenmodes; l_j[0]=%.4f … l_j[-1]=%.4f",
                n_use, eigs[0], eigs[-1])

    # Eigenfunctions on the ξ grid.
    M_arr = np.array([[make_M(eigs[j])[2](x) for x in xax] for j in range(n_use)])

    # Project Gaussian source.  S_j = <S, M_j> / ||M_j||^2
    # Notes emphasize this projection — skip it and Fig. 3 is wrong.
    # Fix: paper Eqs. 7–8, 11 integrate over [0, ξ*] ONLY. Beyond ξ*
    # M_j is analytic continuation (Q_l diverges log-ly at ξ→1); including it
    # inflates α_j by up to ~13% (j=6, Rm=16) and skews the S_j mode mixture.
    inside = xax <= x_star
    S = S0(xax, R_b=R_b, sigma=source_sigma)
    dxi = float(np.diff(xax).mean())
    alpha = np.sum(M_arr[:, inside] ** 2, axis=1) * dxi
    Sj = np.array(
        [np.sum(S[inside] * M_arr[j, inside]) * dxi / alpha[j] for j in range(n_use)]
    )

    vax = np.linspace(0.0, v0, Nv)
    # Avoid v=0 singularity in drag factor: floor tiny v.
    vax = np.maximum(vax, 1e-6 * v0)

    f = np.zeros((Nv, Nxi))
    for j in range(n_use):
        # λ_j = l_j (l_j + 1)  — paper uses λ as the L-eigenvalue.
        # The notebook/walkthrough used u**l with l the Legendre degree
        # in some places and discussed λ=l(l+1) elsewhere.  Paper Eq. 14 uses
        # u^{λ_j} with λ_j the eigenvalue of L, i.e. l(l+1).
        # We follow the PAPER here: λ = l(l+1).
        lamb = eigs[j] * (eigs[j] + 1.0)
        u_v = ufunc(vax, v0=v0, vc=vc, beta=beta) ** lamb
        radial = u_v / (vax**3 + vc**3)  # shape (Nv,)
        f += Sj[j] * np.outer(radial, M_arr[j])

    # Floor negative numerical wiggles from truncated mode sum.
    f = np.clip(f, 0.0, None)
    # Fix: f is only defined on the trapped region; the M_j values at
    # ξ > ξ* are continuation artifacts, not a loss-cone population. Zero them
    # so the midplane contour and any downstream integral are honest.
    f[:, ~inside] = 0.0

    return MidplaneDistribution(
        vax=vax,
        xax=xax,
        f=f,
        eigs=eigs,
        M_arr=M_arr,
        Sj=Sj,
        x_star=x_star,
        Rm=Rm,
        R_b=R_b,
        v0=v0,
        vc=vc,
        beta=beta,
        E_beam_keV=E_beam_keV,
    )


# =============================================================================
# Part C — Analytic B(z) and Eq. 57 density (benchmark)   [paper §4]
# =============================================================================


def Bfunc(
    z: np.ndarray | float,
    w: float = 0.2,
    k: float = 1.3,
    Rm: float = 16.0,
) -> np.ndarray | float:
    """
    Model magnetic field B(z)/B_0 (paper Eq. 56).

    Peaks near |z|=1 (mirror throat). w controls throat sharpness;
    w→0 approaches a square well.
    """
    z = np.asarray(z, dtype=float)
    e1 = np.exp(-(np.abs(np.abs(z) - 1.0) ** k) / w**k)
    e0 = np.exp(-1.0 / w**k)
    return 1.0 + (Rm - 1.0) * (e1 - e0) / (1.0 - e0)


def get_zTP(
    Lambda: np.ndarray,
    B: Callable[[float], float],
    z_max: float = 1.0,
) -> np.ndarray:
    """Turning points: B(z_TP)/B(0) = 1/Λ. Passing → z_max."""
    B0 = float(B(0.0))

    def zTP_single(L: float) -> float:
        try:
            return brentq(lambda z: B(z) / B0 - 1.0 / L, 0.0, z_max)
        except ValueError:
            return z_max

    return np.array([zTP_single(L) for L in Lambda])


def get_tau(Lambda: np.ndarray, B: Callable[[float], float], z_TP: np.ndarray) -> np.ndarray:
    """Bounce integral τ(Λ) ∝ ∫_0^{z_TP} dz / sqrt(1 - Λ B(z))."""
    tau: list[float] = []
    for L, Z in zip(Lambda, z_TP):
        z = np.linspace(0.0, Z, 100, endpoint=False)
        Bz = np.asarray(B(z), dtype=float)
        dz = float(np.diff(z).mean())
        tau.append(float(np.sum(dz / np.sqrt(np.clip(1.0 - L * Bz, 1e-15, None)))))
    return np.asarray(tau)


def get_1overB(
    Lambda: np.ndarray,
    B: Callable[[float], float],
    z_TP: np.ndarray,
    tau: np.ndarray,
) -> np.ndarray:
    """Orbit average ⟨1/B⟩_bounce."""
    out: list[float] = []
    for Z, L, t in zip(z_TP, Lambda, tau):
        z = np.linspace(0.0, Z, 100, endpoint=False)
        dz = float(np.diff(z).mean())
        Bz = np.asarray(B(z), dtype=float)
        integrand = (1.0 / Bz) / np.sqrt(np.clip(1.0 - L * Bz, 1e-15, None))
        out.append(float(dz * np.sum(integrand) / t))
    return np.asarray(out)


def get_M_basis(
    eta_TP: float,
    lmin: float = 0.1,
    lmax: float = 30.1,
    nl: int = 300,
    N_eta: int = 500,
) -> tuple[np.ndarray, np.ndarray]:
    """Eigenbasis M_j(η) with M_j(η_TP)=0 on [0, η_TP]."""
    l_scan = np.linspace(lmin, lmax, nl)
    M_funcs = [make_M(l)[2] for l in l_scan]
    M_TP = np.array([M(eta_TP) for M in M_funcs])
    M_interp = interp1d(l_scan, M_TP)
    indices = np.where(np.diff(np.sign(M_TP)))[0]
    eig_l = np.array([brentq(M_interp, l_scan[i], l_scan[i + 1]) for i in indices])
    M_funcs = [make_M(l)[2] for l in eig_l]
    eta_ax = np.linspace(0.0, eta_TP, N_eta)
    M_basis = np.array([[M(eta) for eta in eta_ax] for M in M_funcs])
    return M_basis, eta_ax


def get_LM(M_basis: np.ndarray, spl_B, l_ax: np.ndarray) -> np.ndarray:
    """Apply bounce-averaged pitch-angle operator L to each basis function."""
    arr: list[np.ndarray] = []
    for M in M_basis:
        spl = make_interp_spline(l_ax[::-1], M[::-1], k=5)
        dM = spl.derivative(nu=1)
        d2M = spl.derivative(nu=2)
        L1 = (4.0 * spl_B(l_ax) - 6.0 * l_ax) * dM(l_ax)
        L2 = 4.0 * l_ax * (spl_B(l_ax) - l_ax) * d2M(l_ax)
        arr.append(L1 + L2)
    return np.asarray(arr)


def get_H(M_basis: np.ndarray, LM: np.ndarray, eta_ax: np.ndarray) -> np.ndarray:
    """Galerkin matrix H_{jk} = ⟨M_j, L M_k⟩ / ||M_j||²."""
    alpha = np.sum(M_basis**2, axis=1)
    cM = M_basis / alpha[:, np.newaxis]
    MLM = cM[:, np.newaxis, :] * LM[np.newaxis, :, :]
    deta = float(np.diff(eta_ax).mean())
    return np.sum(MLM, axis=2) * deta


class Mirror:
    """
    Axisymmetric mirror with model B(z); builds H and Eq. 57 n(z)/⟨n⟩.

    This is the §4 *benchmark* density from I_1 alone — it does NOT
    include NBI injection angle. Use `map_distribution_along_z` for the
    source-aware density / neutron source.
    """

    def __init__(
        self,
        Rm: float = 16.0,
        Nz: int = 1001,
        N_Lambda: int = 200,
        w: float = 0.1,
        k: float = 1.3,
    ) -> None:
        zax = np.linspace(0.0, 1.0, Nz)
        lax = np.linspace(1e-3, 1.0, N_Lambda, endpoint=False)
        dL = float(np.diff(lax).mean())

        Bf = lambda z: Bfunc(z, w=w, k=k, Rm=Rm)
        B = np.asarray(Bf(zax), dtype=float)
        z_TP = get_zTP(lax, Bf)
        tau = get_tau(lax, Bf, z_TP)
        Itau = float(tau.sum() * dL)

        avgB = get_1overB(lax, Bf, z_TP, tau)
        spl_B = make_interp_spline(lax, avgB, k=3)

        eta = 1.0 - tau.cumsum() * dL / Itau
        eta_TP = float(np.interp(1.0 / Rm, lax, eta))

        M_basis, eta_ax = get_M_basis(eta_TP)
        l_of_eta = make_interp_spline(eta[::-1], lax[::-1], k=3)
        l_ax = l_of_eta(eta_ax)
        LM = get_LM(M_basis, spl_B, l_ax)
        H = get_H(M_basis, LM, eta_ax)

        eigval, eigvec = np.linalg.eig(H)
        # eig may return complex dtype with ~0 imag; order on the real
        # part explicitly. Largest (least negative) eigenvalue ↔ smallest λ₁.
        # NOTE: self.eigval carries an extra factor Δη from the discrete α
        # normalization in get_H — fine for ordering/eigvecs, but do NOT read
        # these values as -λ_j without dividing by Δη.
        ik = int(np.argmax(eigval.real))
        Ij = (eigvec[:, ik][:, np.newaxis] * M_basis).sum(axis=0)

        self.Rm = Rm
        self.w = w
        self.k = k
        self.zax = zax
        self.lax = lax
        self.B = B
        self.eta = eta
        self.tau = tau
        self.Itau = Itau
        self.eta_TP = eta_TP
        self.eta_ax = eta_ax
        self.l_ax = np.asarray(l_ax)
        self.H = H
        self.eigval = eigval
        self.M_basis = M_basis
        self.I1 = Ij / Ij[0]
        self.calc_dens()

    def calc_dens(self) -> None:
        """n(z)/⟨n⟩ from bounce-weighted I₁(Λ) — paper Eq. 57."""
        f_tau = make_interp_spline(self.lax, self.tau, k=3)
        f_I1 = make_interp_spline(self.l_ax[::-1], self.I1[::-1], k=3)

        L = np.linspace(1.0 / self.Rm, 1.0, 100)
        dL = float(np.diff(L).mean())
        Inorm = float(np.sum(f_tau(L) * f_I1(L)) * dL)
        norm = self.Itau / Inorm

        ne: list[float] = []
        for B in self.B:
            L_max = max(1.0 / self.Rm + 1e-6, 1.0 / B - 1e-3)
            Lloc = np.linspace(1.0 / self.Rm, L_max, 100)
            dLl = float(np.diff(Lloc).mean())
            radicand = np.clip(1.0 - B * Lloc, 1e-12, None)
            I = B * f_I1(Lloc) / np.sqrt(radicand)
            ne.append(float(I.sum() * dLl * norm))
        self.ne = np.asarray(ne)

    def plot(self, axs: np.ndarray) -> None:
        """Overlay this mirror on a 2x2 axes grid (walkthrough Fig. 10 style)."""
        axs[0, 0].plot(self.zax, self.B, label=rf"$w={self.w}$")
        axs[0, 1].plot(self.lax, self.eta)
        axs[1, 0].plot(self.eta_ax, self.I1)
        axs[1, 1].plot(self.zax, self.ne)


# =============================================================================
# Part D — Map f(v,ξ) along B(z) via Λ   [paper §4 / Fig. 12 idea]
# =============================================================================
#
# f is constant along a bounce orbit when expressed in (E, Λ).
# Midplane cell (v, ξ) ↔ (E, Λ=1-ξ²).  At location z with B̃=B/B0:
#   - cell is ABSENT if Λ ≥ 1/B̃  (turning point before z)
#   - otherwise local pitch satisfies ξ_loc² = 1 - Λ B̃
#
# Density weight generalizes Eq. 57's Jacobian:
#   dn ∝ f(v,Λ) v² dv  *  B̃ / sqrt(1 - B̃ Λ)  dΛ
# (the B̃/sqrt factor is how many particles of that Λ are at this z).


@dataclass
class AxialProfiles:
    """Mapped axial profiles from a midplane distribution."""

    zax: np.ndarray
    B: np.ndarray
    n: np.ndarray              # relative density from mapped f
    n_eq57: np.ndarray         # Eq. 57 benchmark (I1-only), same z grid
    Sn: np.ndarray             # relative beam–beam neutron source
    Sn_norm: np.ndarray        # Sn / ∫ Sn dz  (PDF along z ∈ [0,1])
    midplane: MidplaneDistribution


def _midplane_cells(mid: MidplaneDistribution) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Flatten midplane f into cell centers (v, ξ, Λ, weight).

    weight_i ∝ f * v² Δv Δξ   — relative phase-space content of the cell.
    Only ξ ≥ 0 is stored (distribution is even in ξ for this model).
    """
    vax, xax, f = mid.vax, mid.xax, mid.f
    dv = float(np.diff(vax).mean())
    dxi = float(np.diff(xax).mean())
    V, X = np.meshgrid(vax, xax, indexing="ij")
    W = f * (V**2) * dv * dxi
    Lam = 1.0 - X**2
    # Fix: Λ extent of each ξ cell, for analytic per-cell integration of
    # the axial Jacobian (see map_distribution_along_z). ξ∈[ξ-h, ξ+h] maps to
    # Λ∈[Λ_lo, Λ_hi] with Λ decreasing in ξ.
    h = 0.5 * dxi
    Lam_lo = 1.0 - np.clip(X + h, 0.0, 1.0) ** 2
    Lam_hi = 1.0 - np.clip(X - h, 0.0, 1.0) ** 2
    Wbase = f * (V**2) * dv  # measure-free amplitude; Λ-cell integral supplies dΛ→dξ_local
    # Drop empty / tiny cells for speed.
    mask = W > 1e-14 * W.max()
    return (
        V[mask],
        X[mask],
        Lam[mask],
        W[mask],
        Lam_lo[mask],
        Lam_hi[mask],
        Wbase[mask],
    )


def bosch_hale_dt_sigma_barns(E_cm_keV: np.ndarray) -> np.ndarray:
    """
    D-T fusion cross section σ(E_cm) [barns].

    H.-S. Bosch & G.M. Hale, Nucl. Fusion 32, 611 (1992), Table VII Padé fit.
    E_cm_keV : center-of-mass energy [keV].

    Absolute barns only matter for absolute rates. For a *normalized*
    S_n(z) shape the overall σ scale cancels — but the *energy dependence*
    of σ is what makes S_n(z) differ from [n(z)]² when the local spectrum
    changes along z.
    """
    E = np.asarray(E_cm_keV, dtype=float)
    BG = 34.3827
    A1, A2, A3, A4, A5 = 6.927e4, 7.454e8, 2.050e6, 5.2002e4, 0.0
    B1, B2, B3, B4 = 6.38e1, -9.95e-1, 6.981e-5, 1.728e-4

    out = np.zeros_like(E)
    pos = E > 0.0
    Ep = E[pos]
    Snum = A1 + Ep * (A2 + Ep * (A3 + Ep * (A4 + Ep * A5)))
    Sden = 1.0 + Ep * (B1 + Ep * (B2 + Ep * (B3 + Ep * B4)))
    S = Snum / Sden
    # Bosch–Hale returns σ in millibarns for this S; convert to barns.
    out[pos] = (S / Ep) * np.exp(-BG / np.sqrt(Ep)) * 1e-3
    return out


# --- D/T kinematics ---------------------------------------------------------
# Normalized speed û is defined by ENERGY, E/E_beam = û², mass-agnostic.
# Physical speed differs by species: v_phys = û·sqrt(2 E_beam / m_s).
# In units of the DEUTERON injection speed v0_D = sqrt(2 E_beam / m_D):
#     v̂_phys = û                 (D)
#     v̂_phys = û·sqrt(m_D/m_T)   (T, ≈ 0.817 — same energy, √(2/3) slower)
# CM energy for a D–T pair:  E_cm = ½ μ |Δv_phys|²
#     = [m_T/(m_D+m_T)]·E_beam·|Δv̂|²  = 0.6·E_beam·|Δv̂|²
# with Δv̂ built from the SPECIES-SCALED vectors above. The old MC path fed
# unscaled û for both species (T ~22% too fast) into this same 0.6 prefactor.
_M_D, _M_T = 2.0141, 3.0160
_T_SPEED_SCALE = float(np.sqrt(_M_D / _M_T))


def _beam_beam_rate_deterministic(
    v_par: np.ndarray,
    v_perp: np.ndarray,
    w: np.ndarray,
    E_beam_keV: float,
    n_bins: int = 24,
    n_gyro: int = 16,
) -> float:
    """
    Deterministic ∬ f_D f_T σ(|Δv|)|Δv| from local (|v∥|, v⊥) cell weights.

    Steps
    -----
    1. Bin weights onto an (|v∥|, v⊥) grid (collapses the many midplane cells
       that land on similar local coordinates).
    2. Split each bin 50/50 into D and T; scale T speeds by sqrt(m_D/m_T).
    3. Sum over all D–T bin pairs, averaging |Δv| over relative gyrophase
       (midpoint nodes on [0, π]) and the two relative v∥ signs — both f's are
       gyrotropic and bounce-symmetric, so that average is exact in those
       angles.

    Replaces the old Monte-Carlo macroparticle estimate — no seed, no
    z-to-z noise, no self-pair bias, and D–D / T–T pairs (which the MC wrongly
    counted with the DT cross-section) are excluded by construction.
    Cost per z: O(n_occupied² · 2 · n_gyro) vectorized σ evals; n_bins=24
    keeps this in the few-million range.
    Convergence: n_bins is converged by 16 (<0.5%); the gyro average is only
    ~O(1/n_gyro) because |Δv| has a √-cusp at zero relative speed —
    n_gyro=16 ⇒ ≈3% max on S_n(z), 32 ⇒ <1%.
    """
    if len(w) < 1 or w.sum() <= 0.0:
        return 0.0

    H, p_edges, q_edges = np.histogram2d(
        v_par, v_perp, bins=n_bins, weights=w
    )
    pc = 0.5 * (p_edges[1:] + p_edges[:-1])
    qc = 0.5 * (q_edges[1:] + q_edges[:-1])
    P, Q = np.meshgrid(pc, qc, indexing="ij")
    occ = H > 0.0
    p, q, W = P[occ], Q[occ], H[occ]
    if len(W) < 1:
        return 0.0

    # Species: same normalized-energy distribution, physical-speed scaled.
    pD, qD, wD = p, q, 0.5 * W
    pT, qT, wT = p * _T_SPEED_SCALE, q * _T_SPEED_SCALE, 0.5 * W

    phi = (np.arange(n_gyro) + 0.5) * np.pi / n_gyro
    cphi = np.cos(phi)

    qq = np.outer(qD, qT)
    perp2_base = qD[:, None] ** 2 + qT[None, :] ** 2
    K = np.zeros((len(pD), len(pT)))
    for s in (1.0, -1.0):
        dpar2 = (pD[:, None] - s * pT[None, :]) ** 2
        for c in cphi:
            du2 = dpar2 + perp2_base - 2.0 * qq * c
            du = np.sqrt(np.clip(du2, 0.0, None))
            K += bosch_hale_dt_sigma_barns(0.6 * E_beam_keV * du2) * du
    K /= 2.0 * n_gyro

    # Rate ∝ Σ_ij w_D,i K_ij w_T,j  (each D–T pair counted once).
    return float(wD @ K @ wT)


def map_distribution_along_z(
    mid: MidplaneDistribution,
    Rm: float,
    w_throat: float = 0.1,
    k: float = 1.3,
    Nz: int = 81,
    n_bins: int = 24,
    n_gyro: int = 16,
    include_eq57: bool = True,
) -> AxialProfiles:
    """
    Map midplane f → n(z) and S_n(z) on analytic B(z).

    Parameters
    ----------
    mid : MidplaneDistribution
        From `build_midplane_f`.
    Rm, w_throat, k :
        Analytic mirror field parameters (paper Eq. 56).
    Nz :
        Axial grid resolution on z ∈ [0, 1].
    n_bins, n_gyro :
        (|v∥|, v⊥) binning and gyrophase quadrature for the deterministic
        beam–beam S_n (replaces the old Monte-Carlo n_macro/n_pair_samples;
        no seed — results are exactly reproducible).
    include_eq57 :
        Also compute Mirror/Eq. 57 density on the same z grid for comparison.
    """
    zax = np.linspace(0.0, 1.0, Nz)
    B = np.asarray(Bfunc(zax, w=w_throat, k=k, Rm=Rm), dtype=float)

    v_c, xi_c, Lam_c, w_c, Llo_c, Lhi_c, wb_c = _midplane_cells(mid)
    logger.info("Midplane cells kept: %d (of %d grid points)",
                len(w_c), mid.f.size)

    n_z = np.zeros(Nz)
    Sn = np.zeros(Nz)

    for iz, (z, Bt) in enumerate(zip(zax, B)):
        # Present at this z if any part of the cell's Λ interval is below 1/B̃
        # (and cell is magnetically trapped, Λ ≥ 1/Rm; f is already zeroed in
        # the loss cone, so the second condition is belt-and-suspenders).
        trapped = (Llo_c < 1.0 / Bt) & (Lam_c >= 1.0 / Rm)
        if not np.any(trapped):
            continue

        v_t = v_c[trapped]
        Lam_t = Lam_c[trapped]

        # ---- density: exact per-cell Λ integral of the axial Jacobian ----
        # Fix (two layers):
        # (1) Measure: midplane weights carry dξ, but the axial pile-up is
        #     dξ_local = B̃ dΛ / (2 sqrt(1 - B̃ Λ)); the original code applied
        #     B̃/sqrt(1-B̃Λ) directly to dξ-weights, overweighting deeply
        #     trapped (Λ→1) cells by 1/ξ — a clipped ∫dξ/ξ that grows with Nξ
        #     and artificially midplane-peaks n(z) and S_n(z).
        # (2) Turning points: 1/sqrt(1-B̃Λ) is integrably singular; POINT
        #     evaluation spikes whenever a grid Λ sits ε from 1/B̃ (the old
        #     sawtooth). Integrate analytically over each cell instead:
        #       ∫ B̃ dΛ / (2 sqrt(1-B̃Λ)) = sqrt(1-B̃Λ_lo) - sqrt(1-B̃Λ_ub),
        #     Λ_ub = min(Λ_hi, 1/B̃). Exact, finite, and reduces to dξ at
        #     B̃=1, so mapped n(0) equals the direct midplane integral of f.
        Llo = Llo_c[trapped]
        Lub = np.minimum(Lhi_c[trapped], 1.0 / Bt)
        cellint = np.sqrt(np.clip(1.0 - Bt * Llo, 0.0, None)) - np.sqrt(
            np.clip(1.0 - Bt * Lub, 0.0, None)
        )
        w_z = wb_c[trapped] * cellint
        n_z[iz] = float(np.sum(w_z))

        # ---- beam–beam neutron source: deterministic D–T pair sum ----
        # Local pitch decomposition from magnetic moment conservation:
        # |v∥|/v = sqrt(1 - Λ B̃),  v⊥/v = sqrt(Λ B̃). Signs and gyrophase are
        # averaged exactly inside the kernel (gyrotropic, bounce-symmetric f).
        v_par_loc = v_t * np.sqrt(np.clip(1.0 - Lam_t * Bt, 0.0, None))
        v_perp_loc = v_t * np.sqrt(np.clip(Lam_t * Bt, 0.0, None))
        Sn[iz] = _beam_beam_rate_deterministic(
            v_par_loc, v_perp_loc, w_z, mid.E_beam_keV,
            n_bins=n_bins, n_gyro=n_gyro,
        )

    # Normalize density to mean 1 over z (shape compare with Eq. 57).
    if n_z.mean() > 0:
        n_z = n_z / n_z.mean()

    # Normalized neutron source PDF on z ∈ [0, 1].
    dz = float(np.diff(zax).mean()) if Nz > 1 else 1.0
    integ = float(np.sum(Sn) * dz)
    Sn_norm = Sn / integ if integ > 0 else Sn.copy()

    # Eq. 57 benchmark (I1 only — no injection angle).
    if include_eq57:
        logger.info("Computing Eq. 57 Mirror benchmark (slow: Legendre scan)…")
        mir = Mirror(Rm=Rm, w=w_throat, k=k, Nz=Nz)
        # Mirror uses its own zax; interpolate onto ours if needed.
        n_eq57 = np.interp(zax, mir.zax, mir.ne)
        # Match normalization (mean 1).
        n_eq57 = n_eq57 / n_eq57.mean()
    else:
        n_eq57 = np.full_like(n_z, np.nan)

    return AxialProfiles(
        zax=zax,
        B=B,
        n=n_z,
        n_eq57=n_eq57,
        Sn=Sn,
        Sn_norm=Sn_norm,
        midplane=mid,
    )


def save_axial_csv(profiles: AxialProfiles, path: Path) -> None:
    """Write z, B, n, n_eq57, Sn, Sn_norm to CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    header = "z,B_over_B0,n_mapped,n_eq57,Sn,Sn_norm"
    data = np.column_stack(
        [profiles.zax, profiles.B, profiles.n, profiles.n_eq57, profiles.Sn, profiles.Sn_norm]
    )
    np.savetxt(path, data, delimiter=",", header=header, comments="")
    logger.info("wrote %s", path)
