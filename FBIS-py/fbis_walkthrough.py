#!/usr/bin/env python3
"""
FBIS walkthrough — Python reproduction of the FBIS (2022, NF) analytic model.

Paper: Egedal, Endrizzi, Forest, Fowler — Fast Beam Ion Solver / mirror fusion.
Notes: latex/fbis-notes.tex (Tony Qian).

This script is the cleaned-up form of FBISpy.ipynb. Read it top-to-bottom, or run
individual sections with:

    python fbis_walkthrough.py --section basis
    python fbis_walkthrough.py --section physics
    python fbis_walkthrough.py --section bz
    python fbis_walkthrough.py --section all

For the **mapped neutron source S_n(z)** (multi-mode f → Λ map → beam–beam),
see the heavily commented modules:

    fbis_core.py                  # physics + mapping (feed this to Fable)
    compute_neutron_source.py     # CLI

Plots are written to outputs/ and also shown if a display is available.

Physical picture (short):
  1. Pitch-angle operator L has eigenfunctions M_l(ξ) built from Legendre P/Q.
  2. Loss-cone boundary ξ* = sqrt(1 - 1/R_m) selects discrete eigenvalues l_j.
  3. Source S(ξ) projected onto {M_j} + slowing-down u(v) → distribution f(v,ξ).
  4. For smooth B(z), map ξ → Λ, bounce-average, build matrix H, get n_e(z).
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import matplotlib

# Default to Agg so the script runs headless; override with MPLBACKEND=... for GUIs.
if "MPLBACKEND" not in os.environ:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline

# Shared physics lives in fbis_core (single source of truth for Fable / neutron map).
from fbis_core import (
    OUTPUT_DIR,
    Bfunc,
    Mirror,
    S0,
    find_eigenvalues_from_data,
    get_1overB,
    get_tau,
    get_zTP,
    make_M,
    ufunc,
)

logger = logging.getLogger(__name__)


def section_basis(*, show: bool = True) -> dict:
    """Reproduce notebook §Basis Functions: coefficients, contours, eigenmodes."""
    logger.info("Part 1: Legendre basis and loss-cone eigenvalues")

    # --- sanity check for one l ---
    l_demo = 2.3
    a_demo, b_demo, M_demo = make_M(l_demo)
    logger.info("demo l=%.2f → a=%.6f, b=%.6f, M(0)=%.6f", l_demo, a_demo, b_demo, M_demo(0))

    # --- scan (a,b) vs l ---
    # Half-integer steps in l rotate (a,b) by ~π/2.
    lax_ab = np.linspace(0.11, 12.2, 99)
    abM = [make_M(l) for l in lax_ab]
    a_coef = np.array([t[0] for t in abM])
    b_coef = np.array([t[1] for t in abM])

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    axs[0].plot(lax_ab, a_coef, label="a")
    axs[0].plot(lax_ab, b_coef, label="b")
    axs[0].set_xlabel("l")
    axs[0].legend()
    axs[0].grid(True)
    axs[0].set_title("BC coefficients vs degree l")
    axs[1].plot(a_coef, b_coef)
    axs[1].set_xlabel("a")
    axs[1].set_ylabel("b")
    axs[1].set_title("(a,b) trajectory — π/2 per half-integer l")
    fig.tight_layout()
    _save(fig, "01_ab_coefficients.png", show)

    # --- M(ξ; l) contour: white = zeros = candidate eigenvalues ---
    # Vertical cut at fixed ξ* (mirror ratio) gives the discrete spectrum.
    lax = np.linspace(0.11, 25.3, 159)
    Ms = [make_M(l)[2] for l in lax]
    xax = np.linspace(0.0, 1.0, 100, endpoint=False)
    data = np.array([[m(xi) for xi in xax] for m in Ms])

    fig, ax = plt.subplots()
    cf = ax.contourf(xax, lax, data, levels=50, cmap="RdBu")
    fig.colorbar(cf, ax=ax)
    ax.contour(xax, lax, data, levels=[0], colors="y", linewidths=1)
    ax.set_xlabel(r"$\xi = v_\parallel / v$")
    ax.set_ylabel("l")
    ax.set_title(r"$M_l(\xi)$  (yellow: $M=0$)")
    _save(fig, "02_M_contour.png", show)

    # --- pick a mirror ratio / loss-cone edge and find l_j ---
    x_star = 0.95  # ξ* = sqrt(1 - 1/R_m)  →  R_m ≈ 10.3
    Rm = 1.0 / (1.0 - x_star**2)
    eigs = find_eigenvalues_from_data(x_star, xax, lax, data)
    logger.info("R_m = %.2f at ξ*=%.2f → %d eigenvalues in scan", Rm, x_star, len(eigs))

    fig, ax = plt.subplots()
    cf = ax.contourf(xax, lax, data, levels=50, cmap="RdBu")
    fig.colorbar(cf, ax=ax)
    ax.contour(xax, lax, data, levels=[0], colors="y", linewidths=1)
    ax.scatter([x_star] * len(eigs), eigs, color="r", zorder=5, s=20)
    ax.axvline(x_star, ls="--", color="k", label=rf"$R_m={Rm:.2f}$")
    ax.set_xlabel(r"$\xi$")
    ax.set_ylabel("l")
    ax.legend()
    ax.set_title("eigenvalues (red) at loss-cone boundary")
    _save(fig, "03_eigenvalues.png", show)

    # --- first eigenfunctions M_j(ξ) ---
    n_modes = min(10, len(eigs))
    M_arr = np.array([[make_M(eigs[j])[2](x) for x in xax] for j in range(n_modes)])

    fig, ax = plt.subplots()
    for j in range(n_modes):
        ax.plot(xax, M_arr[j], label=rf"$\lambda_{{{j}}}={eigs[j]:.4f}$")
    ax.axvline(x_star, ls="--", color="k")
    ax.set_xlabel(r"$\xi$")
    ax.set_ylabel(r"$M_j(\xi)$")
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_title("first eigenfunctions")
    _save(fig, "04_eigenfunctions.png", show)

    # mirror ratio ↔ pitch mapping (geometry only)
    Rm_ax = np.linspace(1.1, 100.0, 200)
    xi_of_Rm = np.sqrt(1.0 - 1.0 / Rm_ax)
    fig, ax = plt.subplots()
    ax.plot(Rm_ax, xi_of_Rm)
    ax.set_xlabel(r"$R_m$")
    ax.set_ylabel(r"$\xi^* = \sqrt{1 - 1/R_m}$")
    ax.set_title("loss-cone pitch vs mirror ratio")
    ax.grid(True)
    _save(fig, "05_Rm_vs_xi.png", show)

    return {
        "xax": xax,
        "lax": lax,
        "data": data,
        "x_star": x_star,
        "Rm": Rm,
        "eigs": eigs,
        "M_arr": M_arr,
    }


# ---------------------------------------------------------------------------
# Part 2 — Velocity-space physics (Section II of the paper)
# ---------------------------------------------------------------------------
# ufunc / S0 imported from fbis_core.


def section_physics(basis: dict, *, show: bool = True) -> dict:
    """Reproduce notebook §Physics: u(v), source projection, f1, (v∥,v⊥) map."""
    logger.info("Part 2: velocity-space distribution from projected source")

    xax = basis["xax"]
    x_star = basis["x_star"]
    eigs = basis["eigs"]
    M_arr = basis["M_arr"]
    n_modes = M_arr.shape[0]

    vax = np.linspace(0.0, 60.0, 100)

    # --- u(v) parameter studies ---
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    for vc in [5, 10, 20, 30, 50]:
        axs[0].plot(vax, ufunc(vax, vc=vc), label=rf"$v_c={vc}$")
    for beta in [0.5, 1, 2, 5]:
        axs[1].plot(vax, ufunc(vax, vc=12, beta=beta), label=rf"$\beta={beta}$")
    axs[0].set_title(r"$u(v)$ vs $v_c$")
    axs[1].set_title(r"$u(v)$ vs $\beta$")
    for ax in axs:
        ax.set_xlabel("v")
        ax.legend()
        ax.grid(True)
    fig.tight_layout()
    _save(fig, "06_ufunc.png", show)

    # --- single-mode 2D basis pieces (injection/critical velocity fixed) ---
    vc, beta = 20.0, 5.0
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    for j in range(min(4, n_modes)):
        lamb = eigs[j]
        # shape (Nv, Nξ): M(ξ) / (v^3+vc^3) * u(v)^λ
        mode = np.array(
            [
                [m / (v**3 + vc**3) * ufunc(np.array([v]), vc=vc, beta=beta)[0] ** lamb for m in M_arr[j]]
                for v in vax
            ]
        )
        C = axs[j].contourf(xax, vax, mode, 20)
        fig.colorbar(C, ax=axs[j])
        axs[j].axhline(vc, ls="--", color="w")
        axs[j].axvline(x_star, ls="--", color="w")
        axs[j].set_xlabel(r"$\xi$")
        axs[j].set_ylabel("v")
        axs[j].set_title(rf"mode $j={j}$")
    fig.suptitle("2D eigenbasis pieces for the source expansion")
    fig.tight_layout()
    _save(fig, "07_mode_basis_2d.png", show)

    # --- project Gaussian S(ξ) onto {M_j}:  S_j = <S, M_j> / ||M_j||^2 ---
    S = S0(xax)
    alpha = np.sum(M_arr**2, axis=1)
    dxi = float(np.diff(xax).mean())
    Sj = np.sum([S * M for M in M_arr], axis=1) * dxi / alpha

    fig, axs = plt.subplots(1, 2, figsize=(9, 3))
    axs[0].plot(xax, S)
    axs[0].set_title(r"$S(\xi)$")
    axs[0].set_xlabel(r"$\xi$")
    axs[1].plot(Sj, "o-")
    axs[1].set_title(r"$S_j$ coefficients")
    axs[1].set_xlabel("j")
    fig.tight_layout()
    _save(fig, "08_source_projection.png", show)

    # --- reconstruct f₁(v,ξ) ≈ Σ_j S_j · mode_j   (paper Fig. 3 / Eq. 14) ---
    vc, beta = 22.0, 4.0
    f1 = np.zeros((len(vax), len(xax)))
    for j in range(n_modes):
        lamb = eigs[j]
        f1 += np.array(
            [
                [
                    Sj[j] * m / (v**3 + vc**3) * ufunc(np.array([v]), vc=vc, beta=beta)[0] ** lamb
                    for m in M_arr[j]
                ]
                for v in vax
            ]
        )

    fig, ax = plt.subplots()
    cf = ax.contourf(xax, vax, f1, 20, cmap="plasma")
    fig.colorbar(cf, ax=ax)
    ax.axhline(vc, ls="--", color="w")
    ax.axvline(x_star, ls="--", color="w")
    ax.set_xlabel(r"$\xi$")
    ax.set_ylabel("v")
    ax.set_title(r"$f_1(v,\xi)$ from projected source (Fig. 3)")
    _save(fig, "09_f1_xi_v.png", show)

    # --- map (v, ξ) → (v∥, v⊥) ---
    theta = np.arccos(xax)
    vpar = vax[np.newaxis, :] * np.cos(theta)[:, np.newaxis]
    vperp = vax[np.newaxis, :] * np.sin(theta)[:, np.newaxis]

    fig, ax = plt.subplots()
    cf = ax.contourf(vpar, vperp, f1.T, 20, cmap="plasma")
    fig.colorbar(cf, ax=ax)
    ax.set_xlabel(r"$v_\parallel$")
    ax.set_ylabel(r"$v_\perp$")
    ax.set_title(r"$f_1$ in $(v_\parallel, v_\perp)$")
    _save(fig, "10_f1_vpar_vperp.png", show)

    return {"vax": vax, "Sj": Sj, "f1": f1, "vpar": vpar, "vperp": vperp}


# ---------------------------------------------------------------------------
# Part 3 — Smooth B(z): bounce averages (Section 4)
# ---------------------------------------------------------------------------
#
# Step-function mirrors → continuous B(z). Pitch ξ is replaced by
#     Λ = μ ε^{-1} B_0   (magnetic moment / energy),
# with turning point z_TP(Λ) where B(z_TP)/B_0 = 1/Λ.
# Bounce time τ(Λ) and <1/B> enter the pitch-angle operator and η(Λ).
# Bfunc / get_zTP / get_tau / get_1overB imported from fbis_core.


def section_bz(*, show: bool = True) -> dict:
    """Reproduce notebook §Physics with B(z): profiles, τ, η for several widths w."""
    logger.info("Part 3: B(z) model and bounce averages")

    zax = np.linspace(0.0, 1.0, 1001)
    lax = np.linspace(1e-3, 1.0, 200, endpoint=False)
    dL = float(np.diff(lax).mean())
    Rm = 16.0

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for w in [0.03, 0.1, 0.2, 0.3, 0.5]:
        axs[1].plot(zax, Bfunc(zax, w=w, Rm=Rm), label=rf"$w={w}$")
    axs[1].legend()
    axs[1].set_xlabel("z")
    axs[1].set_title(r"$B(z)$ vs throat width $w$")
    axs[0].plot(zax, Bfunc(zax, w=0.2, Rm=Rm))
    axs[0].set_xlabel("z")
    axs[0].set_title(r"$B(z)$ example ($w=0.2$)")
    fig.tight_layout()
    _save(fig, "11_Bfunc.png", show)

    # single-w diagnostics
    w = 0.1
    Bf = lambda z: Bfunc(z, w=w, Rm=Rm)
    z_TP = get_zTP(lax, Bf)
    tau = get_tau(lax, Bf, z_TP)
    avgB = get_1overB(lax, Bf, z_TP, tau)

    fig, axs = plt.subplots(1, 3, figsize=(11, 4))
    axs[0].plot(zax, Bfunc(zax, w=w, Rm=Rm), "C2", label="B(z)")
    axs[0].set_xlabel("z")
    axs[1].plot(lax, z_TP, label=r"$z_{TP}$")
    axs[1].plot(lax, avgB, label=r"$\langle 1/B \rangle$")
    axs[1].set_xlabel(r"$\Lambda$")
    axs[1].set_ylim(bottom=0)
    axs[2].plot(lax, tau, "C4", label=r"$\tau$")
    axs[2].axhline(np.sum(tau) * dL, ls="--", color="k", lw=0.7, label=r"$\int\tau\,d\Lambda$")
    axs[2].set_xlabel(r"$\Lambda$")
    for ax in axs:
        ax.legend(fontsize=8)
    fig.tight_layout()
    _save(fig, "12_bounce_integrals.png", show)

    # η(Λ) for several throat widths — maps Λ onto a pitch-like coordinate
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
    eta_data: dict[float, tuple[np.ndarray, Callable]] = {}
    for w in [0.03, 0.1, 0.2, 0.3, 0.5]:
        Bf = lambda z, ww=w: Bfunc(z, w=ww, Rm=Rm)
        z_TP = get_zTP(lax, Bf)
        tau = get_tau(lax, Bf, z_TP)
        Itau = float(tau.sum() * dL)
        avgB = get_1overB(lax, Bf, z_TP, tau)
        spl_B = make_interp_spline(lax, avgB, k=3)
        eta = 1.0 - tau.cumsum() * dL / Itau
        eta_TP = float(np.interp(1.0 / Rm, lax, eta))
        eta_data[eta_TP] = (eta, spl_B)

        axs[0].plot(zax, Bfunc(zax, w=w, Rm=Rm), label=rf"$w={w}$")
        axs[1].plot(lax, tau)
        axs[2].plot(lax, eta)

    axs[0].set_title("B")
    axs[1].set_title(r"$\tau$")
    axs[2].set_title(r"$\eta$")
    axs[0].set_xlabel("z")
    axs[1].set_xlabel(r"$\Lambda$")
    axs[2].set_xlabel(r"$\Lambda$")
    axs[2].plot(lax, np.sqrt(1.0 - lax), "k--", lw=0.5, label=r"$\xi=\sqrt{1-\Lambda}$")
    axs[1].axvline(1.0 / Rm, color="k", ls="--", lw=0.7, alpha=0.5)
    axs[2].axvline(1.0 / Rm, color="k", ls="--", lw=0.7, alpha=0.5)
    axs[0].legend(fontsize=8)
    axs[2].legend(fontsize=8)
    fig.tight_layout()
    _save(fig, "13_eta_profiles.png", show)

    logger.info("η_TP values for scanned w: %s", [f"{k:.4f}" for k in eta_data])
    return {"zax": zax, "lax": lax, "Rm": Rm, "eta_data": eta_data}


# ---------------------------------------------------------------------------
# Part 4 — Mirror class: matrix H → I₁(η) → n_e(z)
# ---------------------------------------------------------------------------
# Mirror / get_M_basis / get_LM / get_H imported from fbis_core.


def section_mirror(*, show: bool = True) -> list[Mirror]:
    """Build Mirror objects for several throat widths and plot Fig.-10-style panels."""
    logger.info("Part 4: Mirror class — H eigenmode → n_e(z)")

    mirrors: list[Mirror] = []
    for w in [0.03, 0.1, 0.2, 0.3, 0.5]:
        logger.info("  building Mirror(w=%.2f) ...", w)
        m = Mirror(w=w)
        mirrors.append(m)
        logger.info("    η_TP = %.4f", m.eta_TP)

    fig, axs = plt.subplots(2, 2, figsize=(9, 10))
    for m in mirrors:
        m.plot(axs)

    Rm = mirrors[-1].Rm
    axs[0, 1].axvline(1.0 / Rm, ls="--", color="k")
    axs[0, 0].set_title("B(z)")
    axs[0, 0].set_xlabel("z")
    axs[0, 1].set_title(r"$\eta(\Lambda)$")
    axs[0, 1].set_xlabel(r"$\Lambda$")
    axs[1, 0].set_title(r"$I_1(\eta)$")
    axs[1, 0].set_xlabel(r"$\eta$")
    axs[1, 1].set_title(r"$n_e(z)/\langle n\rangle$")
    axs[1, 1].set_xlabel("z")
    axs[0, 0].legend()
    fig.tight_layout()
    _save(fig, "14_mirror_summary.png", show)

    # H matrix for the last case (symmetry check from notebook)
    H = mirrors[-1].H
    asym = float(np.max(np.abs(H - H.T)))
    logger.info("max |H - Hᵀ| = %.3e (notebook: not exactly symmetric)", asym)

    fig, ax = plt.subplots()
    im = ax.imshow(H)
    fig.colorbar(im, ax=ax)
    ax.set_title("H matrix (last w)")
    _save(fig, "15_H_matrix.png", show)

    return mirrors


# ---------------------------------------------------------------------------
# Utilities / CLI
# ---------------------------------------------------------------------------


def _save(fig: plt.Figure, name: str, show: bool) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    logger.info("wrote %s", path)
    if show:
        plt.show()
    else:
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Walk through the FBIS analytic model (cleaned FBISpy.ipynb)."
    )
    p.add_argument(
        "--section",
        choices=("basis", "physics", "bz", "mirror", "all"),
        default="all",
        help="Which part to run (default: all).",
    )
    p.add_argument(
        "--no-show",
        action="store_true",
        help="Save plots only; do not open interactive windows.",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="DEBUG logging.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s | %(message)s",
    )
    show = not args.no_show
    section = args.section

    basis: dict | None = None
    if section in ("basis", "physics", "all"):
        basis = section_basis(show=show)

    if section in ("physics", "all"):
        assert basis is not None
        section_physics(basis, show=show)

    if section in ("bz", "mirror", "all"):
        # Mirror rebuilds bounce averages; still useful pedagogy before the class.
        if section in ("bz", "all"):
            section_bz(show=show)

    if section in ("mirror", "all"):
        section_mirror(show=show)

    logger.info("done. plots in %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
