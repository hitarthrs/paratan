#!/usr/bin/env python3
"""
Reproduce Egedal et al. 2022 Figure 10 (all six panels).

Figure 10 is §4 (orbit-averaged Lorentz / non-uniform B), not §3.

Panels
------
(a) B(z) for R_M=16, w ∈ {0.03,0.1,0.2,0.3,0.5}
(b) η(Λ) for those fields
(c) I₁(η) first orbit-averaged eigenfunctions
(d) λ₁/λ₁,square vs log10(R_M)  [unweighted bounce average]
(e) n(z)/⟨n⟩ from Eq. (57) with f ∝ I₁
(f) like (d) but with n(z)-weighted scattering average

Example
-------
    python reproduce_fig10.py --no-show
    python reproduce_fig10.py --quick --no-show   # coarser Rm scan
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import matplotlib

if "MPLBACKEND" not in os.environ:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from fbis_core import OUTPUT_DIR, Bfunc
from fbis_orbit_averaged import (
    bounce_averages,
    scan_lambda_ratio_vs_Rm,
    solve_orbit_averaged_I1,
)

logger = logging.getLogger(__name__)

# Paper Fig. 10 color scheme (matches walkthrough / notebook habit)
W_VALUES = [0.03, 0.1, 0.2, 0.3, 0.5]
COLORS = {
    0.03: "black",
    0.1: "magenta",
    0.2: "limegreen",
    0.3: "red",
    0.5: "blue",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Reproduce paper Figure 10")
    p.add_argument("--Rm", type=float, default=16.0, help="R_M for panels a,b,c,e")
    p.add_argument("--k", type=float, default=1.3)
    p.add_argument(
        "--quick",
        action="store_true",
        help="Fewer R_M points / coarser grids for panels d,f",
    )
    p.add_argument("--skip-scans", action="store_true",
                   help="Skip slow panels (d) and (f)")
    p.add_argument("--no-show", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s | %(message)s",
    )

    Rm = args.Rm
    k = args.k
    Nz = 301 if args.quick else 401
    N_Lambda = 100 if args.quick else 160
    N_eta = 200 if args.quick else 300

    # --- panels a–c,e: fixed Rm, scan w ---
    logger.info("Building orbit-averaged I1 for Rm=%.3g, %d widths...", Rm, len(W_VALUES))
    eigs = []
    for w in W_VALUES:
        logger.info("  w=%.3f (unweighted) ...", w)
        eigs.append(
            solve_orbit_averaged_I1(
                Rm=Rm, w=w, k=k, Nz=Nz, N_Lambda=N_Lambda, N_eta=N_eta,
                density_weight=False,
            )
        )

    # Also need bounce η(Λ) on a common Λ grid for panel (b)
    bounces = [
        bounce_averages(Rm=Rm, w=w, k=k, N_Lambda=N_Lambda) for w in W_VALUES
    ]

    # --- panels d,f: Rm scan ---
    scan_d: dict | None = None
    scan_f: dict | None = None
    if not args.skip_scans:
        if args.quick:
            Rm_scan = np.logspace(0.3, 2.5, 8)  # ~2 .. 316
        else:
            Rm_scan = np.logspace(0.3, 3.0, 14)  # ~2 .. 1000
        logger.info("Panel (d): λ1/λ1□ vs Rm (unweighted), %d points × %d w...",
                    len(Rm_scan), len(W_VALUES))
        scan_d = scan_lambda_ratio_vs_Rm(
            W_VALUES, Rm_scan, density_weight=False, k=k,
            Nz=201 if args.quick else 251,
            N_Lambda=80 if args.quick else 100,
            N_eta=150 if args.quick else 200,
        )
        logger.info("Panel (f): λ1/λ1□ vs Rm (n-weighted)...")
        scan_f = scan_lambda_ratio_vs_Rm(
            W_VALUES, Rm_scan, density_weight=True, k=k,
            Nz=201 if args.quick else 251,
            N_Lambda=80 if args.quick else 100,
            N_eta=150 if args.quick else 200,
        )

    # --- figure ---
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))
    ax_a, ax_b, ax_c = axs[0]
    ax_d, ax_e, ax_f = axs[1]

    zax = eigs[0].zax
    for w, eig, bn in zip(W_VALUES, eigs, bounces):
        c = COLORS[w]
        ax_a.plot(zax, eig.B, color=c, lw=2, label=f"w={w:g}")
        ax_b.plot(bn.lax, bn.eta, color=c, lw=2)
        ax_c.plot(eig.eta_ax, eig.I1, color=c, lw=2)
        ax_e.plot(eig.zax, eig.ne, color=c, lw=2)

    # (a)
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, Rm)
    ax_a.set_xlabel(r"$z/L$")
    ax_a.set_ylabel(r"$B/B_0$")
    ax_a.set_title(rf"(a)  $B(z)$  ($R_M={Rm:g}$)")
    ax_a.legend(loc="upper left", fontsize=8, frameon=False)
    ax_a.text(0.05, 0.55, rf"$R_M={int(Rm)}$", transform=ax_a.transAxes)

    # (b)
    ax_b.axvline(1.0 / Rm, color="k", ls="--", lw=1)
    ax_b.text(1.0 / Rm + 0.02, 0.05, r"$1/R_M$", fontsize=10)
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    ax_b.set_xlabel(r"$\Lambda$")
    ax_b.set_ylabel(r"$\eta$")
    ax_b.set_title(r"(b)  $\eta(\Lambda)$")

    # (c)
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1.05)
    ax_c.set_xlabel(r"$\eta$")
    ax_c.set_ylabel(r"$I_1(\eta)$")
    ax_c.set_title(r"(c)  First eigenfunction $I_1$")

    # (d)
    if scan_d is not None:
        for w in W_VALUES:
            logR, ratio = scan_d[w]
            ax_d.plot(logR, ratio, color=COLORS[w], lw=2, label=f"w={w:g}")
        ax_d.axhline(1.0, color="k", ls=":", lw=0.8)
        ax_d.set_xlabel(r"$\log_{10}(R_M)$")
        ax_d.set_ylabel(r"$\lambda_1 / \lambda_{1,\mathrm{square}}$")
        ax_d.set_title(r"(d)  unweighted (λ-calib. WIP)")
        ax_d.legend(fontsize=7, frameon=False)
        ax_d.set_ylim(bottom=0.9)
        ax_d.text(
            0.02, 0.98,
            "paper: $w\\leq0.2\\Rightarrow$ ratio$\\lesssim1.05$\n"
            "ours: absolute λ still high — see notes",
            transform=ax_d.transAxes, va="top", fontsize=7, color="0.3",
        )
    else:
        ax_d.text(0.5, 0.5, "skipped (--skip-scans)", ha="center", va="center",
                  transform=ax_d.transAxes)
        ax_d.set_title("(d)")

    # (e)
    ax_e.set_xlim(0, 1)
    ax_e.set_xlabel(r"$z/L$")
    ax_e.set_ylabel(r"$n / \langle n\rangle$")
    ax_e.set_title(r"(e)  Density profiles (Eq. 57)")
    ax_e.set_ylim(bottom=0)

    # (f)
    if scan_f is not None:
        for w in W_VALUES:
            logR, ratio = scan_f[w]
            ax_f.plot(logR, ratio, color=COLORS[w], lw=2, label=f"w={w:g}")
        ax_f.axhline(1.0, color="k", ls=":", lw=0.8)
        # Paper callout: w=0.2 → ~1.15
        ax_f.axhline(1.15, color="limegreen", ls="--", lw=0.7, alpha=0.7)
        ax_f.set_xlabel(r"$\log_{10}(R_M)$")
        ax_f.set_ylabel(r"$\lambda_1 / \lambda_{1,\mathrm{square}}$")
        ax_f.set_title(r"(f)  $n(z)$-weighted scattering")
        ax_f.legend(fontsize=7, frameon=False)
        ax_f.set_ylim(bottom=0.9)
        # Log ratio at Rm=16, w=0.2 for sanity
        eig_w02 = next(e for e in eigs if abs(e.w - 0.2) < 1e-9)
        logger.info(
            "Sanity @ Rm=16, w=0.2 unweighted ratio=%.3f (paper ~O(1))",
            eig_w02.ratio,
        )
    else:
        ax_f.text(0.5, 0.5, "skipped (--skip-scans)", ha="center", va="center",
                  transform=ax_f.transAxes)
        ax_f.set_title("(f)")

    for ax in axs.ravel():
        ax.grid(True, alpha=0.25)

    fig.suptitle(
        r"Egedal et al. 2022 Fig. 10 — orbit-averaged eigenfunctions (§4)",
        fontsize=12,
    )
    fig.tight_layout()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUTPUT_DIR / "fig10_orbit_averaged.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    logger.info("wrote %s", out)

    # CSV snapshot of panel (e) densities
    csv = OUTPUT_DIR / "fig10_density_Rm16.csv"
    cols = [eigs[0].zax] + [e.ne for e in eigs]
    header = "z," + ",".join(f"n_w{w:g}" for w in W_VALUES)
    np.savetxt(csv, np.column_stack(cols), delimiter=",", header=header, comments="")
    logger.info("wrote %s", csv)

    if args.no_show:
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":
    main()
