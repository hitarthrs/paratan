#!/usr/bin/env python3
"""
compute_neutron_source.py
=========================

CLI: build a normalized axial neutron source S_n(z) from the FBIS §2+§4
pipeline in fbis_core.py.

Example
-------
    # ~45°-class injection (R_b=2), square-ish well (w=0.1)
    python compute_neutron_source.py --Rm 16 --w 0.1 --R-b 2 --no-show

    # Near-perpendicular injection (R_b → 1⁺)
    python compute_neutron_source.py --Rm 16 --w 0.1 --R-b 1.05 --no-show

Outputs (under FBIS-py/outputs/)
--------------------------------
- sn_Rm{Rm}_w{w}_Rb{Rb}.csv
- sn_Rm{Rm}_w{w}_Rb{Rb}.png   (B, n vs Eq.57, S_n, midplane f)

FABLE: read fbis_core.py module docstring + Part D comments first.
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

from fbis_core import (
    OUTPUT_DIR,
    build_midplane_f,
    map_distribution_along_z,
    save_axial_csv,
)

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Map FBIS midplane f(v,ξ) along analytic B(z) and estimate "
            "normalized beam–beam DT neutron source S_n(z)."
        )
    )
    p.add_argument("--Rm", type=float, default=16.0, help="Mirror ratio R_m")
    p.add_argument("--w", type=float, default=0.1, help="Throat width parameter w")
    p.add_argument("--k", type=float, default=1.3, help="Throat shape exponent k")
    p.add_argument(
        "--R-b",
        type=float,
        default=2.0,
        dest="R_b",
        help="Injection pitch: ξ_b=sqrt(1-1/R_b). R_b=2 ~ 45°; R_b=1.05 ~ perp.",
    )
    p.add_argument("--v0", type=float, default=60.0, help="Injection speed (norm. units)")
    p.add_argument("--vc", type=float, default=22.0, help="Critical speed (norm. units)")
    p.add_argument("--beta", type=float, default=4.0, help="β_m mass-ratio factor")
    p.add_argument("--E-beam", type=float, default=60.0, dest="E_beam", help="Beam energy [keV]")
    p.add_argument("--n-modes", type=int, default=10, help="Pitch eigenmodes")
    p.add_argument("--Nz", type=int, default=61, help="Axial grid points")
    p.add_argument("--Nv", type=int, default=40, help="Speed grid points")
    p.add_argument("--Nxi", type=int, default=64, help="Pitch grid points (midplane build)")
    p.add_argument("--n-lam", type=int, default=64, help="Λ grid points for axial map")
    p.add_argument(
        "--max-bins",
        type=int,
        default=100,
        help="Heaviest (v,Λ) bins kept for deterministic S_n pair sum",
    )
    p.add_argument("--n-gyro", type=int, default=4, help="Gyrophase samples per pair")
    p.add_argument(
        "--skip-eq57",
        action="store_true",
        help="Skip Eq. 57 Mirror benchmark (not recommended)",
    )
    p.add_argument("--no-show", action="store_true", help="Save plots only")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args()


def plot_profiles(profiles, out_png: Path, show: bool) -> None:
    """Four-panel summary: B, n vs Eq57, S_n vs [n]^2, midplane f."""
    mid = profiles.midplane
    fig, axs = plt.subplots(2, 2, figsize=(11, 9))

    axs[0, 0].plot(profiles.zax, profiles.B, "C2")
    axs[0, 0].set_xlabel("z")
    axs[0, 0].set_ylabel(r"$B/B_0$")
    axs[0, 0].set_title(rf"$B(z)$  ($R_m={mid.Rm:g}$)")
    axs[0, 0].grid(True, alpha=0.3)

    axs[0, 1].plot(profiles.zax, profiles.n, label=r"$n$ mapped $f(v,\Lambda)$")
    if np.isfinite(profiles.n_eq57).any():
        axs[0, 1].plot(
            profiles.zax,
            profiles.n_eq57,
            "k--",
            label=rf"$n$ Eq.57 (corr={profiles.corr_n_eq57:.2f})",
        )
    axs[0, 1].axvline(profiles.n_peak_z, color="C0", ls=":", lw=0.8, alpha=0.7)
    axs[0, 1].set_xlabel("z")
    axs[0, 1].set_ylabel(r"$n / \langle n\rangle$")
    axs[0, 1].set_title("Density (mapped vs Eq. 57)")
    axs[0, 1].legend(fontsize=8)
    axs[0, 1].grid(True, alpha=0.3)

    # Red = deterministic beam–beam S_n; black dotted = [n]^2 proxy
    axs[1, 0].plot(
        profiles.zax,
        profiles.Sn_norm,
        "C3",
        label=r"$S_n$ deterministic $\sigma v$",
    )
    n2 = profiles.n**2
    n2 = n2 / (np.sum(n2) * np.diff(profiles.zax).mean())
    axs[1, 0].plot(profiles.zax, n2, "k--", lw=0.8, label=r"$[n]^2$ proxy")
    axs[1, 0].axvline(profiles.sn_peak_z, color="C3", ls=":", lw=0.8, alpha=0.7)
    axs[1, 0].set_xlabel("z")
    axs[1, 0].set_ylabel(r"$S_n(z)$ (normalized)")
    axs[1, 0].set_title(rf"Neutron source  ($R_b={mid.R_b:g}$)")
    axs[1, 0].legend(fontsize=8)
    axs[1, 0].grid(True, alpha=0.3)

    cf = axs[1, 1].contourf(mid.xax, mid.vax, mid.f, 30, cmap="plasma")
    fig.colorbar(cf, ax=axs[1, 1])
    axs[1, 1].axvline(mid.x_star, color="w", ls="--", lw=0.8)
    axs[1, 1].set_xlabel(r"$\xi$")
    axs[1, 1].set_ylabel("v (norm.)")
    axs[1, 1].set_title(r"midplane $f(v,\xi)$")

    fig.suptitle(
        r"FBIS-mapped $S_n(z)$ — $\Lambda$-measure map, deterministic pairs "
        r"(no $\Phi$)",
        fontsize=11,
    )
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    logger.info("wrote %s", out_png)
    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s | %(message)s",
    )

    mid = build_midplane_f(
        Rm=args.Rm,
        R_b=args.R_b,
        v0=args.v0,
        vc=args.vc,
        beta=args.beta,
        E_beam_keV=args.E_beam,
        n_modes=args.n_modes,
        Nv=args.Nv,
        Nxi=args.Nxi,
    )

    profiles = map_distribution_along_z(
        mid,
        Rm=args.Rm,
        w_throat=args.w,
        k=args.k,
        Nz=args.Nz,
        n_lam=args.n_lam,
        max_bins=args.max_bins,
        n_gyro=args.n_gyro,
        include_eq57=not args.skip_eq57,
    )

    tag = f"Rm{args.Rm:g}_w{args.w:g}_Rb{args.R_b:g}"
    save_axial_csv(profiles, OUTPUT_DIR / f"sn_{tag}.csv")
    plot_profiles(profiles, OUTPUT_DIR / f"sn_{tag}.png", show=not args.no_show)

    logger.info(
        "done | Rm=%.3f w=%.3f R_b=%.3g | n peak z=%.3f | Sn peak z=%.3f | "
        "n(0)/nmax=%.3f",
        args.Rm,
        args.w,
        args.R_b,
        profiles.n_peak_z,
        profiles.sn_peak_z,
        float(profiles.n[0] / profiles.n.max()) if profiles.n.max() > 0 else 0.0,
    )


if __name__ == "__main__":
    main()
