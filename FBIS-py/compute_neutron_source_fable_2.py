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

FABLE: read fbis_core.py module docstring first — units, Λ mapping, and what
is intentionally NOT included (Φ, absolute Amps, radial profile).
S_n is a deterministic D–T pair sum (no RNG): identical inputs ⇒ identical CSVs.
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
    # Magnetic geometry (paper Eq. 56)
    p.add_argument("--Rm", type=float, default=16.0, help="Mirror ratio R_m")
    p.add_argument("--w", type=float, default=0.1, help="Throat width parameter w")
    p.add_argument("--k", type=float, default=1.3, help="Throat shape exponent k")
    # NBI / distribution
    p.add_argument(
        "--R-b",
        type=float,
        default=2.0,
        dest="R_b",
        help="Injection pitch parameter: ξ_b=sqrt(1-1/R_b). "
             "R_b=2 ~ 45°; R_b=1.05 ~ near-perpendicular (ξ_b small).",
    )
    p.add_argument("--v0", type=float, default=60.0, help="Injection speed (norm. units)")
    p.add_argument(
        "--vc",
        type=float,
        default=22.0,
        help="Critical speed, SAME normalized units as --v0. Physical anchor: "
             "E_c ≈ 20 T_e (Eq. 16) ⇒ vc = v0·sqrt(20·T_e/E_beam). "
             "Default 22 with v0=60, E_beam=60 keV ⇒ T_e ≈ 0.4 keV; "
             "move it if you change --E-beam.",
    )
    p.add_argument(
        "--beta",
        type=float,
        default=0.5,
        help="β_m = z_eff·m_i/(2 m_f); 0.5 for pure DT (paper, below Eq. 3)",
    )
    p.add_argument(
        "--E-beam",
        type=float,
        default=60.0,
        dest="E_beam",
        help="Beam energy [keV] for Bosch–Hale σ(E) conversion",
    )
    p.add_argument("--n-modes", type=int, default=10, help="Number of pitch eigenmodes")
    # Numerics
    p.add_argument("--Nz", type=int, default=61, help="Axial grid points")
    p.add_argument("--Nv", type=int, default=40, help="Speed grid points")
    p.add_argument("--Nxi", type=int, default=64, help="Pitch grid points")
    p.add_argument("--n-bins", type=int, default=24, help="(|v∥|, v⊥) bins per z for S_n")
    p.add_argument("--n-gyro", type=int, default=16, help="Gyrophase quadrature nodes for S_n (16≈3%%, 32<1%%)")
    p.add_argument(
        "--skip-eq57",
        action="store_true",
        help="Skip slow Eq. 57 Mirror benchmark",
    )
    p.add_argument("--no-show", action="store_true", help="Save plots only")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args()


def plot_profiles(profiles, out_png: Path, show: bool) -> None:
    """Four-panel summary: B, n, S_n, midplane f."""
    mid = profiles.midplane
    fig, axs = plt.subplots(2, 2, figsize=(11, 9))

    axs[0, 0].plot(profiles.zax, profiles.B, "C2")
    axs[0, 0].set_xlabel("z")
    axs[0, 0].set_ylabel(r"$B/B_0$")
    axs[0, 0].set_title(rf"$B(z)$  ($R_m={mid.Rm}$, $w$ throat)")
    axs[0, 0].grid(True, alpha=0.3)

    axs[0, 1].plot(profiles.zax, profiles.n, label=r"$n$ from mapped $f$")
    if np.isfinite(profiles.n_eq57).any():
        axs[0, 1].plot(profiles.zax, profiles.n_eq57, "--", label=r"$n$ Eq. 57 ($I_1$ only)")
    axs[0, 1].set_xlabel("z")
    axs[0, 1].set_ylabel(r"$n / \langle n\rangle$")
    axs[0, 1].set_title("Density shape")
    axs[0, 1].legend(fontsize=8)
    axs[0, 1].grid(True, alpha=0.3)

    axs[1, 0].plot(profiles.zax, profiles.Sn_norm, "C3")
    # Also show naive [n]^2 proxy (renormalized) for comparison.
    n2 = profiles.n**2
    n2 = n2 / (np.sum(n2) * np.diff(profiles.zax).mean())
    axs[1, 0].plot(profiles.zax, n2, "k--", lw=0.8, label=r"$[n]^2$ proxy")
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
        "FBIS-mapped axial neutron source "
        r"(no $\Phi$; shape only)",
        fontsize=12,
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

    # Loss-cone Rm for midplane eigenmodes should match the field Rm.
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
        n_bins=args.n_bins,
        n_gyro=args.n_gyro,
        include_eq57=not args.skip_eq57,
    )

    tag = f"Rm{args.Rm:g}_w{args.w:g}_Rb{args.R_b:g}"
    csv_path = OUTPUT_DIR / f"sn_{tag}.csv"
    png_path = OUTPUT_DIR / f"sn_{tag}.png"
    save_axial_csv(profiles, csv_path)
    plot_profiles(profiles, png_path, show=not args.no_show)

    # One-line summary for logs / Fable.
    z_peak = float(profiles.zax[int(np.argmax(profiles.Sn_norm))])
    logger.info(
        "done | Rm=%.3f w=%.3f R_b=%.3g | Sn peak at z=%.3f | "
        "n_mapped range [%.2f, %.2f]",
        args.Rm,
        args.w,
        args.R_b,
        z_peak,
        float(profiles.n.min()),
        float(profiles.n.max()),
    )


if __name__ == "__main__":
    main()
