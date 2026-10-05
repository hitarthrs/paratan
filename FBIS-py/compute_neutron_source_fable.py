#!/usr/bin/env python3
"""Run neutron-source pipeline using fbis_core_fable (cell-integrated Jacobian + MC Sn)."""

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

import fbis_core_fable as core

logger = logging.getLogger(__name__)
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Neutron source with fbis_core_fable")
    p.add_argument("--Rm", type=float, default=16.0)
    p.add_argument("--w", type=float, default=0.1)
    p.add_argument("--k", type=float, default=1.3)
    p.add_argument("--R-b", type=float, default=2.0, dest="R_b")
    p.add_argument("--v0", type=float, default=60.0)
    p.add_argument("--vc", type=float, default=22.0)
    p.add_argument("--beta", type=float, default=4.0)
    p.add_argument("--E-beam", type=float, default=60.0, dest="E_beam")
    p.add_argument("--n-modes", type=int, default=10)
    p.add_argument("--Nz", type=int, default=61)
    p.add_argument("--Nv", type=int, default=40)
    p.add_argument("--Nxi", type=int, default=64)
    p.add_argument("--n-macro", type=int, default=300)
    p.add_argument("--n-pairs", type=int, default=15000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--skip-eq57", action="store_true")
    p.add_argument("--no-show", action="store_true")
    return p.parse_args()


def plot_profiles(profiles: core.AxialProfiles, out_png: Path, show: bool) -> None:
    mid = profiles.midplane
    n_peak_z = float(profiles.zax[int(np.argmax(profiles.n))]) if profiles.n.max() > 0 else 0.0
    sn_peak_z = (
        float(profiles.zax[int(np.argmax(profiles.Sn_norm))])
        if profiles.Sn_norm.max() > 0
        else 0.0
    )
    corr = float("nan")
    if np.isfinite(profiles.n_eq57).any() and np.std(profiles.n) > 0 and np.std(profiles.n_eq57) > 0:
        corr = float(np.corrcoef(profiles.n, profiles.n_eq57)[0, 1])

    fig, axs = plt.subplots(2, 2, figsize=(11, 9))
    axs[0, 0].plot(profiles.zax, profiles.B, "C2")
    axs[0, 0].set_xlabel("z")
    axs[0, 0].set_ylabel(r"$B/B_0$")
    axs[0, 0].set_title(rf"$B(z)$  ($R_m={mid.Rm:g}$) — Fable core")
    axs[0, 0].grid(True, alpha=0.3)

    axs[0, 1].plot(profiles.zax, profiles.n, label=r"$n$ mapped (cell-∫ Jacobian)")
    if np.isfinite(profiles.n_eq57).any():
        axs[0, 1].plot(profiles.zax, profiles.n_eq57, "k--", label=rf"$n$ Eq.57 (corr={corr:.2f})")
    axs[0, 1].axvline(n_peak_z, color="C0", ls=":", lw=0.8, alpha=0.7)
    axs[0, 1].set_xlabel("z")
    axs[0, 1].set_ylabel(r"$n / \langle n\rangle$")
    axs[0, 1].set_title("Density (Fable map vs Eq. 57)")
    axs[0, 1].legend(fontsize=8)
    axs[0, 1].grid(True, alpha=0.3)

    axs[1, 0].plot(profiles.zax, profiles.Sn_norm, "C3", label=r"$S_n$ MC (Fable)")
    n2 = profiles.n**2
    n2 = n2 / (np.sum(n2) * np.diff(profiles.zax).mean())
    axs[1, 0].plot(profiles.zax, n2, "k--", lw=0.8, label=r"$[n]^2$ proxy")
    axs[1, 0].axvline(sn_peak_z, color="C3", ls=":", lw=0.8, alpha=0.7)
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
    axs[1, 1].set_title(r"midplane $f(v,\xi)$ (zeroed past $\xi^*$)")

    fig.suptitle("Fable fbis_core — cell-∫ Jacobian + MC $S_n$", fontsize=11)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    logger.info("wrote %s", out_png)
    if show:
        plt.show()
    else:
        plt.close(fig)

    logger.info(
        "n peak z=%.3f | Sn peak z=%.3f | n(0)/nmax=%.3f | corr vs Eq57=%.3f",
        n_peak_z,
        sn_peak_z,
        float(profiles.n[0] / profiles.n.max()) if profiles.n.max() > 0 else 0.0,
        corr,
    )


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")

    mid = core.build_midplane_f(
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
    profiles = core.map_distribution_along_z(
        mid,
        Rm=args.Rm,
        w_throat=args.w,
        k=args.k,
        Nz=args.Nz,
        n_macro=args.n_macro,
        n_pair_samples=args.n_pairs,
        seed=args.seed,
        include_eq57=not args.skip_eq57,
    )

    tag = f"fable_Rm{args.Rm:g}_w{args.w:g}_Rb{args.R_b:g}"
    core.save_axial_csv(profiles, OUTPUT_DIR / f"sn_{tag}.csv")
    plot_profiles(profiles, OUTPUT_DIR / f"sn_{tag}.png", show=not args.no_show)


if __name__ == "__main__":
    main()
