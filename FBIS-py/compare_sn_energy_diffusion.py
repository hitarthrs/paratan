#!/usr/bin/env python3
"""
Compare mapped axial neutron sources with vs without §5.1 energy diffusion.

Builds midplane f(v,ξ) twice (Eq. 14 drag-only vs Eq. 59 + Rosenbluth), maps
each along analytic B(z), and overlays n(z) and normalized S_n(z).

Default injection: θ = 60° to B  →  R_b = 1/sin²θ = 4/3, ξ_b = cosθ = 0.5.

Example
-------
    python compare_sn_energy_diffusion.py --theta 60 --no-show
    python compare_sn_energy_diffusion.py --theta 60 --Rm 16 -v
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
)

logger = logging.getLogger(__name__)


def rb_from_theta_deg(theta_deg: float) -> float:
    """θ = angle between v and B; R_b = 1/sin²θ."""
    s = float(np.sin(np.deg2rad(theta_deg)))
    if s <= 0.0:
        raise ValueError(f"theta={theta_deg}° has sinθ≤0; cannot form R_b")
    return 1.0 / s**2


def xi_b_from_theta_deg(theta_deg: float) -> float:
    return float(np.cos(np.deg2rad(theta_deg)))


def mean_energy_over_Ebeam(mid) -> float:
    """⟨E⟩/E_beam from midplane f, using E/E_beam = (v/v0)²."""
    v, x, f = mid.vax, mid.xax, mid.f
    V, _ = np.meshgrid(v, x, indexing="ij")
    w = f * (V**2)
    den = float(np.trapezoid(np.trapezoid(w, x, axis=1), v))
    if den <= 0.0:
        return float("nan")
    num = float(np.trapezoid(np.trapezoid(w * (V / mid.v0) ** 2, x, axis=1), v))
    return num / den


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Map S_n(z) with and without energy diffusion; compare."
    )
    p.add_argument(
        "--theta",
        type=float,
        default=60.0,
        help="Injection angle θ to B [deg] (default 60 → R_b=4/3)",
    )
    p.add_argument("--Rm", type=float, default=16.0)
    p.add_argument("--w", type=float, default=0.1)
    p.add_argument("--k", type=float, default=1.3)
    p.add_argument("--v0", type=float, default=60.0)
    p.add_argument("--vc", type=float, default=22.0)
    p.add_argument("--beta", type=float, default=0.5)
    p.add_argument("--E-beam", type=float, default=60.0, dest="E_beam")
    p.add_argument("--n-modes", type=int, default=10, dest="n_modes")
    p.add_argument("--Nz", type=int, default=61)
    p.add_argument("--Nv", type=int, default=48)
    p.add_argument("--Nxi", type=int, default=64)
    p.add_argument("--n-bins", type=int, default=24, dest="n_bins")
    p.add_argument("--n-gyro", type=int, default=16, dest="n_gyro")
    p.add_argument("--ed-n-iter", type=int, default=5, dest="ed_n_iter")
    p.add_argument(
        "--ed-v-tail-mult",
        type=float,
        default=1.6,
        dest="ed_v_tail_mult",
        help="Extend v-grid to this × v0 when energy diffusion is on",
    )
    p.add_argument("--no-show", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s | %(message)s",
    )

    R_b = rb_from_theta_deg(args.theta)
    xi_b = xi_b_from_theta_deg(args.theta)
    logger.info(
        "θ=%.1f° → R_b=%.4f, ξ_b=%.4f | Rm=%.3f w=%.3f",
        args.theta,
        R_b,
        xi_b,
        args.Rm,
        args.w,
    )

    common = dict(
        Rm=args.Rm,
        R_b=R_b,
        v0=args.v0,
        vc=args.vc,
        beta=args.beta,
        E_beam_keV=args.E_beam,
        n_modes=args.n_modes,
        Nv=args.Nv,
        Nxi=args.Nxi,
    )

    logger.info("Building midplane f WITHOUT energy diffusion (Eq. 14)…")
    mid_off = build_midplane_f(**common, energy_diffusion=False)

    logger.info("Building midplane f WITH energy diffusion (Eq. 59)…")
    mid_on = build_midplane_f(
        **common,
        energy_diffusion=True,
        ed_n_iter=args.ed_n_iter,
        ed_v_tail_mult=args.ed_v_tail_mult,
    )

    e_off = mean_energy_over_Ebeam(mid_off)
    e_on = mean_energy_over_Ebeam(mid_on)
    logger.info("⟨E⟩/E_beam: no-diff=%.3f  with-diff=%.3f", e_off, e_on)

    map_kw = dict(
        Rm=args.Rm,
        w_throat=args.w,
        k=args.k,
        Nz=args.Nz,
        n_bins=args.n_bins,
        n_gyro=args.n_gyro,
        include_eq57=False,
    )

    logger.info("Mapping axial profiles (no diffusion)…")
    prof_off = map_distribution_along_z(mid_off, **map_kw)

    logger.info("Mapping axial profiles (with diffusion)…")
    prof_on = map_distribution_along_z(mid_on, **map_kw)

    z = prof_off.zax
    peak_off = float(z[int(np.argmax(prof_off.Sn_norm))])
    peak_on = float(z[int(np.argmax(prof_on.Sn_norm))])
    # L1 shape difference of normalized PDFs
    dz = float(np.diff(z).mean())
    l1 = float(np.sum(np.abs(prof_on.Sn_norm - prof_off.Sn_norm)) * dz)
    logger.info(
        "Sn peak z: no-diff=%.3f  with-diff=%.3f | ∫|ΔSn| dz=%.4f",
        peak_off,
        peak_on,
        l1,
    )

    # --- CSV ---
    tag = (
        f"theta{args.theta:g}_Rm{args.Rm:g}_w{args.w:g}_Eb{args.E_beam:g}"
    )
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_DIR / f"sn_ediff_compare_{tag}.csv"
    header = (
        "z,B_over_B0,"
        "n_no_diff,n_with_diff,"
        "Sn_norm_no_diff,Sn_norm_with_diff"
    )
    data = np.column_stack(
        [
            z,
            prof_off.B,
            prof_off.n,
            prof_on.n,
            prof_off.Sn_norm,
            prof_on.Sn_norm,
        ]
    )
    np.savetxt(csv_path, data, delimiter=",", header=header, comments="")
    logger.info("wrote %s", csv_path)

    # --- plot ---
    fig, axs = plt.subplots(2, 2, figsize=(11, 9))

    axs[0, 0].plot(z, prof_off.B, "C2")
    axs[0, 0].set_xlabel(r"$z/L$")
    axs[0, 0].set_ylabel(r"$B/B_0$")
    axs[0, 0].set_title(rf"$B(z)$  ($R_m={args.Rm:g}$, $w={args.w:g}$)")
    axs[0, 0].grid(True, alpha=0.3)

    axs[0, 1].plot(z, prof_off.n, "C0", label="no diffusion (Eq. 14)")
    axs[0, 1].plot(z, prof_on.n, "C1", ls="--", label="with diffusion (Eq. 59)")
    axs[0, 1].set_xlabel(r"$z/L$")
    axs[0, 1].set_ylabel(r"$n/\langle n\rangle$")
    axs[0, 1].set_title("Mapped density")
    axs[0, 1].legend(fontsize=8)
    axs[0, 1].grid(True, alpha=0.3)

    axs[1, 0].plot(z, prof_off.Sn_norm, "C0", label="no diffusion")
    axs[1, 0].plot(z, prof_on.Sn_norm, "C1", ls="--", label="with diffusion")
    axs[1, 0].axvline(peak_off, color="C0", ls=":", lw=0.8, alpha=0.7)
    axs[1, 0].axvline(peak_on, color="C1", ls=":", lw=0.8, alpha=0.7)
    axs[1, 0].set_xlabel(r"$z/L$")
    axs[1, 0].set_ylabel(r"$S_n(z)$ (normalized)")
    axs[1, 0].set_title(
        rf"Neutron source  ($\theta={args.theta:g}^\circ$, "
        rf"$R_b={R_b:.3f}$, $\xi_b={xi_b:.3f}$)"
    )
    axs[1, 0].legend(fontsize=8)
    axs[1, 0].grid(True, alpha=0.3)

    # Speed-integrated midplane density vs ξ (shape of pitch content)
    def n_of_xi(mid) -> np.ndarray:
        v, f = mid.vax, mid.f
        return np.trapezoid(f * (v[:, None] ** 2), v, axis=0)

    nxi_off = n_of_xi(mid_off)
    nxi_on = n_of_xi(mid_on)
    s_off = float(np.trapezoid(nxi_off, mid_off.xax))
    s_on = float(np.trapezoid(nxi_on, mid_on.xax))
    if s_off > 0.0:
        nxi_off = nxi_off / s_off
    if s_on > 0.0:
        nxi_on = nxi_on / s_on
    axs[1, 1].plot(mid_off.xax, nxi_off, "C0", label="no diffusion")
    axs[1, 1].plot(mid_on.xax, nxi_on, "C1", ls="--", label="with diffusion")
    axs[1, 1].axvline(mid_off.x_star, color="k", ls="--", lw=0.8, label=r"$\xi^*$")
    axs[1, 1].axvline(xi_b, color="0.4", ls=":", lw=0.8, label=rf"$\xi_b={xi_b:.2f}$")
    axs[1, 1].set_xlabel(r"$\xi$")
    axs[1, 1].set_ylabel(r"$\int f\,v^2\,dv$ (norm.)")
    axs[1, 1].set_title(
        rf"midplane pitch  "
        rf"($\langle E\rangle/E_b$: {e_off:.2f} → {e_on:.2f})"
    )
    axs[1, 1].legend(fontsize=7)
    axs[1, 1].grid(True, alpha=0.3)

    fig.suptitle(
        "Energy diffusion vs drag-only — axial neutron source comparison",
        fontsize=12,
    )
    fig.tight_layout()
    png_path = OUTPUT_DIR / f"sn_ediff_compare_{tag}.png"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    logger.info("wrote %s", png_path)
    if args.no_show:
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":
    main()
