#!/usr/bin/env python3
"""
Same layout as the fixed-core injection scan (n(z), Sn(z) curves),
but legend uses beam angle θ to B instead of R_b.

S_n is Monte Carlo in fbis_core_fable — we raise n_macro / n_pairs and
average several seeds to kill jaggedness.
"""

from __future__ import annotations

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

# MC settings (was ~280 / 12k / 1 seed)
N_MACRO = 3200
N_PAIRS = 105000
SEEDS = (0, 1, 2, 3, 4)  # average 5 realizations


def rb_from_theta_deg(theta_deg: float) -> float:
    """θ = angle between v and B; R_b = 1/sin²θ."""
    s = np.sin(np.deg2rad(theta_deg))
    return float(1.0 / s**2)


def map_sn_averaged(
    mid: core.MidplaneDistribution,
    *,
    Rm: float,
    w: float,
    k: float,
    Nz: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Density from one run (deterministic); S_n = mean over SEEDS with fat MC.
    Returns zax, n, Sn_norm.
    """
    sn_stack: list[np.ndarray] = []
    n_ref: np.ndarray | None = None
    zax: np.ndarray | None = None

    for seed in SEEDS:
        prof = core.map_distribution_along_z(
            mid,
            Rm=Rm,
            w_throat=w,
            k=k,
            Nz=Nz,
            n_macro=N_MACRO,
            n_pair_samples=N_PAIRS,
            seed=seed,
            include_eq57=False,
        )
        if n_ref is None:
            n_ref = prof.n.copy()
            zax = prof.zax.copy()
        sn_stack.append(prof.Sn.copy())  # raw before norm; re-norm after mean
        logger.info("  seed=%d done", seed)

    assert zax is not None and n_ref is not None
    sn_mean = np.mean(np.stack(sn_stack, axis=0), axis=0)
    dz = float(np.diff(zax).mean())
    integ = float(np.sum(sn_mean) * dz)
    sn_norm = sn_mean / integ if integ > 0 else sn_mean
    return zax, n_ref, sn_norm


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")

    cases = [
        (77.0, "C0"),
        (45.0, "C1"),
        (30.0, "C2"),
        (89.0, "C3"),
    ]

    Rm, w, k = 16.0, 0.1, 1.3
    Nz = 61
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))

    logger.info(
        "MC: n_macro=%d  n_pairs=%d  seeds=%s (avg)",
        N_MACRO,
        N_PAIRS,
        SEEDS,
    )

    for theta, color in cases:
        Rb = rb_from_theta_deg(theta)
        xi_b = float(np.cos(np.deg2rad(theta)))
        logger.info("θ=%.0f°  ξ_b=%.2f  R_b=%.2f", theta, xi_b, Rb)

        mid = core.build_midplane_f(Rm=Rm, R_b=Rb, Nv=36, Nxi=56, n_modes=10)
        zax, n, sn = map_sn_averaged(mid, Rm=Rm, w=w, k=k, Nz=Nz)

        label = rf"$\theta={theta:.0f}^\circ$  ($\xi_b={xi_b:.2f}$)"
        axs[0].plot(zax, n, color=color, label=label)
        axs[1].plot(zax, sn, color=color, label=label)

    axs[0].set_xlabel(r"$z$")
    axs[0].set_ylabel(r"$n/\langle n\rangle$")
    axs[0].set_title("normalized density")
    axs[0].set_xlim(0, 1)
    axs[0].grid(True, alpha=0.3)
    axs[0].legend(fontsize=9)

    axs[1].set_xlabel(r"$z$")
    axs[1].set_ylabel(r"$S_n$")
    axs[1].set_title(
        rf"neutron source  (MC avg: ${N_MACRO}$ macro × ${N_PAIRS}$ pairs × "
        rf"${len(SEEDS)}$ seeds)"
    )
    axs[1].set_xlim(0, 1)
    axs[1].grid(True, alpha=0.3)
    axs[1].legend(fontsize=9)

    fig.suptitle(
        r"FBIS core: injection-angle scan  ($\theta$ = angle to $B$; $90^\circ$=perp)",
        fontsize=12,
    )
    fig.tight_layout()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUTPUT_DIR / "sn_angle_scan.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    logger.info("wrote %s", out)


if __name__ == "__main__":
    main()
