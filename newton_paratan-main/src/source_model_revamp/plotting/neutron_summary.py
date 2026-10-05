"""Neutron source summary plots"""
from __future__ import annotations
import argparse
from collections.abc import Mapping
from pathlib import Path
from typing import Any
import numpy as np
from source_model_revamp.plotting.better_plots import BLACK, BLUE, GRAY, GREEN, PURPLE, VERMILION, _direction_axis, _event_bank, _plot_domain_shading, _save_better, better_RC
from source_model_revamp.plotting.common import MEV_TO_J, MissingMetadata, _centers_from_edges, _configure_matplotlib, _fusion_component_label
from source_model_revamp.plotting.model_view import build_model_view

def _positive_event_arrays(metadata: Mapping[str, Any]) -> tuple[Any, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    bank = _event_bank(metadata)
    if bank is None:
        load_error = metadata.get("_correlated_event_bank_load_error")
        suffix = "" if load_error is None else f": {load_error}"
        raise MissingMetadata("correlated_neutron_events.npz is required" + suffix)
    normalized_weights = np.asarray(bank.normalized_weights, dtype=float)
    positive = np.isfinite(normalized_weights) & (normalized_weights > 0.0)
    if not np.any(positive):
        raise MissingMetadata("correlated event bank has no positive probability weights")
    energy_MeV = np.asarray(bank.energies_J, dtype=float)[positive] / MEV_TO_J
    direction = np.asarray(bank.directions, dtype=float)[positive]
    reaction = np.asarray(bank.reaction_keys)[positive]
    component = np.asarray(bank.component_labels)[positive]
    physical_weight = normalized_weights[positive] * float(bank.physical_total_rate_s)

    return bank, energy_MeV, direction, reaction, component, physical_weight

def _dt_component_masks(metadata: Mapping[str, Any], reaction: np.ndarray, component: np.ndarray) -> tuple[tuple[str, np.ndarray, str, Any], ...]:
    dt = reaction == "dt_n"
    if not np.any(dt):
        raise MissingMetadata("correlated event bank has no DT neutron events")
    view = build_model_view(metadata)
    styles = ((VERMILION, "--"), (GREEN, "-."), (PURPLE, ":"), (BLUE, (0, (5, 1))))
    values: list[tuple[str, np.ndarray, str, Any]] = [("All DT", dt, BLACK, "-")]
    known = np.zeros(dt.shape, dtype=bool)
    records = tuple(record for record in view.fusion_components if record.reaction == "dt_n")
    for index, record in enumerate(records):
        mask = dt & (component == record.label)
        if not np.any(mask):
            continue
        color, linestyle = styles[index % len(styles)]
        values.append((_fusion_component_label(record), mask, color, linestyle))
        known |= mask
    other = dt & ~known
    if np.any(other):
        values.append(("Unregistered DT components", other, GRAY, "-"))
    return tuple(values)

def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    total = float(np.sum(weights))
    if total <= 0.0:
        return float("nan")

    return float(np.sum(values * weights) / total)

def _weighted_standard_deviation(values: np.ndarray, weights: np.ndarray) -> float:
    mean = _weighted_mean(values, weights)
    total = float(np.sum(weights))
    if total <= 0.0 or not np.isfinite(mean):
        return float("nan")

    return float(np.sqrt(np.sum(weights * (values - mean) ** 2) / total))

def _weighted_pdf(values: np.ndarray, weights: np.ndarray, edges: np.ndarray) -> np.ndarray:
    histogram, _ = np.histogram(values, bins=edges, weights=weights)
    width = np.diff(edges)
    total = float(np.sum(histogram))
    if total <= 0.0:
        return np.zeros_like(histogram, dtype=float)

    return np.divide(histogram, total * width, out=np.zeros_like(histogram, dtype=float), where=width > 0.0)

def _rate_density(values: np.ndarray, weights: np.ndarray, edges: np.ndarray) -> np.ndarray:
    histogram, _ = np.histogram(values, bins=edges, weights=weights)
    width = np.diff(edges)

    return np.divide(histogram, width, out=np.zeros_like(histogram, dtype=float), where=width > 0.0)

def plot_dt_component_spectrum(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    _, energy, _, reaction, component, physical_weight = _positive_event_arrays(metadata)
    components = _dt_component_masks(metadata, reaction, component)
    dt_mask = reaction == "dt_n"
    dt_energy = energy[dt_mask]
    if dt_energy.size == 0:
        raise MissingMetadata("correlated event bank has no DT neutron events")
    span = max(float(np.max(dt_energy) - np.min(dt_energy)), 0.05)
    bin_count = max(80, int(getattr(args, "dt_spectrum_bins", 220)))
    edges = np.linspace(max(0.0, float(np.min(dt_energy)) - 0.05 * span), float(np.max(dt_energy)) + 0.05 * span, bin_count + 1)

    with plt.rc_context(better_RC):
        fig = plt.figure(figsize=(10.2, 7.3), constrained_layout=True)
        grid = fig.add_gridspec(2, 2, height_ratios=(1.0, 0.58))
        ax_rate = fig.add_subplot(grid[0, 0])
        ax_pdf = fig.add_subplot(grid[0, 1])
        ax_table = fig.add_subplot(grid[1, :])

        table_rows: list[list[str]] = []
        for label, mask, color, linestyle in components:
            local_energy = energy[mask]
            local_weight = physical_weight[mask]
            if local_energy.size == 0 or float(np.sum(local_weight)) <= 0.0:
                continue
            rate_density = _rate_density(local_energy, local_weight, edges)
            pdf = _weighted_pdf(local_energy, local_weight, edges)
            ax_rate.stairs(rate_density, edges, label=label, color=color, linestyle=linestyle, linewidth=1.6)
            ax_pdf.stairs(pdf, edges, label=label, color=color, linestyle=linestyle, linewidth=1.6)
            mean = _weighted_mean(local_energy, local_weight)
            standard_deviation = _weighted_standard_deviation(local_energy, local_weight)
            table_rows.append([
                label,
                f"{float(np.sum(local_weight)):.3e}",
                f"{mean:.3f}",
                f"{standard_deviation:.3f}",
            ])

        ax_rate.set_xlabel("Neutron energy [MeV]")
        ax_rate.set_ylabel(r"$dR/dE$ [s$^{-1}$ MeV$^{-1}$]")
        ax_rate.set_title("DT physical neutron spectrum", loc="left")
        ax_rate.grid(True, alpha=0.16)
        ax_rate.legend(loc="best", fontsize=7.5)
        ax_pdf.set_xlabel("Neutron energy [MeV]")
        ax_pdf.set_ylabel("Probability density [MeV$^{-1}$]")
        ax_pdf.set_title("Normalized DT energy PDF", loc="left")
        ax_pdf.grid(True, alpha=0.16)
        ax_pdf.legend(loc="best", fontsize=7.5)
        ax_table.axis("off")
        table = ax_table.table(cellText=table_rows, colLabels=("DT population", "Rate [s$^{-1}$]", "Mean [MeV]", "Std dev [MeV]"), loc="center", cellLoc="center", colLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(7.6)
        table.scale(1.0, 1.45)
        ax_table.set_title("DT spectral width summary", loc="left", pad=5)

        _save_better(fig, output_path)

def _weighted_cdf(values: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(values)
    sorted_values = np.asarray(values, dtype=float)[order]
    sorted_weights = np.asarray(weights, dtype=float)[order]
    total = float(np.sum(sorted_weights))
    if total <= 0.0:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)

    return sorted_values, np.cumsum(sorted_weights) / total

def _anisotropy_metrics(mu: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    total = float(np.sum(weights))
    if total <= 0.0:
        return float("nan"), float("nan")
    mean_mu = float(np.sum(weights * mu) / total)
    mean_mu2 = float(np.sum(weights * mu**2) / total)
    p2 = 0.5 * (3.0 * mean_mu2 - 1.0)

    return mean_mu, p2

def plot_dt_angular_pdf_cdf(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    _, _, direction, reaction, component, physical_weight = _positive_event_arrays(metadata)
    axis_mode = str(getattr(args, "direction_reference", "machine")).lower()
    axis_vector, _, axis_label = _direction_axis(metadata, axis_mode, getattr(args, "direction_beam_id", None))
    mu = direction @ axis_vector
    theta_deg = np.rad2deg(np.arccos(np.clip(mu, -1.0, 1.0)))
    components = _dt_component_masks(metadata, reaction, component)
    bin_count = max(36, int(getattr(args, "dt_angle_bins", 90)))
    edges = np.linspace(0.0, 180.0, bin_count + 1)
    centers = _centers_from_edges(edges)
    isotropic_pdf = 0.5 * np.sin(np.deg2rad(centers)) * np.pi / 180.0
    isotropic_cdf = 0.5 * (1.0 - np.cos(np.deg2rad(centers)))

    with plt.rc_context(better_RC):
        fig = plt.figure(figsize=(10.2, 7.1), constrained_layout=True)
        grid = fig.add_gridspec(2, 2, height_ratios=(1.0, 0.40))
        ax_pdf = fig.add_subplot(grid[0, 0])
        ax_cdf = fig.add_subplot(grid[0, 1])
        ax_table = fig.add_subplot(grid[1, :])

        ax_pdf.plot(centers, isotropic_pdf, color=GRAY, linestyle=":", linewidth=1.4, label="Isotropic reference")
        ax_cdf.plot(centers, isotropic_cdf, color=GRAY, linestyle=":", linewidth=1.4, label="Isotropic reference")
        table_rows: list[list[str]] = []
        for label, mask, color, linestyle in components:
            local_theta = theta_deg[mask]
            local_mu = mu[mask]
            local_weight = physical_weight[mask]
            if local_theta.size == 0 or float(np.sum(local_weight)) <= 0.0:
                continue
            pdf = _weighted_pdf(local_theta, local_weight, edges)
            cdf_theta, cdf = _weighted_cdf(local_theta, local_weight)
            ax_pdf.stairs(pdf, edges, color=color, linestyle=linestyle, linewidth=1.6, label=label)
            ax_cdf.plot(cdf_theta, cdf, color=color, linestyle=linestyle, linewidth=1.6, label=label)
            mean_mu, p2 = _anisotropy_metrics(local_mu, local_weight)
            table_rows.append([label, f"{mean_mu:.4f}", f"{p2:.4f}"])

        ax_pdf.set_xlabel(f"Emission angle relative to the {axis_label} [degree]")
        ax_pdf.set_ylabel("Probability density [degree$^{-1}$]")
        ax_pdf.set_title("DT angular PDF", loc="left")
        ax_pdf.set_xlim(0.0, 180.0)
        ax_pdf.grid(True, alpha=0.16)
        ax_pdf.legend(loc="best", fontsize=7.5)
        ax_cdf.set_xlabel(f"Emission angle relative to the {axis_label} [degree]")
        ax_cdf.set_ylabel("Cumulative probability")
        ax_cdf.set_title("DT angular CDF", loc="left")
        ax_cdf.set_xlim(0.0, 180.0)
        ax_cdf.set_ylim(0.0, 1.0)
        ax_cdf.grid(True, alpha=0.16)
        ax_cdf.legend(loc="best", fontsize=7.5)

        ax_table.axis("off")
        table = ax_table.table(cellText=table_rows, colLabels=("DT population", r"Mean $\mu$", r"$P_2$"), loc="center", cellLoc="center", colLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(7.6)
        table.scale(1.0, 1.45)
        ax_table.set_title("DT direction anisotropy summary", loc="left", pad=5)

        _save_better(fig, output_path)

def _channel_group_rates(metadata: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    view = build_model_view(metadata)
    grid = view.grid("full_device_population")
    dd_profiles = tuple(component.axial_neutron_rate_s for component in view.fusion_components if component.reaction == "dd_n")
    dt_profiles = tuple(component.axial_neutron_rate_s for component in view.fusion_components if component.reaction == "dt_n")
    dd = np.sum(dd_profiles, axis=0) if dd_profiles else np.zeros(grid.centers_m.size, dtype=float)
    dt = np.sum(dt_profiles, axis=0) if dt_profiles else np.zeros(grid.centers_m.size, dtype=float)
    return grid.edges_m, dd, dt, dd + dt

def plot_axial_neutron_birth_histogram(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    z_edges, dd, dt, total = _channel_group_rates(metadata)
    if z_edges.size != total.size + 1:
        raise MissingMetadata("full device neutron grid edges are unavailable")
    if not np.any(total > 0.0):
        raise MissingMetadata("axial neutron birth distribution has no positive rates")

    with plt.rc_context(better_RC):
        fig, ax = plt.subplots(figsize=(9.2, 4.9), constrained_layout=True)
        ax.stairs(total, z_edges, color=BLACK, linewidth=1.7, fill=True, alpha=0.10, label="Total neutron births")
        if np.any(dt > 0.0):
            ax.stairs(dt, z_edges, color=VERMILION, linewidth=1.55, linestyle="--", label="DT neutron births")
        if np.any(dd > 0.0):
            ax.stairs(dd, z_edges, color=BLUE, linewidth=1.55, linestyle=":", label="DD neutron births")
        positive_values = np.concatenate((total[total > 0.0], dt[dt > 0.0], dd[dd > 0.0]))
        if positive_values.size and float(np.max(positive_values) / np.min(positive_values)) > 100.0:
            ax.set_yscale("log")
        _plot_domain_shading(ax, metadata)
        ax.set_xlabel("Axial position z [m]")
        ax.set_ylabel("Neutron birth rate per axial bin [s$^{-1}$]")
        ax.set_title("Axial neutron birth distribution", loc="left")
        ax.grid(True, which="both", alpha=0.16)
        ax.legend(loc="best")
        ax.text(0.015, 0.97, f"Total {float(np.sum(total)):.3e} n/s\nDT {float(np.sum(dt)):.3e} n/s\nDD {float(np.sum(dd)):.3e} n/s", transform=ax.transAxes, va="top", ha="left", fontsize=8.2, bbox={"facecolor": "white", "alpha": 0.84, "edgecolor": "none"})

        _save_better(fig, output_path)