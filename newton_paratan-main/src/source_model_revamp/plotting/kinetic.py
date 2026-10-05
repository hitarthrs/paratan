"""Kinetic and electrostatic diagnostic plots"""
from __future__ import annotations
import argparse
from pathlib import Path
from typing import Any, Mapping
from source_model_revamp.plotting.common import _configure_matplotlib, _z_edges_centers_m, _as_array, MissingMetadata, _save, _cm, _scalar, _log10_positive, _select_z_by_B, _m_to_cm
from source_model_revamp.plotting.model_view import build_model_view
import numpy as np

def plot_quasineutrality_residual(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    _, zc = _z_edges_centers_m(metadata)
    ion = _as_array(metadata, "modal_phi_z_ion_density_profile_m3", ndim=1)
    elec = _as_array(metadata, "modal_phi_z_electron_density_profile_m3", ndim=1)
    if ion.size != zc.size or elec.size != zc.size:
        raise MissingMetadata("missing ion or electron density profiles")
    residual = elec - ion
    scale = np.maximum(np.maximum(np.abs(elec), np.abs(ion)), np.finfo(float).tiny)
    rel = residual / scale
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(9.0, 6.6), sharex=True)
    ax0.plot(_cm(zc), ion, marker="o", label="ΣZ n_i")
    ax0.plot(_cm(zc), elec, marker="s", linestyle="--", label="n_e")
    ax0.set_ylabel("Density [m$^{-3}$]")
    ax0.grid(True, alpha=0.3)
    ax0.legend(loc="best", fontsize=8)
    ax1.plot(_cm(zc), residual, marker="o", label="n_e - ΣZ n_i")
    ax1b = ax1.twinx()
    ax1b.plot(_cm(zc), rel, marker="s", linestyle="--", label="relative residual")
    ax1.set_xlabel("Axial position z [cm]")
    ax1.set_ylabel("Residual [m$^{-3}$]")
    ax1b.set_ylabel("Relative residual")
    ax1.grid(True, alpha=0.3)
    lines = ax1.get_lines() + ax1b.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines], loc="best", fontsize=8)
    fig.suptitle("Quasineutrality residual")
    _save(fig, output_path)

def plot_eta_lambda_and_eigenmodes(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    lam = _as_array(metadata, "modal_eta_lambda_grid", ndim=1)
    eta = _as_array(metadata, "modal_eta_of_lambda_grid", ndim=1)
    eta_grid = _as_array(metadata, "modal_physical_eta_grid", ndim=1)
    eigfunc = _as_array(metadata, "modal_physical_eigenfunctions_eta")
    eigvals = _as_array(metadata, "modal_eigenvalues", ndim=1)
    if lam.size == 0 or eta.size != lam.size or eta_grid.size == 0 or eigfunc.ndim != 2:
        raise MissingMetadata("missing modal eta lambda or eigenfunction arrays")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    axes[0].plot(lam, eta, linewidth=1.8, label="η(Λ)")
    lb = _scalar(metadata, "modal_lambda_boundary")
    eb = _scalar(metadata, "modal_eta_trapped_passing_boundary")
    if lb is not None:
        axes[0].axvline(lb, linestyle=":", label="Λ loss boundary")
    if eb is not None:
        axes[0].axhline(eb, linestyle="--", label="η trapped boundary")
    axes[0].set_xlabel("Λ")
    axes[0].set_ylabel("η")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="best", fontsize=8)
    max_modes = min(eigfunc.shape[0], args.max_modes)
    for j in range(max_modes):
        label = f"j={j + 1}"
        if eigvals.size > j:
            label += f", λ={eigvals[j]:.3g}"
        axes[1].plot(eta_grid, eigfunc[j], label=label)
    axes[1].set_xlabel("η")
    axes[1].set_ylabel("I_j(η)")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(loc="best", fontsize=8)
    fig.suptitle("Shared D/T density weighted modal η(Λ) map and eigenfunctions")
    _save(fig, output_path)

def plot_local_velocity_space(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    view = build_model_view(metadata)
    if not view.fast_species:
        raise MissingMetadata("no active fast ion species")
    B_tilde = _as_array(metadata, "B_tilde_centers", ndim=1)
    confined_grid = view.grid("confined_kinetic")
    if B_tilde.size != confined_grid.centers_m.size:
        raise MissingMetadata("B_tilde_centers does not match the confined kinetic grid")
    indices = _select_z_by_B(metadata, args.local_b_over_b0)
    mirror_ratio = _scalar(metadata, "mirror_ratio") or _scalar(metadata, "fitted_mirror_ratio")
    row_count = len(view.fast_species)
    column_count = len(indices)
    fig = plt.figure(figsize=(4.0 * column_count + 0.8, 3.8 * row_count))
    grid = fig.add_gridspec(row_count, column_count + 1, width_ratios=[1.0] * column_count + [0.045])
    axes = np.empty((row_count, column_count), dtype=object)
    for row, species in enumerate(view.fast_species):
        beams = tuple(beam for beam in view.beams if beam.species_id == species.species_id)
        component_energies = tuple(beam.component_energies_J for beam in beams if beam.component_energies_J.size)
        if component_energies:
            reference_speed = float(np.sqrt(2.0 * max(float(np.max(values)) for values in component_energies) / species.mass_kg))
        else:
            reference_speed = float(np.max(species.local_speed_centers_m_s))
        V, XI = np.meshgrid(species.local_speed_centers_m_s, species.pitch_centers, indexing="ij")
        X = V * XI / reference_speed
        Y = V * np.sqrt(np.maximum(1.0 - XI**2, 0.0)) / reference_speed
        log_values = tuple(_log10_positive(species.local_distribution_z_v_pitch[index]) for index in indices)
        lower = min(float(np.nanmin(values)) for values in log_values)
        upper = max(float(np.nanmax(values)) for values in log_values)
        levels = np.linspace(lower, upper if upper > lower else lower + 1.0, 36)
        contour = None
        for column, (index, values) in enumerate(zip(indices, log_values, strict=True)):
            ax = fig.add_subplot(grid[row, column])
            axes[row, column] = ax
            contour = ax.contourf(X, Y, values, levels=levels)
            Bt = float(B_tilde[index])
            if mirror_ratio is not None and Bt < mirror_ratio:
                xi_abs = np.sqrt(max(0.0, 1.0 - Bt / mirror_ratio))
                vline = np.linspace(0.0, 1.05 * float(np.max(species.local_speed_centers_m_s / reference_speed)), 100)
                ax.plot(vline * xi_abs, vline * np.sqrt(max(0.0, 1.0 - xi_abs**2)), linestyle="--", linewidth=1.0, label="magnetic loss boundary")
                ax.plot(-vline * xi_abs, vline * np.sqrt(max(0.0, 1.0 - xi_abs**2)), linestyle="--", linewidth=1.0)
            ax.set_xlabel(r"$v_\parallel/v_{\mathrm{ref}}$")
            if column == 0:
                ax.set_ylabel(r"$v_\perp/v_{\mathrm{ref}}$")
                ax.legend(loc="best", fontsize=7)
            ax.set_title(f"{species.symbol}, z={_m_to_cm(confined_grid.centers_m[index]):.1f} cm, B/B0={Bt:.2f}")
            ax.set_aspect("equal", adjustable="box")
        if contour is not None:
            colorbar_axis = fig.add_subplot(grid[row, -1])
            colorbar = fig.colorbar(contour, cax=colorbar_axis)
            colorbar.set_label(rf"log10 f, $v_{{\mathrm{{ref}}}}$={reference_speed / 1.0e6:.3g} Mm/s")
    fig.suptitle("Species resolved local fast ion velocity distributions")
    _save(fig, output_path)
