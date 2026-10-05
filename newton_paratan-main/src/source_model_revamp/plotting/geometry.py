"""Geometry and beam path diagnostic plots"""
from __future__ import annotations
import argparse
from pathlib import Path
from typing import Any, Mapping

from source_model_revamp.plotting.common import _beam_path_arrays, _beam_projected_polygon, _configure_matplotlib, _field_arrays, _scalar, _cm, _m_to_cm, _as_array, _z_edges_centers_m, _save, _plasma_radius_profile, MissingMetadata, _vessel_radius_m
from source_model_revamp.plotting.model_view import build_model_view

def plot_magnetic_field_profile(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    z_m, B_T, B_tilde = _field_arrays(metadata)
    B0 = _scalar(metadata, "B0_T")
    fig, ax1 = plt.subplots(figsize=(9.0, 5.2))
    ax1.plot(_cm(z_m), B_T, linewidth=1.8, label="B(z)")
    if B0 is not None:
        ax1.axhline(B0, linestyle="--", linewidth=1.0, label="B0")
    Bmax = _scalar(metadata, "magnetic_field_derived_Bmax_T") or _scalar(metadata, "magnetic_field_effective_Bmax_T")
    if Bmax is not None:
        ax1.axhline(Bmax, linestyle=":", linewidth=1.0, label="Bmax")
    throat = _scalar(metadata, "fitted_mirror_throat_z_m") or _scalar(metadata, "magnetic_field_derived_throat_z_m")
    if throat is not None:
        for sign in (-1.0, 1.0):
            ax1.axvline(sign * abs(_m_to_cm(throat)), linestyle=":", linewidth=0.9)
    source_min = _scalar(metadata, "source_geometry_linkage_source_z_min_m")
    source_max = _scalar(metadata, "source_geometry_linkage_source_z_max_m")
    half_length = _scalar(metadata, "half_length_m")
    if source_min is not None and source_max is not None:
        ax1.axvspan(_m_to_cm(source_min), _m_to_cm(source_max), alpha=0.12, label="source region")
    elif half_length is not None:
        ax1.axvspan(-_m_to_cm(half_length), _m_to_cm(half_length), alpha=0.06, label="source region")
    coil_z = _as_array(metadata, "coil_z_center_m", ndim=1)
    if coil_z.size:
        for i, zc in enumerate(coil_z):
            ax1.axvline(_m_to_cm(zc), alpha=0.20, linewidth=0.8, label="coil z" if i == 0 else None)
    ax1.set_xlabel("Axial position z [cm]")
    ax1.set_ylabel("Magnetic field B [T]")
    ax1.grid(True, alpha=0.3)
    ax2 = ax1.twinx()
    ax2.plot(_cm(z_m), B_tilde, linestyle="--", linewidth=1.0, alpha=0.75, label="B/B0")
    ax2.set_ylabel("B/B0")
    lines = ax1.get_lines() + ax2.get_lines()
    labels = [line.get_label() for line in lines if not line.get_label().startswith("_")]
    lines = [line for line in lines if not line.get_label().startswith("_")]
    ax1.legend(lines, labels, loc="best", fontsize=8)
    ax1.set_title("Magnetic mirror field")
    _save(fig, output_path)

def plot_beam_path_projection(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    view = build_model_view(metadata)
    if not view.beams:
        raise MissingMetadata("no enabled beam records")
    _, zc = _z_edges_centers_m(metadata)
    radius = _plasma_radius_profile(metadata, zc)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.2), sharey=True)
    species_colors = {"deuterium": "C0", "tritium": "C3"}
    for ax, coord, label in ((axes[0], 0, "x"), (axes[1], 1, "y")):
        if radius.size == zc.size:
            ax.fill_between(_cm(zc), -_cm(radius), _cm(radius), alpha=0.10, label="plasma envelope")
            ax.plot(_cm(zc), _cm(radius), linestyle="--", linewidth=1.0)
            ax.plot(_cm(zc), -_cm(radius), linestyle="--", linewidth=1.0)
        vessel = _vessel_radius_m(metadata)
        if vessel is not None:
            ax.axhline(_m_to_cm(vessel), linestyle=":", linewidth=1.0, label="vessel radius")
            ax.axhline(-_m_to_cm(vessel), linestyle=":", linewidth=1.0)
        for beam in view.beams:
            centers, _ = _beam_path_arrays(beam)
            polygon = _beam_projected_polygon(beam, coord)
            color = species_colors.get(beam.species_id)
            if polygon.size:
                ax.fill(_cm(polygon[:, 0]), _cm(polygon[:, 1]), color=color, alpha=0.14)
            species = "D" if beam.species_id == "deuterium" else "T" if beam.species_id == "tritium" else beam.species_id
            details = f"{beam.beam_id} ({species})"
            ax.plot(_cm(centers[:, 2]), _cm(centers[:, coord]), linewidth=1.8, color=color, label=details)
            ax.scatter(_cm([centers[0, 2], centers[-1, 2]]), _cm([centers[0, coord], centers[-1, coord]]), marker="x", s=45, color=color)
        ax.set_xlabel("Axial position z [cm]")
        ax.set_ylabel(f"Transverse {label} position [cm]")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)
    fig.suptitle("Enabled neutral beam paths against plasma and vessel envelopes")
    _save(fig, output_path)
