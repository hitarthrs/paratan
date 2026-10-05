"""Neutron source diagnostic plots"""
from __future__ import annotations
import argparse
from pathlib import Path
from typing import Any, Mapping
import numpy as np
from source_model_revamp.plotting.common import _configure_matplotlib, MissingMetadata, _save, _cm, _fusion_component_label
from source_model_revamp.plotting.model_view import build_model_view

def plot_channel_split_axial_source(metadata: Mapping[str, Any], output_path: Path, args: argparse.Namespace) -> None:
    plt = _configure_matplotlib()
    view = build_model_view(metadata)
    grid = view.grid("full_device_population")
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    made = False
    for component in view.fusion_components:
        if component.reaction not in ("dd_n", "dt_n") or not np.any(component.axial_neutron_rate_s > 0.0):
            continue
        reaction = "DD" if component.reaction == "dd_n" else "DT"
        ax.plot(_cm(grid.centers_m), component.axial_neutron_rate_s, marker="o", label=f"{_fusion_component_label(component)} ({reaction})")
        made = True
    if not made:
        raise MissingMetadata("no positive fusion neutron component profiles")
    ax.set_xlabel("Axial position z [cm]")
    ax.set_ylabel("Neutron source rate per bin [n/s]")
    ax.set_title("Axial neutron source split by reactant populations")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    _save(fig, output_path)