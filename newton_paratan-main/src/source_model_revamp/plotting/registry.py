"""Plot registry, metadata loading, and plot execution"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from typing import Any, Mapping
import numpy as np
from source_model_revamp.export.event_bank import DEFAULT_FILENAME as CORRELATED_EVENT_BANK_FILENAME, read_correlated_neutron_event_bank_npz
from source_model_revamp.plotting.better_plots import plot_better_axial_neutron_source_density, plot_better_axial_profiles, plot_better_current_balance, plot_better_fast_ion_density_profiles, plot_better_neutron_energy_direction, plot_better_radial_neutron_source_end_cell, plot_better_radial_neutron_source_midplane, plot_better_radial_neutron_source_throat, plot_better_source_cloud, plot_full_device_cross_section
from source_model_revamp.plotting.common import METADATA_FILENAME, MissingMetadata, PlotSpec
from source_model_revamp.plotting.geometry import plot_beam_path_projection, plot_magnetic_field_profile
from source_model_revamp.plotting.kinetic import plot_eta_lambda_and_eigenmodes, plot_local_velocity_space, plot_quasineutrality_residual
from source_model_revamp.plotting.neutron_summary import plot_axial_neutron_birth_histogram, plot_dt_angular_pdf_cdf, plot_dt_component_spectrum
from source_model_revamp.plotting.neutrons import plot_channel_split_axial_source

def _plot_specs() -> list[PlotSpec]:
    return [
        PlotSpec("magnetic_field_profile", "geometry", "magnetic_field_profile.png", plot_magnetic_field_profile, "B(z) and B/B0"),
        PlotSpec("full_device_cross_section", "geometry", "full_device_cross_section.png", plot_full_device_cross_section, "Full device shell, flux tube, source region, and beam cross section"),
        PlotSpec("beam_path_projection", "geometry", "beam_path_projection.png", plot_beam_path_projection, "Enabled beam paths against plasma and vessel"),
        PlotSpec("quasineutrality_residual", "axial", "quasineutrality_residual.png", plot_quasineutrality_residual, "Quasineutrality residual"),
        PlotSpec("eta_lambda_eigenmodes", "modal", "eta_lambda_eigenmodes.png", plot_eta_lambda_and_eigenmodes, "Eta lambda map and eigenmodes"),
        PlotSpec("local_velocity_space", "modal", "local_velocity_space.png", plot_local_velocity_space, "Local velocity distributions selected by B/B0"),
        PlotSpec("dt_component_spectrum", "neutron", "dt_component_spectrum.png", plot_dt_component_spectrum, "Component resolved DT neutron spectrum and width summary"),
        PlotSpec("dt_angular_pdf_cdf", "neutron", "dt_angular_pdf_cdf.png", plot_dt_angular_pdf_cdf, "DT angular PDF and CDF relative to the selected axis"),
        PlotSpec("axial_neutron_birth_histogram", "neutron", "axial_neutron_birth_histogram.png", plot_axial_neutron_birth_histogram, "Simplified axial neutron birth histogram"),
        PlotSpec("channel_split_axial_source", "neutron", "channel_split_axial_source.png", plot_channel_split_axial_source, "Fusion channel split source"),
        PlotSpec("better_axial_profiles", "better", "axial_profiles.png", plot_better_axial_profiles, "Full device D/T axial magnetic, electrostatic, density, and source profiles"),
        PlotSpec("better_fast_ion_density_profiles", "better", "fast_ion_density_profiles.png", plot_better_fast_ion_density_profiles, "Full device fast D and T axial density profiles with magnetic turning estimates"),
        PlotSpec("better_axial_neutron_source_density", "better", "axial_neutron_source_density.png", plot_better_axial_neutron_source_density, "Volume normalized total, DD, and DT axial neutron production"),
        PlotSpec("better_radial_neutron_source_midplane", "better", "radial_neutron_source_midplane.png", plot_better_radial_neutron_source_midplane, "Volume normalized radial neutron source density near the magnetic midplane"),
        PlotSpec("better_radial_neutron_source_throat", "better", "radial_neutron_source_throat.png", plot_better_radial_neutron_source_throat, "Volume normalized radial neutron source density near the right mirror throat"),
        PlotSpec("better_radial_neutron_source_end_cell", "better", "radial_neutron_source_end_cell.png", plot_better_radial_neutron_source_end_cell, "Volume normalized radial neutron source density near the right end cell"),
        PlotSpec("better_source_cloud", "better", "neutron_source_cloud.png", plot_better_source_cloud, "Three dimensional neutron source cloud"),
        PlotSpec("better_neutron_energy_direction", "better", "neutron_energy_direction.png", plot_better_neutron_energy_direction, "Correlated neutron energy and direction distributions"),
        PlotSpec("better_current_balance", "better", "current_balance.png", plot_better_current_balance, "Pass 11 represented confined plasma current balance"),
    ]

def _selected_specs(groups: list[str]) -> list[PlotSpec]:
    specs = _plot_specs()
    identifiers = [spec.name for spec in specs]
    filenames = [spec.filename for spec in specs]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("duplicate plot registry identifier")
    if len(filenames) != len(set(filenames)):
        raise ValueError("duplicate plot output filename")
    if any(Path(filename).suffix.lower() != ".png" for filename in filenames):
        raise ValueError("registered plot output filenames must use PNG")
    normalized = {g.strip().lower() for g in groups}
    if not normalized or "all" in normalized:
        return specs
    return [spec for spec in specs if spec.group in normalized or spec.name in normalized]

def _resolve_array_references(value: Any, run_dir: Path, archives: dict[Path, dict[str, np.ndarray]]) -> Any:
    if isinstance(value, Mapping):
        array_file = value.get("array_file")
        dataset = value.get("dataset")
        if isinstance(array_file, str) and isinstance(dataset, str):
            archive_path = run_dir / array_file
            if archive_path not in archives:
                with np.load(archive_path, allow_pickle=False) as archive:
                    archives[archive_path] = {key: np.asarray(archive[key]) for key in archive.files}
            if dataset not in archives[archive_path]:
                raise KeyError(f"dataset {dataset!r} is not present in {archive_path}")
            return archives[archive_path][dataset]
        return {str(key): _resolve_array_references(item, run_dir, archives) for key, item in value.items()}
    if isinstance(value, list):
        return [_resolve_array_references(item, run_dir, archives) for item in value]
    return value

def load_metadata(run_dir: Path) -> dict[str, Any]:
    path = run_dir / METADATA_FILENAME
    if not path.is_file():
        raise FileNotFoundError(f"expected {path}")
    with path.open("r", encoding="utf-8") as handle:
        metadata = _resolve_array_references(json.load(handle), run_dir, {})
    for target_key, source_key in (("modal_phi_z_electron_density_profile_m3", "electron_density_profile_m3"), ("modal_phi_z_ion_density_profile_m3", "total_positive_charge_profile_m3")):
        if metadata.get(target_key) is None and metadata.get(source_key) is not None:
            metadata[target_key] = metadata[source_key]
    metadata["_metadata_run_dir"] = str(run_dir)
    event_bank_path = run_dir / CORRELATED_EVENT_BANK_FILENAME
    if event_bank_path.is_file():
        try:
            metadata["_correlated_event_bank"] = (read_correlated_neutron_event_bank_npz(event_bank_path))
        except Exception as exc:
            metadata["_correlated_event_bank_load_error"] = (f"{type(exc).__name__}: {exc}")
        else:
            metadata["_correlated_event_bank_path"] = str(event_bank_path)
    return metadata

def run_plots(args: argparse.Namespace) -> list[dict[str, Any]]:
    run_dir = Path(args.run_dir)
    metadata = load_metadata(run_dir)
    output_dir = Path(args.output_dir) if args.output_dir else run_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, Any]] = []
    for spec in _selected_specs(args.groups):
        output_path = output_dir / spec.filename
        entry: dict[str, Any] = {"name": spec.name, "group": spec.group, "description": spec.description, "path": str(output_path)}
        try:
            spec.function(metadata, output_path, args)
        except MissingMetadata as exc:
            entry.update({"status": "skipped", "reason": str(exc), "path": None})
        except Exception as exc:
            entry.update({"status": "failed", "reason": f"{type(exc).__name__}: {exc}", "path": None})
            if not args.keep_going:
                manifest.append(entry)
                break
        else:
            entry.update({"status": "made", "reason": None})
        manifest.append(entry)
    manifest_path = output_dir / "plot_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump({"metadata_file": str(run_dir / METADATA_FILENAME), "plots": manifest}, handle, indent=2)
        handle.write("\n")
        
    return manifest
