"""Canonical source model run directory output helpers"""
from __future__ import annotations
from collections.abc import Mapping
from pathlib import Path
from typing import Any
import shutil
import yaml
from source_model_revamp.export.event_bank import DEFAULT_FILENAME as CORRELATED_EVENT_BANK_FILENAME, write_correlated_neutron_event_bank_npz
from source_model_revamp.integration.metadata import write_metadata_json
from source_model_revamp.integration.pipeline_types import SourceModelPipelineResult
from source_model_revamp.integration.config.root import SourceModelRunConfig

SOURCE_MODEL_METADATA_FILENAME = "source_model_metadata.json"
SOURCE_MODEL_ARRAYS_FILENAME = "source_model_arrays.npz"
RUN_SUMMARY_FILENAME = "run_summary.json"
RESOLVED_INPUT_FILENAME = "resolved_input.yaml"
INPUT_USED_FILENAME = "input_used.yaml"

def copy_input_used_yaml(input_path: str | Path, output_dir: str | Path) -> Path:
    """Copy the submitted input YAML into the run directory as input_used.yaml"""
    source = Path(input_path)
    destination = Path(output_dir) / INPUT_USED_FILENAME
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() != destination.resolve():
        shutil.copy2(source, destination)
   
    return destination

def write_resolved_input_yaml(input_data: Mapping[str, Any], run_config: SourceModelRunConfig, output_dir: str | Path) -> Path:
    """Write the input mapping with the fully parsed source_model block as resolved_input.yaml"""
    output_path = Path(output_dir) / RESOLVED_INPUT_FILENAME
    output_path.parent.mkdir(parents=True, exist_ok=True)
    resolved_input = dict(input_data)
    resolved_input["source_model"] = run_config.resolved_mapping(input_data["source_model"])
    with output_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(resolved_input, handle, sort_keys=False)
  
    return output_path

def source_model_output_paths(output_dir: str | Path) -> dict[str, Path]:
    """Return canonical metadata, array archive, and compact summary paths for a run directory"""
    directory = Path(output_dir)

    return {"metadata_json": directory / SOURCE_MODEL_METADATA_FILENAME, "arrays_npz": directory / SOURCE_MODEL_ARRAYS_FILENAME, "run_summary_json": directory / RUN_SUMMARY_FILENAME}

def write_source_model_outputs(pipeline_result: SourceModelPipelineResult, output_dir: str | Path) -> dict[str, str]:
    """
    Write metadata, archived arrays, correlated neutron events, and available OpenMC file source paths for one completed pipeline result
    
    OpenMC files themselves are created by the export stage and are only reported here
    """
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    paths = source_model_output_paths(directory)
    result = {key: str(path) for key, path in paths.items()}
    correlated_event_path = directory / CORRELATED_EVENT_BANK_FILENAME
    neutrons = pipeline_result.neutrons
    correlated_event_bank = None if neutrons is None else neutrons.correlated_event_bank
    if correlated_event_bank is not None:
        write_correlated_neutron_event_bank_npz(correlated_event_bank, correlated_event_path)
        result["correlated_neutron_events_npz"] = str(correlated_event_path)
    metadata = dict(pipeline_result.metadata)
    metadata.update({
        "correlated_neutron_event_bank_written": correlated_event_bank is not None,
        "correlated_neutron_event_bank_path": str(correlated_event_path) if correlated_event_bank is not None else None,
    })
    write_metadata_json(metadata, paths["metadata_json"])
    file_bundle = pipeline_result.openmc_export.file_source_bundle
    if file_bundle is not None:
        result.update({
            "openmc_source_h5": str(file_bundle.source_file_path),
            "openmc_source_validity_json": str(file_bundle.validity_file_path),
            "openmc_source_sampling_metadata_json": str(file_bundle.sampling_metadata_path),
        })
   
    return result