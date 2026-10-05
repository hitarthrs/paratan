"""
ParaTAN <-> source_model_revamp neutron source coupling

This file is the coupling between the paraTAN OpenMC model builder and the source_model_revamp physics blocks
"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any
import re
from source_model_revamp.integration.metadata import metadata_value
from source_model_revamp.integration.pipeline import run_source_model_pipeline
from source_model_revamp.integration.validation import SourceConstructionError
from source_model_revamp.integration.output_writer import copy_input_used_yaml, write_resolved_input_yaml, write_source_model_outputs
from source_model_revamp.integration.pipeline_types import SourceModelPipelineResult
from source_model_revamp.integration.config.root import SourceModelRunConfig, source_model_config_from_root_mapping

DEFAULT_SOURCE_INFORMATION: dict[str, Any] = {
    "settings": {
        "particles_per_batch": 10000,
        "batches": 10,
        "statepoint_frequency": 10,
        "weight_windows": False,
        "photon_transport": False,
        "geometry_plot": False,
        "tallies": True,
    },
}

@dataclass(frozen=True)
class ParaTANSourceModelBuildResult:
    """Result returned by the ParaTAN source hook"""
    sources: tuple[Any, ...]
    pipeline_result: SourceModelPipelineResult
    run_config: SourceModelRunConfig
    source_information: dict[str, Any]
    metadata: dict[str, Any]
    output_files: dict[str, str] = field(default_factory=dict)

    @property
    def total_neutron_rate_s(self) -> float:
        return float(self.pipeline_result.metadata["total_neutron_rate_s"])

    @property
    def metadata_file(self) -> str | None:
        return self.output_files.get("metadata_json")

def _as_mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}

def _deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    result = {str(key): value for key, value in dict(base).items()}
    for key, value in override.items():
        key = str(key)
        if isinstance(value, Mapping) and isinstance(result.get(key), Mapping):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value

    return result

def resolve_source_information(input_data: Mapping[str, Any], source_data: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Resolve ParaTAN/OpenMC run settings for the model builder"""
    resolved = _deep_merge(DEFAULT_SOURCE_INFORMATION, _as_mapping(input_data.get("source_information")))
    if source_data is not None:
        if not isinstance(source_data, Mapping):
            raise ValueError("source_data must be a mapping")
        resolved = _deep_merge(resolved, source_data)
    resolved["settings"] = source_settings_from_source_information(resolved)
    return resolved

def source_settings_from_source_information(source_information: Mapping[str, Any]) -> dict[str, Any]:
    """Return source_information[settings] merged with safe defaults"""
    default_settings = DEFAULT_SOURCE_INFORMATION["settings"]
    requested = _as_mapping(source_information.get("settings"))
    unknown = sorted(set(requested) - set(default_settings))
    if unknown:
        raise ValueError(f"unknown source_information.settings fields: {', '.join(unknown)}")
    settings = _deep_merge(default_settings, requested)
    settings["particles_per_batch"] = int(settings.get("particles_per_batch", 10000))
    settings["batches"] = int(settings.get("batches", 10))
    settings["statepoint_frequency"] = int(settings.get("statepoint_frequency", max(1, settings["batches"])))
    settings["weight_windows"] = bool(settings.get("weight_windows", False))
    settings["photon_transport"] = bool(settings.get("photon_transport", False))
    settings["geometry_plot"] = bool(settings.get("geometry_plot", False))
    settings["tallies"] = bool(settings.get("tallies", True))
    if settings["particles_per_batch"] <= 0:
        raise ValueError("settings.particles_per_batch must be positive")
    if settings["batches"] <= 0:
        raise ValueError("settings.batches must be positive")
    if settings["statepoint_frequency"] <= 0:
        raise ValueError("settings.statepoint_frequency must be positive")
    
    return settings

def sanitize_run_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    cleaned = cleaned.strip("._-")
    return cleaned or "source_model_revamp_run"

def source_model_run_name(input_data: Mapping[str, Any]) -> str:
    run = _as_mapping(input_data.get("run"))
    if run.get("name"):
        return sanitize_run_name(str(run["name"]))
    source_model = _as_mapping(input_data.get("source_model"))
    if source_model.get("run_name"):
        return sanitize_run_name(str(source_model["run_name"]))
    if source_model.get("name"):
        return sanitize_run_name(str(source_model["name"]))
    
    return "source_model_revamp_run"

def default_paratan_run_dir(input_data: Mapping[str, Any], *, base_dir: str | Path = "") -> Path:
    """Default run directory: <base_dir>/runs/<run_name>"""
    return Path(base_dir) / "runs" / source_model_run_name(input_data)

def build_source_model_openmc_source(input_data: Mapping[str, Any], source_data: Mapping[str, Any] | None = None, *, output_dir: str | Path | None = None, write_outputs: bool = True, input_path: str | Path | None = None, progress_callback=None) -> ParaTANSourceModelBuildResult:
    """Build the selected OpenMC source from the source model block"""
    if "source_model" not in input_data or not isinstance(input_data["source_model"], Mapping):
        raise ValueError("ParaTAN source model runs require a source_model block in the main input YAML")
    source_information = resolve_source_information(input_data, source_data)

    out_dir = Path.cwd() if output_dir is None else Path(output_dir)
    run_config = source_model_config_from_root_mapping(input_data, strict=True)
    configured_source_dir = Path(run_config.openmc_export.source_output_directory)
    if not configured_source_dir.is_absolute():
        configured_source_dir = out_dir / configured_source_dir
    run_config = replace(
        run_config,
        openmc_export=replace(
            run_config.openmc_export,
            source_output_directory=str(configured_source_dir),
        ),
    )
    if write_outputs:
        out_dir.mkdir(parents=True, exist_ok=True)
        write_resolved_input_yaml(input_data, run_config, out_dir)
        if input_path is not None:
            copy_input_used_yaml(input_path, out_dir)
    pipeline_result = run_source_model_pipeline(
        run_config,
        write_openmc_source=True,
        progress_callback=progress_callback,
    )
    export_result = pipeline_result.openmc_export
    if export_result is None:
        raise RuntimeError("source model pipeline did not return an OpenMC export stage")

    if export_result.file_source_bundle is None:
        raise SourceConstructionError("missing_correlated_file_source", "ParaTAN requires the correlated HDF5 FileSource bundle", quantity="openmc_file_source_bundle")
    bundle = export_result.file_source_bundle
    sources = (bundle.file_source,)
    metadata = dict(pipeline_result.metadata)
    metadata.update(metadata_value(bundle.metadata.to_dict()))
    metadata.update(
        {
            "paratan_source_model_link_enabled": True,
            "paratan_source_model_link_backend": "source_model_revamp",
            "paratan_geometry_plot_requested": bool(source_information["settings"]["geometry_plot"]),
        }
    )
    pipeline_result = replace(pipeline_result, metadata=metadata)
    files: dict[str, str] = {}
    if write_outputs:
        files = write_source_model_outputs(pipeline_result, out_dir)
        files["resolved_input_yaml"] = str(out_dir / "resolved_input.yaml")
        input_used_file = out_dir / "input_used.yaml"
        if input_used_file.exists():
            files["input_used_yaml"] = str(input_used_file)

    return ParaTANSourceModelBuildResult(sources=sources, pipeline_result=pipeline_result, run_config=run_config, source_information=source_information, metadata=metadata_value(metadata), output_files=files)
