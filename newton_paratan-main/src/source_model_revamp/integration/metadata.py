"""Metadata serialization, array archiving, and compact run summary helpers"""
from __future__ import annotations
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any
import json
import re
import numpy as np

ARRAY_REFERENCE_KEY = "array_file"
ARRAY_REFERENCE_DATASET_KEY = "dataset"

def metadata_value(value: Any) -> Any:
    """Convert dataclasses, NumPy values, mappings, and sequences into JSON compatible values"""
    if is_dataclass(value):
        return metadata_value(asdict(value))
    if isinstance(value, dict):
        return {str(key): metadata_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [metadata_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)

def require_metadata_keys(metadata: dict[str, Any], required_keys: list[str] | tuple[str, ...]) -> None:
    """Require every named canonical run metadata key"""
    missing = [key for key in required_keys if key not in metadata]
    if missing:
        raise ValueError(f"metadata is missing required keys: {missing}")

REQUIRED_RUN_METADATA_KEYS = (
    "source_model_backend",
    "source_model_workflow_model",
    "geometry_model",
    "magnetic_field_model",
    "B0_T",
    "mirror_ratio",
    "plasma_radius_m",
    "half_length_m",
    "beam_source_model",
    "beam_attenuation_model",
    "beam_fast_birth_rate_s",
    "kinetic_model",
    "electron_temperature_mode",
    "electron_temperature_solved_keV",
    "kinetic_convergence_status",
    "fusion_model",
    "active_fusion_channels",
    "total_neutron_rate_s",
    "total_fusion_power_W",
    "neutron_spectrum_model",
    "neutron_angular_model",
    "openmc_export_model",
    "openmc_source_representation",
    "num_openmc_sources",
)

def _dataset_name(prefix: str, arrays: dict[str, np.ndarray]) -> str:
    """Return a sanitized unique NPZ dataset name for one metadata path"""
    base = re.sub(r"[^A-Za-z0-9_]+", "_", prefix.replace(".", "__")).strip("_") or "array"
    if base not in arrays:
        return base
    index = 2
    while f"{base}_{index}" in arrays:
        index += 1

    return f"{base}_{index}"

def _archive_array(array: np.ndarray, *, prefix: str, arrays: dict[str, np.ndarray], archive_name: str) -> dict[str, Any]:
    """Store one array in the pending NPZ archive and return its JSON reference record"""
    value = np.asarray(array)
    dataset = _dataset_name(prefix, arrays)
    arrays[dataset] = value
  
    return {ARRAY_REFERENCE_KEY: archive_name, ARRAY_REFERENCE_DATASET_KEY: dataset, "shape": list(value.shape), "dtype": str(value.dtype)}

def _archive_arrays(value: Any, *, prefix: str, arrays: dict[str, np.ndarray], archive_name: str) -> Any:
    """
    Recursively replace arrays with NPZ references
    
    Numeric lists and tuples with at least sixteen entries are archived as arrays while smaller sequences remain inline JSON values
    """
    if is_dataclass(value):
        return _archive_arrays(asdict(value), prefix=prefix, arrays=arrays, archive_name=archive_name)
    if isinstance(value, dict):
        return {str(key): _archive_arrays(item, prefix=f"{prefix}.{key}" if prefix else str(key), arrays=arrays, archive_name=archive_name) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _archive_array(value, prefix=prefix, arrays=arrays, archive_name=archive_name)
    if isinstance(value, (tuple, list)):
        # Archive long numeric sequences but keep compact records readable in JSON
        if value and all(isinstance(item, (bool, int, float, np.generic)) for item in value):
            array = np.asarray(value)
            if array.ndim > 0 and array.size >= 16:
                return _archive_array(array, prefix=prefix, arrays=arrays, archive_name=archive_name)
        return [_archive_arrays(item, prefix=f"{prefix}.{index}" if prefix else str(index), arrays=arrays, archive_name=archive_name) for index, item in enumerate(value)]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
  
    return str(value)

def compact_run_summary(metadata: dict[str, Any]) -> dict[str, Any]:
    """Return the compact user facing result summary written to run_summary.json"""
    keys = (
        "calculation_valid",
        "numerically_converged",
        "source_constructible",
        "source_generated",
        "source_file_valid",
        "failed_convergence_checks",
        "kinetic_convergence_status",
        "operating_point_converged",
        "B0_T",
        "mirror_ratio",
        "plasma_radius_m",
        "half_length_m",
        "electron_temperature_mode",
        "electron_temperature_solved_keV",
        "electron_temperature_final_qualification_passed",
        "modal_effective_ion_temperature_keV",
        "fast_D_confined_volume_average_density_m3",
        "fast_T_confined_volume_average_density_m3",
        "modal_confined_ion_collisional_confinement_time_s",
        "beam_power_W",
        "beam_energy_keV",
        "beam_injection_angle_deg",
        "beam_fast_birth_rate_s",
        "beam_deposited_birth_power_W",
        "total_fusion_power_W",
        "total_neutron_rate_s",
        "modal_n_lambda_grid",
        "modal_n_eta_grid",
        "modal_n_square_basis_modes",
        "modal_n_physical_modes",
        "openmc_source_representation",
        "num_openmc_sources",
    )
  
    return {key: metadata.get(key) for key in keys if key in metadata}

def _array_archive_path(output_path: Path) -> Path:
    """Return the companion NPZ archive path for one metadata JSON path"""
    if output_path.name == "source_model_metadata.json":
        return output_path.with_name("source_model_arrays.npz")
  
    return output_path.with_name(f"{output_path.stem}_arrays.npz")

def write_metadata_json(metadata: dict[str, Any], path: str | Path, *, require_core_keys: bool = True, write_array_archive: bool = True, write_summary: bool = True) -> None:
    """
    Write canonical metadata JSON, an optional companion NPZ array archive, and an optional compact run summary
    
    Large numerical arrays are referenced from JSON by archive filename and dataset name
    """
    if require_core_keys:
        require_metadata_keys(metadata, REQUIRED_RUN_METADATA_KEYS)
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {}
    archive_path = _array_archive_path(output_path)
    serializable = _archive_arrays(metadata, prefix="", arrays=arrays, archive_name=archive_path.name) if write_array_archive else metadata_value(metadata)
    if arrays:
        np.savez_compressed(archive_path, **arrays)
    elif archive_path.exists():
        archive_path.unlink()
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(serializable, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if write_summary:
        summary_path = output_path.with_name("run_summary.json")
        with summary_path.open("w", encoding="utf-8") as handle:
            json.dump(metadata_value(compact_run_summary(metadata)), handle, indent=2, sort_keys=True)
            handle.write("\n")