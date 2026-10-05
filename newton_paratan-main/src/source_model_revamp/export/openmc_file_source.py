"""
OpenMC HDF5 FileSource export for correlated neutron event banks

The event probability distribution is represented by source particle weights with mean one while the physical neutron rate remains separate metadata used for transport normalization
"""
from __future__ import annotations
from dataclasses import asdict, dataclass, is_dataclass
from enum import Enum
from hashlib import sha256
import importlib
import json
import os
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
from typing import Any, Mapping
import numpy as np
from source_model_revamp.constants import EV_TO_J
from source_model_revamp.neutrons.events import CorrelatedNeutronEventBank
from source_model_revamp.pipeline_errors import OpenMCExportDependencyError, PipelineConservationError

SOURCE_SCHEMA_IDENTIFIER = "openmc_source_bank_correlated_lab_events_v1"
SOURCE_PARTICLE_WEIGHT_MODEL = "exported_event_count_times_normalized_probability"

@dataclass(frozen=True)
class OpenMCFileSourceMetadata:
    """
    Measured metadata and round trip validation results for one correlated OpenMC source bundle
    
    Positions are written in cm energies in eV source particle weights encode normalized event probability and physical_total_rate_s remains the absolute neutron rate in s⁻¹
    """
    source_file_path: str
    validity_file_path: str
    sampling_metadata_path: str
    source_file_sha256: str
    openmc_version: str | None
    physical_total_rate_s: float
    file_source_strength: float
    input_event_count: int
    positive_weight_event_count: int
    excluded_zero_weight_event_count: int
    event_probability_weight_sum: float
    source_particle_weight_sum: float
    source_particle_weight_mean: float
    hdf5_structure_valid: bool
    roundtrip_position_valid: bool
    roundtrip_direction_valid: bool
    roundtrip_energy_valid: bool
    roundtrip_time_valid: bool
    roundtrip_probability_weights_valid: bool
    roundtrip_particle_identity_valid: bool
    roundtrip_weighted_moments_valid: bool
    joint_lab_energy_direction_correlation_preserved: bool
    maximum_position_roundtrip_error_cm: float
    maximum_direction_roundtrip_error: float
    maximum_energy_relative_roundtrip_error: float
    maximum_time_roundtrip_error_s: float
    maximum_particle_weight_relative_roundtrip_error: float
    energy_direction_covariance_before_eV: tuple[float, float, float]
    energy_direction_covariance_after_eV: tuple[float, float, float]
    energy_direction_covariance_difference_eV: tuple[float, float, float]
    source_schema_identifier: str
    source_particle_weight_model_identifier: str
    position_units: str
    energy_units: str
    source_time_s: float
    physical_neutron_rate_interpretation: str
    probability_weight_interpretation: str
    roundtrip_valid: bool
    real_openmc_readback_performed: bool
    real_openmc_roundtrip_valid: bool | None

    def to_dict(self) -> dict[str, Any]:
        """Return the metadata as a plain dictionary"""
        return asdict(self)

@dataclass(frozen=True)
class OpenMCFileSourceBundle:
    """
    Paths OpenMC FileSource object and measured metadata for one correlated source bundle
    
    The bundle contains the HDF5 source file plus JSON validity and sampling metadata sidecars
    """
    source_file_path: Path
    validity_file_path: Path
    sampling_metadata_path: Path
    file_source: Any
    metadata: OpenMCFileSourceMetadata

def _require_h5py():
    """Import h5py or raise the export dependency error used by the pipeline"""
    try:
        return importlib.import_module("h5py")
    except ImportError as exc:
        raise OpenMCExportDependencyError("correlated OpenMC HDF5 export requires the optional openmc_export dependency group containing h5py") from exc

def _json_safe(value: Any) -> Any:
    """
    Convert metadata values to JSON compatible containers and scalar types
    
    Unlike the NPZ sidecar helper this path rejects nonfinite floating metadata instead of replacing it
    """
    if is_dataclass(value) and not isinstance(value, type):
        return _json_safe(asdict(value))
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return _json_safe(value.value)
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError("metadata contains a nonfinite floating point value")
    return value

def _write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write one deterministic JSON sidecar with finite metadata values"""
    path.write_text(json.dumps(_json_safe(payload), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")

def _particle_type_is_neutron(value: Any) -> bool:
    """Return whether a read back OpenMC particle identifier represents a neutron"""
    name = getattr(value, "name", None)
    if isinstance(name, str) and name.strip().lower() == "neutron":
        return True
    try:
        if int(value) == 0:
            return True
    except (TypeError, ValueError):
        pass
    return str(value).strip().lower() in {"0", "neutron", "particletype.neutron"}

def _openmc_neutron_particle(openmc_module) -> Any:
    """Return the OpenMC neutron particle identifier supported by the imported OpenMC module"""
    particle_type = getattr(openmc_module, "ParticleType", None)
    if particle_type is not None and hasattr(particle_type, "NEUTRON"):
        return particle_type.NEUTRON
    return "neutron"

def _source_file_sha256(path: Path) -> str:
    """Return the hexadecimal SHA 256 digest of one written source file"""
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def _hdf5_source_structure_valid(h5py_module, path: Path, expected_count: int) -> bool:
    """
    Check the minimum HDF5 source bank structure expected by this export path
    
    The file must identify itself as a source and contain one source_bank row per exported positive probability event with position direction energy time weight and particle fields
    """
    try:
        with h5py_module.File(path, "r") as handle:
            filetype = handle.attrs.get("filetype")
            if hasattr(filetype, "decode"):
                filetype = filetype.decode("utf-8")
            if str(filetype).strip().lower() != "source":
                return False
            if "source_bank" not in handle:
                return False
            bank = handle["source_bank"]
            names = set(bank.dtype.names or ())
            return bank.shape == (expected_count,) and {"r", "u", "E", "time", "wgt", "particle"}.issubset(names)
    except (OSError, ValueError, TypeError):
        return False

def _readback_arrays(openmc_module, path: Path) -> dict[str, np.ndarray]:
    """
    Read a written source file through OpenMC and return the correlated fields as NumPy arrays
    
    This read back path supplies the measured evidence used by the round trip checks
    """
    if not hasattr(openmc_module, "read_source_file"):
        raise RuntimeError("OpenMC correlated source export requires openmc.read_source_file")
    particles = tuple(openmc_module.read_source_file(path))
    return {
        "position_cm": np.asarray([particle.r for particle in particles], dtype=float),
        "direction": np.asarray([particle.u for particle in particles], dtype=float),
        "energy_eV": np.asarray([particle.E for particle in particles], dtype=float),
        "time_s": np.asarray([particle.time for particle in particles], dtype=float),
        "weight": np.asarray([particle.wgt for particle in particles], dtype=float),
        "particle": np.asarray([particle.particle for particle in particles], dtype=object),
    }

def _maximum_absolute_error(actual: np.ndarray, expected: np.ndarray) -> float:
    """Return the maximum absolute array error or infinity when shapes are incompatible"""
    if actual.shape != expected.shape or actual.size == 0:
        return float("inf")
  
    return float(np.max(np.abs(actual - expected)))

def _maximum_relative_error(actual: np.ndarray, expected: np.ndarray) -> float:
    """Return the maximum elementwise relative array error or infinity when shapes are incompatible"""
    if actual.shape != expected.shape or actual.size == 0:
        return float("inf")
    denominator = np.maximum(np.abs(expected), np.finfo(float).tiny)
  
    return float(np.max(np.abs(actual - expected) / denominator))

def _weighted_energy_direction_covariance(energy_eV: np.ndarray, direction: np.ndarray, probabilities: np.ndarray) -> np.ndarray:
    """
    Return the probability weighted covariance between neutron energy and each lab direction component
    
    The result has shape `(3,)` and is used to verify preservation of the joint lab energy direction state
    """
    mean_energy = float(np.sum(probabilities * energy_eV))
    mean_direction = np.sum(probabilities[:, None] * direction, axis=0)

    return np.sum(probabilities[:, None] * (energy_eV - mean_energy)[:, None] * (direction - mean_direction[None, :]), axis=0)

def _reaction_or_component_summary(labels: np.ndarray, probabilities: np.ndarray, physical_total_rate_s: float) -> list[dict[str, Any]]:
    """
    Summarize event probability and physical rate by reaction key or component label
    
    Each physical subgroup rate is `probability_weight_sum * physical_total_rate_s`
    """
    summary = []
    for label in np.unique(labels):
        selected = labels == label
        probability = float(np.sum(probabilities[selected]))
        summary.append({"label": str(label), "source_particle_count": int(np.count_nonzero(selected)), "probability_weight_sum": probability, "physical_rate_s": probability * physical_total_rate_s})
  
    return summary

def _validate_filenames(source_filename: str, validity_filename: str, sampling_metadata_filename: str) -> None:
    """Validate the three distinct bundle filenames and their required file suffixes"""
    for filename, suffix in ((source_filename, ".h5"), (validity_filename, ".json"), (sampling_metadata_filename, ".json")):
        candidate = Path(filename)
        if (not filename or candidate.name != filename or candidate.parent != Path(".") or candidate.suffix.lower() != suffix):
            raise ValueError(f"output filename must be one {suffix} filename without a directory")
    if len({source_filename, validity_filename, sampling_metadata_filename}) != 3:
        raise ValueError("OpenMC source output filenames must be distinct")

def _atomic_publish_bundle_directory(staging: Path, output: Path, filenames: tuple[str, str, str], *, overwrite_existing: bool) -> None:
    """
    Atomically replace the dedicated output directory with a fully validated staged bundle
    
    When replacement is allowed the previous directory is moved aside until the new bundle has been installed successfully
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    backup = output.with_name(f".{output.name}.previous_bundle")
    if backup.exists():
        raise FileExistsError(f"stale atomic publication backup exists: {backup}")
    if output.exists():
        entries = tuple(output.iterdir())
        unexpected = tuple(entry for entry in entries if entry.name not in filenames)
        if unexpected:
            raise FileExistsError("atomic OpenMC publication requires a dedicated output directory")
        existing = tuple(output / name for name in filenames if (output / name).exists())
        if existing and not overwrite_existing:
            names = ", ".join(path.name for path in existing)
            raise FileExistsError(f"refusing to overwrite existing outputs: {names}")
        os.replace(output, backup)
        try:
            os.replace(staging, output)
        except BaseException:
            os.replace(backup, output)
            raise
        shutil.rmtree(backup)
        return
    os.replace(staging, output)

def _create_bundle(event_bank: CorrelatedNeutronEventBank, *, output_directory: str | Path, source_filename: str, validity_filename: str, sampling_metadata_filename: str, origin_cm: tuple[float, float, float], source_time_s: float, overwrite_existing: bool, openmc_module) -> OpenMCFileSourceBundle:
    """
    Create validate and install one correlated OpenMC HDF5 source bundle
    
    Only positive probability events are exported
    For `N` exported rows with probabilities `p_i` the stored OpenMC particle weights are `w_i = N p_i`, so uniform row selection has mean particle weight one and preserves the event probability measure
    The absolute neutron rate remains `physical_total_rate_s` in the sidecar metadata and is not encoded in FileSource strength
    """
    _validate_filenames(source_filename, validity_filename, sampling_metadata_filename)
    if not np.isfinite(source_time_s) or source_time_s < 0.0:
        raise ValueError("source_time_s must be finite and nonnegative")
    origin = np.asarray(origin_cm, dtype=float)
    if origin.shape != (3,) or np.any(~np.isfinite(origin)):
        raise ValueError("origin_cm must contain three finite values")
    for name in ("SourceParticle", "FileSource", "write_source_file"):
        if not hasattr(openmc_module, name):
            raise RuntimeError(f"OpenMC module is missing required API {name}")
    h5py_module = _require_h5py()
    probabilities_all = np.asarray(event_bank.normalized_weights, dtype=float)
    positive = probabilities_all > 0.0
    positive_count = int(np.count_nonzero(positive))
    if positive_count < 1:
        raise ValueError("event bank contains no positive probability events")
    probabilities = probabilities_all[positive]
    probability_sum = float(np.sum(probabilities))
    if not np.isclose(probability_sum, 1.0, rtol=0.0, atol=2.0e-12):
        raise PipelineConservationError("positive correlated event probabilities do not sum to one", stage="openmc_export", check_name="openmc_file_source_probability_normalization",)
    # Uniform row selection with w_i = N p_i preserves the normalized event probability measure
    stored_weights = positive_count * probabilities
    expected_position_cm = (100.0 * np.asarray(event_bank.positions_m, dtype=float)[positive] + origin[None, :])
    expected_direction = np.asarray(event_bank.directions, dtype=float)[positive]
    expected_energy_eV = np.asarray(event_bank.energies_J, dtype=float)[positive] / EV_TO_J
    expected_time_s = np.full(positive_count, float(source_time_s), dtype=float)
    neutron_particle = _openmc_neutron_particle(openmc_module)
    particles = [openmc_module.SourceParticle(r=tuple(expected_position_cm[index]), u=tuple(expected_direction[index]), E=float(expected_energy_eV[index]), time=float(source_time_s), wgt=float(stored_weights[index]), delayed_group=0, surf_id=0, particle=neutron_particle) for index in range(positive_count)]
    output = Path(output_directory)
    final_source = output / source_filename
    final_validity = output / validity_filename
    final_sampling = output / sampling_metadata_filename
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=f".{output.name}.staging_", dir=output.parent) as temporary:
        staging = Path(temporary)
        staged_source = staging / source_filename
        staged_validity = staging / validity_filename
        staged_sampling = staging / sampling_metadata_filename
        # Validate the staged HDF5 file before replacing any existing complete bundle
        openmc_module.write_source_file(particles, staged_source)
        structure_valid = _hdf5_source_structure_valid(h5py_module, staged_source, positive_count)
        readback = _readback_arrays(openmc_module, staged_source)
        count_valid = readback["energy_eV"].shape == (positive_count,)
        maximum_position_error = _maximum_absolute_error(readback["position_cm"], expected_position_cm)
        maximum_direction_error = _maximum_absolute_error(readback["direction"], expected_direction)
        maximum_energy_relative_error = _maximum_relative_error(readback["energy_eV"], expected_energy_eV)
        maximum_time_error = _maximum_absolute_error(readback["time_s"], expected_time_s)
        maximum_weight_relative_error = _maximum_relative_error(readback["weight"], stored_weights)
        position_valid = bool(count_valid and maximum_position_error <= 1.0e-12)
        direction_valid = bool(count_valid and maximum_direction_error <= 2.0e-12)
        energy_valid = bool(count_valid and maximum_energy_relative_error <= 2.0e-15)
        time_valid = bool(count_valid and maximum_time_error == 0.0)
        weights_valid = bool( count_valid and np.allclose(readback["weight"], stored_weights, rtol=2.0e-15, atol=1.0e-14) and np.isclose(np.sum(readback["weight"]), positive_count, rtol=0.0, atol=2.0e-12 * positive_count))
        particle_valid = bool(count_valid and all(_particle_type_is_neutron(value) for value in readback["particle"]))
        expected_probability = stored_weights / float(np.sum(stored_weights))
        readback_probability = (readback["weight"] / float(np.sum(readback["weight"])) if count_valid and float(np.sum(readback["weight"])) > 0.0 else np.asarray([], dtype=float))
        covariance_before = _weighted_energy_direction_covariance(expected_energy_eV, expected_direction, expected_probability)
        if readback_probability.shape == probabilities.shape:
            covariance_after = _weighted_energy_direction_covariance(readback["energy_eV"], readback["direction"], readback_probability)
        else:
            covariance_after = np.full(3, np.nan)
        covariance_difference = covariance_after - covariance_before
        weighted_moments_valid = bool(count_valid  and readback_probability.shape == probabilities.shape  and np.allclose(covariance_after, covariance_before, rtol=2.0e-15, atol=2.0e-12))
        joint_state_valid = bool(count_valid and direction_valid and energy_valid and np.allclose(np.column_stack((readback["energy_eV"], readback["direction"])), np.column_stack((expected_energy_eV, expected_direction)), rtol=2.0e-15, atol=2.0e-12))
        roundtrip_valid = bool(structure_valid and position_valid and direction_valid and energy_valid and time_valid and weights_valid and particle_valid and weighted_moments_valid and joint_state_valid)
        if not roundtrip_valid:
            raise PipelineConservationError(f"OpenMC correlated source file failed measured readback validation: structure={structure_valid}, count={count_valid}, position={position_valid}, direction={direction_valid}, energy={energy_valid}, time={time_valid}, weights={weights_valid}, particle={particle_valid}, weighted_moments={weighted_moments_valid}, joint_state={joint_state_valid}, covariance_difference_eV={covariance_difference.tolist()}", stage="openmc_export", check_name="openmc_file_source_roundtrip_validation")
        openmc_module.FileSource(str(staged_source), strength=1.0)
        source_hash = _source_file_sha256(staged_source)
        metadata = OpenMCFileSourceMetadata(
            source_file_path=str(final_source),
            validity_file_path=str(final_validity),
            sampling_metadata_path=str(final_sampling),
            source_file_sha256=source_hash,
            openmc_version=str(getattr(openmc_module, "__version__", "unknown")),
            physical_total_rate_s=float(event_bank.physical_total_rate_s),
            file_source_strength=1.0,
            input_event_count=event_bank.event_count,
            positive_weight_event_count=positive_count,
            excluded_zero_weight_event_count=event_bank.event_count - positive_count,
            event_probability_weight_sum=probability_sum,
            source_particle_weight_sum=float(np.sum(stored_weights)),
            source_particle_weight_mean=float(np.mean(stored_weights)),
            hdf5_structure_valid=structure_valid,
            roundtrip_position_valid=position_valid,
            roundtrip_direction_valid=direction_valid,
            roundtrip_energy_valid=energy_valid,
            roundtrip_time_valid=time_valid,
            roundtrip_probability_weights_valid=weights_valid,
            roundtrip_particle_identity_valid=particle_valid,
            roundtrip_weighted_moments_valid=weighted_moments_valid,
            joint_lab_energy_direction_correlation_preserved=joint_state_valid,
            maximum_position_roundtrip_error_cm=maximum_position_error,
            maximum_direction_roundtrip_error=maximum_direction_error,
            maximum_energy_relative_roundtrip_error=maximum_energy_relative_error,
            maximum_time_roundtrip_error_s=maximum_time_error,
            maximum_particle_weight_relative_roundtrip_error=maximum_weight_relative_error,
            energy_direction_covariance_before_eV=tuple(float(value) for value in covariance_before),
            energy_direction_covariance_after_eV=tuple(float(value) for value in covariance_after),
            energy_direction_covariance_difference_eV=tuple(float(value) for value in covariance_difference),
            source_schema_identifier=SOURCE_SCHEMA_IDENTIFIER,
            source_particle_weight_model_identifier=SOURCE_PARTICLE_WEIGHT_MODEL,
            position_units="cm",
            energy_units="eV",
            source_time_s=float(source_time_s),
            physical_neutron_rate_interpretation=("separate_neutrons_per_second_metadata_not_probability_weight"),
            probability_weight_interpretation=("uniform_row_selection_with_mean_one_source_particle_weights"),
            roundtrip_valid=roundtrip_valid,
            real_openmc_readback_performed=True,
            real_openmc_roundtrip_valid=roundtrip_valid,
        )
        validity_payload = metadata.to_dict()
        validity_payload["validation_status"] = "measured_and_passed_with_real_openmc"
        sampling_payload = {
            "source_schema_identifier": SOURCE_SCHEMA_IDENTIFIER,
            "source_particle_weight_model_identifier": SOURCE_PARTICLE_WEIGHT_MODEL,
            "position_units": "cm",
            "energy_units": "eV",
            "source_time_s": float(source_time_s),
            "event_model": event_bank.metadata.get("event_model"),
            "physical_total_rate_s": float(event_bank.physical_total_rate_s),
            "physical_rate_units": "neutrons_per_second",
            "physical_neutron_rate_interpretation": metadata.physical_neutron_rate_interpretation,
            "probability_weight_interpretation": metadata.probability_weight_interpretation,
            "source_probability_weight_sum": probability_sum,
            "source_particle_weight_convention": SOURCE_PARTICLE_WEIGHT_MODEL,
            "source_particle_weight_sum": float(np.sum(stored_weights)),
            "source_particle_weight_mean": float(np.mean(stored_weights)),
            "input_event_count": event_bank.event_count,
            "exported_positive_weight_event_count": positive_count,
            "excluded_zero_weight_event_count": event_bank.event_count - positive_count,
            "reaction_summary": _reaction_or_component_summary(np.asarray(event_bank.reaction_keys)[positive], probabilities, float(event_bank.physical_total_rate_s)),
            "component_summary": _reaction_or_component_summary(np.asarray(event_bank.component_labels)[positive], probabilities, float(event_bank.physical_total_rate_s)),
            "event_bank_metadata": event_bank.metadata,
        }
        _write_json(staged_validity, validity_payload)
        _write_json(staged_sampling, sampling_payload)
        _atomic_publish_bundle_directory(staging, output, (source_filename, validity_filename, sampling_metadata_filename), overwrite_existing=overwrite_existing)
    file_source = openmc_module.FileSource(str(final_source), strength=1.0)

    return OpenMCFileSourceBundle(source_file_path=final_source, validity_file_path=final_validity, sampling_metadata_path=final_sampling, file_source=file_source, metadata=metadata)

def create_openmc_correlated_file_source(event_bank: CorrelatedNeutronEventBank, *, output_directory: str | Path, source_filename: str = "source.h5", validity_filename: str = "source_file_validity.json", sampling_metadata_filename: str = "source_sampling_metadata.json", origin_cm: tuple[float, float, float] = (0.0, 0.0, 0.0), source_time_s: float = 0.0, overwrite_existing: bool = False) -> OpenMCFileSourceBundle:
    """
    Create one correlated OpenMC FileSource bundle using the installed OpenMC Python module
    
    The returned FileSource has strength one and the physical neutron rate is retained separately in the bundle metadata
    """
    try:
        openmc_module = importlib.import_module("openmc")
    except ImportError as exc:
        raise OpenMCExportDependencyError("correlated source export requires the OpenMC Python module") from exc
    
    return _create_bundle(
        event_bank,
        output_directory=output_directory,
        source_filename=source_filename,
        validity_filename=validity_filename,
        sampling_metadata_filename=sampling_metadata_filename,
        origin_cm=origin_cm,
        source_time_s=source_time_s,
        overwrite_existing=overwrite_existing,
        openmc_module=openmc_module,
    )

__all__ = [
    "OpenMCFileSourceBundle",
    "OpenMCFileSourceMetadata",
    "create_openmc_correlated_file_source",
]