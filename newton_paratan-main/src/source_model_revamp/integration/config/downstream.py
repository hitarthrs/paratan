"""Fusion, neutron, and OpenMC source export configuration"""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from source_model_revamp.integration.config.common import _bool, _check_unknown, _integer, _number, canonical_model_name
from source_model_revamp.integration.config.constants import *
from source_model_revamp.radial_profiles import KOTELNIKOV_PARABOLIC_FLUX_K1, canonical_radial_profile_model

@dataclass(frozen=True)
class FusionConfig:
    """Fusion pair integration controls used by the downstream fusion stage"""
    include_fast_fast: bool = True
    num_gyroangle_points: int = 3
    max_pair_states_per_batch: int = 262144
    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, strict: bool = False) -> "FusionConfig":
        """Parse fusion pair selection, gyroangle quadrature size, and pair batching controls"""
        allowed = {'include_fast_fast', 'num_gyroangle_points', 'max_pair_states_per_batch'}
        _check_unknown(data, allowed, "source_model.fusion", strict)

        return cls(include_fast_fast=_bool(data.get('include_fast_fast'), True), num_gyroangle_points=_integer(data.get('num_gyroangle_points'), 'fusion.num_gyroangle_points', minimum=1, default=3), max_pair_states_per_batch=_integer(data.get('max_pair_states_per_batch'), 'fusion.max_pair_states_per_batch', minimum=1, default=262144))

@dataclass(frozen=True)
class NeutronConfig:
    """
    Neutron spectrum and correlated event controls
    
    Energy bounds are in MeV and the selected angular model is applied in the center of momentum frame by the neutron event path
    """
    model: str = "distribution_kinematics_spectrum"
    energy_min_MeV: float = 0.0
    energy_max_MeV: float = 20.0
    energy_bins: int = 80
    angular_model: str = "endf_b_viii1_evaluated_cm"
    num_gyroangle_points: int = 2
    num_emission_directions: int = 4
    allow_out_of_range: bool = False
    max_pair_states_per_batch: int = 131072
    correlated_event_count: int = 10000
    correlated_event_seed: int = 12345
    radial_source_profile_model: str = KOTELNIKOV_PARABOLIC_FLUX_K1

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, strict: bool = False) -> "NeutronConfig":
        """Parse neutron energy, angular, batching, event sampling, and radial source controls"""
        allowed = {"model", "energy_min_MeV", "energy_max_MeV", "energy_bins", "angular_model", "num_gyroangle_points", "num_emission_directions", "allow_out_of_range", "max_pair_states_per_batch", "correlated_event_count", "correlated_event_seed", "radial_source_profile_model"}
        _check_unknown(data, allowed, "source_model.neutrons", strict)
        energy_min = _number(data.get("energy_min_MeV"), "neutrons.energy_min_MeV", nonnegative=True, default=0.0,)
        energy_max = _number(data.get("energy_max_MeV"), "neutrons.energy_max_MeV", positive=True, default=20.0)
        if energy_max <= energy_min:
            raise ValueError("neutrons.energy_max_MeV must be greater than neutrons.energy_min_MeV")
        angular_model = str(data.get("angular_model", "endf_b_viii1_evaluated_cm")).strip().lower()
        if angular_model not in {"isotropic_cm","endf_b_viii1_evaluated_cm",}:
            raise ValueError("neutrons.angular_model must be isotropic_cm or endf_b_viii1_evaluated_cm")

        return cls(
            model=canonical_model_name(data.get("model"), NEUTRON_MODEL_ALIASES, "distribution_kinematics_spectrum",),
            energy_min_MeV=energy_min,
            energy_max_MeV=energy_max,
            energy_bins=_integer(data.get("energy_bins"), "neutrons.energy_bins", minimum=1, default=80),
            angular_model=angular_model,
            num_gyroangle_points=_integer(data.get("num_gyroangle_points"), "neutrons.num_gyroangle_points", minimum=1, default=2),
            num_emission_directions=_integer(data.get("num_emission_directions"), "neutrons.num_emission_directions", minimum=1, default=4),
            allow_out_of_range=_bool(data.get("allow_out_of_range"), False),
            max_pair_states_per_batch=_integer(data.get("max_pair_states_per_batch"), "neutrons.max_pair_states_per_batch", minimum=1, default=131072),
            correlated_event_count=_integer(data.get("correlated_event_count"), "neutrons.correlated_event_count", minimum=1, default=10000),
            correlated_event_seed=_integer(data.get("correlated_event_seed"), "neutrons.correlated_event_seed", minimum=0, default=12345),
            radial_source_profile_model=canonical_radial_profile_model(data.get("radial_source_profile_model", KOTELNIKOV_PARABOLIC_FLUX_K1)),
        )

def _filename_only(value: Any, field_name: str, suffix: str) -> str:
    """Validate one output filename with no directory component and the required suffix"""
    text = str(value).strip()
    candidate = Path(text)
    if not text or candidate.name != text or candidate.parent != Path("."):
        raise ValueError(f"{field_name} must be one filename without a directory")
    if candidate.suffix.lower() != suffix:
        raise ValueError(f"{field_name} must use the {suffix} suffix")
  
    return text

@dataclass(frozen=True)
class OpenMCExportConfig:
    """
    Correlated OpenMC file source output controls
    
    `source_time_s` is written into exported source particles and the three configured output filenames must remain distinct
    """
    source_output_directory: str = "openmc_source"
    source_filename: str = "source.h5"
    validity_filename: str = "source_file_validity.json"
    sampling_metadata_filename: str = "source_sampling_metadata.json"
    source_time_s: float = 0.0
    overwrite_existing: bool = False

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, strict: bool = False) -> "OpenMCExportConfig":
        """Parse output directory, filenames, source time, and overwrite behavior"""
        allowed = {"source_output_directory", "source_filename", "validity_filename", "sampling_metadata_filename", "source_time_s", "overwrite_existing"}
        _check_unknown(data, allowed, "source_model.openmc_export", strict)
        output_directory = str(data.get("source_output_directory", "openmc_source")).strip()
        if not output_directory:
            raise ValueError("openmc_export.source_output_directory must not be empty")
        source_filename = _filename_only(data.get("source_filename", "source.h5"), "openmc_export.source_filename", ".h5")
        validity_filename = _filename_only(data.get("validity_filename", "source_file_validity.json"), "openmc_export.validity_filename", ".json")
        sampling_metadata_filename = _filename_only(data.get("sampling_metadata_filename", "source_sampling_metadata.json"), "openmc_export.sampling_metadata_filename", ".json")
        if len({source_filename, validity_filename, sampling_metadata_filename}) != 3:
            raise ValueError("openmc_export output filenames must be distinct")
        source_time_s = _number(data.get("source_time_s"), "openmc_export.source_time_s", nonnegative=True, default=0.0)
       
        return cls(
            source_output_directory=output_directory,
            source_filename=source_filename,
            validity_filename=validity_filename,
            sampling_metadata_filename=sampling_metadata_filename,
            source_time_s=source_time_s,
            overwrite_existing=_bool(data.get("overwrite_existing"), False),
        )
