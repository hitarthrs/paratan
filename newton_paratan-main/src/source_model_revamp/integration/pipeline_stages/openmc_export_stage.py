"""Correlated OpenMC HDF5 file source export stage"""
from __future__ import annotations
from time import perf_counter
from pathlib import Path
import numpy as np
from source_model_revamp.export.openmc_file_source import create_openmc_correlated_file_source
from source_model_revamp.integration.pipeline_types import NeutronStageResult, OpenMCExportStageResult
from source_model_revamp.integration.config.root import SourceModelRunConfig
from source_model_revamp.pipeline_errors import PipelineConservationError
from source_model_revamp.runtime_profile import finalize_runtime_profile, new_runtime_profile, record_runtime

def build_openmc_export_stage(config: SourceModelRunConfig, neutrons: NeutronStageResult, *, write_source: bool = True) -> OpenMCExportStageResult:
    """
    Validate and optionally write the single correlated OpenMC file source bundle
    
    Positive event probabilities must sum to one and the event bank physical rate must equal the neutron stage total before any files are written
    """
    openmc_runtime_started = perf_counter()
    runtime_profile = new_runtime_profile()
    event_bank = neutrons.correlated_event_bank
   
    if event_bank is None:
        raise ValueError("correlated OpenMC export requires a correlated neutron event bank")
  
    probabilities = np.asarray(event_bank.normalized_weights, dtype=float)
    positive = probabilities > 0.0
   
    # Sampling probabilities are normalized independently of the retained physical neutron rate
    probability_sum = float(np.sum(probabilities[positive]))
    physical_rate_s = float(event_bank.physical_total_rate_s)
   
    if not np.isclose(probability_sum, 1.0, rtol=0.0, atol=2.0e-12):
        raise PipelineConservationError("correlated source probabilities do not sum to one", stage="openmc_export", check_name="openmc_file_source_probability_normalization")
    if not np.isclose(physical_rate_s, float(neutrons.total_neutron_rate_s), rtol=1.0e-12, atol=0.0):
        raise PipelineConservationError("correlated source does not conserve physical rate", stage="openmc_export", check_name="openmc_export_source_conservation")
  
    file_source_bundle = None
  
    if write_source:
        arguments = {
            "output_directory": Path(config.openmc_export.source_output_directory),
            "source_filename": config.openmc_export.source_filename,
            "validity_filename": config.openmc_export.validity_filename,
            "sampling_metadata_filename": config.openmc_export.sampling_metadata_filename,
            "source_time_s": config.openmc_export.source_time_s,
            "overwrite_existing": config.openmc_export.overwrite_existing,
        }
        file_source_started = perf_counter()
        file_source_bundle = create_openmc_correlated_file_source(event_bank, **arguments)
      
        record_runtime(runtime_profile, "openmc_file_source_creation", perf_counter() - file_source_started)
   
    metadata = {
        "openmc_export_model": "correlated_openmc_hdf5_file_source",
        "openmc_source_representation": "openmc_hdf5_correlated_file_source" if write_source else "correlated_file_source_not_written",
        "num_openmc_sources": 1,
        "openmc_file_source_path": None if file_source_bundle is None else str(file_source_bundle.source_file_path),
        "openmc_file_source_validity_path": None if file_source_bundle is None else str(file_source_bundle.validity_file_path),
    }
    metadata = {**metadata, "runtime_openmc_export_profile": finalize_runtime_profile(runtime_profile, total_s=perf_counter() - openmc_runtime_started)}
   
    return OpenMCExportStageResult(file_source_bundle=file_source_bundle, metadata=metadata)

__all__ = ["build_openmc_export_stage"]
