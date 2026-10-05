"""Integration layer for the source model"""
from source_model_revamp.integration.config.root import SourceModelRunConfig, source_model_config_from_root_mapping
from source_model_revamp.integration.yaml_loader import load_source_model_run_config, load_yaml_mapping
from source_model_revamp.integration.assessment import AssessmentRecord, RunAssessment
from source_model_revamp.integration.validation import SourceConstructionError
from source_model_revamp.pipeline_errors import InvalidPipelineStateError, PipelineConservationError
from source_model_revamp.integration.metadata import write_metadata_json
from source_model_revamp.integration.progress import OperatingPointConsoleProgressReporter

def __getattr__(name: str):
    """Load the orchestration entry point without creating import cycles"""
    if name == "run_source_model_pipeline":
        from source_model_revamp.integration.pipeline import run_source_model_pipeline
        return run_source_model_pipeline
    raise AttributeError(name)

__all__ = [
    "SourceModelRunConfig",
    "source_model_config_from_root_mapping",
    "load_source_model_run_config",
    "load_yaml_mapping",
    "run_source_model_pipeline",
    "write_metadata_json",
    "AssessmentRecord",
    "RunAssessment",
    "SourceConstructionError",
    "InvalidPipelineStateError",
    "PipelineConservationError",
    "OperatingPointConsoleProgressReporter",
]
