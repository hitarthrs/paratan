"""YAML loading and typed source model configuration helpers"""
from __future__ import annotations
from pathlib import Path
from collections.abc import Mapping
from typing import Any
import yaml
from source_model_revamp.integration.config.root import SourceModelRunConfig, source_model_config_from_root_mapping

def load_yaml_mapping(path: str | Path) -> dict[str, Any]:
    """Load one YAML file and require its top level value to be a mapping"""
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if payload is None:
        payload = {}
    if not isinstance(payload, Mapping):
        raise ValueError("YAML input must contain a top level mapping")
    
    return dict(payload)

def load_source_model_run_config(path: str | Path, *, strict: bool | None = None) -> SourceModelRunConfig:
    """Load one YAML file and parse its source_model block into SourceModelRunConfig"""
    root = load_yaml_mapping(path)

    return source_model_config_from_root_mapping(root, strict=strict)