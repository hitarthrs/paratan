"""
Top-level runner API for Paratan 2.

Typical usage
-------------
    from paratan2.runner import build, run
    from paratan2.config.models import SimpleMachineConfig

    config = SimpleMachineConfig.from_yaml("my_machine.yaml")
    model = build(config, output_dir="output")      # geometry + XML only
    run(config, output_dir="output", source_yaml="source_information.yaml")
"""

from __future__ import annotations
import openmc

from paratan2.config.models import SimpleMachineConfig
from paratan2.models.simple_machine_builder import build_simple_model


def build(
    config: SimpleMachineConfig | str,
    output_dir: str = "output",
    materials_path: str | None = None,
) -> openmc.Model:
    """
    Build the model and export XML files without running the simulation.

    Parameters
    ----------
    config : SimpleMachineConfig or path to YAML
    output_dir : str
    materials_path : str, optional

    Returns
    -------
    openmc.Model
    """
    if isinstance(config, str):
        config = SimpleMachineConfig.from_yaml(config)
    return build_simple_model(
        config,
        output_dir=output_dir,
        source_yaml=None,
        materials_path=materials_path,
    )


def run(
    config: SimpleMachineConfig | str,
    source_yaml: str,
    output_dir: str = "output",
    materials_path: str | None = None,
) -> openmc.Model:
    """
    Build the model and run the simulation.

    Parameters
    ----------
    config : SimpleMachineConfig or path to YAML
    source_yaml : str — path to source_information.yaml
    output_dir : str
    materials_path : str, optional

    Returns
    -------
    openmc.Model (already run)
    """
    if isinstance(config, str):
        config = SimpleMachineConfig.from_yaml(config)
    model = build_simple_model(
        config,
        output_dir=output_dir,
        source_yaml=source_yaml,
        materials_path=materials_path,
    )
    model.run(cwd=output_dir)
    return model
