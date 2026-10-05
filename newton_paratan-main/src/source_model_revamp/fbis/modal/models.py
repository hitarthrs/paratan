"""Canonical identifiers and validators for selectable modal physics models"""
from __future__ import annotations

COLD_ION_EQ14 = "cold_ion_eq14"
HOT_ION_ROSENBLUTH_EQ59 = "hot_ion_rosenbluth_eq59"
MAGNETIC_ONLY_FBIS = "magnetic_only_fbis"
ZERO_POTENTIAL_MAGNETIC_REFERENCE = "zero_potential_magnetic_reference"
EGEDAL_2022_PUBLISHED_ELECTROSTATIC_FBIS_APPROXIMATION = "egedal_2022_published_electrostatic_fbis_approximation"
VELOCITY_SOLUTION_MODELS = (COLD_ION_EQ14, HOT_ION_ROSENBLUTH_EQ59)
ELECTROSTATIC_FEEDBACK_MODELS = (MAGNETIC_ONLY_FBIS, ZERO_POTENTIAL_MAGNETIC_REFERENCE, EGEDAL_2022_PUBLISHED_ELECTROSTATIC_FBIS_APPROXIMATION)

def velocity_solution_model(value: str) -> str:
    """Return a canonical supported modal velocity solution model identifier"""
    model = str(value).strip().lower()
    if model not in VELOCITY_SOLUTION_MODELS:
        raise ValueError(f"Unsupported modal velocity solution model {value!r}")
  
    return model

def electrostatic_feedback_model(value: str) -> str:
    """Return a canonical supported electrostatic feedback model identifier"""
    model = str(value).strip().lower()
    if model not in ELECTROSTATIC_FEEDBACK_MODELS:
        raise ValueError(f"Unsupported electrostatic feedback model {value!r}")
   
    return model

__all__ = [
    "COLD_ION_EQ14",
    "HOT_ION_ROSENBLUTH_EQ59",
    "MAGNETIC_ONLY_FBIS",
    "ZERO_POTENTIAL_MAGNETIC_REFERENCE",
    "EGEDAL_2022_PUBLISHED_ELECTROSTATIC_FBIS_APPROXIMATION",
    "VELOCITY_SOLUTION_MODELS",
    "ELECTROSTATIC_FEEDBACK_MODELS",
    "velocity_solution_model",
    "electrostatic_feedback_model",
]
