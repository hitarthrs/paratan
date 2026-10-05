"""Neutral beam component, path, and source configuration"""
from __future__ import annotations
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any
import numpy as np
from source_model_revamp.integration.config.common import canonical_model_name, _mapping, _check_unknown, _number, _bool, _point3
from source_model_revamp.integration.config.constants import *
from source_model_revamp.beam.hydrogenic_dt_atomic_data import HYDROGENIC_DT_GROUND_STATE_MODEL, SUPPORTED_COMPONENT_ENERGY_MAX_KEV, SUPPORTED_COMPONENT_ENERGY_MIN_KEV

@dataclass(frozen=True)
class BeamComponentConfig:
    """
    One beam energy component
    
    `power_fraction` is the fraction of parent beam power in the component and `energy_multiplier` multiplies the parent `energy_keV`
    """
    id: str = "full"
    power_fraction: float = 1.0
    energy_multiplier: float = 1.0

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, beam_id: str, strict: bool = False) -> "BeamComponentConfig":
        """Parse one beam energy component and validate its identifier, power fraction, and energy multiplier"""
        allowed = {"id", "power_fraction", "energy_multiplier"}
        _check_unknown(data, allowed, f"source_model.beams[{beam_id}].components", strict)
        component_id = str(data.get("id", "")).strip()
        if not component_id:
            raise ValueError(f"source_model.beams[{beam_id}].components.id must be nonempty")

        return cls(
            id=component_id,
            power_fraction=_number(data.get("power_fraction"), f"beams[{beam_id}].components[{component_id}].power_fraction", nonnegative=True),
            energy_multiplier=_number(data.get("energy_multiplier"), f"beams[{beam_id}].components[{component_id}].energy_multiplier", positive=True),
        )

@dataclass(frozen=True)
class BeamlineConfig:
    """
    Straight neutral beam path and finite beam envelope
    
    Endpoints and radius are in m while `divergence_half_angle_deg` is in degrees
    """
    start_m: tuple[float, float, float]
    end_m: tuple[float, float, float]
    radius_m: float = 0.05
    divergence_half_angle_deg: float = 0.0
    radial_overlap_model: str = "hard_edge_overlap"

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, beam_id: str = "beam") -> "BeamlineConfig":
        """Parse one beamline and require distinct start and end coordinates"""
        start_m = _point3(data.get("start_m"), f"beams[{beam_id}].beamline.start_m")
        end_m = _point3(data.get("end_m"), f"beams[{beam_id}].beamline.end_m")
        if float(np.linalg.norm(np.asarray(end_m, dtype=float) - np.asarray(start_m, dtype=float))) <= 0.0:
            raise ValueError(f"source_model.beams[{beam_id}].beamline must have nonzero path length")

        return cls(
            start_m=start_m,
            end_m=end_m,
            radius_m=_number(data.get("radius_m"), f"beams[{beam_id}].beamline.radius_m", positive=True),
            divergence_half_angle_deg=_number(data.get("divergence_half_angle_deg"), f"beams[{beam_id}].beamline.divergence_half_angle_deg", nonnegative=True, default=0.0),
            radial_overlap_model=str(data.get("radial_overlap_model", "hard_edge_overlap")),
        )

@dataclass(frozen=True)
class BeamConfig:
    """
    One deuterium or tritium neutral beam definition
    
    The configuration combines injected power and energy, fractional energy components, pitch integration controls, atomic data selection, and the physical beamline
    """
    id: str = "d_nbi_1"
    enabled: bool = True
    model: str = "geometry_linked_attenuation"
    species: str = "deuterium"
    power_W: float = 1.0e6
    energy_keV: float = 100.0
    injection_angle_deg: float = 45.0
    injection_pitch_source: str = "beamline_endpoints"
    injection_angle_tolerance_deg: float = 1.0e-6
    components: tuple[BeamComponentConfig, ...] = field(default_factory=lambda: (BeamComponentConfig(),))
    atomic_data_model: str = HYDROGENIC_DT_GROUND_STATE_MODEL
    channel: str = "fast_ion_birth"
    pitch_distribution_model: str = "geometry_linked_finite_width"
    pitch_quadrature_points: int = 7
    beamline: BeamlineConfig = field(default_factory=lambda: BeamlineConfig(start_m=(-0.75, 0.0, -0.75), end_m=(0.75, 0.0, 0.75)))

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, strict: bool = False) -> "BeamConfig":
        """
        Parse one beam and validate component normalization, species, pitch controls, and supported component energies
        
        Component energies are `energy_keV * energy_multiplier` and component power fractions must sum to one
        """
        allowed = {"id", "enabled", "model", "species", "power_W", "energy_keV", "injection_angle_deg", "injection_pitch_source", "injection_angle_tolerance_deg", "components", "atomic_data_model", "channel", "pitch_distribution_model", "pitch_quadrature_points", "beamline"}
        _check_unknown(data, allowed, "source_model.beams", strict)
        beam_id = str(data.get("id", "")).strip()
        if not beam_id:
            raise ValueError("source_model.beams.id must be nonempty")
        beamline_data = _mapping(data.get("beamline"), f"source_model.beams[{beam_id}].beamline")
        component_data = data.get("components")
        if not isinstance(component_data, Sequence) or isinstance(component_data, (str, bytes)) or not component_data:
            raise ValueError(f"source_model.beams[{beam_id}].components must be a nonempty sequence")
        components = tuple(BeamComponentConfig.from_mapping(_mapping(item, f"source_model.beams[{beam_id}].components"), beam_id=beam_id, strict=strict) for item in component_data)
        component_ids = tuple(component.id for component in components)
        if len(set(component_ids)) != len(component_ids):
            raise ValueError(f"source_model.beams[{beam_id}].components IDs must be unique")
        if not np.isclose(sum(component.power_fraction for component in components), 1.0):
            raise ValueError(f"source_model.beams[{beam_id}].components power fractions must sum to 1")
        model = canonical_model_name(data.get("model"), BEAM_MODEL_ALIASES, "geometry_linked_attenuation")
        species_input = str(data.get("species", "deuterium")).strip().lower()
        species_by_input = {"d": "deuterium", "deuterium": "deuterium", "t": "tritium", "tritium": "tritium"}
        if species_input not in species_by_input:
            raise ValueError(f"source_model.beams[{beam_id}].species must be deuterium or tritium")
        species = species_by_input[species_input]
        enabled = _bool(data.get("enabled"), True)
        power_W = _number(data.get("power_W"), f"beams[{beam_id}].power_W", nonnegative=True, default=1.0e6)
        if enabled and power_W <= 0.0:
            raise ValueError(f"source_model.beams[{beam_id}].power_W must be positive when enabled")
        injection_angle_deg = _number(data.get("injection_angle_deg"), f"beams[{beam_id}].injection_angle_deg", nonnegative=True, default=45.0)
        if injection_angle_deg > 90.0:
            raise ValueError(f"source_model.beams[{beam_id}].injection_angle_deg must lie between 0 and 90 degrees")
        pitch_distribution_model = str(data.get("pitch_distribution_model", "geometry_linked_finite_width")).strip().lower().replace("-", "_")
        if pitch_distribution_model not in {"geometry_linked_finite_width","mono_pitch_benchmark",}:
            raise ValueError(f"source_model.beams[{beam_id}].pitch_distribution_model must be geometry_linked_finite_width or mono_pitch_benchmark")
        pitch_quadrature_points = int(_number(data.get("pitch_quadrature_points"), f"beams[{beam_id}].pitch_quadrature_points", positive=True, default=7))
        if pitch_quadrature_points < 1:
            raise ValueError(f"source_model.beams[{beam_id}].pitch_quadrature_points must be at least one")
        if pitch_distribution_model == "geometry_linked_finite_width" and pitch_quadrature_points < 2:
            raise ValueError("finite width beam pitch integration requires at least two quadrature points")
        channel = str(data.get("channel", "fast_ion_birth")).strip().lower()
        if not channel:
            raise ValueError(f"source_model.beams[{beam_id}].channel must be nonempty")
        energy_keV = _number(data.get("energy_keV"), f"beams[{beam_id}].energy_keV", positive=True, default=100.0)
        # Validate the physical energy reached by every configured fractional component
        component_energies_keV = energy_keV * np.asarray([component.energy_multiplier for component in components], dtype=float)
        if enabled and (np.any(component_energies_keV < SUPPORTED_COMPONENT_ENERGY_MIN_KEV) or np.any(component_energies_keV > SUPPORTED_COMPONENT_ENERGY_MAX_KEV)):
            raise ValueError(f"source_model.beams[{beam_id}] component energies must lie from {SUPPORTED_COMPONENT_ENERGY_MIN_KEV:g} to {SUPPORTED_COMPONENT_ENERGY_MAX_KEV:g} keV")
        atomic_data_model = str(data.get("atomic_data_model", HYDROGENIC_DT_GROUND_STATE_MODEL)).strip().lower().replace("-", "_")
        if atomic_data_model != HYDROGENIC_DT_GROUND_STATE_MODEL:
            raise ValueError(f"source_model.beams[{beam_id}].atomic_data_model must be {HYDROGENIC_DT_GROUND_STATE_MODEL}")

        return cls(
            id=beam_id,
            enabled=enabled,
            model=model,
            species=species,
            power_W=power_W,
            energy_keV=energy_keV,
            injection_angle_deg=injection_angle_deg,
            injection_pitch_source=str(data.get("injection_pitch_source", "beamline_endpoints")),
            injection_angle_tolerance_deg=_number(data.get("injection_angle_tolerance_deg"), f"beams[{beam_id}].injection_angle_tolerance_deg", nonnegative=True, default=1.0e-6,),
            components=components,
            atomic_data_model=atomic_data_model,
            channel=channel,
            pitch_distribution_model=pitch_distribution_model,
            pitch_quadrature_points=pitch_quadrature_points,
            beamline=BeamlineConfig.from_mapping(beamline_data, beam_id=beam_id),
        )
