"""Root source model configuration assembly and cross block validation"""
from __future__ import annotations
from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import asdict, dataclass, field, replace
from typing import Any
from source_model_revamp.integration.config.common import canonical_model_name, _mapping, _check_unknown
from source_model_revamp.integration.config.constants import *
from source_model_revamp.integration.config.geometry import GeometryConfig, PlasmaBoundaryConfig
from source_model_revamp.integration.config.plasma import PlasmaClosureConfig
from source_model_revamp.integration.config.beam import BeamConfig
from source_model_revamp.integration.config.kinetic import KineticElectrostaticConfig
from source_model_revamp.integration.config.operating_point import BeamDensityCouplingConfig
from source_model_revamp.integration.config.power import PowerBalanceConfig
from source_model_revamp.integration.config.downstream import FusionConfig, NeutronConfig, OpenMCExportConfig
from source_model_revamp.integration.config.output import OutputConfig

@dataclass(frozen=True)
class SourceModelRunConfig:
    """
    Resolved configuration consumed by the source model pipeline
    
    The object gathers geometry, plasma, beams, modal kinetics, operating point coupling, energy balance, fusion, neutron, export, and reporting controls
    """
    model: str = PRODUCTION_WORKFLOW_MODEL
    strict: bool = True
    geometry: GeometryConfig = field(default_factory=GeometryConfig)
    plasma_boundaries: PlasmaBoundaryConfig = field(default_factory=PlasmaBoundaryConfig)
    plasma_closure: PlasmaClosureConfig = field(default_factory=PlasmaClosureConfig)
    beams: tuple[BeamConfig, ...] = field(default_factory=lambda: (BeamConfig(),))
    kinetic_electrostatic: KineticElectrostaticConfig = field(default_factory=KineticElectrostaticConfig)
    beam_density_coupling: BeamDensityCouplingConfig = field(default_factory=BeamDensityCouplingConfig)
    power_balance: PowerBalanceConfig = field(default_factory=PowerBalanceConfig)
    fusion: FusionConfig = field(default_factory=FusionConfig)
    neutrons: NeutronConfig = field(default_factory=NeutronConfig)
    openmc_export: OpenMCExportConfig = field(default_factory=OpenMCExportConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

    @property
    def enabled_beams(self) -> tuple[BeamConfig, ...]:
        """Return enabled beams in their configured order"""
        return tuple(beam for beam in self.beams if beam.enabled)

    def resolved_mapping(self, source_model: Mapping[str, Any]) -> dict[str, Any]:
        """
        Return a copy of the source model mapping with canonical parsed controls written back
        
        Derived beam structures and numerical settings are serialized while unrelated user fields such as run labels are preserved
        """
        data = deepcopy(dict(source_model))
        data["model"] = self.model
        data["beams"] = [asdict(beam) for beam in self.beams]
        data["kinetic_electrostatic"] = asdict(self.kinetic_electrostatic)
        data["beam_density_coupling"] = asdict(self.beam_density_coupling)
        power_balance = asdict(self.power_balance)
        if power_balance.get("temperature_search_initial_trial_keV") is None:
            power_balance.pop("temperature_search_initial_trial_keV", None)
        data["power_balance"] = power_balance
        openmc_export = dict(data.get("openmc_export", {}))
        openmc_export["source_output_directory"] = self.openmc_export.source_output_directory
        data["openmc_export"] = openmc_export
        data["output"] = asdict(self.output)

        return data

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any], *, root: Mapping[str, Any] | None = None, strict: bool | None = None) -> "SourceModelRunConfig":
        """
        Parse the complete `source_model` block and apply cross block consistency rules
        
        The plasma closure controls the modal density closure model and startup seed closures default to the self consistent electron temperature mode when no mode is explicitly supplied
        """
        root = root or {}
        sm = _mapping(data, "source_model")
        local_strict = True if strict is None else bool(strict)
        if "beam" in sm:
            raise ValueError("source_model.beam is unsupported; use source_model.beams")
        allowed = {"model", "name", "run_name", "geometry", "plasma_boundaries", "plasma_closure", "beams", "kinetic_electrostatic", "beam_density_coupling", "power_balance", "fusion", "neutrons", "openmc_export", "output"}
        _check_unknown(sm, allowed, "source_model", local_strict)
        model = canonical_model_name(sm.get("model"), MODEL_ALIASES, PRODUCTION_WORKFLOW_MODEL)
        geometry_data = dict(_mapping(sm.get("geometry"), "source_model.geometry"))
        plasma_data = _mapping(sm.get("plasma_closure"), "source_model.plasma_closure")
        boundary_data = _mapping(sm.get("plasma_boundaries"), "source_model.plasma_boundaries")
        beams_data = sm.get("beams")
        if not isinstance(beams_data, Sequence) or isinstance(beams_data, (str, bytes)) or not beams_data:
            raise ValueError("source_model.beams must be a nonempty sequence")
        beams = tuple(BeamConfig.from_mapping(_mapping(item, "source_model.beams"), strict=local_strict) for item in beams_data)
        beam_ids = tuple(beam.id for beam in beams)
        if len(set(beam_ids)) != len(beam_ids):
            raise ValueError("source_model.beams IDs must be unique")
        enabled_beams = tuple(beam for beam in beams if beam.enabled)
        if not enabled_beams:
            raise ValueError("source_model.beams must contain at least one enabled deuterium or tritium beam")
        if len({beam.channel for beam in enabled_beams}) != 1:
            raise ValueError("enabled source_model.beams must use one common source channel")
        kinetic_data = _mapping(sm.get("kinetic_electrostatic"), "source_model.kinetic_electrostatic")
        beam_density_coupling_data = _mapping(sm.get("beam_density_coupling"), "source_model.beam_density_coupling")
        power_balance_data = _mapping(sm.get("power_balance"), "source_model.power_balance")
        fusion_data = _mapping(sm.get("fusion"), "source_model.fusion")
        neutron_data = dict(_mapping(sm.get("neutrons"), "source_model.neutrons"))
        export_data = _mapping(sm.get("openmc_export"), "source_model.openmc_export")
        output_data = _mapping(sm.get("output"), "source_model.output")
        plasma_closure = PlasmaClosureConfig.from_mapping(plasma_data, strict=local_strict)
        kinetic_electrostatic = KineticElectrostaticConfig.from_mapping(kinetic_data, strict=local_strict)
        # Keep one density closure authority across the plasma and modal configuration blocks
        if ("modal_density_closure_model" in kinetic_data and kinetic_electrostatic.modal_density_closure_model != plasma_closure.model):
            raise ValueError("kinetic_electrostatic.modal_density_closure_model disagrees with plasma_closure.model")
        kinetic_electrostatic = replace(kinetic_electrostatic, modal_density_closure_model=plasma_closure.model)
        beam_density_coupling = BeamDensityCouplingConfig.from_mapping(beam_density_coupling_data, strict=local_strict)
        power_balance = PowerBalanceConfig.from_mapping(power_balance_data, strict=local_strict)
        # Startup seed closures solve the electron temperature unless the user selects another mode explicitly
        if plasma_closure.uses_startup_seed_only and "electron_temperature_mode" not in power_balance_data:
            power_balance = replace(power_balance, electron_temperature_mode="self_consistent_electron_energy")

        return cls(
            model=model,
            strict=local_strict,
            geometry=GeometryConfig.from_mapping(geometry_data, root=root, strict=local_strict),
            plasma_boundaries=PlasmaBoundaryConfig.from_mapping(boundary_data, strict=local_strict),
            plasma_closure=plasma_closure,
            beams=beams,
            kinetic_electrostatic=kinetic_electrostatic,
            beam_density_coupling=beam_density_coupling,
            power_balance=power_balance,
            fusion=FusionConfig.from_mapping(fusion_data, strict=local_strict),
            neutrons=NeutronConfig.from_mapping(neutron_data, strict=local_strict),
            openmc_export=OpenMCExportConfig.from_mapping(export_data, strict=local_strict),
            output=OutputConfig.from_mapping(output_data, strict=local_strict),
        )

def source_model_config_from_root_mapping(root: Mapping[str, Any], *, strict: bool | None = None) -> SourceModelRunConfig:
    """Extract and validate the `source_model` block from the full input mapping"""
    if "source_model" not in root:
        raise ValueError("input mapping must contain a source_model block")
    return SourceModelRunConfig.from_mapping(_mapping(root["source_model"], "source_model"), root=root, strict=strict)