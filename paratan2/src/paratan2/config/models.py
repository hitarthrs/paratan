"""
Pydantic config models for all simple mirror machine components.

These replace raw YAML dicts and provide type checking + validation
before anything gets handed to OpenMC.
"""

from __future__ import annotations
from typing import Any
from pydantic import BaseModel, Field, model_validator


# ---------------------------------------------------------------------------
# Shared primitives
# ---------------------------------------------------------------------------

class MaterialLayer(BaseModel):
    thickness: float = Field(gt=0, description="Layer thickness (cm)")
    material: str = Field(description="Material name — must exist in the materials module")


class TallyEntry(BaseModel):
    scores: list[str]
    filters: list[str] = []
    nuclides: list[str] = []


class MeshTallyEntry(BaseModel):
    scores: list[str]
    filters: list[str] = []
    nuclides: list[str] = []
    dimensions: list[int] = Field(default=[10, 10, 10], min_length=3, max_length=3)


class CellTallyConfig(BaseModel):
    cell_tallies: list[TallyEntry] = []
    mesh_tallies: list[MeshTallyEntry] = []


# ---------------------------------------------------------------------------
# Vacuum vessel
# ---------------------------------------------------------------------------

class VacuumVesselStructureLayer(BaseModel):
    thickness: float = Field(gt=0)
    material: str


class VacuumVesselConfig(BaseModel):
    outer_axial_length: float = Field(gt=0, description="Axial length of outer (bottleneck) cylinders (cm)")
    central_axial_length: float = Field(gt=0, description="Axial length of central cylinder (cm)")
    central_radius: float = Field(gt=0, description="Radius of central cylinder (cm)")
    bottleneck_radius: float = Field(gt=0, description="Radius of bottleneck cylinders (cm)")
    left_bottleneck_length: float = Field(gt=0, description="Left bottleneck axial length (cm)")
    right_bottleneck_length: float = Field(gt=0, description="Right bottleneck axial length (cm)")
    axial_midplane: float = Field(default=0.0, description="Z-coordinate of midplane (cm)")
    # Ordered dict of layer_name -> properties
    structure: dict[str, VacuumVesselStructureLayer] = Field(
        default={},
        description="Ordered structural layers from inner to outer (e.g. first_wall, blanket, shield)"
    )

    @model_validator(mode="after")
    def check_radii(self) -> VacuumVesselConfig:
        if self.central_radius <= self.bottleneck_radius:
            raise ValueError(
                f"central_radius ({self.central_radius}) must be greater than "
                f"bottleneck_radius ({self.bottleneck_radius})"
            )
        if self.central_axial_length <= self.outer_axial_length:
            raise ValueError(
                f"central_axial_length ({self.central_axial_length}) must be greater than "
                f"outer_axial_length ({self.outer_axial_length})"
            )
        return self


# ---------------------------------------------------------------------------
# Central cell
# ---------------------------------------------------------------------------

class CentralCellBreederTallyConfig(BaseModel):
    cell_tallies: list[TallyEntry] = []
    mesh_tallies: list[MeshTallyEntry] = []


class CentralCellLayerTallyEntry(BaseModel):
    position: int = Field(ge=0)
    cell_tallies: list[TallyEntry] = []
    mesh_tallies: list[MeshTallyEntry] = []


class CentralCellTallyConfig(BaseModel):
    breeder: CentralCellBreederTallyConfig = Field(default_factory=CentralCellBreederTallyConfig)
    layer_tallies: list[CentralCellLayerTallyEntry] = []


class CentralCellConfig(BaseModel):
    axial_length: float = Field(gt=0, description="Total axial length of the central cell (cm)")
    layers: list[MaterialLayer] = Field(min_length=1, description="Radial layers from inner to outer")
    tallies: CentralCellTallyConfig = Field(default_factory=CentralCellTallyConfig)


# ---------------------------------------------------------------------------
# LF coil
# ---------------------------------------------------------------------------

class LFCoilShellThicknesses(BaseModel):
    front: float = Field(gt=0, description="Radial shell thickness on inner (plasma) side (cm)")
    back: float = Field(gt=0, description="Radial shell thickness on outer side (cm)")
    axial: float = Field(gt=0, description="Axial shell thickness (cm)")


class LFCoilInnerDimensions(BaseModel):
    axial_length: float = Field(gt=0, description="Axial length of coil winding pack (cm)")
    radial_thickness: float = Field(gt=0, description="Radial thickness of coil winding pack (cm)")


class LFCoilMaterials(BaseModel):
    shield: str = Field(description="Shield/casing material name")
    magnet: str = Field(description="Winding pack material name")


class LFCoilTallyConfig(BaseModel):
    cell_tallies: list[TallyEntry] = []
    mesh_tallies: list[MeshTallyEntry] = []


class LFCoilConfig(BaseModel):
    shell_thicknesses: LFCoilShellThicknesses
    inner_dimensions: LFCoilInnerDimensions
    positions: list[float] = Field(min_length=1, description="Z-positions of coil centres (cm)")
    materials: LFCoilMaterials
    tallies: LFCoilTallyConfig = Field(default_factory=LFCoilTallyConfig)


# ---------------------------------------------------------------------------
# HF coil
# ---------------------------------------------------------------------------

class HFCoilMagnet(BaseModel):
    bore_radius: float = Field(gt=0, description="Inner bore radius (cm)")
    radial_thickness: float = Field(gt=0, description="Radial thickness of winding pack (cm)")
    axial_thickness: float = Field(gt=0, description="Axial thickness of winding pack (cm)")
    material: str


class HFCoilShield(BaseModel):
    radial_gap_before_casing: float = Field(ge=0, description="Radial gap between VV and casing (cm)")
    shield_central_cell_gap: float = Field(ge=0, description="Axial gap between CC and shield (cm)")
    radial_thickness: list[float] = Field(min_length=2, max_length=2,
                                           description="[inner side, outer side] radial thickness (cm)")
    axial_thickness: list[float] = Field(min_length=2, max_length=2,
                                          description="[toward midplane, away] axial thickness (cm)")
    material: str


class HFCoilTallyConfig(BaseModel):
    cell_tallies: list[TallyEntry] = []
    mesh_tallies: list[MeshTallyEntry] = []


class HFCoilConfig(BaseModel):
    magnet: HFCoilMagnet
    casing_layers: list[MaterialLayer] = Field(min_length=1)
    shield: HFCoilShield
    tallies: HFCoilTallyConfig = Field(default_factory=HFCoilTallyConfig)


# ---------------------------------------------------------------------------
# End cell
# ---------------------------------------------------------------------------

class EndCellTallyConfig(BaseModel):
    cell_tallies: list[TallyEntry] = []
    mesh_tallies: list[MeshTallyEntry] = []


class EndCellConfig(BaseModel):
    axial_length: float = Field(gt=0, description="Total axial length of each end cell (cm)")
    shell_thickness: float = Field(gt=0, description="Shell thickness (cm)")
    diameter: float = Field(gt=0, description="Outer diameter of end cell (cm)")
    shell_material: str = Field(default="stainless")
    inner_material: str = Field(default="vacuum")
    tallies: EndCellTallyConfig = Field(default_factory=EndCellTallyConfig)


# ---------------------------------------------------------------------------
# Room (bounding box)
# ---------------------------------------------------------------------------

class RoomConfig(BaseModel):
    xmin: float = -400.0
    xmax: float = 400.0
    ymin: float = -400.0
    ymax: float = 400.0
    zmin: float = -1350.0
    zmax: float = 1350.0


# ---------------------------------------------------------------------------
# Top-level simple machine config
# ---------------------------------------------------------------------------

class SimpleMachineConfig(BaseModel):
    vacuum_vessel: VacuumVesselConfig
    central_cell: CentralCellConfig
    lf_coil: LFCoilConfig
    hf_coil: HFCoilConfig
    end_cell: EndCellConfig = Field(default_factory=lambda: EndCellConfig(
        axial_length=50, shell_thickness=2, diameter=100
    ))
    room: RoomConfig = Field(default_factory=RoomConfig)

    @classmethod
    def from_yaml(cls, path: str) -> SimpleMachineConfig:
        import yaml
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls.model_validate(data)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> SimpleMachineConfig:
        return cls.model_validate(data)
