from __future__ import annotations
import numpy as np
import openmc

from paratan2.models.components.base import MachineComponent
from paratan2.models.registry import register_component
from paratan2.config.models import VacuumVesselConfig
from paratan2.geometry.core import single_vacuum_vessel_region


@register_component("vacuum_vessel")
class SimpleVacuumVessel(MachineComponent):
    """Builds the vacuum vessel structure for a simple mirror machine."""

    def __init__(self, config: VacuumVesselConfig, material_ns):
        self.config = config
        self.material_ns = material_ns
        self.vacuum_section_regions: list[openmc.Region] = []
        self.vacuum_section_cells: list[openmc.Cell] = []

    def build(self, context) -> list[openmc.Cell]:
        vv_main_region, _ = single_vacuum_vessel_region(
            self.config.outer_axial_length,
            self.config.central_axial_length,
            self.config.central_radius,
            self.config.bottleneck_radius,
            self.config.left_bottleneck_length,
            self.config.right_bottleneck_length,
            self.config.axial_midplane,
        )

        vv_main_cell = openmc.Cell(
            1000,
            region=vv_main_region,
            fill=self.material_ns.vacuum,
        )

        self.vacuum_section_regions.append(vv_main_region)
        self.vacuum_section_cells.append(vv_main_cell)

        if self.config.structure:
            structural_thicknesses = [
                layer.thickness for layer in self.config.structure.values()
            ]
            structural_materials = [
                getattr(self.material_ns, layer.material)
                for layer in self.config.structure.values()
            ]

            central_radii = np.cumsum(
                [self.config.central_radius] + structural_thicknesses
            )
            bottleneck_radii = np.cumsum(
                [self.config.bottleneck_radius] + structural_thicknesses
            )

            for i, (thickness, material) in enumerate(
                zip(structural_thicknesses, structural_materials)
            ):
                vv_layer_region, _ = single_vacuum_vessel_region(
                    self.config.outer_axial_length,
                    self.config.central_axial_length,
                    central_radii[i + 1],
                    bottleneck_radii[i + 1],
                    self.config.left_bottleneck_length,
                    self.config.right_bottleneck_length,
                    self.config.axial_midplane,
                )

                if i == 0:
                    structural_region = vv_layer_region & ~self.vacuum_section_regions[0]
                else:
                    structural_region = vv_layer_region & ~self.vacuum_section_regions[-1]

                structural_cell = openmc.Cell(
                    1000 + i + 1,
                    region=structural_region,
                    fill=material,
                )

                self.vacuum_section_regions.append(vv_layer_region)
                self.vacuum_section_cells.append(structural_cell)

        # Write to context
        context.vv_outermost_radius = self.outermost_radius
        context.vv_outermost_bottleneck_radius = self.outermost_bottleneck_radius
        context.vv_regions = self.vacuum_section_regions
        context.subtract_from_room(self.vacuum_section_regions)

        return self.vacuum_section_cells

    def get_tally_descriptors(self) -> list[dict]:
        return []

    @property
    def outermost_region(self) -> openmc.Region:
        return self.vacuum_section_regions[-1] if self.vacuum_section_regions else None

    @property
    def outermost_radius(self) -> float:
        if not self.config.structure:
            return self.config.central_radius
        total = sum(layer.thickness for layer in self.config.structure.values())
        return self.config.central_radius + total

    @property
    def outermost_bottleneck_radius(self) -> float:
        if not self.config.structure:
            return self.config.bottleneck_radius
        total = sum(layer.thickness for layer in self.config.structure.values())
        return self.config.bottleneck_radius + total
