from __future__ import annotations
import numpy as np
import openmc

from paratan2.models.components.base import MachineComponent
from paratan2.models.registry import register_component
from paratan2.config.models import CentralCellConfig


@register_component("central_cell")
class SimpleCentralCell(MachineComponent):
    """Builds the central cell (breeding blanket + structural layers)."""

    def __init__(self, config: CentralCellConfig, material_ns):
        self.config = config
        self.material_ns = material_ns
        self.central_cell_cells: list[openmc.Cell] = []
        self.central_cell_cylinders: list[openmc.ZCylinder] = []
        self.central_cell_cylinders_radii: np.ndarray = np.array([])
        self._tally_descriptors: list[dict] = []

    def build(self, context) -> list[openmc.Cell]:
        inner_radius = context.vv_outermost_radius
        self.central_cell_cylinders_radii = inner_radius + np.cumsum(
            [layer.thickness for layer in self.config.layers]
        )

        for radius in self.central_cell_cylinders_radii:
            self.central_cell_cylinders.append(openmc.ZCylinder(r=radius))

        half_length = self.config.axial_length / 2
        left_plane = openmc.ZPlane(-half_length)
        right_plane = openmc.ZPlane(half_length)

        vv_regions = context.vv_regions

        # First (innermost) layer
        first_region = -right_plane & +left_plane & -self.central_cell_cylinders[0]
        for vv_region in vv_regions:
            first_region = first_region & ~vv_region

        first_cell = openmc.Cell(
            2000,
            region=first_region,
            fill=self.material_ns.__dict__[self.config.layers[0].material],
        )
        self.central_cell_cells.append(first_cell)

        # Subsequent layers
        for i in range(1, len(self.central_cell_cylinders_radii)):
            region = (
                -right_plane
                & +left_plane
                & -self.central_cell_cylinders[i]
                & +self.central_cell_cylinders[i - 1]
            )
            for vv_region in vv_regions:
                region = region & ~vv_region

            cell = openmc.Cell(
                2000 + i,
                region=region,
                fill=self.material_ns.__dict__[self.config.layers[i].material],
            )
            self.central_cell_cells.append(cell)

        # Build tally descriptors
        tally_data = self.config.tallies
        if tally_data.breeder.cell_tallies or tally_data.breeder.mesh_tallies:
            if self.central_cell_cells:
                self._tally_descriptors.append({
                    "type": "central_cell",
                    "location": "breeder",
                    "description": "breeding_blanket",
                    "cell": self.central_cell_cells[0],
                    "cell_tallies": [t.model_dump() for t in tally_data.breeder.cell_tallies],
                    "mesh_tallies": [t.model_dump() for t in tally_data.breeder.mesh_tallies],
                })

        for layer_entry in tally_data.layer_tallies:
            pos = layer_entry.position
            if pos < len(self.central_cell_cells):
                layer_name = self.config.layers[pos].material
                self._tally_descriptors.append({
                    "type": "central_cell",
                    "location": f"layer_{pos}",
                    "description": layer_name,
                    "cell": self.central_cell_cells[pos],
                    "cell_tallies": [t.model_dump() for t in layer_entry.cell_tallies],
                    "mesh_tallies": [t.model_dump() for t in layer_entry.mesh_tallies],
                })

        # Write to context
        context.cc_outermost_radius = self.outermost_radius
        context.cc_half_length = half_length
        context.subtract_from_room([cell.region for cell in self.central_cell_cells])

        return self.central_cell_cells

    def get_tally_descriptors(self) -> list[dict]:
        return self._tally_descriptors

    @property
    def outermost_region(self) -> openmc.Region | None:
        if self.central_cell_cells:
            return self.central_cell_cells[-1].region
        return None

    @property
    def outermost_radius(self) -> float:
        if len(self.central_cell_cylinders_radii) > 0:
            return float(self.central_cell_cylinders_radii[-1])
        return 0.0
