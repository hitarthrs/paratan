from __future__ import annotations
import openmc

from paratan2.models.components.base import MachineComponent
from paratan2.models.registry import register_component
from paratan2.config.models import LFCoilConfig
from paratan2.geometry.core import hollow_cylinder_with_shell


@register_component("lf_coil")
class SimpleLFCoilBuilder(MachineComponent):
    """Builds the LF (low-field / solenoid) coils."""

    def __init__(self, config: LFCoilConfig, material_ns):
        self.config = config
        self.material_ns = material_ns
        self.lf_coil_cells: list[openmc.Cell] = []
        self.lf_coil_regions: list[openmc.Region] = []
        self.lf_coil_shield_cells: list[openmc.Cell] = []
        self.lf_coil_shield_regions: list[openmc.Region] = []
        self._tally_descriptors: list[dict] = []

    def build(self, context) -> list[openmc.Cell]:
        outer_radius = context.cc_outermost_radius
        shield_mat = getattr(self.material_ns, self.config.materials.shield)
        magnet_mat = getattr(self.material_ns, self.config.materials.magnet)

        for i, center_z in enumerate(self.config.positions):
            shell_region, inner_region = hollow_cylinder_with_shell(
                center_z,
                outer_radius,
                self.config.inner_dimensions.radial_thickness,
                self.config.inner_dimensions.axial_length,
                self.config.shell_thicknesses.front,
                self.config.shell_thicknesses.back,
                self.config.shell_thicknesses.axial,
            )

            shield_cell = openmc.Cell(4000 + i, region=shell_region, fill=shield_mat)
            inner_cell = openmc.Cell(4100 + i, region=inner_region, fill=magnet_mat)

            self.lf_coil_shield_regions.append(shell_region)
            self.lf_coil_shield_cells.append(shield_cell)
            self.lf_coil_regions.append(inner_region)
            self.lf_coil_cells.append(inner_cell)

            if self.config.tallies.cell_tallies or self.config.tallies.mesh_tallies:
                self._tally_descriptors.append({
                    "type": "lf_coil",
                    "location": f"coil_{i}",
                    "description": "magnet",
                    "cell": inner_cell,
                    "cell_tallies": [t.model_dump() for t in self.config.tallies.cell_tallies],
                    "mesh_tallies": [t.model_dump() for t in self.config.tallies.mesh_tallies],
                })

        all_cells = self.lf_coil_cells + self.lf_coil_shield_cells
        context.subtract_from_room([c.region for c in all_cells])
        return all_cells

    def get_tally_descriptors(self) -> list[dict]:
        return self._tally_descriptors
