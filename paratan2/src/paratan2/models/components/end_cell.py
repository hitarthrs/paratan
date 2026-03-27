from __future__ import annotations
import openmc

from paratan2.models.components.base import MachineComponent
from paratan2.models.registry import register_component
from paratan2.config.models import EndCellConfig
from paratan2.geometry.core import cylinder_with_shell


@register_component("end_cell")
class SimpleEndCellBuilder(MachineComponent):
    """Builds the end cells (beyond the HF coils at each axial end)."""

    def __init__(self, config: EndCellConfig, material_ns):
        self.config = config
        self.material_ns = material_ns
        self.end_cell_cells: list[openmc.Cell] = []
        self.end_cell_regions: list[openmc.Region] = []
        self._tally_descriptors: list[dict] = []

    def build(self, context) -> list[openmc.Cell]:
        cfg = self.config
        shell_mat = getattr(self.material_ns, cfg.shell_material)
        inner_mat = getattr(self.material_ns, cfg.inner_material)

        radial_thickness = cfg.diameter / 2

        right_z0 = (
            context.hf_coil_center_z0
            + context.hf_shield_axial_toward_midplane
            + context.hf_casing_total_thickness
            + context.hf_magnet_axial_half
            + cfg.shell_thickness
            + cfg.axial_length / 2
        )
        left_z0 = -right_z0

        left_shell, left_inner = cylinder_with_shell(
            left_z0, radial_thickness, cfg.axial_length, cfg.shell_thickness
        )
        right_shell, right_inner = cylinder_with_shell(
            right_z0, radial_thickness, cfg.axial_length, cfg.shell_thickness
        )

        # Subtract all VV regions to avoid overlap
        for region in context.vv_regions:
            left_shell &= ~region
            right_shell &= ~region
            left_inner &= ~region
            right_inner &= ~region

        left_shell_cell = openmc.Cell(5001, region=left_shell, fill=shell_mat)
        left_inner_cell = openmc.Cell(5002, region=left_inner, fill=inner_mat)
        right_shell_cell = openmc.Cell(5003, region=right_shell, fill=shell_mat)
        right_inner_cell = openmc.Cell(5004, region=right_inner, fill=inner_mat)

        self.end_cell_cells.extend([
            left_shell_cell, left_inner_cell,
            right_shell_cell, right_inner_cell,
        ])
        self.end_cell_regions.extend([left_shell, left_inner, right_shell, right_inner])

        if cfg.tallies.cell_tallies or cfg.tallies.mesh_tallies:
            for loc, cell in [("left", left_shell_cell), ("right", right_shell_cell)]:
                self._tally_descriptors.append({
                    "type": "end_cell",
                    "location": loc,
                    "description": "shell",
                    "cell": cell,
                    "cell_tallies": [t.model_dump() for t in cfg.tallies.cell_tallies],
                    "mesh_tallies": [t.model_dump() for t in cfg.tallies.mesh_tallies],
                })

        context.subtract_from_room([c.region for c in self.end_cell_cells])
        return self.end_cell_cells

    def get_tally_descriptors(self) -> list[dict]:
        return self._tally_descriptors
