from __future__ import annotations
import numpy as np
import openmc

from paratan2.models.components.base import MachineComponent
from paratan2.models.registry import register_component
from paratan2.config.models import HFCoilConfig
from paratan2.geometry.core import nested_cylindrical_shells


@register_component("hf_coil")
class SimpleHFCoilBuilder(MachineComponent):
    """Builds the HF (high-field / mirror) coils at the axial ends."""

    def __init__(self, config: HFCoilConfig, material_ns):
        self.config = config
        self.material_ns = material_ns
        self.hf_coil_cells: list[openmc.Cell] = []
        self.hf_coil_regions: list[openmc.Region] = []
        self._tally_descriptors: list[dict] = []
        self._left_magnet_cell: openmc.Cell | None = None
        self._right_magnet_cell: openmc.Cell | None = None

    def build(self, context) -> list[openmc.Cell]:
        cfg = self.config
        casing_thicknesses = np.array([layer.thickness for layer in cfg.casing_layers])

        # Offset: if VV central length > CC axial length, shift HF coil outward
        vv_half = context.vv_half_length if hasattr(context, "vv_half_length") else context.cc_half_length
        hf_z0_offset = max(0.0, vv_half - context.cc_half_length)

        hf_coil_center_z0 = (
            context.cc_half_length
            + cfg.shield.shield_central_cell_gap
            + cfg.shield.axial_thickness[0]
            + np.sum(casing_thicknesses)
            + cfg.magnet.axial_thickness / 2
            + hf_z0_offset
        )

        vv_outermost_bottleneck = context.vv_outermost_bottleneck_radius
        inner_start_radius = vv_outermost_bottleneck + cfg.shield.radial_gap_before_casing

        all_layers_front = np.concatenate([casing_thicknesses, [cfg.shield.radial_thickness[0]]])
        all_layers_back = np.concatenate([casing_thicknesses, [cfg.shield.radial_thickness[1]]])
        all_layers_axial = np.concatenate([casing_thicknesses, [cfg.shield.axial_thickness[0]]])

        shield_regions_right = nested_cylindrical_shells(
            hf_coil_center_z0,
            inner_start_radius,
            cfg.magnet.radial_thickness,
            cfg.magnet.axial_thickness,
            all_layers_front,
            all_layers_back,
            all_layers_axial,
        )
        shield_regions_left = nested_cylindrical_shells(
            -hf_coil_center_z0,
            inner_start_radius,
            cfg.magnet.radial_thickness,
            cfg.magnet.axial_thickness,
            all_layers_front,
            all_layers_back,
            all_layers_axial,
        )

        magnet_mat = getattr(self.material_ns, cfg.magnet.material)

        self._left_magnet_cell = openmc.Cell(
            6101, region=shield_regions_left[0], fill=magnet_mat
        )
        self._right_magnet_cell = openmc.Cell(
            6102, region=shield_regions_right[0], fill=magnet_mat
        )

        self.hf_coil_cells.extend([self._left_magnet_cell, self._right_magnet_cell])
        self.hf_coil_regions.extend([shield_regions_left[0], shield_regions_right[0]])

        casing_mats = [
            getattr(self.material_ns, layer.material) for layer in cfg.casing_layers
        ]
        shield_mat = getattr(self.material_ns, cfg.shield.material)

        for i in range(1, len(shield_regions_left)):
            combined = shield_regions_left[i] | shield_regions_right[i]
            if i <= len(casing_mats):
                material = casing_mats[i - 1]
                cell_id = 6500 + i
            else:
                material = shield_mat
                cell_id = 6200 + (i - len(casing_mats))
            cell = openmc.Cell(cell_id, region=combined, fill=material)
            self.hf_coil_cells.append(cell)
            self.hf_coil_regions.append(combined)

        if cfg.tallies.cell_tallies or cfg.tallies.mesh_tallies:
            for loc, cell in [("right", self._right_magnet_cell), ("left", self._left_magnet_cell)]:
                self._tally_descriptors.append({
                    "type": "hf_coil",
                    "location": loc,
                    "description": "magnet",
                    "cell": cell,
                    "cell_tallies": [t.model_dump() for t in cfg.tallies.cell_tallies],
                    "mesh_tallies": [t.model_dump() for t in cfg.tallies.mesh_tallies],
                })

        # Write to context for EndCell
        context.hf_coil_center_z0 = hf_coil_center_z0
        context.hf_casing_total_thickness = float(np.sum(casing_thicknesses))
        context.hf_magnet_axial_half = cfg.magnet.axial_thickness / 2
        context.hf_shield_axial_toward_midplane = cfg.shield.axial_thickness[0]

        context.subtract_from_room([c.region for c in self.hf_coil_cells])
        return self.hf_coil_cells

    def get_tally_descriptors(self) -> list[dict]:
        return self._tally_descriptors

    def get_left_magnet_cell(self) -> openmc.Cell:
        return self._left_magnet_cell

    def get_right_magnet_cell(self) -> openmc.Cell:
        return self._right_magnet_cell
