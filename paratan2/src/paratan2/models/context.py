"""
BuildContext — shared state passed between component builders during assembly.

Each component reads what it needs from context and writes back what
downstream components will need. This replaces the tangled parameter-passing
in the original SimpleMachineBuilder.
"""

from __future__ import annotations
import openmc


class BuildContext:
    """
    Mutable bag of state that flows through the build sequence.

    Components update it as they are built so the next component can
    consume accurate geometry information without hard-coded assumptions.
    """

    def __init__(self, room_region: openmc.Region, material_ns):
        # The room region — gets carved up as components are added
        self.room_region: openmc.Region = room_region

        # Materials module / namespace
        self.material_ns = material_ns

        # Set by VacuumVessel after build
        self.vv_outermost_radius: float = 0.0
        self.vv_outermost_bottleneck_radius: float = 0.0
        self.vv_regions: list[openmc.Region] = []

        # Set by CentralCell after build
        self.cc_outermost_radius: float = 0.0
        self.cc_half_length: float = 0.0  # half of axial_length

        # Set by HFCoil after build — used by EndCell for positioning
        self.hf_coil_center_z0: float = 0.0
        self.hf_casing_total_thickness: float = 0.0
        self.hf_magnet_axial_half: float = 0.0
        self.hf_shield_axial_toward_midplane: float = 0.0

    def subtract_from_room(self, regions: list[openmc.Region]) -> None:
        for region in regions:
            self.room_region &= ~region
