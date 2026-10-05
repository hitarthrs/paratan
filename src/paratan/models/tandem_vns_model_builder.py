"""Tandem VNS geometry with the central blanket reserved for test modules.

All dimensions follow the tandem input convention (centimetres). This first
stage builds the existing tandem machine and a vacuum-filled central test
envelope. It does not construct individual test modules.
"""

from copy import deepcopy

import openmc

from src.paratan.materials import material as default_materials
from src.paratan.models.tandem_model_builder import (
    TandemMachineBuilder,
    parse_machine_input,
)


def parse_vns_machine_input(input_data, material_ns=default_materials):
    """Parse tandem inputs, replacing only the central blanket with an envelope.

    ``central_cell.test_region`` defines the axial length and radial thickness
    measured outward from the central first wall. ``test_modules.count`` records
    the intended number of modules; their internal geometry is a later step.
    """
    central = input_data["central_cell"]
    if "blanket" in central:
        raise ValueError("VNS input uses central_cell.test_region, not central_cell.blanket")

    test_region = central["test_region"]
    for key in ("axial_length", "radial_thickness"):
        value = test_region[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
            raise ValueError(f"central_cell.test_region.{key} must be positive")

    module_count = central["test_modules"]["count"]
    if isinstance(module_count, bool) or not isinstance(module_count, int) or module_count < 1:
        raise ValueError("central_cell.test_modules.count must be a positive integer")

    # The shared parser expects a central blanket layer. A single vacuum layer
    # preserves its radius and the LF-coil stand-off without introducing any
    # blanket material, shield, back wall, or breeder tally in the central cell.
    tandem_input = deepcopy(input_data)
    tandem_input["central_cell"]["blanket"] = {
        "axial_length": test_region["axial_length"],
        "layers": [{"thickness": test_region["radial_thickness"], "material": "vacuum"}],
    }
    return parse_machine_input(tandem_input, material_ns)


class TandemVNSMachineBuilder(TandemMachineBuilder):
    """Shared tandem machine with a replaceable central test-volume cell."""

    def __init__(self, *args, module_count, **kwargs):
        super().__init__(*args, **kwargs)
        self.module_count = module_count
        self._universe.name = "TandemVNSMachine"
        self.test_region_cell = None

    def subtract_from_room(self, regions):
        super().subtract_from_room(regions)
        # The parent stores the initial room region in the Cell. Its later
        # region reassignment must also reach that Cell to avoid overlaps.
        self._universe.cells[69].region = self._room_region

    def build_central_cylinders(self):
        super().build_central_cylinders()
        self.test_region_cell = self.ccyl_builder.get_cells_by_region()["central"][0]
        self.test_region_cell.name = "central_test_region_placeholder"
        self._regions["test_region"] = self.ccyl_builder.get_regions_by_region()["central"][0]

    def set_test_module_universe(self, universe):
        """Replace the placeholder fill with a module universe when available.

        The module universe must tile its local space, including any gaps
        between modules. No module geometry is created or checked here.
        """
        if self.test_region_cell is None:
            raise RuntimeError("Build central cylinders before installing modules")
        if not isinstance(universe, openmc.Universe):
            raise TypeError("Expected an OpenMC Universe")
        self.test_region_cell.fill = universe


def build_tandem_vns_universe(input_data, material_ns=default_materials):
    """Build the shared tandem geometry and return its builder.

    The returned builder exposes ``get_universe()``, ``get_all_tallies()``, and
    ``test_region_cell`` for later module construction and scoring.
    """
    params = parse_vns_machine_input(input_data, material_ns)
    builder = TandemVNSMachineBuilder(
        *params, material_ns=material_ns,
        module_count=input_data["central_cell"]["test_modules"]["count"],
    )
    builder.build_vacuum_vessel()
    builder.build_first_wall()
    builder.build_central_cylinders()
    builder.build_lf_coils()
    builder.build_hf_coils()
    builder.build_end_cells()
    return builder
