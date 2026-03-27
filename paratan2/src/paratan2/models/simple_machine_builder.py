"""
SimpleMachineBuilder — assembles a complete simple mirror machine model.

The build sequence is:
  1. Room (bounding box)
  2. VacuumVessel
  3. CentralCell
  4. LFCoils
  5. HFCoils
  6. EndCells

Each component is retrieved from the registry, built in order, and updates
the shared BuildContext. To add a new component type, register it with
@register_component and add it to the config — no changes here needed.
"""

from __future__ import annotations
import contextlib
import os

import openmc
import matplotlib.pyplot as plt

from paratan2.config.models import SimpleMachineConfig
from paratan2.models.context import BuildContext
from paratan2.models.components.vacuum_vessel import SimpleVacuumVessel
from paratan2.models.components.central_cell import SimpleCentralCell
from paratan2.models.components.lf_coil import SimpleLFCoilBuilder
from paratan2.models.components.hf_coil import SimpleHFCoilBuilder
from paratan2.models.components.end_cell import SimpleEndCellBuilder
from paratan2.tallies.builder import TallyBuilder
from paratan2.tallies.base import hollow_mesh_from_domain


def _load_materials(materials_path: str | None = None):
    """Load materials module. Override via MATERIALS_PATH env var or argument."""
    path = materials_path or os.environ.get("MATERIALS_PATH")
    if path:
        import importlib.util, sys
        spec = importlib.util.spec_from_file_location("custom_materials", path)
        m = importlib.util.module_from_spec(spec)
        sys.modules["custom_materials"] = m
        spec.loader.exec_module(m)
        return m
    from paratan2.materials import material as m
    return m


@contextlib.contextmanager
def _cd(path: str):
    origin = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(origin)


class SimpleMachineBuilder:
    """
    Orchestrates the full build sequence for a simple mirror machine.

    Parameters
    ----------
    config : SimpleMachineConfig
        Validated Pydantic config object.
    material_ns : module, optional
        Materials module. If None, loads the default paratan2 materials.
    """

    def __init__(self, config: SimpleMachineConfig, material_ns=None):
        self.config = config
        self.material_ns = material_ns or _load_materials()

        self.universe = openmc.Universe(786)
        self.all_cells: list[openmc.Cell] = []
        self.context: BuildContext | None = None

        # Component builders — instantiated here so they can be accessed
        # externally (e.g. for post-build queries like get_left_magnet_cell)
        self.vacuum_vessel = SimpleVacuumVessel(config.vacuum_vessel, self.material_ns)
        self.central_cell = SimpleCentralCell(config.central_cell, self.material_ns)
        self.lf_coils = SimpleLFCoilBuilder(config.lf_coil, self.material_ns)
        self.hf_coils = SimpleHFCoilBuilder(config.hf_coil, self.material_ns)
        self.end_cells = SimpleEndCellBuilder(config.end_cell, self.material_ns)

        self.tally_builder = TallyBuilder()

    def _add_cells(self, cells: list[openmc.Cell]) -> None:
        for cell in cells:
            self.universe.add_cells([cell])
            self.all_cells.append(cell)

    def build(self) -> None:
        """Run the full build sequence."""
        rc = self.config.room

        # --- Room ---
        bounding_surface = openmc.model.RectangularParallelepiped(
            xmin=rc.xmin, xmax=rc.xmax,
            ymin=rc.ymin, ymax=rc.ymax,
            zmin=rc.zmin, zmax=rc.zmax,
            boundary_type="vacuum",
        )
        room_region = -bounding_surface
        room_cell = openmc.Cell(100, region=room_region, fill=self.material_ns.air)
        self._add_cells([room_cell])

        self.context = BuildContext(room_region, self.material_ns)

        # Store VV half-lengths so HF coil can compute its offset
        vv_cfg = self.config.vacuum_vessel
        self.context.vv_half_length = vv_cfg.central_axial_length / 2

        # --- Build in sequence ---
        self._add_cells(self.vacuum_vessel.build(self.context))
        self._add_cells(self.central_cell.build(self.context))
        self._add_cells(self.lf_coils.build(self.context))
        self._add_cells(self.hf_coils.build(self.context))
        self._add_cells(self.end_cells.build(self.context))

        # Spacing sanity check
        vv_total = (
            vv_cfg.left_bottleneck_length
            + vv_cfg.central_axial_length
            + vv_cfg.right_bottleneck_length
        )
        cc_length = self.config.central_cell.axial_length
        available = (vv_total - cc_length) / 2
        if available < 250:
            print(
                f"WARNING: available axial space per side ({available:.1f} cm) "
                f"is less than 250 cm. Check coil positioning."
            )
        else:
            print(f"Spacing OK: {available:.1f} cm available per side")

    def get_geometry(self) -> openmc.Geometry:
        geometry = openmc.Geometry([openmc.Cell(fill=self.universe)])
        geometry.merge_surfaces = True
        return geometry

    def get_tallies(self) -> list[openmc.Tally]:
        all_descriptors = (
            self.central_cell.get_tally_descriptors()
            + self.lf_coils.get_tally_descriptors()
            + self.hf_coils.get_tally_descriptors()
            + self.end_cells.get_tally_descriptors()
        )
        self.tally_builder.add_descriptors(all_descriptors)
        return self.tally_builder.get_tallies()

    def get_weight_window_generators(
        self, hf_mesh_dims: list[int] | None = None
    ) -> list[openmc.WeightWindowGenerator]:
        dims = hf_mesh_dims or [20, 1, 20]
        left_cell = self.hf_coils.get_left_magnet_cell()
        right_cell = self.hf_coils.get_right_magnet_cell()
        wwgs = []
        for cell in [left_cell, right_cell]:
            mesh = hollow_mesh_from_domain(cell, dims)
            wwgs.append(
                openmc.WeightWindowGenerator(
                    mesh,
                    energy_bounds=[0, 14e6],
                    particle_type="neutron",
                    method="magic",
                    max_realizations=25,
                    update_interval=3,
                    on_the_fly=True,
                )
            )
        return wwgs


def build_simple_model(
    config: SimpleMachineConfig,
    output_dir: str = "output",
    source_yaml: str | None = None,
    materials_path: str | None = None,
) -> openmc.Model:
    """
    Build a complete OpenMC model for a simple mirror machine.

    Parameters
    ----------
    config : SimpleMachineConfig
        Validated machine configuration.
    output_dir : str
        Directory to write XML files and plots into.
    source_yaml : str, optional
        Path to source_information.yaml. If None the model is built without
        a source (useful for geometry-only checks).
    materials_path : str, optional
        Override path to a custom materials module.

    Returns
    -------
    openmc.Model
    """
    import yaml as _yaml

    os.makedirs(output_dir, exist_ok=True)

    m = _load_materials(materials_path)
    builder = SimpleMachineBuilder(config, m)

    with _cd(output_dir):
        builder.build()

        geometry = builder.get_geometry()

        geometry.root_universe.plot(
            basis="xz",
            width=(1000, 2800),
            pixels=(2400, 4000),
            color_by="material",
        )
        plt.savefig("simple_mirror_cross_section.png", bbox_inches="tight")
        plt.close()

        geometry.export_to_xml("geometry.xml")

        tallies = openmc.Tallies(builder.get_tallies())
        tallies.export_to_xml("tallies.xml")

        materials = m.materials

        settings = openmc.Settings()
        settings.run_mode = "fixed source"

        if source_yaml is not None:
            with open(source_yaml) as f:
                src_data = _yaml.safe_load(f)

            from paratan2.source.volumetric import load_source_from_yaml
            source = load_source_from_yaml(source_yaml)
            settings.source = source

            sim = src_data.get("settings", {})
            settings.particles = int(sim.get("particles_per_batch", 1000))
            settings.batches = int(sim.get("batches", 10))
            settings.output = {"tallies": False}

            freq = sim.get("statepoint_frequency", settings.batches)
            settings.statepoint = {
                "batches": [1]
                + list(range(freq, settings.batches, freq))
                + [settings.batches]
            }

            if sim.get("weight_windows", False):
                wwgs = builder.get_weight_window_generators()
                settings.weight_window_generators = wwgs
                settings.weight_windows_on = True
                settings.weight_window_checkpoints = {"collision": True, "surface": True}

            settings.photon_transport = sim.get("photon_transport", False)

        model = openmc.Model(geometry, materials, settings, tallies)
        model.export_to_xml("model_xml_files")

    return model
