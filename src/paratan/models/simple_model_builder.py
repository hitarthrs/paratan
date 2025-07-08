import openmc
import numpy as np
import openmc.lib
import matplotlib.pyplot as plt
from src.paratan.geometry.core import *
from src.paratan.tallies.tandem_tallies import TallyBuilder
from pathlib import Path
import src.paratan.materials.material as m
import yaml

from types import SimpleNamespace

def parse_simple_machine_input(input_data):
    """
    Parse structured input for simple mirror model.

    Parameters
    ----------
    input_data : dict
        User-provided YAML input.

    Returns
    -------
    Tuple of structured parameter dictionaries.
    """
    # --------- 1. Vacuum Vessel Parameters ---------
    vv_info = input_data["vacuum_vessel"]

    vv_params = {
        "machine_length_from_midplane": vv_info["machine_length_from_midplane"],
        "first_vv_plane_from_midplane": vv_info["first_vv_plane_from_midplane"],
        "vacuum_chamber": vv_info["vacuum_chamber"],
        "bottleneck_cylinder": vv_info["bottleneck_cylinder"],
        "cone_angle": vv_info["cone_angle"],
        "structure": vv_info["structure"]
    }

    # --------- 2. Central Cell Parameters ---------
    cc_info = input_data["central_cell"]

    central_cell_params = {
        "axial_length": cc_info["axial_length"],
        "layers": cc_info["layers"]
    }

    # --------- 3. LF Coil Parameters ---------
    lf_info = input_data["lf_coil"]
    lf_tallies = input_data.get("lf_coil_tallies", {})

    lf_coil_params = {
        "shell_thicknesses": lf_info["shell_thicknesses"],
        "inner_dimensions": lf_info["inner_dimensions"],
        "positions": lf_info["positions"],
        "materials": lf_info["materials"],
        "tallies": lf_tallies
    }

    # --------- 4. HF Coil Parameters ---------
    hf_info = input_data["hf_coil"]

    hf_coil_params = {
        "magnet": hf_info["magnet"],
        "casing_layers": hf_info["casing_layers"],
        "shield": hf_info["shield"],
        "tallies": hf_info.get("tallies", {})
    }

    # --------- 5. End Cell Parameters ---------
    end_cell_params = input_data.get("end_cell", {})

    return vv_params, central_cell_params, lf_coil_params, hf_coil_params, end_cell_params

class SimpleVacuumVessel:
    """
    Class to construct OpenMC cells and expose regions for a simple mirror vacuum vessel:
    - Single axisymmetric innermost vacuum region
    """

    def __init__(self, machine_length_from_midplane, first_vv_plane_from_midplane, cone_angle, vacuum_chamber, bottleneck_cylinder, vacuum_material):
        self.machine_length_from_midplane = machine_length_from_midplane
        self.first_vv_plane_from_midplane = first_vv_plane_from_midplane
        self.cone_angle = cone_angle
        self.vacuum_chamber_radius = vacuum_chamber["radius"]
        self.bottleneck_cylinder_radius = bottleneck_cylinder["cylinder_radius"]
        self.bottleneck_plane_distance = bottleneck_cylinder["plane_distance"]
        self.vacuum_material = vacuum_material

        self._region = None
        self._cell = None

    def build_cell(self):

        self._region = vacuum_vessel_region(
            self.first_vv_plane_from_midplane,
            self.machine_length_from_midplane,
            self.vacuum_chamber_radius,
            self.bottleneck_cylinder_radius,
            self.cone_angle
        )

        self._cell = openmc.Cell(1000, region=self._region, fill=self.vacuum_material)
        return self._cell

    def get_region(self):
        return self._region

    def get_cell(self):
        return self._cell
    

class SimpleFWStructure:
    """
    Class to construct concentric cylindrical structural layers around the simple vacuum vessel.
    """

    def __init__(self, structure_layers, vacuum_chamber_radius, bottleneck_cylinder_radius, vacuum_chamber_material, first_plane, machine_half_length, cone_angle):
        self.structure_layers = structure_layers
        self.vacuum_chamber_radius = vacuum_chamber_radius
        self.bottleneck_cylinder_radius = bottleneck_cylinder_radius
        self.vacuum_chamber_material = vacuum_chamber_material
        self.first_plane = first_plane
        self.machine_half_length = machine_half_length
        self.cone_angle = cone_angle

        self._cells = []
        self._regions = []
        self._full_region_list = []

    def build_cells(self, mat_ns, room_region, add_cell_callback):

        structural_thicknesses = [layer.get("thickness") for _, layer in self.structure_layers.items()]
        structural_materials = [getattr(mat_ns, layer.get("material")) for _, layer in self.structure_layers.items()]

        all_radii = np.cumsum([self.vacuum_chamber_radius] + structural_thicknesses)
        all_bottlenecks = np.cumsum([self.bottleneck_cylinder_radius] + structural_thicknesses)

        all_materials = [self.vacuum_chamber_material] + structural_materials

        regions = []

        for i in range(len(all_radii)):
            region = vacuum_vessel_region(
                self.first_plane,
                self.machine_half_length,
                all_radii[i],
                all_bottlenecks[i],
                self.cone_angle
            )
            regions.append(region)

        self._regions = regions

        for i in range(1, len(regions)):
            shell = regions[i] & ~regions[i-1]
            room_region &= ~shell
            cell = openmc.Cell(1000 + i, region=shell, fill=all_materials[i])
            self._cells.append(cell)
            self._full_region_list.append(shell)
            add_cell_callback(cell)

        return self._cells

    def get_cells(self):
        return self._cells

    def get_regions(self):
        return self._regions

    def get_full_region_list(self):
        return self._full_region_list

class SimpleCentralCell:
    """
    Class to construct the central cell region of a simple mirror machine.
    """
    def __init__(self, axial_length, layers):
        self.axial_length = axial_length
        self.layers = layers
        self._cells = []
        self._regions = []
        self._full_region_list = []

    def build_cells(self, vacuum_section_regions, room_region, add_cell_callback):
        central_cell_length_from_midplane = self.axial_length / 2
        
        # Define axial boundary planes
        left_plane = openmc.ZPlane(-central_cell_length_from_midplane)
        right_plane = openmc.ZPlane(central_cell_length_from_midplane)
        
        # Compute radii for each layer
        central_cell_cylinders_radii = np.cumsum([layer["thickness"] for layer in self.layers])
        
        # Create cells for each layer
        for i in range(len(central_cell_cylinders_radii)):
            if i == 0:
                region = (-right_plane & +left_plane & ~vacuum_section_regions[-1] & 
                         -openmc.ZCylinder(r=central_cell_cylinders_radii[i]))
            else:
                region = (-right_plane & +left_plane & 
                         -openmc.ZCylinder(r=central_cell_cylinders_radii[i]) & 
                         +openmc.ZCylinder(r=central_cell_cylinders_radii[i-1]))
            
            cell = openmc.Cell(2000 + i, region=region, fill=getattr(m, self.layers[i]["material"]))
            self._cells.append(cell)
            self._regions.append(region)
            self._full_region_list.append(region)
            room_region &= ~region
            add_cell_callback(cell)
            
        return self._cells

    def get_cells(self):
        return self._cells

    def get_regions(self):
        return self._regions

    def get_full_region_list(self):
        return self._full_region_list

class SimpleLFCoil:
    """
    Class to construct low-field coils for a simple mirror machine.
    """
    def __init__(self, shell_thicknesses, inner_dimensions, positions, materials):
        self.shell_thicknesses = shell_thicknesses
        self.inner_dimensions = inner_dimensions
        self.positions = positions
        self.materials = materials
        self._cells = []
        self._regions = []
        self._full_region_list = []

    def build_cells(self, central_cell_cylinders_radii, room_region, add_cell_callback):
        for i, center in enumerate(self.positions):
            # Generate cylindrical shell and inner region
            lf_coil_shell, lf_coil_inner = hollow_cylinder_with_shell(
                center,
                central_cell_cylinders_radii[-1],
                self.inner_dimensions["radial_thickness"],
                self.inner_dimensions["axial_length"],
                self.shell_thicknesses["front"],
                self.shell_thicknesses["back"],
                self.shell_thicknesses["axial"]
            )

            # Create cells for shell and inner region
            shell_cell = openmc.Cell(
                4000 + i,
                region=lf_coil_shell,
                fill=getattr(m, self.materials["shield"])
            )
            
            inner_cell = openmc.Cell(
                4100 + i,
                region=lf_coil_inner,
                fill=getattr(m, self.materials["magnet"])
            )

            self._cells.extend([shell_cell, inner_cell])
            self._regions.extend([lf_coil_shell, lf_coil_inner])
            self._full_region_list.extend([lf_coil_shell, lf_coil_inner])
            
            room_region &= ~lf_coil_shell & ~lf_coil_inner
            add_cell_callback(shell_cell)
            add_cell_callback(inner_cell)

        return self._cells

    def get_cells(self):
        return self._cells

    def get_regions(self):
        return self._regions

    def get_full_region_list(self):
        return self._full_region_list

class SimpleHFCoil:
    """
    Class to construct high-field coils for a simple mirror machine.
    """
    def __init__(self, magnet, casing_layers, shield):
        self.magnet = magnet
        self.casing_layers = casing_layers
        self.shield = shield
        self._cells = []
        self._regions = []
        self._full_region_list = []

    def build_cells(self, bottleneck_cylinders_radii, central_cell_length_from_midplane, room_region, add_cell_callback):
        # Compute HF coil center position
        casing_layers_thickness = np.array([layer["thickness"] for layer in self.casing_layers])
        hf_coil_center_z0 = (
            central_cell_length_from_midplane
            + self.shield["shield_central_cell_gap"]
            + self.shield["axial_thickness"][0]
            + np.sum(casing_layers_thickness)
            + self.magnet["axial_thickness"] / 2
        )

        # Compute inner radius for magnet and casings
        inner_radius_magnet_casings = (
            bottleneck_cylinders_radii[-1]
            + self.shield["radial_thickness"][0]
            + self.shield["radial_gap_before_casing"]
        )

        # Build left and right HF coils
        for z0, cell_id_base in [(hf_coil_center_z0, 6101), (-hf_coil_center_z0, 6102)]:
            # Build shield regions
            shield_regions = nested_cylindrical_shells(
                z0,
                inner_radius_magnet_casings,
                self.magnet["radial_thickness"],
                self.magnet["axial_thickness"],
                casing_layers_thickness,
                casing_layers_thickness,
                casing_layers_thickness
            )

            # Create magnet cell
            magnet_cell = openmc.Cell(
                cell_id_base,
                region=shield_regions[0],
                fill=getattr(m, self.magnet["material"])
            )
            self._cells.append(magnet_cell)
            self._regions.append(shield_regions[0])
            self._full_region_list.append(shield_regions[0])
            room_region &= ~shield_regions[0]
            add_cell_callback(magnet_cell)

            # Create casing cells
            for i in range(1, len(shield_regions)):
                casing_cell = openmc.Cell(
                    6500 + i,
                    region=shield_regions[i],
                    fill=getattr(m, self.casing_layers[i-1]["material"])
                )
                self._cells.append(casing_cell)
                self._regions.append(shield_regions[i])
                self._full_region_list.append(shield_regions[i])
                room_region &= ~shield_regions[i]
                add_cell_callback(casing_cell)

        return self._cells

    def get_cells(self):
        return self._cells

    def get_regions(self):
        return self._regions

    def get_full_region_list(self):
        return self._full_region_list

class SimpleEndCell:
    """
    Class to construct end cells for a simple mirror machine.
    """
    def __init__(self, axial_length, shell_thickness, diameter, shell_material, inner_material):
        self.axial_length = axial_length
        self.shell_thickness = shell_thickness
        self.outer_radius = diameter / 2
        self.shell_material = shell_material
        self.inner_material = inner_material
        self._cells = []
        self._regions = []
        self._full_region_list = []

    def build_cells(self, hf_center_z0, hf_magnet, hf_shield, casing_thicknesses, vacuum_section_regions, room_region, add_cell_callback):
        for z0, cell_id_base in [(hf_center_z0, 5001), (-hf_center_z0, 5003)]:
            # Calculate end cell position
            z0 = z0 + (
                hf_shield["axial_thickness"][0]
                + np.sum(casing_thicknesses)
                + hf_magnet["axial_thickness"] / 2
                + self.shell_thickness
                + self.axial_length / 2
            )

            # Create shell and inner regions
            shell_region, inner_region = cylinder_with_shell(
                z0,
                self.outer_radius - self.shell_thickness,
                self.axial_length,
                self.shell_thickness
            )

            # Exclude vacuum section regions
            shell_region &= ~vacuum_section_regions[-1]
            inner_region &= ~vacuum_section_regions[-1]

            # Create cells
            shell_cell = openmc.Cell(
                cell_id_base,
                region=shell_region,
                fill=getattr(m, self.shell_material)
            )
            
            inner_cell = openmc.Cell(
                cell_id_base + 1,
                region=inner_region,
                fill=getattr(m, self.inner_material)
            )

            self._cells.extend([shell_cell, inner_cell])
            self._regions.extend([shell_region, inner_region])
            self._full_region_list.extend([shell_region, inner_region])
            
            room_region &= ~shell_region & ~inner_region
            add_cell_callback(shell_cell)
            add_cell_callback(inner_cell)

        return self._cells

    def get_cells(self):
        return self._cells

    def get_regions(self):
        return self._regions

    def get_full_region_list(self):
        return self._full_region_list

class SimpleTallyBuilder:
    """
    Collects tally descriptors and generates OpenMC Tally objects for simple mirror model.
    Supports both cell and mesh tallies.
    """
    def __init__(self):
        self._tallies = []
        self._tally_descriptors = []

    def add_descriptors(self, descriptors):
        for desc in descriptors:
            self._add_cell_tallies(desc)
            self._add_mesh_tallies(desc)

    def _add_cell_tallies(self, desc):
        for i, entry in enumerate(desc.get("cell_tallies", [])):
            filters = [openmc.CellFilter(desc["cell"])] + strings_to_openmc_filters(entry.get("filters", []))
            tally = openmc.Tally(name=f"{desc['type']}_{desc['location']}_{desc['description']}_cell_tally_{i+1}")
            tally.filters = filters
            tally.scores = entry.get("scores", [])
            if "nuclides" in entry:
                tally.nuclides = entry["nuclides"]
            self._tallies.append(tally)

    def _add_mesh_tallies(self, desc):
        for i, entry in enumerate(desc.get("mesh_tallies", [])):
            mesh = hollow_mesh_from_domain(desc["cell"], entry["dimensions"])
            filters = [openmc.MeshFilter(mesh)] + strings_to_openmc_filters(entry.get("filters", []))
            tally = openmc.Tally(name=f"{desc['type']}_{desc['location']}_{desc['description']}_mesh_tally_{i+1}")
            tally.filters = filters
            tally.scores = entry.get("scores", [])
            if "nuclides" in entry:
                tally.nuclides = entry["nuclides"]
            self._tallies.append(tally)

    def get_tallies(self):
        return self._tallies

    def add_tally_descriptor(self, descriptor):
        self._tally_descriptors.append(descriptor)

class SimpleMachineBuilder:
    """
    Modular builder for simple mirror fusion machines.
    Components: VV, FW, Central Cell, LF/HF coils, End cells, Room, Tallies.
    """
    def __init__(self, vv_params, central_cell_params, lf_coil_params, hf_coil_params, end_cell_params, material_ns):
        self.vv_params = vv_params
        self.central_cell_params = central_cell_params
        self.lf_coil_params = lf_coil_params
        self.hf_coil_params = hf_coil_params
        self.end_cell_params = end_cell_params
        self.material_ns = material_ns

        # Initialize bounding box and room
        self._bounding_surface = openmc.model.RectangularParallelepiped(
            xmin=-400, xmax=400, ymin=-400, ymax=400,
            zmin=-1350, zmax=1350, boundary_type='vacuum'
        )
        self._room_region = -self._bounding_surface

        # Initialize component builders
        self.vv_builder = None
        self.fw_builder = None
        self.central_cell_builder = None
        self.lf_coil_builder = None
        self.hf_coil_builder = None
        self.end_cell_builder = None

        # Initialize universe and storage
        self._universe = openmc.Universe(786)
        self._all_cells = []
        self._regions = {
            "vv": None,
            "fw": None,
            "central_cell": None,
            "lf_coils": None,
            "hf_coils": None,
            "end_cell": None
        }

        # Initialize tally builder
        self.tally_builder = SimpleTallyBuilder()

        # Create room cell
        room_cell = openmc.Cell(100, region=self._room_region, fill=getattr(self.material_ns, "air"))
        self._universe.add_cell(room_cell)

    def _add_cell(self, cell):
        self._universe.add_cell(cell)
        self._all_cells.append(cell)

    def build_vacuum_vessel(self):
        self.vv_builder = SimpleVacuumVessel(
            self.vv_params["machine_length_from_midplane"],
            self.vv_params["first_vv_plane_from_midplane"],
            self.vv_params["cone_angle"],
            self.vv_params["vacuum_chamber"],
            self.vv_params["bottleneck_cylinder"],
            getattr(self.material_ns, "vacuum")
        )
        vv_cell = self.vv_builder.build_cell()
        self._add_cell(vv_cell)
        self._regions["vv"] = self.vv_builder.get_region()
        self._room_region &= ~self.vv_builder.get_region()

    def build_first_wall(self):
        self.fw_builder = SimpleFWStructure(
            self.vv_params["structure"],
            self.vv_params["vacuum_chamber"]["radius"],
            self.vv_params["bottleneck_cylinder"]["cylinder_radius"],
            getattr(self.material_ns, self.vv_params["vacuum_chamber"]["material"]),
            self.vv_params["first_vv_plane_from_midplane"],
            self.vv_params["machine_length_from_midplane"],
            self.vv_params["cone_angle"]
        )
        self.fw_builder.build_cells(
            self.material_ns,
            self._room_region,
            self._add_cell
        )
        self._regions["fw"] = self.fw_builder.get_regions()

    def build_central_cell(self):
        self.central_cell_builder = SimpleCentralCell(
            self.central_cell_params["axial_length"],
            self.central_cell_params["layers"]
        )
        self.central_cell_builder.build_cells(
            self.fw_builder.get_regions(),
            self._room_region,
            self._add_cell
        )
        self._regions["central_cell"] = self.central_cell_builder.get_regions()

    def build_lf_coils(self):
        self.lf_coil_builder = SimpleLFCoil(
            self.lf_coil_params["shell_thicknesses"],
            self.lf_coil_params["inner_dimensions"],
            self.lf_coil_params["positions"],
            self.lf_coil_params["materials"]
        )
        self.lf_coil_builder.build_cells(
            [r.r for r in self.fw_builder.get_regions()],
            self._room_region,
            self._add_cell
        )
        self._regions["lf_coils"] = self.lf_coil_builder.get_regions()

        # Add LF coil tallies
        if "tallies" in self.lf_coil_params:
            for i, cell in enumerate(self.lf_coil_builder.get_cells()):
                if cell.id >= 4100:  # LF coil inner cells
                    self.tally_builder.add_tally_descriptor({
                        "type": "lf_coil",
                        "location": f"coil_{i}",
                        "description": "magnet",
                        "cell": cell,
                        "cell_tallies": self.lf_coil_params["tallies"].get("cell_tallies", []),
                        "mesh_tallies": self.lf_coil_params["tallies"].get("mesh_tallies", [])
                    })

    def build_hf_coils(self):
        self.hf_coil_builder = SimpleHFCoil(
            self.hf_coil_params["magnet"],
            self.hf_coil_params["casing_layers"],
            self.hf_coil_params["shield"]
        )
        self.hf_coil_builder.build_cells(
            [r.r for r in self.fw_builder.get_regions()],
            self.central_cell_params["axial_length"] / 2,
            self._room_region,
            self._add_cell
        )
        self._regions["hf_coils"] = self.hf_coil_builder.get_regions()

        # Add HF coil tallies
        if "tallies" in self.hf_coil_params:
            for cell in self.hf_coil_builder.get_cells():
                if cell.id == 6101:  # Right HF coil
                    self.tally_builder.add_tally_descriptor({
                        "type": "hf_coil",
                        "location": "right",
                        "description": "magnet",
                        "cell": cell,
                        "cell_tallies": self.hf_coil_params["tallies"].get("cell_tallies", []),
                        "mesh_tallies": self.hf_coil_params["tallies"].get("mesh_tallies", [])
                    })

    def build_end_cells(self):
        self.end_cell_builder = SimpleEndCell(
            self.end_cell_params["axial_length"],
            self.end_cell_params["shell_thickness"],
            self.end_cell_params["diameter"],
            self.end_cell_params["shell_material"],
            self.end_cell_params["inner_material"]
        )
        self.end_cell_builder.build_cells(
            self.central_cell_params["axial_length"] / 2,
            self.hf_coil_params["magnet"],
            self.hf_coil_params["shield"],
            [layer["thickness"] for layer in self.hf_coil_params["casing_layers"]],
            self.fw_builder.get_regions(),
            self._room_region,
            self._add_cell
        )
        self._regions["end_cell"] = self.end_cell_builder.get_regions()

    def get_universe(self):
        return self._universe

    def get_all_cells(self):
        return self._all_cells

    def get_all_regions(self):
        return self._regions

    def get_all_tallies(self):
        """
        Get all tallies configured for the model.
        """
        self.tally_builder.add_descriptors(self.tally_builder._tally_descriptors)
        return self.tally_builder.get_tallies()

def build_simple_model_from_input(input_data, output_dir="."):
    """
    Build a complete simple mirror model from input parameters.
    """
    # Parse input parameters
    vv_params, central_cell_params, lf_coil_params, hf_coil_params, end_cell_params = parse_simple_machine_input(input_data)

    # Create builder
    builder = SimpleMachineBuilder(
        vv_params, central_cell_params,
        lf_coil_params, hf_coil_params, end_cell_params,
        material_ns=m
    )

    # Build components
    builder.build_vacuum_vessel()
    builder.build_first_wall()
    builder.build_central_cell()
    builder.build_lf_coils()
    builder.build_hf_coils()
    builder.build_end_cells()

    # Create geometry
    universe_machine = builder.get_universe()
    geometry = openmc.Geometry([openmc.Cell(fill=universe_machine)])
    geometry.merge_surfaces = True

    # Generate cross-sectional plot
    geometry.root_universe.plot(
        basis='xz',
        width=(1000, 2800),
        pixels=(700, 700),
        color_by='material'
    )
    plt.savefig('simple_mirror_cross_section.png', bbox_inches="tight")

    # Export geometry
    geometry.export_to_xml("geometry.xml")

    # Create and export tallies
    tallies = openmc.Tallies(builder.get_all_tallies())
    tallies.export_to_xml("tallies.xml")

    return geometry, tallies

        def get_neutron_wall_loading_model(self):
        """
        Get neutron wall loading model for the simple mirror machine.
        """
        # Create materials
        tungsten = openmc.Material(31, name="tungsten")
        tungsten.set_density("g/cm3", 19.25)
        tungsten.add_element("W", 100)

        my_materials = openmc.Materials([tungsten])
        my_materials.export_to_xml(path='nwl/materials.xml')

        # Get vacuum vessel geometric extents
        vv_geometric_extents = self.vv_builder.get_vv_components()
        
        # Extract key dimensions
        left_end_vv = vv_geometric_extents['z_planes']['left_end']
        right_end_vv = vv_geometric_extents['z_planes']['right_end']
        central_radius_vv = vv_geometric_extents['central_cylinder_radius']
        bottleneck_radius_vv = vv_geometric_extents['bottleneck_radius']

        # Define model boundaries
        model_right = openmc.ZPlane(z0=right_end_vv.z0 + 10, boundary_type='vacuum')
        model_left = openmc.ZPlane(z0=left_end_vv.z0 - 10, boundary_type='vacuum')
        model_side_outer_radius = openmc.ZCylinder(r=central_radius_vv + 25, boundary_type='vacuum', surface_id=123)

        # Define first wall armor planes
        centralfwarmor_right = vv_geometric_extents['z_planes']['central_right']
        centralfwarmor_left = vv_geometric_extents['z_planes']['central_left']
        
        rightendfwarmor_right = right_end_vv
        rightendfwarmor_left = vv_geometric_extents['z_planes']['right_cone_plane']
        
        leftendfwarmor_right = vv_geometric_extents['z_planes']['left_cone_plane']
        leftendfwarmor_left = left_end_vv

        # Set transmission boundaries
        leftendfwarmor_left.boundary_type = 'transmission'
        leftendfwarmor_right.boundary_type = 'transmission'
        rightendfwarmor_right.boundary_type = 'transmission'

        # Define plasma region
        centralplasma_outer_radius = openmc.ZCylinder(r=central_radius_vv)

        return {
            'materials': my_materials,
            'surfaces': {
                'model_boundaries': {
                    'right': model_right,
                    'left': model_left,
                    'outer_radius': model_side_outer_radius
                },
                'first_wall': {
                    'central': {
                        'right': centralfwarmor_right,
                        'left': centralfwarmor_left
                    },
                    'right': {
                        'right': rightendfwarmor_right,
                        'left': rightendfwarmor_left
                    },
                    'left': {
                        'right': leftendfwarmor_right,
                        'left': leftendfwarmor_left
                    }
                },
                'plasma': {
                    'outer_radius': centralplasma_outer_radius
                }
            }
        }
