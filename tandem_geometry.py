import openmc
import numpy as np
import openmc.lib
import matplotlib.pyplot as plt
from geometry_lib import *

from types import SimpleNamespace

def parse_machine_input(input_data, material_ns):
    vv_info = input_data["vacuum_vessel"]
    cc_vv = vv_info["central_cell"]
    ep_vv = vv_info["end_plug"]
    sep = vv_info["central_cell_end_plug_separation_distance"]
    bottleneck = vv_info["bottleneck radius"]
    end_axial = vv_info["end_axial_distance"]

    # --- Midplanes ---
    right_midplane = cc_vv["central_axis_length"]/2 + sep + ep_vv["central_axis_length"]/2
    left_midplane = -right_midplane

    # --------- 1. Vacuum Vessel Params ---------
    vv_params = {
        "central_cell_vv_parameters": SimpleNamespace(radius=cc_vv["central_radius"],
                                              central_axis_length=cc_vv["central_axis_length"],
                                              outer_length=cc_vv["outer_axial_length"]),
        "end_plug_vv_parameters": SimpleNamespace(radius=ep_vv["central_radius"],
                                          central_axis_length=ep_vv["central_axis_length"],
                                          outer_length=ep_vv["outer_axial_length"]),

        "central_cell_end_plug_separation_distance": sep,
        "bottleneck_radius": bottleneck,
        "end_axial_distance": end_axial,
        "vacuum_material": getattr(material_ns, "vacuum"),
    }

    # --------- 2. First Wall Params ---------
    fw_info = input_data["first_wall"]
    cc_fw_layers = fw_info["central_cell"]["layers"]
    ep_fw_layers = fw_info["end_plug"]["layers"]

    fw_params = {
        "left_fw": ep_fw_layers,
        "central_fw": cc_fw_layers,
        "right_fw": ep_fw_layers,  # reuse
        "param_dict": {
            "central": (cc_vv["outer_axial_length"], cc_vv["central_axis_length"], sep/2, sep/2, 0.0),
            "left": (ep_vv["outer_axial_length"], ep_vv["central_axis_length"], end_axial, sep/2, left_midplane),
            "right": (ep_vv["outer_axial_length"], ep_vv["central_axis_length"], sep/2, end_axial, right_midplane),
        },
        "base_radii": {
            "central": cc_vv["central_radius"],
            "left": ep_vv["central_radius"],
            "right": ep_vv["central_radius"],
        },
        "bottleneck_radii": {
            "central": bottleneck,
            "left": bottleneck,
            "right": bottleneck,
        },
    }

    # --------- 3. Central Cylinder Params ---------
    cc_cyl_layers = input_data["central_cell"]["blanket"]["layers"]
    ep_cyl_layers = input_data["end_plug"]["central_cylinder"]["layers"]

    cyl_params = {
        "left": ep_cyl_layers,
        "central": cc_cyl_layers,
        "right": ep_cyl_layers,
        "midplanes": {
            "central": 0.0,
            "left": left_midplane,
            "right": right_midplane
        },
        "axial_lengths": {
            "central": input_data["central_cell"]["blanket"]["axial_length"],
            "left": input_data["end_plug"]["central_cylinder"]["axial_length"],
            "right": input_data["end_plug"]["central_cylinder"]["axial_length"],
        }
    }

    return vv_params, fw_params, cyl_params



class TandemVacuumVessel:
    """
    Class to construct OpenMC cells and expose regions for a tandem mirror vacuum vessel:
    - Central cell vacuum vessel
    - Two end plug vacuum vessels (left and right)
    """

    def __init__(self, end_plug_vv_parameters, central_cell_vv_parameters, central_cell_end_plug_separation_distance, bottleneck_radius, end_axial_distance, vacuum_material):
        self.end_plug_vv_parameters = end_plug_vv_parameters
        self.central_cell_vv_parameters = central_cell_vv_parameters
        self.central_cell_end_plug_separation_distance = central_cell_end_plug_separation_distance
        self.bottleneck_radius = bottleneck_radius
        self.end_axial_distance = end_axial_distance
        self.vacuum_material = vacuum_material
        self._regions = {}
        self._cells = {}
        self._z_extents = {}

    def build_cells(self):
        end_plug_vv_central_radius = self.end_plug_vv_parameters.radius
        end_plug_vv_central_axis_length = self.end_plug_vv_parameters.central_axis_length
        end_plug_vv_outer_length = self.end_plug_vv_parameters.outer_length

        central_cell_central_radius = self.central_cell_vv_parameters.radius
        central_cell_vv_central_axis_length = self.central_cell_vv_parameters.central_axis_length
        central_cell_vv_outer_length = self.central_cell_vv_parameters.outer_length

        right_midplane = (
            central_cell_vv_central_axis_length / 2
            + self.central_cell_end_plug_separation_distance
            + end_plug_vv_central_axis_length / 2
        )
        left_midplane = -right_midplane

        self._z_extents["central"] = (-central_cell_vv_outer_length/2, central_cell_vv_outer_length/2)
        self._z_extents["right"] = (
            right_midplane - end_plug_vv_outer_length/2,
            right_midplane + end_plug_vv_outer_length/2
        )
        self._z_extents["left"] = (
            left_midplane - end_plug_vv_outer_length/2,
            left_midplane + end_plug_vv_outer_length/2
        )

        self._regions["central"] = redefined_vacuum_vessel_region(
            central_cell_vv_outer_length,
            central_cell_vv_central_axis_length,
            central_cell_central_radius,
            self.bottleneck_radius,
            self.central_cell_end_plug_separation_distance / 2,
            self.central_cell_end_plug_separation_distance / 2,
            axial_midplane=0.0,
        )

        self._regions["right"] = redefined_vacuum_vessel_region(
            end_plug_vv_outer_length,
            end_plug_vv_central_axis_length,
            end_plug_vv_central_radius,
            self.bottleneck_radius,
            self.central_cell_end_plug_separation_distance / 2,
            self.end_axial_distance,
            axial_midplane=right_midplane,
        )

        self._regions["left"] = redefined_vacuum_vessel_region(
            end_plug_vv_outer_length,
            end_plug_vv_central_axis_length,
            end_plug_vv_central_radius,
            self.bottleneck_radius,
            self.end_axial_distance,
            self.central_cell_end_plug_separation_distance / 2,
            axial_midplane=left_midplane,
        )

        self._cells = {
            key: openmc.Cell(name=f"{key}_vv_cell", region=region, fill=self.vacuum_material)
            for key, region in self._regions.items()
        }

        return self._cells

    def get_regions(self):
        return self._regions

    def get_cells(self):
        return self._cells

    def get_z_extents(self):
        return self._z_extents
    

class TandemFWStructure:
    """
    Class to construct concentric cylindrical structural layers around vacuum vessel sections
    in a tandem mirror (left plug, central cell, right plug).
    """

    def __init__(self, left_fw, central_fw, right_fw=None):
        self.left_fw = left_fw
        self.central_fw = central_fw
        self.right_fw = right_fw if right_fw else left_fw
        self._cells_by_region = {"left": [], "central": [], "right": []}
        self._regions_by_region = {"left": [], "central": [], "right": []}
        self._z_extents = {}

    def _process_layers(self, layer_dict, base_radius, bottleneck_radius, mat_ns):
        thicknesses = [layer["thickness"] for layer in layer_dict]
        materials = [getattr(mat_ns, layer["material"]) for layer in layer_dict]
        radii = np.cumsum([base_radius] + thicknesses)
        bottlenecks = np.cumsum([bottleneck_radius] + thicknesses)
        return radii, bottlenecks, materials

    def build_all_sections(self, region_dict, param_dict, mat_ns, base_radii, bottleneck_radii, add_cell_callback):
        for key, fw_layers in zip(["left", "central", "right"], [self.left_fw, self.central_fw, self.right_fw]):
            outer_len, central_len, left_len, right_len, midplane = param_dict[key]
            base_region = region_dict[key]
            base_r = base_radii[key]
            bott_r = bottleneck_radii[key]

            z_min = midplane - central_len / 2
            z_max = midplane + central_len / 2
            self._z_extents[key] = (z_min, z_max)

            vac_radii, bott_radii, materials = self._process_layers(fw_layers, base_r, bott_r, mat_ns)

            regions = [base_region]
            cells = []

            for i in range(1, len(vac_radii)):
                new_region = redefined_vacuum_vessel_region(
                    outer_len,
                    central_len,
                    vac_radii[i],
                    bott_radii[i],
                    left_len,
                    right_len,
                    axial_midplane=midplane,
                )
                shell = new_region & ~regions[i - 1]
                cell = openmc.Cell(region=shell, fill=materials[i-1])
                cells.append(cell)
                add_cell_callback(cell)
                regions.append(new_region)

            self._cells_by_region[key] = cells
            self._regions_by_region[key] = regions[1:]  # exclude base region

        return self._cells_by_region

    def get_cells_by_region(self):
        return self._cells_by_region

    def get_regions_by_region(self):
        return self._regions_by_region

    def get_z_extents(self):
        return self._z_extents
    
class TandemCentralCylinders:
    """
    Builds concentric cylindrical axial layers for left plug, central cell, and right plug regions.
    Uses outermost FW region as base to define first shell.
    """

    def __init__(self, left_layers, central_layers, right_layers=None):
        self.left_layers = left_layers
        self.central_layers = central_layers
        self.right_layers = right_layers if right_layers else left_layers
        self._cells_by_region = {"left": [], "central": [], "right": []}
        self._regions_by_region = {"left": [], "central": [], "right": []}
        self._z_extents = {}

    def build_all_sections(self, base_radii, midplanes, axial_lengths, base_regions, mat_ns, add_cell_callback):
        for key, layers in zip(["left", "central", "right"], [self.left_layers, self.central_layers, self.right_layers]):
            z0 = midplanes[key]
            half_len = axial_lengths[key] / 2
            zmin, zmax = z0 - half_len, z0 + half_len
            self._z_extents[key] = (zmin, zmax)

            # Radial shell construction
            thicknesses = [layer["thickness"] for layer in layers]
            materials = [getattr(mat_ns, layer["material"]) for layer in layers]
            radii = base_radii[key] + np.cumsum(thicknesses)

            print(radii)

            zplanes = -openmc.ZPlane(z0=zmax) & +openmc.ZPlane(z0=zmin)
            prev_cyl = openmc.ZCylinder(r=radii[0])

            first_region = zplanes & ~base_regions[key] & -prev_cyl
            self._regions_by_region[key].append(first_region)

            cells = [openmc.Cell(region=first_region, fill=materials[0])]
            add_cell_callback(cells[-1])

            for i in range(1, len(radii)):
                next_cyl = openmc.ZCylinder(r=radii[i])
                region = zplanes & -next_cyl & +prev_cyl
                self._regions_by_region[key].append(region)
                cell = openmc.Cell(region=region, fill=materials[i])
                cells.append(cell)
                add_cell_callback(cell)
                prev_cyl = next_cyl

            self._cells_by_region[key] = cells

    def get_cells_by_region(self):
        return self._cells_by_region

    def get_regions_by_region(self):
        return self._regions_by_region

    def get_z_extents(self):
        return self._z_extents


class TandemMachineBuilder:
    """
    Assembles vacuum vessel, first wall, and central cylinder structures into a complete tandem mirror model.
    """

    def __init__(self, vv_params, fw_params, central_cyl_params, material_ns):
        self.vv_params = vv_params
        self.fw_params = fw_params
        self.central_cyl_params = central_cyl_params
        self.material_ns = material_ns

        self.vv_builder = None
        self.fw_builder = None
        self.ccyl_builder = None

        self._universe = openmc.Universe(name="TandemMachine")
        self._all_cells = []
        self._z_extents = {}
        self._regions = {"vv": {}, "fw": {}, "cyl": {}}

    def _add_cell(self, cell):
        self._universe.add_cell(cell)
        self._all_cells.append(cell)

    def build(self):
        # ----- 1. Vacuum Vessel -----
        self.vv_builder = TandemVacuumVessel(**self.vv_params)
        vv_cells = self.vv_builder.build_cells()

        for cell in vv_cells.values():
            self._add_cell(cell)

        self._regions["vv"] = self.vv_builder.get_regions()
        self._z_extents["vv"] = self.vv_builder.get_z_extents()

        # ----- 2. First Wall -----
        self.fw_builder = TandemFWStructure(
            self.fw_params["left_fw"],
            self.fw_params["central_fw"],
            self.fw_params.get("right_fw", None)
        )

        fw_cells = self.fw_builder.build_all_sections(
            region_dict=self._regions["vv"],
            param_dict=self.fw_params["param_dict"],
            mat_ns=self.material_ns,
            base_radii=self.fw_params["base_radii"],
            bottleneck_radii=self.fw_params["bottleneck_radii"],
            add_cell_callback=self._add_cell
        )

        self._regions["fw"] = self.fw_builder.get_regions_by_region()
        self._z_extents["fw"] = self.fw_builder.get_z_extents()

        # ----- 3. Central Cylinders -----
        base_regions = {
            key: self._regions["fw"][key][-1]
            for key in ["left", "central", "right"]
        }

        base_radii = {
            key: self.fw_params["base_radii"][key] + sum(
                layer["thickness"] for layer in layers
            ) for key, layers in zip(
                ["left", "central", "right"],
                [self.fw_params["left_fw"], self.fw_params["central_fw"], self.fw_params.get("right_fw", self.fw_params["left_fw"])]
            )
        }

        self.ccyl_builder = TandemCentralCylinders(
            self.central_cyl_params["left"],
            self.central_cyl_params["central"],
            self.central_cyl_params.get("right", None)
        )

        self.ccyl_builder.build_all_sections(
            base_radii=base_radii,
            midplanes=self.central_cyl_params["midplanes"],
            axial_lengths=self.central_cyl_params["axial_lengths"],
            base_regions=base_regions,
            mat_ns=self.material_ns,
            add_cell_callback=self._add_cell
        )

        self._regions["cyl"] = self.ccyl_builder.get_regions_by_region()
        self._z_extents["cyl"] = self.ccyl_builder.get_z_extents()

    # ---------- Accessors ----------
    def get_universe(self):
        return self._universe

    def get_all_cells(self):
        return self._all_cells

    def get_all_regions(self):
        return self._regions

    def get_z_extents(self):
        return self._z_extents
