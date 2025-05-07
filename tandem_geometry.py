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
            "central": (cc_vv["outer_axial_length"], cc_vv["central_axis_length"], sep/10, sep/10, 0.0),
            "left": (ep_vv["outer_axial_length"], ep_vv["central_axis_length"], end_axial, 9*sep/10, left_midplane),
            "right": (ep_vv["outer_axial_length"], ep_vv["central_axis_length"], 9*sep/10, end_axial, right_midplane),
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

    # --------- 4. LF Coil Params ---------
    lf_coil_params = {
        "central": input_data.get("central_cell", {}).get("lf_coil", {}),
        "end": input_data.get("end_plug", {}).get("lf_coil", {})
    }

    # --------- 5. HF Coil Params ---------
    hf_coil_data = input_data.get("end_plug", {}).get("hf_coil", {})

    hf_coil_params = {
        key: {
            "magnet": hf_coil_data.get(key, {}).get("magnet", {}),
            "casing_layers": {
                "thicknesses": [layer["thickness"] for layer in hf_coil_data.get(key, {}).get("casing_layers", [])],
                "materials": [layer["material"] for layer in hf_coil_data.get(key, {}).get("casing_layers", [])]
            },
            "shield": hf_coil_data.get(key, {}).get("shield", {})
        } for key in ["left", "right"]}

    # ---------- 6. End Cell Params --------
    end_cell_params = input_data.get("end_cell", {})

    return vv_params, fw_params, cyl_params, lf_coil_params, hf_coil_params, end_cell_params



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
        self._region_list = []
        self._components = {}
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

        self._regions["central"], self._components["central"]  = redefined_vacuum_vessel_region(
            central_cell_vv_outer_length,
            central_cell_vv_central_axis_length,
            central_cell_central_radius,
            self.bottleneck_radius,
            self.central_cell_end_plug_separation_distance / 10,
            self.central_cell_end_plug_separation_distance / 10,
            axial_midplane=0.0,
        )


        self._regions["right"], self._components["right"]  = redefined_vacuum_vessel_region(
            end_plug_vv_outer_length,
            end_plug_vv_central_axis_length,
            end_plug_vv_central_radius,
            self.bottleneck_radius,
            9 * self.central_cell_end_plug_separation_distance / 10,
            self.end_axial_distance,
            axial_midplane=right_midplane,
        )

        self._regions["left"],  self._components["left"] = redefined_vacuum_vessel_region(
            end_plug_vv_outer_length,
            end_plug_vv_central_axis_length,
            end_plug_vv_central_radius,
            self.bottleneck_radius,
            self.end_axial_distance,
            9 * self.central_cell_end_plug_separation_distance / 10,
            axial_midplane=left_midplane,
        )

        self._region_list.append(self._regions["central"])
        self._region_list.append(self._regions["left"])
        self._region_list.append(self._regions["right"])

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
    
    def get_full_region_list(self):
        return self._region_list
    

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
        self._fw_radii = {"left": [], "central": [], "right": []}
        self._region_list = []

    def _process_layers(self, layer_dict, base_radius, bottleneck_radius, mat_ns):
        thicknesses = [layer["thickness"] for layer in layer_dict]
        materials = [getattr(mat_ns, layer["material"]) for layer in layer_dict]
        radii = np.cumsum([base_radius] + thicknesses)
        bottlenecks = np.cumsum([bottleneck_radius] + thicknesses)
        return radii, bottlenecks, materials

    def build_all_sections(self, region_dict, param_dict, mat_ns, base_radii, bottleneck_radii, add_cell_callback):
        for key, fw_layers in zip(["left", "central", "right"], [self.left_fw, self.central_fw, self.right_fw]):
            outer_len, central_len, left_len, right_len, midplane = param_dict[key]
            
            # Flags
            # print(f"For {[key]} section the left legnth is {left_len}")
            # print(f"For {[key]} section the right legnth is {right_len}")
            
            base_region = region_dict[key]
            base_r = base_radii[key]
            bott_r = bottleneck_radii[key]

            z_min = midplane - central_len / 2
            z_max = midplane + central_len / 2
            self._z_extents[key] = (z_min, z_max)

            vac_radii, bott_radii, materials = self._process_layers(fw_layers, base_r, bott_r, mat_ns)

            self._fw_radii[key] = bott_radii

            regions = [base_region]
            self._region_list.append(base_region)

            cells = []

            for i in range(1, len(vac_radii)):
                new_region, parts = redefined_vacuum_vessel_region(
                    outer_len,
                    central_len,
                    vac_radii[i],
                    bott_radii[i],
                    left_len,
                    right_len,
                    axial_midplane=midplane,
                )
                shell = new_region & ~regions[i - 1]
                self._region_list.append(shell)
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
    
    def fw_radii(self):
        return self._fw_radii
    
    def get_full_region_list(self):
        return self._region_list
        

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
        self._region_list = []
        self._radii_by_region = {"left": [], "central": [], "right": []}

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

            self._radii_by_region[key] = radii

            zplanes = -openmc.ZPlane(z0=zmax) & +openmc.ZPlane(z0=zmin)
            prev_cyl = openmc.ZCylinder(r=radii[0])

            first_region = zplanes & ~base_regions[key] & -prev_cyl
            self._regions_by_region[key].append(first_region)
            self._region_list.append(first_region)

            cells = [openmc.Cell(region=first_region, fill=materials[0])]
            add_cell_callback(cells[-1])

            for i in range(1, len(radii)):
                next_cyl = openmc.ZCylinder(r=radii[i])
                region = zplanes & -next_cyl & +prev_cyl
                self._regions_by_region[key].append(region)
                self._region_list.append(region)
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
    
    def get_radii_by_region(self):
        return self._radii_by_region
    
    def get_full_region_list(self):
        return self._region_list
    
class LFCoilBuilder:
    """
    Builds LF Coils over the central cylinder of the central cell, the left end plug, and the right end plug.
    Uses outermost radii of central cylinder regions in each section and the midplane of each section.
    """

    def __init__(self, lf_coil_params):
        """
        Parameters
        ----------
        lf_coil_params : dict
            Dictionary with 'central' and 'end' keys from parse_machine_input.
        """

        self.coil_data = {
            "central": lf_coil_params.get("central", {}),
            "end": lf_coil_params.get("end", {})
        }

        self._cells_by_region = {region: {"coil": [], "shield": []} for region in ["left", "central", "right"]}
        self._regions_by_region = {region: {"coil": [], "shield": []} for region in ["left", "central", "right"]}
        self._region_list = []

    def build_all_sections(self, outer_radii_by_region, midplanes, mat_ns, add_cell_callback):
        """
        Builds coils for all three sections: left, central, and right.

        Parameters
        ----------
        outer_radii_by_region : dict
            Keys: 'left', 'central', 'right' → outermost radius from each central cylinder stack.
        midplanes : dict
            Keys: 'left', 'central', 'right' → axial center of each region.
        mat_ns : module
            Material namespace module.
        """
        i = 0
        for key in ["left", "central", "right"]:
            data_key = "central" if key == "central" else "end"
            coil_data = self.coil_data.get(data_key, {})

            positions = coil_data.get("positions", [])
            inner_dims = coil_data.get("inner_dimensions", {})
            shell_dims = coil_data.get("shell_thicknesses", {})
            materials = coil_data.get("materials", {})

            coil_regions = []
            shield_regions = []
            coil_cells = []
            shield_cells = []

            for z_offset in positions:
                z_center = midplanes[key] + z_offset

                shell_region, coil_region = hollow_cylinder_with_shell(
                    z_center,
                    outer_radii_by_region[key],
                    inner_dims["radial_thickness"],
                    inner_dims["axial_length"],
                    shell_dims.get("front", 0),
                    shell_dims.get("back", 0),
                    shell_dims.get("axial", 0))
                
                shield_cell = openmc.Cell(
                    cell_id=4000 + i,
                    region=shell_region,
                    fill=getattr(mat_ns, materials["shield"])
                )
                coil_cell = openmc.Cell(
                    cell_id=4100 + i,
                    region=coil_region,
                    fill=getattr(mat_ns, materials["magnet"])
                )

                add_cell_callback(shield_cell)
                add_cell_callback(coil_cell)

                self._regions_by_region[key]["shield"].append(shell_region)
                self._regions_by_region[key]["coil"].append(coil_region)
                self._region_list.append(shell_region)
                self._region_list.append(coil_region)
                self._cells_by_region[key]["shield"].append(shield_cell)
                self._cells_by_region[key]["coil"].append(coil_cell)
                i += 1

    def get_cells_by_region(self):
        return self._cells_by_region

    def get_regions_by_region(self):
        return self._regions_by_region
    
    def get_full_region_list(self):
        return self._region_list

class HFCoilBuilder:
    """
    Builds HF Coils around the left and right end plugs.
    Verifies that the bore radius matches the combined bottleneck + shield + casing geometry.
    """

    def __init__(self, hf_coil_params):
        """
        Parameters
        ----------
        hf_coil_params : dict
            Output from parse_machine_input()[4]
        """
        self.coil_data = hf_coil_params
        self._cells_by_region = {region: {"coil": [], "shield": []} for region in ["left", "right"]}
        self._regions_by_region = {region: {"coil": [], "shield": []} for region in ["left", "right"]}
        self._outermost_coil_z0 = {"left": [] , "right": [] }
        self._region_list = []

    def build_all_sections(self, bottleneck_radii, cc_half_lengths, midplanes, mat_ns, add_cell_callback):
        i = 0
        for key in ["left", "right"]:
            hf_data = self.coil_data[key]
            magnet = hf_data["magnet"]
            casing = hf_data["casing_layers"]
            shield = hf_data["shield"]

            casing_thicknesses = np.array(casing["thicknesses"])
            casing_materials = [getattr(mat_ns, mat) for mat in casing["materials"]]
            bore_radius = magnet["bore_radius"]
            bottleneck_r = bottleneck_radii[key]

            # Check bore radius
            computed_bore = bottleneck_r + shield["radial_thickness"][0] + shield["radial_gap_before_casing"] + np.sum(casing_thicknesses)
            if not np.isclose(computed_bore, bore_radius, atol=1e-4):
                raise ValueError(f"HF {key} bore radius mismatch: expected {bore_radius}, got {computed_bore:.4f}")

            # Compute base radius before casing begins
            r_inner = bottleneck_r + shield["radial_thickness"][0] + shield["radial_gap_before_casing"]

            for direction in ["inward", "outward"]:
                if direction == "inward":
                    axial_sign = 1 if key == "left" else -1
                    axial_thickness = shield["axial_thickness"][0]
                else:
                    axial_sign = -1 if key == "left" else 1
                    axial_thickness = shield["axial_thickness"][0]
                z0 = (
                    midplanes[key]
                    + axial_sign * (
                        cc_half_lengths[key]
                        + shield["shield_central_cell_gap"]
                        + axial_thickness
                        + np.sum(casing_thicknesses)
                        + magnet["axial_thickness"] / 2
                    )
                )

                if direction == "outward":
                    self._outermost_coil_z0[key] = z0


                print(f"The z0 for the {direction} magnet on the {key} side is {z0}.")

                # Magnet and casing
                casing_regions = nested_cylindrical_shells(
                    z0=z0,
                    innermost_radius=r_inner,
                    inner_radial_thickness=magnet["radial_thickness"],
                    inner_axial_thickness=magnet["axial_thickness"],
                    layer_front_thickness=casing_thicknesses,
                    layer_back_thickness=casing_thicknesses,
                    layer_axial_thickness=casing_thicknesses
                )

                magnet_region = casing_regions[0]
                magnet_cell = openmc.Cell(
                    cell_id=6100 + i,
                    region=magnet_region,
                    fill=getattr(mat_ns, magnet["material"])
                )
                self._region_list.append(magnet_region)
                self._regions_by_region[key]["coil"].append(magnet_region)
                self._cells_by_region[key]["coil"].append(magnet_cell)
                add_cell_callback(magnet_cell)

                for j, region in enumerate(casing_regions[1:]):
                    casing_cell = openmc.Cell(
                        cell_id=6500 + i * 10 + j,
                        region=region,
                        fill=casing_materials[j]
                    )
                    self._region_list.append(region)
                    self._regions_by_region[key]["coil"].append(region)
                    self._cells_by_region[key]["coil"].append(casing_cell)
                    add_cell_callback(casing_cell)

                # Shield
                outer_r_thick = magnet["radial_thickness"] + 2 * np.sum(casing_thicknesses)
                outer_z_thick = magnet["axial_thickness"] + 2 * np.sum(casing_thicknesses)
                r_shield_base = bottleneck_r + shield["radial_gap_before_casing"]

                shield_region = hollow_cylinder_with_shell(
                    z0,
                    r_shield_base,
                    outer_r_thick,
                    outer_z_thick,
                    shield["radial_thickness"][0],
                    shield["radial_thickness"][1],
                    axial_thickness
                )[0]

                shield_cell = openmc.Cell(
                    cell_id=6200 + i,
                    region=shield_region,
                    fill=getattr(mat_ns, shield["material"])
                )
                self._region_list.append(shield_region)
                self._regions_by_region[key]["shield"].append(shield_region)
                self._cells_by_region[key]["shield"].append(shield_cell)
                add_cell_callback(shield_cell)

                i += 1

    def get_cells_by_region(self):
        return self._cells_by_region

    def get_regions_by_region(self):
        return self._regions_by_region
    
    def outermost_coil_z0(self):
        return self._outermost_coil_z0
    
    def get_full_region_list(self):
        return self._region_list
    
class EndCellBuilder:
    """
    Builds the left and right end cells that cap the HF coils in a tandem mirror system.
    """

    def __init__(self, end_cell_params):
        self.params = end_cell_params
        self._cells_by_side = {"left": [], "right": []}
        self._regions_by_side = {"left": {}, "right": {}}
        self._get_outer_limits = {"left": [], "right": []}
        self._region_list = []

    def build(self, hf_center_z0_dict, hf_coil_params_dict, vacuum_exclusion_region_dict, mat_ns, add_cell_callback):
        for key in ["left", "right"]:

            axial_sign = 1 if key == "right" else -1

            axial_length = self.params["axial_length"]
            shell_thickness = self.params["shell_thickness"]
            outer_radius = self.params["diameter"] / 2

            hf_center_z0 = hf_center_z0_dict[key]
            hf_magnet = hf_coil_params_dict[key]["magnet"]
            hf_shield = hf_coil_params_dict[key]["shield"]
            casing_thicknesses = np.array(hf_coil_params_dict[key]["casing_layers"]["thicknesses"])
            casing_sum = np.sum(casing_thicknesses)

            shield_toward = hf_shield["axial_thickness"][0]
            magnet_half = hf_magnet["axial_thickness"] / 2

            z0 = hf_center_z0 + axial_sign * (
                + shield_toward
                + casing_sum
                + magnet_half
                + shell_thickness
                + axial_length / 2
            )

            self._get_outer_limits[key] = z0 + axial_length / 2 + shell_thickness + 75


            shell_region, inner_region = cylinder_with_shell(
                z0,
                outer_radius - shell_thickness,
                axial_length,
                shell_thickness
            )

            shell_region &= ~vacuum_exclusion_region_dict[key]

            shell_cell = openmc.Cell(
                cell_id=5001 if key == "left" else 5003,
                region=shell_region,
                fill=getattr(mat_ns, self.params["shell_material"])
            )

            inner_cell = openmc.Cell(
                cell_id=5002 if key == "left" else 5004,
                region=inner_region,
                fill=getattr(mat_ns, self.params["inner_material"])
            )

            self._region_list.append(shell_region)
            self._region_list.append(inner_region)

            add_cell_callback(shell_cell)
            add_cell_callback(inner_cell)

            self._cells_by_side[key] = [shell_cell, inner_cell]
            self._regions_by_side[key] = {"shell": shell_region, "inner": inner_region}

    def get_cells_by_side(self):
        return self._cells_by_side

    def get_regions_by_side(self):
        return self._regions_by_side
    
    def get_outer_limits(self):
        return self._get_outer_limits
    
    def get_full_region_list(self):
        return self._region_list
            

class TandemMachineBuilder:
    """
    Assembles vacuum vessel, first wall, and central cylinder structures into a complete tandem mirror model.
    """

    def __init__(self, vv_params, fw_params, central_cyl_params, lf_coil_params, hf_coil_params, end_cell_params, material_ns):
        self.vv_params = vv_params
        self.fw_params = fw_params
        self.central_cyl_params = central_cyl_params
        self.lf_coil_params = lf_coil_params
        self.hf_coil_params = hf_coil_params
        self.end_cell_params = end_cell_params
        self.material_ns = material_ns

        self.vv_builder = None
        self.fw_builder = None
        self.ccyl_builder = None
        self.lf_coil_builder = None
        self.hf_coil_builder = None
        self.end_cell_builder = None

        self._universe = openmc.Universe(name="TandemMachine")
        self._all_cells = []
        self._z_extents = {}
        self._regions = {"vv": {}, "fw": {}, "cyl": {}, "lf_coils": {}, "hf_coils": {}, "end_cell": {}}

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

        self.fw_builder.build_all_sections(
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

        # ----- LF Coils ------
        
        self.lf_coil_builder = LFCoilBuilder(self.lf_coil_params)

        outer_radii = {key: self.ccyl_builder.get_radii_by_region()[key][-1] for key in ["left", "central", "right"]}

        self.lf_coil_builder.build_all_sections(
            outer_radii_by_region=outer_radii,
            midplanes=self.central_cyl_params["midplanes"],
            mat_ns=self.material_ns,
            add_cell_callback=self._add_cell
        )

        self._regions["lf_coils"] = self.lf_coil_builder.get_regions_by_region()

        # ----- 5. HF Coils -----
        
        self.hf_coil_builder = HFCoilBuilder(self.hf_coil_params)
        
        #print(self.fw_builder.fw_radii()["left"])
        
        bottleneck_radii = {
            "left": self.fw_builder.fw_radii()["left"][-1],
            "right": self.fw_builder.fw_radii()["right"][-1]
        }
        
        cc_half_lengths = {
            key: self.central_cyl_params["axial_lengths"][key] / 2
            for key in ["left", "right"]
        }
        self.hf_coil_builder.build_all_sections(
            bottleneck_radii=bottleneck_radii,
            cc_half_lengths=cc_half_lengths,
            midplanes={
                "left": self.central_cyl_params["midplanes"]["left"],
                "right": self.central_cyl_params["midplanes"]["right"]
            },
            mat_ns=self.material_ns,
            add_cell_callback=self._add_cell
        )
        self._regions["hf_coils"] = self.hf_coil_builder.get_regions_by_region()

    # ----- 6. End Cells
        
        self.end_cell_builder = EndCellBuilder(self.end_cell_params)

        hf_coils_params_dict = {
            key: {
                "magnet": self.hf_coil_params[key]["magnet"],
                "shield": self.hf_coil_params[key]["shield"],
                "casing_layers": self.hf_coil_params[key]["casing_layers"]
            } for key in ["left", "right"]
        }

        hf_center_z0_dict = self.hf_coil_builder.outermost_coil_z0()
        vacuum_outermost_regions_dict = {key: self.fw_builder.get_regions_by_region()[key][-1] for key in ["left", "right"]}

        self.end_cell_builder.build(
            hf_center_z0_dict=hf_center_z0_dict,
            hf_coil_params_dict=hf_coils_params_dict,
            vacuum_exclusion_region_dict=vacuum_outermost_regions_dict,
            mat_ns=self.material_ns,
            add_cell_callback=self._add_cell
        )

        self._regions["end_cell"] = self.end_cell_builder.get_regions_by_side()

    # ------ Add outside room ----------
        
    ################# Defining the room inside which the machine exists #####################
    

        z_min = self.end_cell_builder.get_outer_limits()["left"]
        z_max = self.end_cell_builder.get_outer_limits()["right"]
        bounding_surface = openmc.model.RectangularParallelepiped(xmin=-450, xmax=450, ymin=-450, ymax=450, zmin=-2000, zmax=2000, boundary_type = 'vacuum')
        room_region = -bounding_surface

        machine_compononet_region_list = self.vv_builder.get_full_region_list() + self.fw_builder.get_full_region_list() + self.ccyl_builder.get_full_region_list() + self.lf_coil_builder.get_full_region_list() + self.hf_coil_builder.get_full_region_list() + self.end_cell_builder.get_full_region_list()

        for i in range(len(machine_compononet_region_list)):
            room_region &= ~machine_compononet_region_list[i]

        room_cell = openmc.Cell(name="Room", region=room_region, fill=getattr(self.material_ns, "air"))
        self._add_cell(room_cell)

    # ---------- Accessors ----------
    def get_universe(self):
        return self._universe

    def get_all_cells(self):
        return self._all_cells

    def get_all_regions(self):
        return self._regions

    def get_z_extents(self):
        return self._z_extents
    
    def get_nwl(self):
        return 
