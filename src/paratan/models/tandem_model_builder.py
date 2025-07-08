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
            "left": (ep_vv["outer_axial_length"], ep_vv["central_axis_length"], end_axial, 9 * sep/10, left_midplane),
            "right": (ep_vv["outer_axial_length"], ep_vv["central_axis_length"], 9 * sep/10, end_axial, right_midplane),
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

    blanket_info_cc = input_data.get("central_cell", {}).get("blanket", {}).get("tallies", {})
    blanket_info_ep = input_data.get("end_plug", {}).get("central_cylinder", {}).get("tallies", {})

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
        },
        "tallies": {
            "central": {
                "breeder": blanket_info_cc.get("breeder", {}),
                "layers": {
                    entry["position"]: {
                        "cell_tallies": entry.get("cell_tallies", []),
                        "mesh_tallies": entry.get("mesh_tallies", [])
                    } for entry in blanket_info_cc.get("layer_tallies", [])
                }
            },
            "left": {
                "breeder": blanket_info_ep.get("breeder", {}),
                "layers": {
                    entry["position"]: {
                        "cell_tallies": entry.get("cell_tallies", []),
                        "mesh_tallies": entry.get("mesh_tallies", [])
                    } for entry in blanket_info_ep.get("layer_tallies", [])
                }
            },
            "right": {
                "breeder": blanket_info_ep.get("breeder", {}),
                "layers": {
                    entry["position"]: {
                        "cell_tallies": entry.get("cell_tallies", []),
                        "mesh_tallies": entry.get("mesh_tallies", [])
                    } for entry in blanket_info_ep.get("layer_tallies", [])
                }
            }
        }
    }

    # --------- 4. LF Coil Params ---------
    lf_coil_params = {
        "central": input_data.get("central_cell", {}).get("lf_coil", {}),
        "end": input_data.get("end_plug", {}).get("lf_coil", {}),
        "central_tallies": input_data.get("central_cell", {}).get("lf_coil", {}).get("lf_coil_tallies", {}),
        "end_tallies": input_data.get("end_plug", {}).get("lf_coil", {}).get("lf_coil_tallies", {})
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
            "shield": hf_coil_data.get(key, {}).get("shield", {}),
            "tallies": hf_coil_data.get("hf_coil_tallies", {})
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

        self._z_extents["central"] = (-self.central_cell_end_plug_separation_distance/10 -central_cell_vv_central_axis_length/2, central_cell_vv_central_axis_length/2 + self.central_cell_end_plug_separation_distance/10)
        self._z_extents["right"] = (
            right_midplane - end_plug_vv_central_axis_length/2 - 0.9 * self.central_cell_end_plug_separation_distance,
            right_midplane + end_plug_vv_central_axis_length/2 + self.end_axial_distance
        )
        self._z_extents["left"] = (
            left_midplane - end_plug_vv_central_axis_length/2 -  self.end_axial_distance,
            left_midplane + end_plug_vv_central_axis_length/2 + 0.9 * self.central_cell_end_plug_separation_distance
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
        

        cell_id = {"left": 1500, "central": 1501, "right": 1502}
        
        self._cells = {
            key: openmc.Cell(cell_id[key],name=f"{key}_vv_cell", region=region, fill=self.vacuum_material)
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
    
    def get_vv_components(self):
        return self._components
    

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
            
            if key == "left":
                k = 1
            elif key == "right":
                k = 2
            else:
                k = 3
            
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
                cell = openmc.Cell(1000 + k*100 + i,region=shell, fill=materials[i-1])
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

    def __init__(self, left_layers, central_layers, right_layers=None, tally_data=None):
        self.left_layers = left_layers
        self.central_layers = central_layers
        self.right_layers = right_layers if right_layers else left_layers

        self.tally_data = tally_data or {"left": {}, "central": {}, "right": {}}

        self._cells_by_region = {"left": [], "central": [], "right": []}
        self._regions_by_region = {"left": [], "central": [], "right": []}
        self._z_extents = {}
        self._region_list = []
        self._radii_by_region = {"left": [], "central": [], "right": []}
        self._tally_descriptors = []

    def build_all_sections(self, base_radii, midplanes, axial_lengths, base_regions, mat_ns, add_cell_callback):
        for key, layers in zip(["left", "central", "right"], [self.left_layers, self.central_layers, self.right_layers]):
            z0 = midplanes[key]
            half_len = axial_lengths[key] / 2
            zmin, zmax = z0 - half_len, z0 + half_len
            self._z_extents[key] = (zmin, zmax)

            thicknesses = [layer["thickness"] for layer in layers]
            materials = [getattr(mat_ns, layer["material"]) for layer in layers]
            radii = base_radii[key] + np.cumsum(thicknesses)
            self._radii_by_region[key] = radii

            zplanes = -openmc.ZPlane(z0=zmax) & +openmc.ZPlane(z0=zmin)
            prev_cyl = openmc.ZCylinder(r=radii[0])

            for i in range(len(radii)):
                if i == 0:
                    region = zplanes & ~base_regions[key] & -prev_cyl
                else:
                    next_cyl = openmc.ZCylinder(r=radii[i])
                    region = zplanes & -next_cyl & +prev_cyl
                    prev_cyl = next_cyl

                cell = openmc.Cell(region=region, fill=materials[i])
                add_cell_callback(cell)

                self._regions_by_region[key].append(region)
                self._region_list.append(region)
                self._cells_by_region[key].append(cell)

                desc = "breeder" if i == 0 else f"layer_{i}"
                tally = self.tally_data.get(key, {})
                layer_tallies = tally.get("breeder", {}) if i == 0 else tally.get("layers", {}).get(i, {})

                self._tally_descriptors.append({
                    "type": "blanket",
                    "location": key,
                    "description": desc,
                    "cell": cell,
                    "region": region,
                    "cell_tallies": layer_tallies.get("cell_tallies", []),
                    "mesh_tallies": layer_tallies.get("mesh_tallies", [])
                })

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

    def get_tally_descriptors(self):
        return self._tally_descriptors
    
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
        
        self.tally_data = {
            "central": lf_coil_params.get("central_tallies", {}),
            "end": lf_coil_params.get("end_tallies", {})
        }

        self._cells_by_region = {region: {"coil": [], "shield": []} for region in ["left", "central", "right"]}
        self._regions_by_region = {region: {"coil": [], "shield": []} for region in ["left", "central", "right"]}
        self._region_list = []
        self._tally_descriptors = []

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
            tally_data = self.tally_data.get(data_key, {})

            positions = coil_data.get("positions", [])
            inner_dims = coil_data.get("inner_dimensions", {})
            shell_dims = coil_data.get("shell_thicknesses", {})
            materials = coil_data.get("materials", {})

            coil_regions = []
            shield_regions = []
            coil_cells = []
            shield_cells = []

            pos_i = 1
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
                
                # Add tally descriptor for coil (not shield)
                self._tally_descriptors.append({
                    "type": "lf_coil",
                    "location": key,
                    "description": pos_i,
                    "cell": coil_cell,
                    "region": coil_region,
                    "cell_tallies": tally_data.get("cell_tallies", []),
                    "mesh_tallies": tally_data.get("mesh_tallies", [])
                })
                
                i += 1

    def get_cells_by_region(self):
        return self._cells_by_region

    def get_regions_by_region(self):
        return self._regions_by_region
    
    def get_full_region_list(self):
        return self._region_list
    
    def get_tally_descriptors(self):
        return self._tally_descriptors

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
        self._tally_descriptors = []

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
                
                self._tally_descriptors.append({
                    "type": "hf_coil",
                    "location": key,
                    "description": direction,
                    "cell": magnet_cell,
                    "region": magnet_region,
                    "cell_tallies": self.coil_data.get(key, {}).get("tallies", {}).get("cell_tallies", []),
                    "mesh_tallies": self.coil_data.get(key, {}).get("tallies", {}).get("mesh_tallies", [])
                })

                for j, region in enumerate(casing_regions[1:-1]):
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
    
    def get_all_coil_cells(self):
        return [cell for side in self._cells_by_region.values() for cell in side["coil"]]

    def get_tally_descriptors(self):
        return self._tally_descriptors
    
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
            inner_region &= ~vacuum_exclusion_region_dict[key]

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
    Modular builder for tandem mirror fusion machines.
    Components: VV, FW, Cylinders, LF/HF coils, End cells, Room, Tallies.
    """

    def __init__(self, vv_params, fw_params, central_cyl_params, lf_coil_params, hf_coil_params, end_cell_params, material_ns):
        self.vv_params = vv_params
        self.fw_params = fw_params
        self.central_cyl_params = central_cyl_params
        self.lf_coil_params = lf_coil_params
        self.hf_coil_params = hf_coil_params
        self.end_cell_params = end_cell_params
        self.material_ns = material_ns

        self._bounding_surface = openmc.model.RectangularParallelepiped(
            xmin=-450, xmax=450, ymin=-450, ymax=450,
            zmin=-5000, zmax=5000, boundary_type='vacuum'
        )
        self._room_region = -self._bounding_surface  # Start with full volume

        self.vv_builder = None
        self.fw_builder = None
        self.ccyl_builder = None
        self.lf_coil_builder = None
        self.hf_coil_builder = None
        self.end_cell_builder = None
        self.tally_builder = TallyBuilder()

        self._universe = openmc.Universe(name="TandemMachine")
        self._all_cells = []
        self._z_extents = {}
        self._regions = {"vv": {}, "fw": {}, "cyl": {}, "lf_coils": {}, "hf_coils": {}, "end_cell": {}}

        room_cell = openmc.Cell(69, name="Room", region=self._room_region, fill=getattr(self.material_ns, "air"))
        self._universe.add_cell(room_cell)

    def _add_cell(self, cell):
        self._universe.add_cell(cell)
        self._all_cells.append(cell)

    def subtract_from_room(self, regions):
        for region in regions:
            self._room_region &= ~region

    def build_vacuum_vessel(self):
        self.vv_builder = TandemVacuumVessel(**self.vv_params)
        vv_cells = self.vv_builder.build_cells()
        for cell in vv_cells.values():
            self._add_cell(cell)
        self._regions["vv"] = self.vv_builder.get_regions()
        self._z_extents["vv"] = self.vv_builder.get_z_extents()
        self.subtract_from_room(self.vv_builder.get_full_region_list())

    def build_first_wall(self):
        self.fw_builder = TandemFWStructure(
            self.fw_params["left_fw"], self.fw_params["central_fw"], self.fw_params.get("right_fw", None)
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
        self.subtract_from_room(self.fw_builder.get_full_region_list())

    def build_central_cylinders(self):
        combined_fw = self._regions["fw"]["left"][-1] | self._regions["fw"]["central"][-1] | self._regions["fw"]["right"][-1]
        base_regions = {key: combined_fw for key in ["left", "central", "right"]}
        base_radii = {
            key: self.fw_params["base_radii"][key] + sum(layer["thickness"] for layer in layers)
            for key, layers in zip(["left", "central", "right"],
                                   [self.fw_params["left_fw"], self.fw_params["central_fw"],
                                    self.fw_params.get("right_fw", self.fw_params["left_fw"])])
        }
        self.ccyl_builder = TandemCentralCylinders(
            self.central_cyl_params["left"],
            self.central_cyl_params["central"],
            self.central_cyl_params.get("right", None),
            self.central_cyl_params.get("tallies", None)
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
        self.tally_builder.add_descriptors(self.ccyl_builder.get_tally_descriptors())
        self.subtract_from_room(self.ccyl_builder.get_full_region_list())

    def build_lf_coils(self):
        self.lf_coil_builder = LFCoilBuilder(self.lf_coil_params)
        outer_radii = {key: self.ccyl_builder.get_radii_by_region()[key][-1] for key in ["left", "central", "right"]}
        self.lf_coil_builder.build_all_sections(
            outer_radii_by_region=outer_radii,
            midplanes=self.central_cyl_params["midplanes"],
            mat_ns=self.material_ns,
            add_cell_callback=self._add_cell
        )
        self._regions["lf_coils"] = self.lf_coil_builder.get_regions_by_region()
        self.tally_builder.add_descriptors(self.lf_coil_builder.get_tally_descriptors())
        self.subtract_from_room(self.lf_coil_builder.get_full_region_list())

    def build_hf_coils(self):
        self.hf_coil_builder = HFCoilBuilder(self.hf_coil_params)
        bottleneck_radii = {
            "left": self.fw_builder.fw_radii()["left"][-1],
            "right": self.fw_builder.fw_radii()["right"][-1]
        }
        cc_half_lengths = {key: self.central_cyl_params["axial_lengths"][key] / 2 for key in ["left", "right"]}
        self.hf_coil_builder.build_all_sections(
            bottleneck_radii=bottleneck_radii,
            cc_half_lengths=cc_half_lengths,
            midplanes=self.central_cyl_params["midplanes"],
            mat_ns=self.material_ns,
            add_cell_callback=self._add_cell
        )
        self._regions["hf_coils"] = self.hf_coil_builder.get_regions_by_region()
        self.tally_builder.add_descriptors(self.hf_coil_builder.get_tally_descriptors())
        self.subtract_from_room(self.hf_coil_builder.get_full_region_list())

    def build_end_cells(self):
        self.end_cell_builder = EndCellBuilder(self.end_cell_params)
        hf_coils_params_dict = {
            key: {
                "magnet": self.hf_coil_params[key]["magnet"],
                "shield": self.hf_coil_params[key]["shield"],
                "casing_layers": self.hf_coil_params[key]["casing_layers"]
            } for key in ["left", "right"]
        }
        hf_center_z0_dict = self.hf_coil_builder.outermost_coil_z0()
        vacuum_outermost_regions_dict = {
            key: self._regions["fw"][key][-1] | self._regions["fw"][key][-2] | self._regions["fw"]["central"][-1]
            for key in ["left", "right"]
        }
        self.end_cell_builder.build(
            hf_center_z0_dict=hf_center_z0_dict,
            hf_coil_params_dict=hf_coils_params_dict,
            vacuum_exclusion_region_dict=vacuum_outermost_regions_dict,
            mat_ns=self.material_ns,
            add_cell_callback=self._add_cell
        )
        self._regions["end_cell"] = self.end_cell_builder.get_regions_by_side()
        self.subtract_from_room(self.end_cell_builder.get_full_region_list())

    def get_universe(self):
        return self._universe

    def get_all_cells(self):
        return self._all_cells

    def get_all_regions(self):
        return self._regions

    def get_z_extents(self):
        return self._z_extents

    def get_all_tallies(self):
        return self.tally_builder.get_tallies()
    
    def get_neutron_wall_loading_model(self):
        
        
        Path("nwl").mkdir(parents=True, exist_ok=True)
        
        tungsten = openmc.Material(31, name="tungsten")
        tungsten.set_density("g/cm3", 19.25)
        tungsten.add_element("W", 100)

        my_materials = openmc.Materials([tungsten])


        my_materials.export_to_xml(path='nwl/materials.xml') # write the materials.xml file

        
        vv_geometric_extents = {}
        
        vv_geometric_extents["central"] = self.vv_builder.get_vv_components()["central"]
        
        left_end_vv = vv_geometric_extents["central"]['z_planes']['left_end']
        right_end_vv = vv_geometric_extents["central"]['z_planes']['right_end']
        central_radius_vv = vv_geometric_extents["central"]['central_cylinder_radius']
        bottleneck_radius_vv = vv_geometric_extents["central"]['bottleneck_radius']
        
        
        # surfaces
        model_right=openmc.ZPlane(z0=right_end_vv.z0 + 10, boundary_type='vacuum') 
        model_left=openmc.ZPlane(z0= left_end_vv.z0 - 10, boundary_type='vacuum') 
        model_side_outer_radius=openmc.ZCylinder(r=central_radius_vv + 25, boundary_type='vacuum', surface_id = 123) 
        
        #
        centralfwarmor_right=vv_geometric_extents["central"]['z_planes']['central_right']
        centralfwarmor_left=vv_geometric_extents["central"]['z_planes']['central_left']
        #
        rightendfwarmor_right=right_end_vv
        rightendfwarmor_left=vv_geometric_extents["central"]['z_planes']['right_cone_plane']
        
        #
        leftendfwarmor_right=vv_geometric_extents["central"]['z_planes']['left_cone_plane']
        leftendfwarmor_left=left_end_vv
        leftendfwarmor_left.boundary_type = 'transmission'
        leftendfwarmor_right.boundary_type = 'transmission'

        rightendfwarmor_right.boundary_type = 'transmission'

        
        #
        centralplasma_outer_radius=openmc.ZCylinder(r=central_radius_vv)
        centralfwarmor_outer_radius=openmc.ZCylinder(r=central_radius_vv+0.2)
        
        #
        endplasma_outer_radius=openmc.ZCylinder(r=bottleneck_radius_vv, surface_id = 2012)
        endfwarmor_outer_radius=openmc.ZCylinder(r=bottleneck_radius_vv+0.2)
        
        #
        leftconeplasma_outer_radius = vv_geometric_extents["central"]['surfaces']['cone_left']
        rightconeplasma_outer_radius = vv_geometric_extents["central"]['surfaces']['cone_right']
        
        centralplasma_outer_radius.id = 2010
        endplasma_outer_radius.id = 2012
        leftconeplasma_outer_radius.cone.id = 2014
        rightconeplasma_outer_radius.cone.id = 2018
        #
        leftconefwarmor_outer_radius = vv_geometric_extents["central"]['surfaces']['cone_fw_left']
        rightconefwarmor_outer_radius = vv_geometric_extents["central"]['surfaces']['cone_fw_right']
        
        # regions
        # plasma regions
        leftendplasma_region=-endplasma_outer_radius & +leftendfwarmor_left & -leftendfwarmor_right
        leftconeplasma_region=-leftconeplasma_outer_radius & -centralfwarmor_left & +leftendfwarmor_right
        centralplasma_region=-centralplasma_outer_radius & -centralfwarmor_right & +centralfwarmor_left
        rightconeplasma_region=-rightconeplasma_outer_radius & +centralfwarmor_right & -rightendfwarmor_left
        rightendplasma_region=-endplasma_outer_radius & +rightendfwarmor_left & -rightendfwarmor_right
        # # fwarmor regions
        leftendfwarmor_region=+endplasma_outer_radius & -endfwarmor_outer_radius & +leftendfwarmor_left & -leftendfwarmor_right
        lowerconefwarmor_region=+leftconeplasma_outer_radius & -leftconefwarmor_outer_radius & -centralfwarmor_left & +leftendfwarmor_right
        centralfwarmor_region=+centralplasma_outer_radius & -centralfwarmor_outer_radius & -centralfwarmor_right & +centralfwarmor_left
        upperconefwarmor_region=+rightconeplasma_outer_radius & -rightconefwarmor_outer_radius & +centralfwarmor_right & -rightendfwarmor_left
        rightendfwarmor_region=+endplasma_outer_radius & -endfwarmor_outer_radius & +rightendfwarmor_left & -rightendfwarmor_right
        # # outside regions
        leftdisk_region=+model_left & -leftendfwarmor_left & -model_side_outer_radius
        leftendoutside_region=+endfwarmor_outer_radius & -model_side_outer_radius & +leftendfwarmor_left & -leftendfwarmor_right
        leftconeoutside_region=+leftconefwarmor_outer_radius & -model_side_outer_radius & -centralfwarmor_left & +leftendfwarmor_right
        centraloutside_region=+centralfwarmor_outer_radius & -model_side_outer_radius & -centralfwarmor_right & +centralfwarmor_left
        rightconeoutside_region=+rightconefwarmor_outer_radius & -model_side_outer_radius & +centralfwarmor_right & -rightendfwarmor_left
        rightendoutside_region=+endfwarmor_outer_radius & -model_side_outer_radius & +rightendfwarmor_left & -rightendfwarmor_right
        rightdisk_region=-model_right & +rightendfwarmor_right & -model_side_outer_radius
        
        print(centralfwarmor_left.z0, centralfwarmor_right.z0)
        
        # create cell and assign material
        leftendplasma_cell=openmc.Cell(101, region=leftendplasma_region,fill=None)
        leftconeplasma_cell=openmc.Cell(102, region=leftconeplasma_region,fill=None)
        centralplasma_cell=openmc.Cell(103, region=centralplasma_region,fill=None)
        rightconeplasma_cell=openmc.Cell(104, region=rightconeplasma_region,fill=None)
        rightendplasma_cell=openmc.Cell(105, region=rightendplasma_region,fill=None)
        #
        leftendfwarmor_cell=openmc.Cell(201, region=leftendfwarmor_region,fill=None)
        leftconefwarmor_cell=openmc.Cell(202, region=lowerconefwarmor_region,fill=None)
        centralfwarmor_cell=openmc.Cell(203, region=centralfwarmor_region,fill=None)
        rightconefwarmor_cell=openmc.Cell(204, region=upperconefwarmor_region,fill=None)
        rightendfwarmor_cell=openmc.Cell(205, region=rightendfwarmor_region,fill=None)
        
        # #
        leftdisk_cell=openmc.Cell(50, region=leftdisk_region,fill=None)
        leftendoutside_cell=openmc.Cell(301, region=leftendoutside_region,fill=None)
        leftconeoutside_cell=openmc.Cell(302, region=leftconeoutside_region,fill=None)
        centraloutside_cell=openmc.Cell(303, region=centraloutside_region,fill=None)
        rightconeoutside_cell=openmc.Cell(304, region=rightconeoutside_region,fill=None)
        rightendoutside_cell=openmc.Cell(305, region=rightendoutside_region,fill=None)
        rightdisk_cell=openmc.Cell(51, region=rightdisk_region,fill=None)
        #
        # make a universe to contain all the cells
        geometry = openmc.Geometry([leftconeplasma_cell, centralplasma_cell, rightconeplasma_cell, leftendplasma_cell, rightendplasma_cell,
                                    leftconefwarmor_cell, centralfwarmor_cell, rightconefwarmor_cell, rightendfwarmor_cell, leftendfwarmor_cell,
                                    leftconeoutside_cell,centraloutside_cell,rightconeoutside_cell, leftdisk_cell, rightdisk_cell, leftendoutside_cell, rightendoutside_cell])
        
        # geometry = openmc.Geometry([centralplasma_cell])

        # Get all surfaces from the geometry
        all_surfaces = geometry.get_all_surfaces() 
        print("\n📋 Surfaces with boundary_type='vacuum':")
        for sid, surface in all_surfaces.items():
            if surface.boundary_type == 'vacuum':
                print(f"Surface ID {sid}: Type {type(surface).__name__}, z0/r = {getattr(surface, 'z0', getattr(surface, 'r', 'N/A'))}")
        #
        geometry.export_to_xml(path='nwl/geometry.xml') # writes the geometry.xml file.
        
        
        #
        # plot section
        #
        # generate some plot slices
        plot_xz_overall= geometry.plot(basis='xz', origin=(0,0,0), width=(1200,1200), pixels=(1000,1000))
        plot_xz_overall.figure.savefig('nwl/xz-cell_overall' + 'central' + '.png')
        #
        plot_xz_lower = geometry.plot(basis='xz', origin=(0,0,-300), width=(100,100), pixels=(1000,1000))
        plot_xz_lower.figure.savefig('nwl/xz-cell_lower'+ 'central' + '.png')
        #
        plot_xz_upper = geometry.plot(basis='xz', origin=(0,0,300), width=(100,100), pixels=(1000,1000))
        plot_xz_upper.figure.savefig('nwl/xz-cell_upper'+ 'central'+ '.png')
        
        #
        #
        settings = openmc.Settings()
        settings.run_mode = 'fixed source'
        settings.batches = 10
        settings.particles = 100000

        # Cylindrical source: radius ~ 0, z from -125 to 125 cm
        source_space = openmc.stats.CylindricalIndependent(
            r=openmc.stats.Discrete([0.0], [1.0]),                  # point source at r = 0
            phi=openmc.stats.Uniform(0.0, 2*np.pi),                 # angle uniform (won't matter at r=0)
            z=openmc.stats.Uniform(-230,100)                   
        )

        source_energy = openmc.stats.Discrete([14.1e6], [1.0])
        
        pointsource = openmc.IndependentSource()

        settings.source = [openmc.IndependentSource(space=source_space, energy=source_energy)]
        
        #settings.source = [source_space]

        settings.export_to_xml(path='nwl/settings.xml')
        #
        settings.surf_source_write = {
        'surface_ids': [2012, 2014,2018, 2010],
        'max_particles': 1000000
        } # test writing lower end (actually writes lower and upper ends), lower cone, central, upper cone 
        
        # tallies section
        cell_filter = openmc.CellFilter([102, 103, 104])
        mytally = openmc.Tally(14) # set tally number to 14 (arbitrary)
        mytally.filters = [cell_filter]
        mytally.scores = ['flux']
        #
        tallies = openmc.Tallies([mytally])
        tallies.export_to_xml(path='nwl/tallies.xml')

        #
        # now run
        nwl_model = openmc.model.Model(geometry=geometry, materials=my_materials, settings=settings, tallies=tallies)

        # Creating a dictionary of relevant FW surface z-positions
        fwarmor_zplanes = {
            "leftendfwarmor_left_z0":     leftendfwarmor_left.z0,
            "leftendfwarmor_right_z0":    leftendfwarmor_right.z0,
            "centralfwarmor_left_z0":     centralfwarmor_left.z0,
            "centralfwarmor_right_z0":    centralfwarmor_right.z0,
            "rightendfwarmor_left_z0":    rightendfwarmor_left.z0,
            "rightendfwarmor_right_z0":   rightendfwarmor_right.z0
        }
                
        return nwl_model, fwarmor_zplanes
    

# def build_tandem_model_from_input(input_data):
#     # Parse geometry params
#     vv_params, fw_params, cyl_params, lf_coil_params, hf_coil_params, end_cell_params = parse_machine_input(input_data, m)

#     # Build machine
#     builder = TandemMachineBuilder(
#         vv_params, fw_params, cyl_params,
#         lf_coil_params, hf_coil_params, end_cell_params,
#         material_ns=m
#     )
#     builder.build_vacuum_vessel()
#     builder.build_first_wall()
#     builder.build_central_cylinders()
#     builder.build_lf_coils()
#     builder.build_hf_coils()
#     builder.build_end_cells()

#     # Geometry and universe
#     universe_machine = builder.get_universe()
#     geometry = openmc.Geometry([openmc.Cell(fill=universe_machine)])
#     geometry.merge_surfaces = True

#     # Optional: Plot geometry
#     geometry.root_universe.plot(
#         basis='xz',
#         width=(1000, 5100),
#         pixels=(1000, 1200),
#         color_by='material'
#     )
#     plt.savefig('modular_tandem_mirror_cross_section.png', bbox_inches="tight")

#     # Tallies
#     tallies = openmc.Tallies(builder.get_all_tallies())

#     # Materials
#     materials = m.materials

#     # Load source config
#     with open('input_files/source_information.yaml', 'r') as f:
#         source_data = yaml.safe_load(f)

#     # Source
#     source = openmc.Source()
#     source.space = openmc.stats.CylindricalIndependent(
#         r=openmc.stats.Discrete([0], [1.0]),
#         phi=openmc.stats.Uniform(a=0.0, b=2 * np.pi),
#         z=openmc.stats.Uniform(a=-600.0, b=600.0),
#     )
#     source.angle = openmc.stats.Isotropic()
#     source.energy = openmc.stats.Discrete([14e6], [1.0])

#     # Settings
#     settings = openmc.Settings()
#     settings.run_mode = "fixed source"
#     settings.particles = int(source_data['settings']['particles_per_batch'])
#     settings.batches = source_data['settings']['batches']
#     settings.output = {'tallies': False}

#     freq = source_data['settings']['statepoint_frequency']
#     settings.statepoint = {
#         'batches': [1] + list(range(freq, settings.batches, freq)) + [settings.batches]
#     }

#     settings.weight_windows_on = source_data['settings']['weight_windows']
#     settings.weight_window_checkpoints = {'collision': True, 'surface': True}

#     wwg = openmc.WeightWindowGenerator(
#         openmc.RegularMesh(),
#         energy_bounds=[0, 14e6],
#         particle_type='neutron',
#         method='magic',
#         max_realizations=25,
#         update_interval=1,
#         on_the_fly=True
#     )
#     settings.weight_windows_generator = [wwg]

#     settings.photon_transport = source_data['settings']['photon_transport']
#     settings.source = source

#     # Final model
#     model = openmc.Model(geometry, materials, settings, tallies)
#     model.export_to_xml('model_xml_files')

#     return model

import os
import contextlib

@contextlib.contextmanager
def change_dir(destination):
    current = os.getcwd()
    os.makedirs(destination, exist_ok=True)
    os.chdir(destination)
    try:
        yield
    finally:
        os.chdir(current)

def build_tandem_model_from_input(input_data, output_dir="."):
    with change_dir(output_dir):
        # Parse geometry params
        vv_params, fw_params, cyl_params, lf_coil_params, hf_coil_params, end_cell_params = parse_machine_input(input_data, m)

        # Build machine
        builder = TandemMachineBuilder(
            vv_params, fw_params, cyl_params,
            lf_coil_params, hf_coil_params, end_cell_params,
            material_ns=m
        )
        builder.build_vacuum_vessel()
        builder.build_first_wall()
        builder.build_central_cylinders()
        builder.build_lf_coils()
        builder.build_hf_coils()
        # builder.build_end_cells()

        # Geometry and universe
        universe_machine = builder.get_universe()
        geometry = openmc.Geometry([openmc.Cell(fill=universe_machine)])
        geometry.merge_surfaces = True

        # Optional: Plot geometry
        geometry.root_universe.plot(
            basis='xz',
            width=(1000, 5100),
            pixels=(3600, 4200),
            color_by='material'
        )
        plt.savefig('hf_coils_modular_tandem_mirror_cross_section.png', bbox_inches="tight")

        # Tallies
        tallies = openmc.Tallies(builder.get_all_tallies())

        # Materials
        materials = m.materials

        # Load source config
        with open('source_information.yaml', 'r') as f:
            source_data = yaml.safe_load(f)

        source_z_start = -source_data["source"]["uniform"]["length"]/2
        source_z_end = source_data["source"]["uniform"]["length"]/2

        source_r = source_data["source"]["uniform"]["radius"]

        # Source
        source = openmc.Source()
        source.space = openmc.stats.CylindricalIndependent(
            r=openmc.stats.PowerLaw(a = 0.0, b = source_r, n = 2),
            phi=openmc.stats.Uniform(a=0.0, b=2 * np.pi),
            z=openmc.stats.Uniform(a=source_z_start, b=source_z_end),
        )
        source.angle = openmc.stats.Isotropic()
        source.energy = openmc.stats.Discrete([14e6], [1.0])

        # Settings
        settings = openmc.Settings()
        settings.run_mode = "fixed source"
        settings.particles = int(source_data['settings']['particles_per_batch'])
        settings.batches = source_data['settings']['batches']
        settings.output = {'tallies': False}

        freq = source_data['settings']['statepoint_frequency']
        settings.statepoint = {
            'batches': [1] + list(range(freq, settings.batches, freq)) + [settings.batches]
        }

        settings.weight_windows_on = source_data['settings']['weight_windows']
        settings.weight_window_checkpoints = {'collision': True, 'surface': True}

        wwg = openmc.WeightWindowGenerator(
            openmc.RegularMesh(),
            energy_bounds=[0, 14e6],
            particle_type='neutron',
            method='magic',
            max_realizations=25,
            update_interval=1,
            on_the_fly=True
        )
        settings.weight_windows_generator = [wwg]

        settings.photon_transport = source_data['settings']['photon_transport']
        settings.source = source

        # Final model
        model = openmc.Model(geometry, materials, settings, tallies)
        model.export_to_xml('model_xml_files')

    return model

