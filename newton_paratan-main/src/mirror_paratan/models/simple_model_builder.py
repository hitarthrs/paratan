import openmc
import numpy as np
import hashlib
import shutil
import xml.etree.ElementTree as ET

# Materials library that is used to build the model. To change materials, use --materials flag in CLI
# Materials will be loaded dynamically in the build function

from mirror_paratan.geometry.core import *
from mirror_paratan.tallies.base_tallies import TallyBuilder
import os
import contextlib
from pathlib import Path


def load_materials():
    """Load materials module based on environment variable or default."""
    import os
    if 'MATERIALS_PATH' in os.environ:
        # Import custom materials module specified by environment variable
        import importlib.util
        import sys
        materials_path = os.environ['MATERIALS_PATH']
        spec = importlib.util.spec_from_file_location("custom_materials", materials_path)
        m = importlib.util.module_from_spec(spec)
        sys.modules["custom_materials"] = m
        spec.loader.exec_module(m)
        return m
    else:
        # Use default materials
        from mirror_paratan.materials import material as m
        return m


def parse_simple_machine_input(input_data, material_ns):
    """
    Parse input data for simple mirror machine and return structured parameters.
    """
    # Vacuum vessel parameters
    vv_info = input_data["vacuum_vessel"]
    vv_params = {
        "outer_axial_length": vv_info["outer_axial_length"],
        "central_axial_length": vv_info["central_axial_length"],
        "central_radius": vv_info["central_radius"],
        "bottleneck_radius": vv_info["bottleneck_radius"],
        "left_bottleneck_length": vv_info["left_bottleneck_length"],
        "right_bottleneck_length": vv_info["right_bottleneck_length"],
        "axial_midplane": vv_info.get("axial_midplane", 0.0),
        "structure": vv_info.get("structure", {}),
        "vacuum_material": getattr(material_ns, "vacuum"),
        "geometry_style": vv_info.get("geometry_style", "axisymmetric"),
    }

    # Central cell parameters
    cc_info = input_data["central_cell"]
    cc_params = {
        "axial_length": cc_info["axial_length"],
        "layers": cc_info["layers"],
        "materials": [getattr(material_ns, layer["material"]) for layer in cc_info["layers"]],
        "tallies": cc_info.get("tallies", {})  # Updated to use nested tallies
    }

    # LF coil parameters
    lf_info = input_data["lf_coil"]
    lf_params = {
        "shell_thicknesses": lf_info["shell_thicknesses"],
        "inner_dimensions": lf_info["inner_dimensions"],
        "positions": lf_info["positions"],
        "materials": {
            "shield": getattr(material_ns, lf_info["materials"]["shield"]),
            "magnet": getattr(material_ns, lf_info["materials"]["magnet"])
        },
        "tallies": lf_info.get("lf_coil_tallies", {})  # Updated to use nested tallies
    }

    # HF coil parameters
    hf_info = input_data["hf_coil"]
    hf_params = {
        "magnet": hf_info["magnet"],
        "casing_layers": hf_info["casing_layers"],
        "casing_materials": [getattr(material_ns, layer["material"]) for layer in hf_info["casing_layers"]],
        "shield": hf_info["shield"],
        "shield_material": getattr(material_ns, hf_info["shield"]["material"]),
        "axial_offset_cm": hf_info.get(
            "axial_offset_cm",
            hf_info.get(
                "z0_offset_cm",
                hf_info.get(
                    "center_offset_cm",
                    hf_info.get("outward_shift_cm", 0.0),
                ),
            ),
        ),
        "tallies": hf_info.get("hf_coil_tallies", {})  # Updated to use nested tallies
    }

    # Share one pure-data layout derivation with the source-model geometry
    from source_model_revamp.geometry.device_domains import derive_paratan_layout
    source_layout = derive_paratan_layout(input_data)

    # End cell parameters
    end_info = input_data.get("end_cell", {})
    end_params = {
        "axial_length": end_info["axial_length"],
        "shell_thickness": end_info["shell_thickness"],
        "diameter": end_info["diameter"],
        "hf_to_end_cell_gap_cm": end_info.get("hf_to_end_cell_gap_cm", 10.0),
        "derived_center_offset_cm": 100.0 * (source_layout.right_end_cell_center_m - source_layout.midplane_z_m),
        "derived_full_z_min_cm": 100.0 * source_layout.full_device_domain.z_min_m,
        "derived_full_z_max_cm": 100.0 * source_layout.full_device_domain.z_max_m,
        "shell_material": getattr(material_ns, end_info.get("shell_material", "stainless")),
        "inner_material": getattr(material_ns, end_info.get("inner_material", "vacuum")),
        "tallies": end_info.get("end_cell_tallies", {})  # Updated to use nested tallies
    }

    return vv_params, cc_params, lf_params, hf_params, end_params


class SimpleVacuumVessel:
    """Builds the vacuum vessel structure for simple mirror machine."""
    
    def __init__(self, vv_params, material_ns):
        self.vv_params = vv_params
        self.material_ns = material_ns
        self.vacuum_section_regions = []
        self.vacuum_section_cells = []
        self.room_region = None
        
    def build_vacuum_vessel(self, room_region):
        """Build the vacuum vessel (axisymmetric hourglass or perpendicular straight section)."""
        self.room_region = room_region
        
        # Create main vacuum vessel region (see vacuum_vessel.geometry_style in YAML)
        vv_main_region, _vv_components = simple_mirror_vacuum_vessel_layer_region(
            self.vv_params["geometry_style"],
            self.vv_params["outer_axial_length"],
            self.vv_params["central_axial_length"],
            self.vv_params["central_radius"],
            self.vv_params["bottleneck_radius"],
            self.vv_params["left_bottleneck_length"],
            self.vv_params["right_bottleneck_length"],
            self.vv_params["axial_midplane"],
        )
        
        # Create main vacuum vessel cell
        vv_main_cell = openmc.Cell(
            1000, 
            region=vv_main_region, 
            fill=self.vv_params["vacuum_material"]
        )
        
        self.vacuum_section_regions.append(vv_main_region)
        self.vacuum_section_cells.append(vv_main_cell)
        
        # Build structural layers as nested composite shapes (following run_simple_mirror_model logic)
        if "structure" in self.vv_params:
            structural_thicknesses = [properties.get("thickness") for layer_name, properties in self.vv_params["structure"].items()]
            structural_materials = [getattr(self.material_ns, properties.get("material")) for layer_name, properties in self.vv_params["structure"].items()]
            
            # Calculate cumulative radii for BOTH central and bottleneck sections
            # This follows the exact logic from run_simple_mirror_model
            central_radii = np.cumsum([self.vv_params["central_radius"]] + structural_thicknesses)
            bottleneck_radii = np.cumsum([self.vv_params["bottleneck_radius"]] + structural_thicknesses)
            
            # Build structural layers as nested composite shapes
            for i, (thickness, material) in enumerate(zip(structural_thicknesses, structural_materials)):
                # Create the full vacuum vessel region at this structural layer's outer radius
                # Use BOTH updated radii: central_radii[i+1] and bottleneck_radii[i+1]
                vv_layer_region, _ = simple_mirror_vacuum_vessel_layer_region(
                    self.vv_params["geometry_style"],
                    self.vv_params["outer_axial_length"],
                    self.vv_params["central_axial_length"],
                    central_radii[i + 1],
                    bottleneck_radii[i + 1],
                    self.vv_params["left_bottleneck_length"],
                    self.vv_params["right_bottleneck_length"],
                    self.vv_params["axial_midplane"],
                )
                
                # The structural layer region is the full region minus the previous layer's region
                # This creates the cylindrical shell that follows the complex vacuum vessel geometry
                if i == 0:
                    # First structural layer: full region minus main vacuum vessel region
                    structural_region = vv_layer_region & ~self.vacuum_section_regions[0]
                else:
                    # Subsequent layers: full region minus previous structural layer region
                    structural_region = vv_layer_region & ~self.vacuum_section_regions[-1]
                
                # Create the structural layer cell
                structural_cell = openmc.Cell(
                    1000 + i + 1,
                    region=structural_region,
                    fill=material
                )
                
                self.vacuum_section_regions.append(vv_layer_region)  # Store the full region
                self.vacuum_section_cells.append(structural_cell)
        
        return self.vacuum_section_cells, self.vacuum_section_regions
    
    def get_outermost_radius(self):
        """Get the outermost radius of the vacuum vessel."""
        base_radius = self.vv_params["central_radius"]
        if "structure" in self.vv_params:
            structural_thickness = sum(
                [properties.get("thickness", 0) for layer_name, properties in self.vv_params["structure"].items()]
            )
            return base_radius + structural_thickness
        return base_radius
    def get_outermost_bottleneck_radius(self):
        """Get the outermost bottleneck radius of the vacuum vessel."""
        bottleneck_radius = self.vv_params["bottleneck_radius"]
        if "structure" in self.vv_params:
            structural_thickness = sum(
                [properties.get("thickness", 0) for layer_name, properties in self.vv_params["structure"].items()]
            )
            return bottleneck_radius + structural_thickness
        return bottleneck_radius


class SimpleCentralCell:
    """Builds the central cell structure for simple mirror machine."""
    
    def __init__(self, cc_params, material_ns):
        self.cc_params = cc_params
        self.material_ns = material_ns
        self.central_cell_cells = []
        self.central_cell_cylinders = []
        self.central_cell_cylinders_radii = []
        self.tally_descriptors = []
        
    def build_central_cell(self, outer_radius, room_region, vacuum_vessel_regions=None, axial_midplane=0.0):
        """Build the central cell structure."""
        # The central cell should start from the outermost radius of the vacuum vessel
        # and build outward, following the same logic as the tandem model
        inner_radius = outer_radius  # Start from the vacuum vessel's outer radius
        self.central_cell_cylinders_radii = inner_radius + np.cumsum([layer["thickness"] for layer in self.cc_params["layers"]])
        
        # Create ZCylinders for each layer
        for radius in self.central_cell_cylinders_radii:
            cylinder = openmc.ZCylinder(r=radius)
            self.central_cell_cylinders.append(cylinder)
        
        # Define axial boundary planes
        central_cell_length_from_midplane = self.cc_params["axial_length"] / 2
        left_plane_central_cell = openmc.ZPlane(float(axial_midplane) - central_cell_length_from_midplane)
        right_plane_central_cell = openmc.ZPlane(float(axial_midplane) + central_cell_length_from_midplane)
        
        # Create the first cell (innermost region)
        # IMPORTANT: Subtract the vacuum vessel regions to avoid overlap
        if vacuum_vessel_regions:
            first_cell_region = (-right_plane_central_cell & +left_plane_central_cell & -self.central_cell_cylinders[0])
            # Subtract all vacuum vessel regions
            for vv_region in vacuum_vessel_regions:
                first_cell_region = first_cell_region & ~vv_region
        else:
            first_cell_region = (-right_plane_central_cell & +left_plane_central_cell & -self.central_cell_cylinders[0])
        
        first_cell = openmc.Cell(
            2000,
            region=first_cell_region,
            fill=self.cc_params["materials"][0]
        )
        self.central_cell_cells.append(first_cell)
        
        # Create subsequent layers
        for i in range(1, len(self.central_cell_cylinders_radii)):
            cyl_region = (-right_plane_central_cell
                         & +left_plane_central_cell
                         & -self.central_cell_cylinders[i]
                         & +self.central_cell_cylinders[i - 1])
            
            # Subtract vacuum vessel regions from subsequent layers too
            if vacuum_vessel_regions:
                for vv_region in vacuum_vessel_regions:
                    cyl_region = cyl_region & ~vv_region
            
            new_cell = openmc.Cell(2000 + i, region=cyl_region, fill=self.cc_params["materials"][i])
            self.central_cell_cells.append(new_cell)
        
        # Add tally descriptors for central cell
        if "tallies" in self.cc_params:
            tally_data = self.cc_params["tallies"]
            
            # Handle breeder tallies (for the first layer - breeding blanket)
            if "breeder" in tally_data:
                breeder_tallies = tally_data["breeder"]
                if len(self.central_cell_cells) > 0:
                    self.tally_descriptors.append({
                        "type": "central_cell",
                        "location": "breeder",
                        "description": "breeding_blanket",
                        "cell": self.central_cell_cells[0],
                        "cell_tallies": breeder_tallies.get("cell_tallies", []),
                        "mesh_tallies": breeder_tallies.get("mesh_tallies", [])
                    })
            
            # Handle layer tallies (for specific layers)
            if "layer_tallies" in tally_data:
                for layer_entry in tally_data["layer_tallies"]:
                    position = layer_entry["position"]
                    if position < len(self.central_cell_cells):
                        layer_name = self.cc_params["layers"][position]["material"] if position < len(self.cc_params["layers"]) else "unknown"
                        self.tally_descriptors.append({
                            "type": "central_cell",
                            "location": f"layer_{position}",
                            "description": layer_name,
                            "cell": self.central_cell_cells[position],
                            "cell_tallies": layer_entry.get("cell_tallies", []),
                            "mesh_tallies": layer_entry.get("mesh_tallies", [])
                        })
        
        return self.central_cell_cells
    
    def get_outermost_radius(self):
        """Get the outermost radius of the central cell."""
        if len(self.central_cell_cylinders_radii) > 0:
            return self.central_cell_cylinders_radii[-1]
        return 0
    
    def get_tally_descriptors(self):
        """Get tally descriptors for central cell."""
        return self.tally_descriptors


class SimpleLFCoilBuilder:
    """Builds the LF (Low Field) coils for simple mirror machine."""
    
    def __init__(self, lf_params, material_ns):
        self.lf_params = lf_params
        self.material_ns = material_ns
        self.lf_coil_regions = []
        self.lf_coil_cells = []
        self.lf_coil_shield_regions = []
        self.lf_coil_shield_cells = []
        self.tally_descriptors = []
        
    def build_lf_coils(self, outer_radius, room_region, axial_midplane=0.0):
        """Build the LF coils."""
        for i, center_position in enumerate(self.lf_params["positions"]):
            # Generate cylindrical shell and inner region for the LF coil
            lf_coil_shell, lf_coil_inner = hollow_cylinder_with_shell(
                float(axial_midplane) + center_position,
                outer_radius,
                self.lf_params["inner_dimensions"]["radial_thickness"],
                self.lf_params["inner_dimensions"]["axial_length"],
                self.lf_params["shell_thicknesses"]["front"],
                self.lf_params["shell_thicknesses"]["back"],
                self.lf_params["shell_thicknesses"]["axial"]
            )
            
            # Create coil shell cell
            lf_coil_shell_cell = openmc.Cell(
                4000 + i,
                region=lf_coil_shell,
                fill=self.lf_params["materials"]["shield"]
            )
            
            self.lf_coil_shield_regions.append(lf_coil_shell)
            self.lf_coil_shield_cells.append(lf_coil_shell_cell)
            
            # Create coil inner cell
            lf_coil_inner_cell = openmc.Cell(
                4100 + i,
                region=lf_coil_inner,
                fill=self.lf_params["materials"]["magnet"]
            )
            
            self.lf_coil_regions.append(lf_coil_inner)
            self.lf_coil_cells.append(lf_coil_inner_cell)
            
            # Add tally descriptors for LF coils
            if "tallies" in self.lf_params:
                # Add tally for magnet cell
                self.tally_descriptors.append({
                    "type": "lf_coil",
                    "location": f"coil_{i}",
                    "description": "magnet",
                    "cell": lf_coil_inner_cell,
                    "cell_tallies": self.lf_params["tallies"].get("cell_tallies", []),
                    "mesh_tallies": self.lf_params["tallies"].get("mesh_tallies", [])
                })
        
        return self.lf_coil_cells + self.lf_coil_shield_cells
    
    def get_tally_descriptors(self):
        """Get tally descriptors for LF coils."""
        return self.tally_descriptors


class SimpleHFCoilBuilder:
    """Builds the HF (High Field) coils for simple mirror machine."""
    
    def __init__(self, hf_params, material_ns):
        self.hf_params = hf_params
        self.material_ns = material_ns
        self.hf_coil_cells = []
        self.hf_coil_regions = []
        self.tally_descriptors = []
        
    def build_hf_coils(self, central_cell_length_from_midplane, bottleneck_cylinders_radii, room_region, hf_z0_offset=0, axial_midplane=0.0):
        """Build the HF coils."""
        # Calculate HF coil center position
        casing_layers_thickness = np.array([layer["thickness"] for layer in self.hf_params["casing_layers"]])
        
        hf_coil_center_z0 = (
            central_cell_length_from_midplane
            + self.hf_params["shield"]["shield_central_cell_gap"]
            + self.hf_params["shield"]["axial_thickness"][0]
            + np.sum(casing_layers_thickness)
            + self.hf_params["magnet"]["axial_thickness"] / 2
            + hf_z0_offset
        )
        
        # Compute inner radius for magnet and casings
        # HF coils should start from the outermost radius of the vacuum vessel structural layers
        # NOT from the bottleneck radius to avoid overlap with vacuum vessel
        inner_radius_magnet_casings = (
            bottleneck_cylinders_radii[-1]  # This should be the VV outer radius, not bottleneck
            + self.hf_params["shield"]["radial_thickness"][0]
            + self.hf_params["shield"]["radial_gap_before_casing"]
        )
        
        # Define regions for left and right HF coils
        # Start from the vacuum vessel's outermost radius
        vv_outermost_radius = bottleneck_cylinders_radii[0]  # This is the VV outer radius
        
        # Create arrays that include casing layers + shield
        # Sequence: magnet -> casing layers -> shield
        all_layers_front_thickness = np.concatenate([casing_layers_thickness, [self.hf_params["shield"]["radial_thickness"][0]]])
        all_layers_back_thickness = np.concatenate([casing_layers_thickness, [self.hf_params["shield"]["radial_thickness"][1]]])
        symmetric_casing_axial = np.column_stack((casing_layers_thickness, casing_layers_thickness))
        shield_axial = np.asarray(self.hf_params["shield"]["axial_thickness"], dtype=float)
        if shield_axial.shape != (2,):
            raise ValueError("hf_coil.shield.axial_thickness must contain central-facing and outward-facing values")
        right_axial_thickness = np.vstack((symmetric_casing_axial, shield_axial))
        left_axial_thickness = np.vstack((symmetric_casing_axial, shield_axial[::-1]))
        
        shield_regions_right = nested_cylindrical_shells(
            float(axial_midplane) + hf_coil_center_z0,
            vv_outermost_radius,  # Start from the VV outer radius
            self.hf_params["magnet"]["radial_thickness"],
            self.hf_params["magnet"]["axial_thickness"],
            all_layers_front_thickness, all_layers_back_thickness, right_axial_thickness
        )
        
        shield_regions_left = nested_cylindrical_shells(
            float(axial_midplane) - hf_coil_center_z0,
            vv_outermost_radius,  # Start from the VV outer radius
            self.hf_params["magnet"]["radial_thickness"],
            self.hf_params["magnet"]["axial_thickness"],
            all_layers_front_thickness, all_layers_back_thickness, left_axial_thickness
        )
        
        # Create magnet cells
        hf_magnet_region_left = shield_regions_left[0]
        hf_magnet_region_right = shield_regions_right[0]
        
        hf_magnet_left_cell = openmc.Cell(
            6101, region=hf_magnet_region_left, 
            fill=getattr(self.material_ns, self.hf_params["magnet"]["material"])
        )
        hf_magnet_right_cell = openmc.Cell(
            6102, region=hf_magnet_region_right, 
            fill=getattr(self.material_ns, self.hf_params["magnet"]["material"])
        )
        
        self.hf_coil_cells.extend([hf_magnet_left_cell, hf_magnet_right_cell])
        self.hf_coil_regions.extend([hf_magnet_region_left, hf_magnet_region_right])
        
        # Create casing cells and shield
        casing_layers_materials = [getattr(self.material_ns, layer["material"]) for layer in self.hf_params["casing_layers"]]
        shield_material = getattr(self.material_ns, self.hf_params["shield"]["material"])
        
        for i in range(1, len(shield_regions_left)):
            left_region = shield_regions_left[i]
            right_region = shield_regions_right[i]

            # Determine material: casing layers first, then shield
            if i <= len(casing_layers_materials):
                material = casing_layers_materials[i - 1]
                left_cell_id = 6501 + 10 * (i - 1)
                right_cell_id = 6502 + 10 * (i - 1)
                description = f"casing layer {i}"
            else:
                material = shield_material
                shield_index = i - len(casing_layers_materials)
                left_cell_id = 6201 + 10 * (shield_index - 1)
                right_cell_id = 6202 + 10 * (shield_index - 1)
                description = f"shield layer {shield_index}"

            left_region_cell = openmc.Cell(
                left_cell_id, name=f"HF left {description}", region=left_region, 
                fill=material
            )
            right_region_cell = openmc.Cell(
                right_cell_id, name=f"HF right {description}", region=right_region, 
                fill=material
            )
            self.hf_coil_cells.extend([left_region_cell, right_region_cell])
            self.hf_coil_regions.extend([left_region, right_region])
        

        
        # Add tally descriptors for HF coils
        if "tallies" in self.hf_params:
            # Add tallies for magnet cells
            self.tally_descriptors.append({
                "type": "hf_coil",
                "location": "right",
                "description": "magnet",
                "cell": hf_magnet_right_cell,
                "cell_tallies": self.hf_params["tallies"].get("cell_tallies", []),
                "mesh_tallies": self.hf_params["tallies"].get("mesh_tallies", [])
            })
            self.tally_descriptors.append({
                "type": "hf_coil",
                "location": "left",
                "description": "magnet",
                "cell": hf_magnet_left_cell,
                "cell_tallies": self.hf_params["tallies"].get("cell_tallies", []),
                "mesh_tallies": self.hf_params["tallies"].get("mesh_tallies", [])
            })
        
        return self.hf_coil_cells
    
    def get_tally_descriptors(self):
        """Get tally descriptors for HF coils."""
        return self.tally_descriptors


class SimpleEndCellBuilder:
    """Builds the end cells for simple mirror machine."""
    
    def __init__(self, end_params, material_ns):
        self.end_params = end_params
        self.material_ns = material_ns
        self.end_cell_cells = []
        self.end_cell_regions = []
        self.tally_descriptors = []
        
    def build_end_cells(self, end_cell_z0_offset, hf_coil_shield, casing_layers_thickness, hf_coil_magnet, vacuum_section_regions, axial_midplane=0.0):
        """Build the end cells."""
        end_cell_radial_thickness = self.end_params["diameter"] / 2
        
        # Compute Z-Positions for End Cells following run_simple_mirror_model logic
        # The right end cell is positioned beyond the HF coil, accounting for:
        # 1. The axial thickness of the HF coil shield
        # 2. The total casing thickness  
        # 3. Half of the HF coil's axial thickness
        # 4. The end cell's shell thickness
        # 5. Half of the end cell's axial length
        
        right_end_cell_z0 = float(axial_midplane) + end_cell_z0_offset
        
        # The left end cell is symmetrically positioned at the negative Z-coordinate
        left_end_cell_z0 = float(axial_midplane) - end_cell_z0_offset
        
        # Define End Cell Shells and Inner Volumes
        # Each end cell consists of an outer shell and an inner vacuum region
        # Use cylinder_with_shell function as in the reference implementation
        
        left_end_cell_shell, left_end_cell_inner = cylinder_with_shell(
            left_end_cell_z0, 
            end_cell_radial_thickness,  # Use the already calculated radial thickness
            self.end_params["axial_length"], 
            self.end_params["shell_thickness"]
        )
        
        right_end_cell_shell, right_end_cell_inner = cylinder_with_shell(
            right_end_cell_z0, 
            end_cell_radial_thickness,  # Use the already calculated radial thickness
            self.end_params["axial_length"], 
            self.end_params["shell_thickness"]
        )
        
        # Exclude all regions of the vacuum section from the inner volumes and shells
        for region in vacuum_section_regions:
            left_end_cell_shell &= ~region
            right_end_cell_shell &= ~region
            left_end_cell_inner &= ~region
            right_end_cell_inner &= ~region
        
        # Create cells
        end_cell_left_shell_cell = openmc.Cell(
            5001, region=left_end_cell_shell,
            fill=self.end_params["shell_material"]
        )
        
        end_cell_left_inner_cell = openmc.Cell(
            5002, region=left_end_cell_inner,
            fill=self.end_params["inner_material"]
        )
        
        end_cell_right_shell_cell = openmc.Cell(
            5003, region=right_end_cell_shell,
            fill=self.end_params["shell_material"]
        )
        
        end_cell_right_inner_cell = openmc.Cell(
            5004, region=right_end_cell_inner,
            fill=self.end_params["inner_material"]
        )
        
        self.end_cell_cells.extend([
            end_cell_left_shell_cell, end_cell_left_inner_cell,
            end_cell_right_shell_cell, end_cell_right_inner_cell
        ])
        
        self.end_cell_regions.extend([
            left_end_cell_shell, left_end_cell_inner,
            right_end_cell_shell, right_end_cell_inner
        ])
        
        # Add tally descriptors for end cells
        if "tallies" in self.end_params:
            self.tally_descriptors.append({
                "type": "end_cell",
                "location": "left",
                "description": "shell",
                "cell": end_cell_left_shell_cell,
                "cell_tallies": self.end_params["tallies"].get("cell_tallies", []),
                "mesh_tallies": self.end_params["tallies"].get("mesh_tallies", [])
            })
            self.tally_descriptors.append({
                "type": "end_cell",
                "location": "right",
                "description": "shell",
                "cell": end_cell_right_shell_cell,
                "cell_tallies": self.end_params["tallies"].get("cell_tallies", []),
                "mesh_tallies": self.end_params["tallies"].get("mesh_tallies", [])
            })
        
        return self.end_cell_cells
    
    def get_tally_descriptors(self):
        """Get tally descriptors for end cells."""
        return self.tally_descriptors


class SimpleMachineBuilder:
    """Main builder class for simple mirror machine."""
    
    def __init__(self, vv_params, cc_params, lf_params, hf_params, end_params, material_ns):
        self.vv_params = vv_params
        self.cc_params = cc_params
        self.lf_params = lf_params
        self.hf_params = hf_params
        self.end_params = end_params
        self.material_ns = material_ns
        
        self.universe_machine = openmc.Universe(786)
        self.room_region = None
        self.all_cells = []
        
        # Initialize component builders
        self.vacuum_vessel = SimpleVacuumVessel(vv_params, material_ns)
        self.central_cell = SimpleCentralCell(cc_params, material_ns)
        self.lf_coils = SimpleLFCoilBuilder(lf_params, material_ns)
        self.hf_coils = SimpleHFCoilBuilder(hf_params, material_ns)
        self.end_cells = SimpleEndCellBuilder(end_params, material_ns)
        
        # Initialize tally builder
        self.tally_builder = TallyBuilder()
        
    def _add_cell(self, cell):
        """Add a cell to the universe."""
        self.universe_machine.add_cells([cell])
        self.all_cells.append(cell)
        
    def _subtract_from_room(self, regions):
        """Subtract regions from the room."""
        for region in regions:
            self.room_region &= ~region

    def _model_bounds_cm(self):
        """Return geometry-derived computational enclosure bounds.

        The enclosure is an OpenMC model boundary, not a charged-particle wall.
        A dimensionless padding factor avoids coincident surfaces without
        embedding a machine length or coordinate.
        """
        mid = float(self.vv_params.get("axial_midplane", 0.0))
        structure = sum(float(item.get("thickness", 0.0)) for item in self.vv_params.get("structure", {}).values())
        casing = sum(float(item["thickness"]) for item in self.hf_params["casing_layers"])
        shield_radial = [float(value) for value in self.hf_params["shield"]["radial_thickness"]]
        hf_outer_radius = (
            float(self.vv_params["bottleneck_radius"])
            + structure
            + casing
            + shield_radial[0]
            + float(self.hf_params["magnet"]["radial_thickness"])
            + casing
            + shield_radial[1]
        )
        radial_extent = max(
            float(self.vv_params["central_radius"]) + structure + sum(float(item["thickness"]) for item in self.cc_params["layers"]),
            hf_outer_radius,
            0.5 * float(self.end_params["diameter"]) + float(self.end_params["shell_thickness"]),
        )
        z_min = float(self.end_params["derived_full_z_min_cm"])
        z_max = float(self.end_params["derived_full_z_max_cm"])
        padding = 1.05
        return (-padding * radial_extent, padding * radial_extent, mid + padding * (z_min - mid), mid + padding * (z_max - mid))
            
    def build_vacuum_vessel(self):
        """Build the vacuum vessel."""
        # Create a geometry-derived computational room
        r_min, r_max, z_min, z_max = self._model_bounds_cm()
        bounding_surface = openmc.model.RectangularParallelepiped(
            xmin=r_min, xmax=r_max, ymin=r_min, ymax=r_max,
            zmin=z_min, zmax=z_max, boundary_type='vacuum'
        )
        self.room_region = -bounding_surface
        room_cell = openmc.Cell(100, region=self.room_region, fill=self.material_ns.air)
        self._add_cell(room_cell)
        
        # Build vacuum vessel
        vv_cells, vv_regions = self.vacuum_vessel.build_vacuum_vessel(self.room_region)
        for cell in vv_cells:
            self._add_cell(cell)
        self._subtract_from_room(vv_regions)
        
    def build_central_cell(self):
        """Build the central cell."""
        outer_radius = self.vacuum_vessel.get_outermost_radius()
        # Pass vacuum vessel regions to avoid overlap
        vacuum_vessel_regions = self.vacuum_vessel.vacuum_section_regions
        cc_cells = self.central_cell.build_central_cell(
            outer_radius,
            self.room_region,
            vacuum_vessel_regions,
            axial_midplane=self.vv_params.get("axial_midplane", 0.0),
        )
        for cell in cc_cells:
            self._add_cell(cell)
        self._subtract_from_room([cell.region for cell in cc_cells])
        
    def build_lf_coils(self):
        """Build the LF coils."""
        outer_radius = self.central_cell.get_outermost_radius()
        lf_cells = self.lf_coils.build_lf_coils(
            outer_radius,
            self.room_region,
            axial_midplane=self.vv_params.get("axial_midplane", 0.0),
        )
        for cell in lf_cells:
            self._add_cell(cell)
        self._subtract_from_room([cell.region for cell in lf_cells])
        
    def build_hf_coils(self):
        """Build the HF coils."""
        central_cell_length_from_midplane = self.cc_params["axial_length"] / 2
        
        # Calculate HF coil offset based on axial length relationship
        hf_z0_offset = (self.vv_params["central_axial_length"] >= self.cc_params["axial_length"]) * (self.vv_params["central_axial_length"] / 2 - self.cc_params["axial_length"] / 2)
        hf_z0_offset += float(self.hf_params.get("axial_offset_cm", 0.0))
        
        # HF coils should start from the outermost radius of the vacuum vessel structural layers
        # NOT from the bottleneck radius to avoid overlap with vacuum vessel
        vv_outermost_radius = self.vacuum_vessel.get_outermost_bottleneck_radius()
        hf_starting_radii = [vv_outermost_radius]
        
        hf_cells = self.hf_coils.build_hf_coils(
            central_cell_length_from_midplane, 
            hf_starting_radii, 
            self.room_region,
            hf_z0_offset,
            axial_midplane=self.vv_params.get("axial_midplane", 0.0),
        )
        for cell in hf_cells:
            self._add_cell(cell)
        self._subtract_from_room([cell.region for cell in hf_cells])
        
    def build_end_cells(self):
        """Build the end cells."""
        # The pure-data layout is the sole source of the end-cell center
        casing_layers_thickness = np.array([layer["thickness"] for layer in self.hf_params["casing_layers"]])
        end_cell_z0_offset = float(self.end_params["derived_center_offset_cm"])

        # # Combine all the regions of the vacuum vessel
        # vacuum_section_regions = self.vacuum_vessel.vacuum_section_regions[0]
        # # Remove the first region from the list
        # for region in self.vacuum_vessel.vacuum_section_regions[1:]:
        #     vacuum_section_regions |= region

        end_cells = self.end_cells.build_end_cells(
            end_cell_z0_offset,  # Position end cells beyond VV boundaries
            self.hf_params["shield"],
            casing_layers_thickness,
            self.hf_params["magnet"],
            self.vacuum_vessel.vacuum_section_regions,
            axial_midplane=self.vv_params.get("axial_midplane", 0.0),
        )
        for cell in end_cells:
            self._add_cell(cell)
        self._subtract_from_room([cell.region for cell in end_cells])
        
    def get_universe(self):
        """Get the machine universe."""
        return self.universe_machine
        
    def get_all_cells(self):
        """Get all cells in the machine."""
        return self.all_cells
        
    def get_all_regions(self):
        """Get all regions in the machine."""
        return [cell.region for cell in self.all_cells]
    
    def get_all_tallies(self):
        """Get all tallies configured for the model."""
        # Collect tally descriptors from all components
        all_descriptors = []
        all_descriptors.extend(self.central_cell.get_tally_descriptors())
        all_descriptors.extend(self.lf_coils.get_tally_descriptors())
        all_descriptors.extend(self.hf_coils.get_tally_descriptors())
        all_descriptors.extend(self.end_cells.get_tally_descriptors())
        
        # Add descriptors to tally builder
        self.tally_builder.add_descriptors(all_descriptors)
        
        return self.tally_builder.get_tallies()


@contextlib.contextmanager
def change_dir(destination):
    """Context manager to change directory temporarily."""
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    original_dir = os.getcwd()
    os.chdir(destination)
    try:
        yield
    finally:
        os.chdir(original_dir)


def _sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _package_validated_file_source_for_model(source_build, model_xml_dir):
    bundle = source_build.pipeline_result.openmc_export.file_source_bundle
    if bundle is None:
        from source_model_revamp.integration.validation import SourceConstructionError
        raise SourceConstructionError(
            "missing_correlated_file_source",
            "ParaTAN requires the correlated HDF5 FileSource bundle",
            quantity="openmc_file_source_bundle",
        )

    configured_source_filename = str(source_build.run_config.openmc_export.source_filename)
    source_filename = Path(configured_source_filename).name
    if not source_filename:
        raise ValueError("configured OpenMC source filename has no basename")
    authoritative_source = Path(bundle.source_file_path)
    if not authoritative_source.is_file():
        raise FileNotFoundError(
            f"validated OpenMC source file does not exist: {authoritative_source}"
        )

    model_xml_dir = Path(model_xml_dir)
    model_xml_dir.mkdir(parents=True, exist_ok=True)
    packaged_source = model_xml_dir / source_filename

    if authoritative_source.resolve() != packaged_source.resolve():
        shutil.copy2(authoritative_source, packaged_source)

    packaged_hash = _sha256_file(packaged_source)
    validated_hash = str(bundle.metadata.source_file_sha256)
    if packaged_hash != validated_hash:
        raise OSError(
            "packaged OpenMC source hash does not match the validated source"
        )

    with openmc.config.patch("resolve_paths", False):
        file_source = openmc.FileSource(source_filename, strength=1.0)
    return file_source, packaged_source


def _validate_packaged_settings_source(settings_path, source_filename):
    source = ET.parse(settings_path).find("source")
    if source is None:
        raise OSError("exported settings.xml is missing the OpenMC source element")
    if source.attrib.get("type") != "file":
        raise OSError("exported settings.xml source is not a file source")
    if source.attrib.get("file") != source_filename:
        raise OSError(
            "exported settings.xml does not reference the packaged source by filename"
        )
    try:
        strength = float(source.attrib.get("strength", "nan"))
    except ValueError as exc:
        raise OSError("exported settings.xml source strength is invalid") from exc
    if strength != 1.0:
        raise OSError("exported settings.xml source strength must equal 1.0")


def build_simple_model_from_input(input_data, output_dir=None, *, source_information=None, source_model_progress_callback=None):
    """
    Build a complete simple mirror model from input parameters.
    Returns a complete OpenMC model ready to run.
    """
    # Load materials dynamically based on environment variable
    m = load_materials()
    from mirror_paratan.source.mirror_source_model import default_paratan_run_dir

    original_dir = Path.cwd()
    if output_dir is None or str(output_dir) in {"", "."}:
        output_dir = default_paratan_run_dir(input_data, base_dir=original_dir)

    with change_dir(output_dir):
        # Parse input parameters
        vv_params, cc_params, lf_params, hf_params, end_params = parse_simple_machine_input(input_data, m)

        # Create builder
        builder = SimpleMachineBuilder(
            vv_params, cc_params, lf_params, hf_params, end_params, m
        )

        # Build components - add central cell, LF coils, HF coils, and end cells
        builder.build_vacuum_vessel()
        builder.build_central_cell()
        builder.build_lf_coils()
        builder.build_hf_coils()
        builder.build_end_cells()
        
        # Create geometry
        universe_machine = builder.get_universe()
        geometry = openmc.Geometry([openmc.Cell(fill=universe_machine)])
        geometry.merge_surfaces = True

        # Materials
        materials = m.materials

        # Resolve OpenMC run settings and create neutron sources from source_model_revamp.
        # ParaTAN remains responsible for geometry/materials/tallies/settings XML export;
        # source_model_revamp is responsible for the magnetic-mirror neutron source.
        from mirror_paratan.source.mirror_source_model import build_source_model_openmc_source

        source_build = build_source_model_openmc_source(
            input_data,
            source_information,
            output_dir=Path.cwd(),
            write_outputs=False,
            progress_callback=source_model_progress_callback,
        )
        source = list(source_build.sources)
        if len(source) != 1 or source_build.pipeline_result.openmc_export.file_source_bundle is None:
            from source_model_revamp.integration.validation import SourceConstructionError
            raise SourceConstructionError("invalid_paratan_source_bundle", "ParaTAN requires exactly one correlated FileSource", quantity="settings.source", value=len(source))
        from dataclasses import replace
        from source_model_revamp.integration.output_writer import write_resolved_input_yaml, write_source_model_outputs

        output_files = write_source_model_outputs(source_build.pipeline_result, Path.cwd())
        output_files["resolved_input_yaml"] = str(write_resolved_input_yaml(input_data, source_build.run_config, Path.cwd()))
        input_used_file = Path.cwd() / "input_used.yaml"
        if input_used_file.exists():
            output_files["input_used_yaml"] = str(input_used_file)
        source_build = replace(source_build, output_files=output_files)

        geometry_plot_requested = bool(source_build.source_information["settings"].get("geometry_plot", False))
        source_build.metadata["paratan_geometry_plot_generated"] = geometry_plot_requested
        if geometry_plot_requested:
            import matplotlib.pyplot as plt
            r_min, r_max, z_min, z_max = builder._model_bounds_cm()
            geometry.root_universe.plot(
                basis='xz',
                width=(r_max - r_min, z_max - z_min),
                pixels=(2400, 4000),
                color_by='material'
            )
            plt.savefig('simple_mirror_modular_cross_section.png', bbox_inches="tight")
            plt.close()

        # Create settings
        source_settings = source_build.source_information["settings"]
        tallies = openmc.Tallies(builder.get_all_tallies()) if source_settings['tallies'] else openmc.Tallies()
        settings = openmc.Settings()
        settings.run_mode = "fixed source"
        settings.particles = int(source_settings['particles_per_batch'])
        settings.batches = int(source_settings['batches'])
        settings.output = {'tallies': bool(source_settings['tallies'])}

        freq = int(source_settings['statepoint_frequency'])
        settings.statepoint = {
            'batches': [1] + list(range(freq, settings.batches, freq)) + [settings.batches]
        }

        settings.weight_windows_on = bool(source_settings['weight_windows'])
        if settings.weight_windows_on:
            settings.weight_window_checkpoints = {'collision': True, 'surface': True}
            settings.weight_windows_generator = [openmc.WeightWindowGenerator(
                openmc.RegularMesh(),
                energy_bounds=[0, 14e6],
                particle_type='neutron',
                method='magic',
                max_realizations=25,
                update_interval=1,
                on_the_fly=True
            )]

        settings.photon_transport = bool(source_settings['photon_transport'])

        model_xml_dir = Path("model_xml_files")
        source_filename = Path(
            str(source_build.run_config.openmc_export.source_filename)
        ).name
        packaged_file_source, packaged_source_path = _package_validated_file_source_for_model(
            source_build,
            model_xml_dir,
        )
        settings.source = packaged_file_source

        model = openmc.Model(geometry, materials, settings, tallies)
        model.export_to_xml(model_xml_dir)

        settings_path = model_xml_dir / "settings.xml"
        _validate_packaged_settings_source(settings_path, source_filename)

        packaged_output_files = dict(source_build.output_files)
        packaged_output_files["openmc_model_settings_xml"] = str(settings_path)
        packaged_output_files["openmc_model_source_h5"] = str(packaged_source_path)
        source_build.metadata["paratan_openmc_source_packaged_with_settings"] = True
        source_build.metadata["paratan_openmc_source_reference"] = source_filename
        source_build = replace(source_build, output_files=packaged_output_files)

        model.source_model_build_result = source_build
        model.source_model_metadata = dict(source_build.metadata)

    return model
