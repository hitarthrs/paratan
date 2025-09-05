import openmc
import numpy as np
import matplotlib.pyplot as plt

# Materials library that is used to build the model. To channge materials, import other directory as m

from src.paratan.materials import material as m

from src.paratan.geometry.core import *
from src.paratan.source.core import *
from src.paratan.tallies.base_tallies import hollow_mesh_from_domain, strings_to_openmc_filters
from src.paratan.tallies.tandem_tallies import TallyBuilder
import yaml
import os
import contextlib
from types import SimpleNamespace


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
        "vacuum_material": getattr(material_ns, "vacuum")
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
        "tallies": hf_info.get("hf_coil_tallies", {})  # Updated to use nested tallies
    }

    # End cell parameters
    end_info = input_data.get("end_cell", {})
    end_params = {
        "axial_length": end_info.get("axial_length", 50),
        "shell_thickness": end_info.get("shell_thickness", 2),
        "diameter": end_info.get("diameter", 100),
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
        """Build the vacuum vessel structure using single_vacuum_vessel_region."""
        self.room_region = room_region
        
        # Create main vacuum vessel region using the new function
        vv_main_region, vv_components = single_vacuum_vessel_region(
            self.vv_params["outer_axial_length"],
            self.vv_params["central_axial_length"], 
            self.vv_params["central_radius"],
            self.vv_params["bottleneck_radius"],
            self.vv_params["left_bottleneck_length"],
            self.vv_params["right_bottleneck_length"],
            self.vv_params["axial_midplane"]
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
                vv_layer_region, _ = single_vacuum_vessel_region(
                    self.vv_params["outer_axial_length"],
                    self.vv_params["central_axial_length"], 
                    central_radii[i + 1],      # Updated central radius
                    bottleneck_radii[i + 1],    # Updated bottleneck radius
                    self.vv_params["left_bottleneck_length"],
                    self.vv_params["right_bottleneck_length"],
                    self.vv_params["axial_midplane"]
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
        
    def build_central_cell(self, outer_radius, room_region, vacuum_vessel_regions=None):
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
        left_plane_central_cell = openmc.ZPlane(-central_cell_length_from_midplane)
        right_plane_central_cell = openmc.ZPlane(central_cell_length_from_midplane)
        
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
        
    def build_lf_coils(self, outer_radius, room_region):
        """Build the LF coils."""
        for i, center_position in enumerate(self.lf_params["positions"]):
            # Generate cylindrical shell and inner region for the LF coil
            lf_coil_shell, lf_coil_inner = hollow_cylinder_with_shell(
                center_position,
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
        
    def build_hf_coils(self, central_cell_length_from_midplane, bottleneck_cylinders_radii, room_region, hf_z0_offset=0):
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
        all_layers_axial_thickness = np.concatenate([casing_layers_thickness, [self.hf_params["shield"]["axial_thickness"][0]]])
        
        shield_regions_right = nested_cylindrical_shells(
            hf_coil_center_z0,
            vv_outermost_radius,  # Start from the VV outer radius
            self.hf_params["magnet"]["radial_thickness"],
            self.hf_params["magnet"]["axial_thickness"],
            all_layers_front_thickness, all_layers_back_thickness, all_layers_axial_thickness
        )
        
        shield_regions_left = nested_cylindrical_shells(
            -hf_coil_center_z0,
            vv_outermost_radius,  # Start from the VV outer radius
            self.hf_params["magnet"]["radial_thickness"],
            self.hf_params["magnet"]["axial_thickness"],
            all_layers_front_thickness, all_layers_back_thickness, all_layers_axial_thickness
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
            combined_region = shield_regions_left[i] | shield_regions_right[i]
            
            # Determine material: casing layers first, then shield
            if i <= len(casing_layers_materials):
                material = casing_layers_materials[i - 1]
                cell_id = 6500 + i
            else:
                material = shield_material
                cell_id = 6200 + (i - len(casing_layers_materials))
            
            combined_region_cell = openmc.Cell(
                cell_id, region=combined_region, 
                fill=material
            )
            self.hf_coil_cells.append(combined_region_cell)
            self.hf_coil_regions.append(combined_region)
        

        
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
        
    def build_end_cells(self, end_cell_z0_offset, hf_coil_shield, casing_layers_thickness, hf_coil_magnet, vacuum_section_regions):
        """Build the end cells."""
        end_cell_radial_thickness = self.end_params["diameter"] / 2
        
        # Compute Z-Positions for End Cells following run_simple_mirror_model logic
        # The right end cell is positioned beyond the HF coil, accounting for:
        # 1. The axial thickness of the HF coil shield
        # 2. The total casing thickness  
        # 3. Half of the HF coil's axial thickness
        # 4. The end cell's shell thickness
        # 5. Half of the end cell's axial length
        
        right_end_cell_z0 = (
            end_cell_z0_offset  # This is already calculated  Half of the end cell's axial length
        )
        
        # The left end cell is symmetrically positioned at the negative Z-coordinate
        left_end_cell_z0 = -right_end_cell_z0
        
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
            
    def build_vacuum_vessel(self):
        """Build the vacuum vessel."""
        # Create room
        bounding_surface = openmc.model.RectangularParallelepiped(
            xmin=-400, xmax=400, ymin=-400, ymax=400, 
            zmin=-1350, zmax=1350, boundary_type='vacuum'
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
        cc_cells = self.central_cell.build_central_cell(outer_radius, self.room_region, vacuum_vessel_regions)
        for cell in cc_cells:
            self._add_cell(cell)
        self._subtract_from_room([cell.region for cell in cc_cells])
        
    def build_lf_coils(self):
        """Build the LF coils."""
        outer_radius = self.central_cell.get_outermost_radius()
        lf_cells = self.lf_coils.build_lf_coils(outer_radius, self.room_region)
        for cell in lf_cells:
            self._add_cell(cell)
        self._subtract_from_room([cell.region for cell in lf_cells])
        
    def build_hf_coils(self):
        """Build the HF coils."""
        central_cell_length_from_midplane = self.cc_params["axial_length"] / 2
        
        # Calculate HF coil offset based on axial length relationship
        hf_z0_offset = (self.vv_params["central_axial_length"] >= self.cc_params["axial_length"]) * (self.vv_params["central_axial_length"] / 2 - self.cc_params["axial_length"] / 2)
        
        # HF coils should start from the outermost radius of the vacuum vessel structural layers
        # NOT from the bottleneck radius to avoid overlap with vacuum vessel
        vv_outermost_radius = self.vacuum_vessel.get_outermost_bottleneck_radius()
        hf_starting_radii = [vv_outermost_radius]
        
        hf_cells = self.hf_coils.build_hf_coils(
            central_cell_length_from_midplane, 
            hf_starting_radii, 
            self.room_region,
            hf_z0_offset
        )
        for cell in hf_cells:
            self._add_cell(cell)
        self._subtract_from_room([cell.region for cell in hf_cells])
        
    def build_end_cells(self):
        """Build the end cells."""
        # Calculate parameters needed for end cells
        central_cell_length_from_midplane = self.cc_params["axial_length"] / 2
        casing_layers_thickness = np.array([layer["thickness"] for layer in self.hf_params["casing_layers"]])
        
        # Calculate HF coil offset based on axial length relationship
        hf_z0_offset = (self.vv_params["central_axial_length"] >= self.cc_params["axial_length"]) * (self.vv_params["central_axial_length"] / 2 - self.cc_params["axial_length"] / 2) + 10

        # Calculate HF coil center position
        hf_coil_center_z0 = (
            central_cell_length_from_midplane
            + self.hf_params["shield"]["shield_central_cell_gap"]
            + self.hf_params["shield"]["axial_thickness"][0]
            + np.sum(casing_layers_thickness)
            + self.hf_params["magnet"]["axial_thickness"] / 2
            + hf_z0_offset
        )
        
        # Position end cells following the exact logic from run_simple_mirror_model
        # End cells are positioned right after the HF coils, not way out beyond the VV
        # We need to add the shield thickness, casing thickness, and magnet thickness
        # to get to the edge of the HF coil, then add end cell dimensions
        end_cell_z0_offset = (
            hf_coil_center_z0 
            + self.hf_params["shield"]["axial_thickness"][0]  # Axial thickness of the shield (toward midplane)
            + np.sum(casing_layers_thickness)  # Total casing thickness
            + self.hf_params["magnet"]["axial_thickness"] / 2  # Half of the magnet's axial thickness
            + self.end_params["shell_thickness"]  # Shell thickness of the end cell
            + self.end_params["axial_length"] / 2  # Half of the end cell's axial length
        )

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
            self.vacuum_vessel.vacuum_section_regions
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
    original_dir = os.getcwd()
    os.chdir(destination)
    try:
        yield
    finally:
        os.chdir(original_dir)


def build_simple_model_from_input(input_data, output_dir="."):
    """
    Build a complete simple mirror model from input parameters.
    Returns a complete OpenMC model ready to run.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Calculate the source file path BEFORE changing directories
    # The output_dir is now the config directory itself, so source file is in the same directory
    source_file_path = os.path.abspath(os.path.join(output_dir, 'source_information.yaml'))
    
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
        
        # Check spacing for other components
        total_vv_length = (builder.vv_params["left_bottleneck_length"] + 
                          builder.vv_params["central_axial_length"] + 
                          builder.vv_params["right_bottleneck_length"])
        central_cell_length = builder.cc_params.get("axial_length", 0)
        available_space_per_side = (total_vv_length - central_cell_length) / 2
        
        if available_space_per_side < 250:
            print(f"WARNING: Available space per side ({available_space_per_side:.1f} cm) is less than 250 cm!")
            print(f"Total VV length: {total_vv_length} cm, Central cell: {central_cell_length} cm")
        else:
            print(f"✓ Spacing OK: {available_space_per_side:.1f} cm available on each side")

        # Create geometry
        universe_machine = builder.get_universe()
        geometry = openmc.Geometry([openmc.Cell(fill=universe_machine)])
        geometry.merge_surfaces = True

        # Generate cross-sectional plot
        geometry.root_universe.plot(
            basis='xz',
            width=(1000, 2800),
            pixels=(2400, 4000),
            color_by='material'
        )
        plt.savefig('simple_mirror_modular_cross_section.png', bbox_inches="tight")

        # Export geometry
        geometry.export_to_xml("geometry.xml")

        # Create and export tallies
        tallies = openmc.Tallies(builder.get_all_tallies())
        tallies.export_to_xml("tallies.xml")

        # Materials
        materials = m.materials

        # Load source configuration
        # Use the source_file_path calculated before changing directories
        
        with open(source_file_path, 'r') as f:
            source_data = yaml.safe_load(f)

        # Create source based on type
        source_type = source_data['source']['type']
        power_output = source_data['source']['power_output']
        
        if source_type == "Volumetric":
            # Volumetric source using vacuum vessel parameters
            from src.paratan.source.core import VolumetricSource
            source = VolumetricSource(
                power_output=power_output,
                vacuum_vessel_axial_length=builder.vv_params["central_axial_length"],
                vacuum_vessel_outer_axial_length=builder.vv_params["outer_axial_length"],
                vacuum_vessel_central_radius=builder.vv_params["central_radius"],
                throat_radius=builder.vv_params["bottleneck_radius"],
                z_origin=builder.vv_params["axial_midplane"],
                conical_sources=5
            ).create_openmc_source()
            
        elif source_type == "Uniform":
            # Uniform cylindrical source
            from src.paratan.source.core import UniformSource
            length = source_data['source']['uniform']['length']
            radius = source_data['source']['uniform']['radius']
            source = UniformSource(power_output, length, radius).create_openmc_source()
            
        elif source_type == "1D_Varying":
            # 1D varying source from file
            from src.paratan.source.core import Source1D
            file_name = source_data['source']['source_1D']['file_name']
            radius = source_data['source']['source_1D']['radius']
            source = Source1D(power_output, radius, file_name).create_openmc_source()
            
        elif source_type == "2D_Varying":
            # 2D varying source from file
            from src.paratan.source.core import Source2D
            file_name = source_data['source']['source_2D']['file_name']
            source = Source2D(power_output, file_name).create_openmc_source()
            
        elif source_type == "Custom":
            # Custom source - handle later
            print("Custom source type not yet implemented")
            # Fallback to default source
            source = openmc.Source()
            source.space = openmc.stats.CylindricalIndependent(
                r=openmc.stats.Discrete([0], [1.0]),
                phi=openmc.stats.Uniform(a=0.0, b=2 * np.pi),
                z=openmc.stats.Uniform(a=-600.0, b=600.0),
            )
            source.angle = openmc.stats.Isotropic()
            source.energy = openmc.stats.Discrete([14e6], [1.0])
        else:
            raise ValueError(f"Unsupported source type: {source_type}")

        # Create settings
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

        # Weight window generator
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

        # Create and export final model
        model = openmc.Model(geometry, materials, settings, tallies)
        model.export_to_xml('model_xml_files')

    return model
