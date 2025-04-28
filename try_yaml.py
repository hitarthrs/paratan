import yaml
import openmc
from source_lib import *
from geometry_lib import *
import material as m
import matplotlib.pyplot as plt
# ________________________________________ Get the source from Settings information file ________________________________________ #

# with open('input_files/source_information.yaml', 'r') as f:
#     source_data = yaml.safe_load(f)

# openmc_source = load_source_from_yaml('input_files/source_information.yaml')

# settings = openmc.Settings()
# settings.run_mode = "fixed source"
# settings.particles = int(source_data['settings']['particles_per_batch'])
# settings.batches = source_data['settings']['batches']
# settings.output = {'tallies': False}
# settings.statepoint = {
#     'batches': [1] + list(range(source_data['settings']['statepoint_frequency'], settings.batches, source_data['settings']['statepoint_frequency'])) + [settings.batches]
# }
# settings.weight_windows_on = source_data['settings']['weight_windows']
# settings.weight_window_checkpoints = {'collision': True, 'surface': True}
# wwg = openmc.WeightWindowGenerator(openmc.RegularMesh(), 
#                                    [0, 14e6], 
#                                    'neutron', 
#                                    'magic', 
#                                    max_realizations=25, 
#                                    update_interval=1, 
#                                    on_the_fly=True)
# settings.weight_windows_generator = [wwg]
# settings.photon_transport = source_data['settings']['photon_transport']
# settings.source = openmc_source



# lf_coil_shell, lf_coil_inner = hollow_cylinder_with_shell(
#         2,      # Coil center position
#         10,  # Outer reference radius
#         15,    # Radial thickness of the coil
#         15,  # Inner axial length
#         3,  # Shell front thickness
#         3,   # Shell back thickness
#         3   # Shell axial thickness
#     )



# def hollow_mesh_from_domain(region, dimensions= [10, 10, 10], phi_grid_bounds=(0.0, 2 * np.pi)):
#     """
#     Generate a cylindrical mesh overs a hollow region defined by an OpenMC region.
    
#     Parameters:
#         region (openmc.Region): The region to bound and mesh (not necessarily hollow).
#         dimensions (tuple): Number of divisions in (r, phi, z), i.e., (nr, nphi, nz).
#         phi_grid_bounds (tuple): Angular bounds in radians for phi. Default is (0, 2π).
    
#     Returns:
#         openmc.CylindricalMesh: A cylindrical mesh over the hollow region.
#     """
#     # Get the bounding box of the region
#     bounding_box = region.bounding_box
    
#     # Determine max radial extent from bounding box corners
#     max_radius = max(
#         bounding_box[0][0],  # x-min
#         bounding_box[0][1],  # y-min
#         bounding_box[1][0],  # x-max
#         bounding_box[1][1]   # y-max
#     )
    
#     # Create outer bounding cylindrical surfaces
#     outer_cylinder = openmc.ZCylinder(r=max_radius)
#     lower_z = openmc.ZPlane(bounding_box[0][2])
#     upper_z = openmc.ZPlane(bounding_box[1][2])
    
#     outer_region = -outer_cylinder & +lower_z & -upper_z
    
#     # Subtract the original region to define hollow space
#     hollow_region = outer_region & ~region
    
#     # Extract all surfaces in the resulting region
#     surfaces = hollow_region.get_surfaces()
    
#     # Find all z-cylindrical surfaces and collect their radii
#     radii = [
#         surface.coefficients['r']
#         for surface in surfaces.values()
#         if surface.type == 'z-cylinder'
#     ]
    
#     # Set inner radius based on smallest detected cylindrical surface
#     if radii:
#         min_radius = min(radii)
#     else:
#         min_radius = 0.0  # fallback if no cylinders are found
    
#     # Build the r, phi, z grids
#     r_grid = np.linspace(min_radius, max_radius, num=dimensions[0] + 1)
#     phi_grid = np.linspace(phi_grid_bounds[0], phi_grid_bounds[1], num=dimensions[1] + 1)
#     z_grid = np.linspace(bounding_box[0][2], bounding_box[1][2], num=dimensions[2] + 1)


#     origin = (bounding_box.center[0], bounding_box.center[1], z_grid[0])

#     z_grid -= origin[2]

#     # Construct and return the cylindrical mesh

#     cyl_mesh = openmc.CylindricalMesh(r_grid=r_grid, phi_grid=phi_grid, z_grid=z_grid, origin=origin)
    
#     return cyl_mesh



# with open('input_files/parametric_input.yaml', 'r') as f:
#     input_data = yaml.safe_load(f)

# vacuum_vessel = input_data.get("vacuum_vessel",{})

# first_vv_plane_from_midplane = vacuum_vessel.get("first_vv_plane_from_midplane")
# machine_length_from_midplane = vacuum_vessel.get("machine_length_from_midplane")

# vacuum_chamber = vacuum_vessel.get("vacuum_chamber", {})
# vacuum_chamber_radius = vacuum_chamber.get("radius")

# structure = vacuum_vessel.get("structure", {})

# for layer_name, props in structure.items():
#     thickness = props["thickness"]
#     material_name = props["material"]
#     material_obj = getattr(m, material_name)

#     print(f"{layer_name}: thickness={thickness}, material={material_obj}")


# central_cell = input_data.get("central_cell", {})

# layers_cc = [np.array([layer["thickness"], getattr(m, layer["material"]).name]) for layer in central_cell.get("layers", [])]
# central_cell_layers = np.array(layers_cc)

# print(central_cell_layers)

# lf_coil = input_data.get("lf_coil", [])

# lf_coil_centers = lf_coil.get("positions",[])
# lf_coil_inner_dimensions = lf_coil.get("inner_dimensions", [])

# print(lf_coil_inner_dimensions["radial_thickness"])


# hf_coil = input_data.get("hf_coil", [])

# casing_thickness  = [layer["thickness"] for layer in hf_coil.get("casing_layers", [])]
# print(f"Casing layer thicknesses {casing_thickness}")

import xml.etree.ElementTree as ET

# Load the XML file
tree = ET.parse("your_file.xml")