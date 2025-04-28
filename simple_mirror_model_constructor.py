import openmc
import numpy as np
import matplotlib.pyplot as plt
import material as m
from geometry_lib import *
from source_lib import *
import yaml


with open('input_files/parametric_input.yaml', 'r') as f:
    input_data = yaml.safe_load(f)

################### Universe that contains the ANVIL machine #####################

universe_machine = openmc.Universe(786)

################# Defining the room inside which the machine exists #####################
bounding_surface = openmc.model.RectangularParallelepiped(xmin=-400, xmax=400, ymin=-400, ymax=400, zmin=-1350, zmax=1350, boundary_type = 'vacuum')
room_region = -bounding_surface
room_cell = openmc.Cell(100, region=room_region, fill = m.air)

universe_machine.add_cells([room_cell])

#________________________________________VACUUM VESSEL  CONSTRUCTION ________________________________________#

#--------------------- Innermost vacuum region ---------------------#

vacuum_vessel = input_data.get("vacuum_vessel",{})

first_vv_plane_from_midplane = vacuum_vessel.get("first_vv_plane_from_midplane")
machine_length_from_midplane = vacuum_vessel.get("machine_length_from_midplane")
vv_conical_angle = vacuum_vessel.get("cone_angle")

vacuum_chamber = vacuum_vessel.get("vacuum_chamber", {})
vacuum_chamber_radius = vacuum_chamber.get("radius")
vacuum_chamber_material = vacuum_chamber.get("material")


bottleneck_cylinder = vacuum_vessel.get("bottleneck_cylinder", {})
bottleneck_cylinder_radius = bottleneck_cylinder.get("cylinder_radius")
bottleneck_plane_distance = bottleneck_cylinder.get("plane_distance")




vv_layer_innermost_region = vacuum_vessel_region(first_vv_plane_from_midplane, machine_length_from_midplane, vacuum_chamber_radius, bottleneck_cylinder_radius, vv_conical_angle) # & ~ slanted_ports
vv_layer_innermost_cell = openmc.Cell(1000, region=vv_layer_innermost_region, fill= m.vacuum)

universe_machine.add_cells([vv_layer_innermost_cell])

room_region &= ~vv_layer_innermost_region

# -------------------- Computed Radii ------------------------#

structure = vacuum_vessel.get("structure", {})

structural_thicknesses = [properties.get("thickness") for layer_name, properties in structure.items()]
structural_materials = [getattr(m, properties.get("material")) for layer_name, properties in structure.items()]

vacuum_vessel_cylinders_radii = np.cumsum([vacuum_chamber_radius] + structural_thicknesses)
bottleneck_cylinders_radii = np.cumsum([bottleneck_cylinder_radius] + structural_thicknesses)

vacuum_vessel_materials = [vacuum_chamber_material] + structural_materials

vacuum_section_regions = [vv_layer_innermost_region]

vacuum_section_cells = []

for i in range(1, len(vacuum_vessel_cylinders_radii)):
    vv_layer = vacuum_vessel_region(first_vv_plane_from_midplane, machine_length_from_midplane, vacuum_vessel_cylinders_radii[i], bottleneck_cylinders_radii[i], vv_conical_angle)
    vacuum_section_regions.append(vv_layer)
    add_region = vv_layer & ~vacuum_section_regions[i-1] # & ~slanted_ports
    room_region &= ~add_region
    add_cell = openmc.Cell(1000+i, region=add_region, fill = vacuum_vessel_materials[i])
    vacuum_section_cells.append(add_cell)
    
universe_machine.add_cells(vacuum_section_cells)

############### Central Cell Construction ####################

# ________________________________________ CENTRAL CELL CONSTRUCTION ________________________________________ #

# Define the axial length of the central cell (half-length from the midplane)

central_cell = input_data.get("central_cell", {})

central_cell_axial_length = central_cell.get("axial_length")

central_cell_length_from_midplane = central_cell_axial_length / 2

# Load predefined central cell layers (thickness, material) from parameters

layers_cc = [np.array([layer["thickness"], getattr(m, layer["material"])]) for layer in central_cell.get("layers", [])]
central_cell_layers = np.array(layers_cc)

# Define axial boundary planes for the central cell
left_plane_central_cell = openmc.ZPlane(-central_cell_length_from_midplane)
right_plane_central_cell = openmc.ZPlane(central_cell_length_from_midplane)

# Initialize list to store radial cylinders of the central cell
central_cell_cylinders = []

# Compute the radii of central cell cylinders by accumulating layer thicknesses
central_cell_cylinders_radii = vacuum_vessel_cylinders_radii[-1] + np.cumsum([central_cell_layers[:, 0]])

# Create ZCylinders for each layer in the central cell
for i in range(len(central_cell_cylinders_radii)):
    cylinder = openmc.ZCylinder(r=central_cell_cylinders_radii[i])
    central_cell_cylinders.append(cylinder)

# Create the first cell for the central cell (innermost region)
central_cell_cells = [
    openmc.Cell(
        2000,
        region=(-right_plane_central_cell & +left_plane_central_cell & ~vacuum_section_regions[-1] & -central_cell_cylinders[0]),
                #& ~vertical_port_region & ~slanted_ports,
        fill=central_cell_layers[0][1]  # Assign material from layer data
    )
]

# Update room region to exclude the central cell
room_region &= ~(-right_plane_central_cell & +left_plane_central_cell & ~vacuum_section_regions[-1] & -central_cell_cylinders[0])

# Iterate through layers to define cylindrical shell regions and their corresponding cells
for i in range(1, len(central_cell_cylinders_radii)):
    cyl_region = (-right_plane_central_cell
                  & +left_plane_central_cell
                  & -central_cell_cylinders[i]
                  & +central_cell_cylinders[i - 1]) #& ~vertical_port_region & ~slanted_ports
    
    new_cell = openmc.Cell(2000 + i, region=cyl_region, fill=central_cell_layers[i][1])
    
    # Update room region to exclude the newly added cylindrical region
    room_region &= ~cyl_region
    
    central_cell_cells.append(new_cell)

# Add all central cell regions to the machine universe
universe_machine.add_cells(central_cell_cells)


# ________________________________________ LF COILS ________________________________________ #


lf_coil = input_data.get("lf_coil", [])

lf_coil_centers = lf_coil.get("positions",[])
lf_coil_inner_dimensions = lf_coil.get("inner_dimensions", [])

lf_coil_shell_thicknesses = lf_coil.get("shell_thicknesses", {})


lf_coil_regions = []
lf_coil_cells = []
lf_coil_shield_regions = []
lf_coil_shield_cells = []

for i in range(len(lf_coil_centers)):
    
    # Generate the cylindrical shell and inner region for the LF coil
    # Parameters include center position, outer radius reference, thickness, and axial dimensions
    lf_coil_shell, lf_coil_inner = hollow_cylinder_with_shell(
        lf_coil_centers[i],      # Coil center position
        central_cell_cylinders_radii[-1],  # Outer reference radius
        lf_coil_inner_dimensions["radial_thickness"],    # Radial thickness of the coil
        lf_coil_inner_dimensions["axial_length"],  # Inner axial length
        lf_coil_shell_thicknesses["front"],  # Shell front thickness
        lf_coil_shell_thicknesses["back"],   # Shell back thickness
        lf_coil_shell_thicknesses["axial"]   # Shell axial thickness
    )

    # Define OpenMC cells for the coil shell (structural) and coil winding pack (inner)
    lf_coil_shell_cell = openmc.Cell(
        4000 + i, 
        region=lf_coil_shell, 
        fill=getattr(m, lf_coil["materials"]["shield"])  # Material for the coil shell
    )

    lf_coil_shield_regions.append(lf_coil_shell)
    lf_coil_shield_cells.append(lf_coil_shell_cell)

    lf_coil_inner_cell = openmc.Cell(
        4100 + i, 
        region=lf_coil_inner, 
        fill=getattr(m, lf_coil["materials"]["magnet"])  # Material for the coil winding pack
    )

    lf_coil_regions.append(lf_coil_inner)
    lf_coil_cells.append(lf_coil_inner_cell)

    # Add the coil shell and inner region to the universe
    universe_machine.add_cells([lf_coil_shell_cell, lf_coil_inner_cell])
    
    # Update room region to exclude the space occupied by the coil shell and inner region
    room_region &= ~lf_coil_shell & ~lf_coil_inner

# ________________________________________ HF COILS ________________________________________ #

# --- Determine the HF Coil Center Position ---
# The HF coil is positioned along the Z-axis, starting from the central cell length,
# with additional spacing for the shield, casing, and half the coil thickness.
    
hf_coil = input_data.get("hf_coil", [])

hf_coil_shield = hf_coil.get("shield",[])
hf_coil_magnet = hf_coil.get("magnet", [])

casing_layers_thickness  = np.array([layer["thickness"] for layer in hf_coil.get("casing_layers", [])])
casing_layers_materials = np.array([getattr(m, layer["material"]) for layer in hf_coil.get("casing_layers",[])])

hf_coil_center_z0 = (
    central_cell_length_from_midplane  
    + hf_coil_shield["shield_central_cell_gap"]           # Gap between HF shield and central cell
    + hf_coil_shield["axial_thickness"][0]    # Axial thickness of the shield (toward midplane)
    + np.sum(casing_layers_thickness)   # Total casing thickness
    + hf_coil_magnet["axial_thickness"] / 2          # Half of the magnet's axial thickness
)

# --- Compute Inner Radius for Magnet and Casing ---
# The inner radius is determined by the bottleneck cylinder radius, the shield thickness,
# and an additional radial gap before the casing.
inner_radius_magnet_casings = (
    bottleneck_cylinders_radii[-1]               # Outermost bottleneck cylinder radius
    + hf_coil_shield["radial_thickness"][0]   # Shield thickness (toward central axis)
    + hf_coil_shield["radial_gap_before_casing"]            # Extra spacing before the casing begins
)

# --- Define Regions for Left and Right HF Coils ---
# The `nested_cylindrical_shells` function generates multiple concentric shell regions,
# representing the magnet and its casing layers.
shield_regions_right = nested_cylindrical_shells(
    hf_coil_center_z0, 
    inner_radius_magnet_casings, 
    hf_coil_magnet["radial_thickness"], 
    hf_coil_magnet["axial_thickness"], 
    casing_layers_thickness, casing_layers_thickness, casing_layers_thickness
)

shield_regions_left = nested_cylindrical_shells(
    -hf_coil_center_z0, 
    inner_radius_magnet_casings, 
    hf_coil_magnet["radial_thickness"], 
    hf_coil_magnet["axial_thickness"], 
    casing_layers_thickness, casing_layers_thickness, casing_layers_thickness
)

# --- Define the HF Magnet Region ---
# The magnet region consists of the **innermost layer** of the left and right coil casings.
hf_magnet_region_left = shield_regions_left[0]
hf_magnet_region_right = shield_regions_right[0]
hf_magnet_left_cell = openmc.Cell(6101, region=hf_magnet_region_left, fill=getattr(m, hf_coil_magnet["material"]))
hf_magnet_right_cell = openmc.Cell(6102, region=hf_magnet_region_right, fill=getattr(m, hf_coil_magnet["material"]))

# Add the HF magnet to the simulation
universe_machine.add_cells([hf_magnet_left_cell, hf_magnet_right_cell])
room_region &= ~ (hf_magnet_region_left | hf_magnet_region_right)  # Exclude magnet region from the room

# --- Define the Casing Regions (Excluding Shield for Now) ---
# Each subsequent layer in `shield_regions_left` and `shield_regions_right` represents casing materials.
# We iterate over these layers, assigning them the corresponding casing material.
for i in range(1, len(shield_regions_left)):  
    combined_region = shield_regions_left[i] | shield_regions_right[i]
    combined_region_cell = openmc.Cell(6500 + i, region=combined_region, fill=casing_layers_materials[i - 1])

    universe_machine.add_cells([combined_region_cell])
    room_region &= ~combined_region  # Exclude each casing layer from the room

# --- Compute the Inner and Outer Dimensions of the Shield ---
# The outermost shell must encapsulate the coil, including its casing layers.

# Total radial thickness of the HF outermost shell (coil + casing layers)
inner_radial_thickness_hf_outermost_shell = (
    hf_coil_magnet["radial_thickness"] 
    + 2 * np.sum(casing_layers_thickness)  # Accounts for both sides of the casing
)

# Total axial thickness of the HF outermost shell (coil + casing layers)
inner_axial_thickness_hf_outermost_shell = (
    hf_coil_magnet["axial_thickness"] 
    + 2 * np.sum(casing_layers_thickness)  # Accounts for both sides of the casing
)

# Compute the inner radius of the HF main shield (before applying the shell thickness)
inner_radius_hf_main_shield = (
    bottleneck_cylinders_radii[-1]  # Outermost bottleneck radius
    + hf_coil_shield["radial_gap_before_casing"]  # Gap before casing begins
)

# --- Define the HF Main Shield Region ---
# The shield consists of two hollow cylindrical shells, one for each HF coil.
# The function `hollow_cylinder_with_shell` defines a cylindrical shell with an inner region.

hf_main_shield_right_region = hollow_cylinder_with_shell(
        hf_coil_center_z0,                      # Z-position (right side)
        inner_radius_hf_main_shield,            # Inner radius
        inner_radial_thickness_hf_outermost_shell,  # Radial thickness of the outermost shell
        inner_axial_thickness_hf_outermost_shell,  # Axial thickness of the outermost shell
        hf_coil_shield["radial_thickness"][0],  # Inner shell thickness (towards axis)
        hf_coil_shield["radial_thickness"][1],  # Outer shell thickness (away from axis)
        hf_coil_shield["axial_thickness"][0]    # Axial thickness of the shield (towards midplane)
    )[0] 
hf_main_shield_left_region = hollow_cylinder_with_shell(
        -hf_coil_center_z0,                     # Z-position (left side)
        inner_radius_hf_main_shield,            # Inner radius
        inner_radial_thickness_hf_outermost_shell,  # Radial thickness of the outermost shell
        inner_axial_thickness_hf_outermost_shell,  # Axial thickness of the outermost shell
        hf_coil_shield["radial_thickness"][0],  # Inner shell thickness (towards axis)
        hf_coil_shield["radial_thickness"][1],  # Outer shell thickness (away from axis)
        hf_coil_shield["axial_thickness"][0]    # Axial thickness of the shield (towards midplane)
    )[0]


# --- Create the HF Main Shield Cell ---
hf_main_shield_right_cell = openmc.Cell(
    6201, 
    region=hf_main_shield_right_region, 
    fill=getattr(m, hf_coil_shield["material"])   # Assign material to the shield
)

hf_main_shield_left_cell = openmc.Cell(
    6202, 
    region=hf_main_shield_left_region, 
    fill=getattr(m, hf_coil_shield["material"])  # Assign material to the shield
)

# Add the HF main shield cell to the simulation
universe_machine.add_cells([hf_main_shield_left_cell, hf_main_shield_right_cell])
room_region &= ~ (hf_main_shield_left_region | hf_main_shield_right_region)  # Exclude main shield region from the room

# ________________________________________ END CELLS ________________________________________ #

# --- Define End Cell Dimensions ---
# These parameters define the structure of the end cells, which cap the HF coils.

# Access the `end_cell` section
end_cell = input_data.get("end_cell", {})



end_cell_axial_length = end_cell.get("axial_length")  # Total axial length of the end cell
end_cell_shell_thickness = end_cell.get("shell_thickness")  # Shell thickness of the end cell
end_cell_diameter = end_cell.get("diameter")  # Outer diameter of the end cell
end_cell_radial_thickness = end_cell_diameter / 2  # Convert diameter to radius

# --- Compute Z-Positions for End Cells ---
# The right end cell is positioned beyond the HF coil, accounting for:
# 1. The axial thickness of the HF coil shield
# 2. The total casing thickness
# 3. Half of the HF coil's axial thickness
# 4. The end cell's shell thickness
# 5. Half of the end cell's axial length

right_end_cell_z0 = (
    hf_coil_center_z0 
    + hf_coil_shield["axial_thickness"][0]  # Axial thickness of the shield (toward midplane)
    + np.sum(casing_layers_thickness)  # Total casing thickness
    + hf_coil_magnet["axial_thickness"] / 2  # Half of the magnet's axial thickness
    + end_cell_shell_thickness  # Shell thickness of the end cell
    + end_cell_axial_length / 2  # Half of the end cell's axial length
)

# The left end cell is symmetrically positioned at the negative Z-coordinate.
left_end_cell_z0 = -right_end_cell_z0

# --- Define End Cell Shells and Inner Volumes ---
# Each end cell consists of an outer shell and an inner vacuum region.
# The `cylinder_with_shell` function generates both the shell and inner volume.

left_end_cell_shell, left_end_cell_inner = cylinder_with_shell(
    left_end_cell_z0, 
    end_cell_radial_thickness, 
    end_cell_axial_length, 
    end_cell_shell_thickness
)

right_end_cell_shell, right_end_cell_inner = cylinder_with_shell(
    right_end_cell_z0, 
    end_cell_radial_thickness, 
    end_cell_axial_length, 
    end_cell_shell_thickness
)

# --- Exclude End Cell Shells from the Vacuum Section ---
# This ensures that the end cell shells do not overlap with the last vacuum section.
left_end_cell_shell &= ~vacuum_section_regions[-1]
right_end_cell_shell &= ~vacuum_section_regions[-1]

# --- Create OpenMC Cells for End Cells (Left Side) ---
end_cell_left_shell_cell = openmc.Cell(
    5001, 
    region=left_end_cell_shell, 
    fill=getattr(m, end_cell["shell_material"])  # Material for the shell
)

end_cell_left_inner_cell = openmc.Cell(
    5002, 
    region=left_end_cell_inner, 
    fill=getattr(m, end_cell["inner_material"])  # Material for inside the end cell
)

# Add left end cell components to the simulation
universe_machine.add_cells([end_cell_left_shell_cell, end_cell_left_inner_cell])

# --- Create OpenMC Cells for End Cells (Right Side) ---
end_cell_right_shell_cell = openmc.Cell(
    5003, 
    region=right_end_cell_shell, 
    fill=getattr(m, end_cell["shell_material"])  # Material for the shell
)

end_cell_right_inner_cell = openmc.Cell(
    5004, 
    region=right_end_cell_inner, 
    fill=getattr(m, end_cell["inner_material"])  # Material for inside the end cell
)

# Add right end cell components to the simulation
universe_machine.add_cells([end_cell_right_shell_cell, end_cell_right_inner_cell])

# --- Update Room Region to Exclude End Cells ---
# This ensures that the room region does not include the end cell volumes.
room_region &= ~left_end_cell_shell & ~left_end_cell_inner
room_region &= ~right_end_cell_shell & ~right_end_cell_inner


# ________________________________________ FINALIZE MODEL CROSS-SECTION ________________________________________ #

# Create the OpenMC geometry with the machine universe as the root
geometry = openmc.Geometry([openmc.Cell(fill=universe_machine)])

# Generate and save the cross-sectional plot of the model
geometry.root_universe.plot(
    basis='xz', 
    width=(1000, 2800), 
    pixels=(700, 700), 
    color_by='material', 
    #openmc_exec='/opt/openmc/bin/openmc'  # Specify the OpenMC executable path
)

# Save the generated plot to the results directory
plt.savefig('simple_mirror_cross_section.png', bbox_inches="tight")

# ________________________________________ EXPORT GEOMETRY TO XML ________________________________________ #

# Export the finalized geometry to an OpenMC XML file
geometry.export_to_xml("geometry.xml")


################### Settings #####################

# ________________________________________ Get the source from Settings information file ________________________________________ #

with open('input_files/source_information.yaml', 'r') as f:
    source_data = yaml.safe_load(f)

openmc_source = load_source_from_yaml('input_files/source_information.yaml')

settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.particles = int(source_data['settings']['particles_per_batch'])
settings.batches = source_data['settings']['batches']
settings.output = {'tallies': False}
settings.statepoint = {
    'batches': [1] + list(range(source_data['settings']['statepoint_frequency'], settings.batches, source_data['settings']['statepoint_frequency'])) + [settings.batches]
}
settings.weight_windows_on = source_data['settings']['weight_windows']
settings.weight_window_checkpoints = {'collision': True, 'surface': True}
wwg = openmc.WeightWindowGenerator(openmc.RegularMesh(), 
                                   [0, 14e6], 
                                   'neutron', 
                                   'magic', 
                                   max_realizations=25, 
                                   update_interval=1, 
                                   on_the_fly=True)
settings.weight_windows_generator = [wwg]
settings.photon_transport = source_data['settings']['photon_transport']
settings.source = openmc_source

# Export the finalized geometry to an OpenMC XML file
settings.export_to_xml("settings.xml")

# ################### Tallies ##################################


# hf_coil_cell = hf_magnet_region_right






