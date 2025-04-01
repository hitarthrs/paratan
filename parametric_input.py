import numpy as np
import material as m

def add_layer(array, layer_thickness, layer_material):
    """
    Adds a new layer to the array.

    Parameters:
    array (list): The array to which the layer will be added.
    layer_thickness (float): The thickness of the new layer.
    layer_material (OpenMC material): The material of the new layer.

    Returns:
    list: The updated array with the new layer added.
    """
    new_layer = np.array([layer_thickness, layer_material])
    array.append(new_layer)
    return array



# ________________________________________ VACUUM VESSEL INFORMATION ________________________________________ #

# --- Vessel Dimensions ---
machine_length_from_midplane = 597.5      # Distance from center plane to vessel end
first_vv_plane_from_midplane = 127.9675   # Distance of central cell end from mid-plane

# --- Central Vacuum Chamber ---
vacuum_chamber_radius = 75
vacuum_chamber_material = m.vacuum

# --- Bottleneck Geometry ---
bottleneck_plane_distance = 325           # Distance of bottleneck plane from mid-plane
bottleneck_cylinder_radius = 22.3         # Radius of bottleneck cylinder

# --- Cone Properties ---
angle_vv = 15                             # Cone angle (degrees)

# --- Structural Thicknesses ---
fw_vacuum_structure_thickness = 0.2       # First-wall thickness
fw_vacuum_structure_material = m.tungsten

vacuum_structure_thickness = 2.5          # Outer structural thickness
vacuum_structure_material = m.FW_FNSF


# ________________________________________ CENTRAL CELL INFORMATION ________________________________________ #

# Length of the central cell

central_cell_axial_length = 530

# List to store central cell layers
central_cell_layers = []

# --- Breeding Blanket ---
breeding_blanket_thickness = 77
breeding_blanket_material = m.vacuum
add_layer(central_cell_layers, breeding_blanket_thickness, breeding_blanket_material)

# --- Back Wall ---
back_wall_thickness = 2
back_wall_material = m.BW_FNSF
add_layer(central_cell_layers, back_wall_thickness, back_wall_material)

# --- Helium Manifold ---
helium_manifold_thickness = 6
helium_manifold_material = m.he_manifold
add_layer(central_cell_layers, helium_manifold_thickness, helium_manifold_material)

# --- High-Temperature Shield ---
HT_shield_face_plate_thickness = 2.5
HT_shield_filler_thickness = 30
HT_shield_back_plate_thickness = 2.5

HT_shield_filler_material = m.HT_Shield_filler
HT_shield_material = m.water_cooled_wc

add_layer(central_cell_layers, HT_shield_face_plate_thickness, HT_shield_material)
add_layer(central_cell_layers, HT_shield_filler_thickness, HT_shield_filler_material)
add_layer(central_cell_layers, HT_shield_back_plate_thickness, HT_shield_material)

# --- Additional Layers ---
# If you need to add more layers, use the `add_layer()` function.
# Example:
# additional_layer_thickness = 5  # Thickness in cm
# additional_layer_material = m.some_material  # Replace with desired material
# add_layer(central_cell_layers, additional_layer_thickness, additional_layer_material)

# ________________________________________ LF COIL PARAMETERS ________________________________________ #

# --- Structural Thicknesses ---
lf_coil_shell_front_thickness = 8.3  # Front shell thickness
lf_coil_shell_back_thickness = 4     # Back shell thickness
lf_coil_shell_axial_thickness = 8.3  # Axial thickness of the shell

# --- Coil Inner Dimensions ---
lf_coil_inner_axial_length = 25  # Inner axial length
lf_coil_radial_thickness = 20    # Radial thickness of the coil

# --- Coil Positions ---
magnet_coil_centers = [-95, 95]  # Positions of LF coil centers along the Z-axis

# --- Material Definitions ---
lf_coil_shield_material = m.stainless            # Material for LF coil shell (shield)
lf_coil_magnet_material = m.Magnet_Winding_Pack_2  # Material for coil winding pack

# ________________________________________ LF COIL TALLY FLAGS ________________________________________ #

lf_coil_tally_flags = {
    "neutron_flux": True,  # Total neutron flux
    "photon_flux": True,  # Total photon flux
    "neutron_fast_flux": True,     # Neutron flux above a certain energy (fast flux)
    "neutron_heating": False       # Nuclear heating in the LF coil
}

# ________________________________________ HF COIL PARAMETERS ________________________________________ #

# --- #--- HF Magnet Coil dimensions---
hf_coil_bore_radius = 60        # Bore radius of the HF coil
hf_coil_radial_thickness = 170   # Radial thickness of the HF coil
hf_coil_axial_thickness = 250    # Axial thickness of the HF coil
hf_coil_material = m.Magnet_Winding_Pack_2 # Material used for the HF coil

# --- # --- HF Magnet Casing and Shield---

# --- HF Magnet Casing and Shield ---
hf_coil_casing_layers = {
    "vacuum_1": (1.4, m.vacuum),
    "aluminum_layer": (0.6, m.aluminum_1050),
    "vacuum_2": (1.4, m.vacuum),
    "stainless_steel": (0.6, m.ss316ln)
}

# Convert casing dictionary to NumPy arrays (Move this to the model making file)
hf_coil_casing_thicknesses = np.array([layer[0] for layer in hf_coil_casing_layers.values()])
hf_coil_casing_materials = np.array([layer[1] for layer in hf_coil_casing_layers.values()])


# --- # --- HF Coil Main Shield ---

# Radial gap between the
radial_gap_before_casing = 1

hf_shield_central_cell_gap = 5

hf_coil_shield_radial_thickness = np.array([30, 30])  # First value is thickness of shield towards the central axis, second is away
hf_coil_shield_axial_thickness = np.array([20, 10])   # First value is thickness of shield towards the midplane, second is away
hf_coil_main_shield_material = m.cooled_tungsten_boride

# __________________________________________END CELL PARAMETERS ________________________________________#

# End Cells

end_cell_axial_length = 600
end_cell_shell_thickness = 2.5
end_cell_diameter = 495.0

end_cell_shell_material = m.stainless












