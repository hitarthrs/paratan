#!/usr/bin/python
#
import openmc
import numpy as np
import csv
import matplotlib.pyplot as plt

# create minimum functions needed to create 2D source from anvil_2D_V3_input.py

def csv_columns_to_arrays(file_name):
    # Initialize empty lists to store the first and second columns
    column1 = []
    column2 = []

    # Open and read the CSV file
    with open(file_name, 'r') as file:
        reader = csv.reader(file)
        
        # Iterate through each row in the CSV
        for row in reader:
            # Append the first and second column values to respective lists
            column1.append(row[0])
            column2.append(row[1])

    # Convert lists to arrays
    array1 = np.array(column1).astype(np.float64)
    array2 = np.array(column2).astype(np.float64)

    return array1, array2

def extract_coordinates_and_data(file_name):
    with open(file_name, 'r') as file:
        reader = csv.reader(file)
        data = list(reader)
    
    # Convert to NumPy array for easier slicing
    data_array = np.array(data)

    # Extract z-coordinates (first row, excluding the first element)
    z_coordinates = data_array[0, 1:].astype(np.float64)

    # Extract r-coordinates (first column, excluding the first element)
    r_coordinates = data_array[1:, 0].astype(np.float64)

    # Extract the main data matrix (excluding the first row and first column)
    data_matrix = data_array[1:, 1:].astype(np.float64)

    return z_coordinates*100, r_coordinates*100, data_matrix

# Extracting csv information
z_coords, r_coords, data_matrix = extract_coordinates_and_data('Anvil2Dsourcesheet1.csv')


def extend_and_plot_r_vs_z(r_coords, z_coords, data_matrix):
    # Mirror the z-coordinates
    z_negative = -z_coords[::-1]  # Create negative z-coordinates
    extended_z_coords = np.concatenate((z_negative, z_coords))  # Combine negative and positive z-coordinates
    
    # Mirror the data matrix
    mirrored_data_matrix = data_matrix[:, ::-1]  # Mirror the data matrix along the z-axis
    combined_data_matrix = np.concatenate((mirrored_data_matrix, data_matrix), axis=1)  # Combine the mirrored and original data matrices
    extended_data_matrix = combined_data_matrix/np.sum(combined_data_matrix)

    print(f"Probability: {np.sum(extended_data_matrix)}")

    # Plot the extended 2D data matrix
    plt.figure(figsize=(12, 6))
    plt.pcolormesh(extended_z_coords, r_coords, extended_data_matrix, shading='auto', cmap='viridis')
    plt.colorbar(label='Source strength')
    
    plt.xlabel('Z Coordinates (cm)')
    plt.ylabel('R Coordinates (cm)')
    plt.title('Extended 2D Source Strength Probability: R vs Z (Case: R3M2C2NB)')
    
    plt.savefig('plot2d-source_ext')
    plt.close()

    return r_coords, extended_z_coords, extended_data_matrix

####################################################################################

r_coords_2d, z_coords_2d, source_strength_2d = extend_and_plot_r_vs_z(r_coords, z_coords, data_matrix)

print(f'Total 2D source strength: {np.sum(source_strength_2d)}')

uniform_phi  = openmc.stats.Uniform(0, 2*np.pi)
sources_2d = []
for j in range(len(r_coords_2d)):
    for i in range(len(z_coords_2d)):
        source_term = openmc.IndependentSource()
        source_term.particle = 'neutron'
        source_term.energy = openmc.stats.Discrete([14.1e6], [1])
        source_term.angle = openmc.stats.Isotropic()
        r_discrete = openmc.stats.Discrete([r_coords_2d[j]], [1.0])
        z_discrete = openmc.stats.Discrete([z_coords_2d[i]], [1.0])
        source_term.space = openmc.stats.CylindricalIndependent(r_discrete, uniform_phi,z_discrete)
        source_term.strength = source_strength_2d[j][i]
        sources_2d.append(source_term)
#
print("\n All done! \n")
