import openmc
import yaml
import csv
import numpy as np
import pandas as pd

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
    plt.title('Extended 2D Source Strength Probability: R vs Z')
    
    plt.savefig('/res_dir/2d-source_ext')

    return r_coords, extended_z_coords, extended_data_matrix 

# Parent class representing a generic neutron source
class Source:
    def __init__(self, power_output):
        # Calculate neutron emission strength from given thermal power (MW)
        self.strength = power_output * 0.8 * 1e6 / (14.1 * 1e6 * 1.602e-19)

    def describe(self):
        raise NotImplementedError("This method should be implemented by subclasses.")


# Subclass for a uniform cylindrical neutron source
class UniformSource(Source):
    def __init__(self, power_output, length, radius):
        super().__init__(power_output)  # Initialize the parent class (calculates strength)
        self.length = length            # Axial length of the cylindrical source
        self.radius = radius            # Radial extent of the cylindrical source

    def create_openmc_source(self):
        """
        Creates and returns an OpenMC neutron source object with a uniform cylindrical spatial distribution
        and a monoenergetic neutron energy (14.1 MeV).
        """
        # Uniform axial (z) distribution along cylinder length
        z_distribution = openmc.stats.Uniform(-self.length / 2, self.length / 2)

        # Radial distribution following a power-law to ensure uniform distribution in the cylinder
        r_distribution = openmc.stats.PowerLaw(0, self.radius, -2)

        # Uniform angular distribution from 0 to 2π radians
        phi_distribution = openmc.stats.Uniform(0, 2 * np.pi)

        # Monoenergetic neutron energy distribution (typical fusion neutron energy)
        energy_distribution = openmc.stats.Discrete([14.1e6], [1])

        # Initialize OpenMC independent neutron source
        source = openmc.IndependentSource()
        source.particle = 'neutron'
        source.angle = openmc.stats.Isotropic()  # Isotropic angular emission
        source.energy = energy_distribution

        # Cylindrical spatial distribution (r, θ, z)
        source.space = openmc.stats.CylindricalIndependent(
            r_distribution,
            phi_distribution,
            z_distribution
        )

        return source
    
    def describe(self):
        print("Uniform Cylindrical Source:")
        print(f" Strength (neutrons/sec): {self.strength:.3e}")
        print(f" Length: {self.length} cm, Radius: {self.radius} cm")

    

    
# Subclass representing a 1D varying neutron source, based on distribution provided from a file
class Source1D(Source):
    def __init__(self, power_output, radius, file_name):
        super().__init__(power_output)  # Initialize parent class (calculates neutron strength)
        self.radius = radius            # Radial extent of the cylindrical source
        self.file_name = file_name      # Path to file containing source distribution data

    def create_openmc_source(self):
        """
        Creates a list of OpenMC neutron sources with spatial distribution based on a 1D profile
        provided in an external CSV file.
        """
        # Load data from CSV (user-defined function assumed to exist)
        source_z_position_m, one_d_source_distribution_pos = csv_columns_to_arrays(self.file_name)

        # Mirror data to generate symmetric distribution around z=0
        one_d_source_distribution_neg = one_d_source_distribution_pos[1:][::-1]
        source_z_position_m_neg = -source_z_position_m[1:][::-1]

        # Combine negative and positive sides
        full_distribution = np.concatenate((one_d_source_distribution_neg, one_d_source_distribution_pos))
        full_z_positions_cm = np.concatenate((source_z_position_m_neg, source_z_position_m)) * 100  # convert to cm

        # Normalize distribution to form probabilities
        source_probabilities = full_distribution / np.sum(full_distribution)

        # Monoenergetic neutron energy (14.1 MeV)
        energy_distribution = openmc.stats.Discrete([14.1e6], [1.0])

        # Uniform angular distribution around cylindrical axis
        phi_distribution = openmc.stats.Uniform(0, 2 * np.pi)

        # Radial distribution ensuring uniform radial emission
        r_distribution = openmc.stats.PowerLaw(0, self.radius, -2)

        sources_1d = []
        # Create individual discrete OpenMC sources for each axial position
        for z_pos, probability in zip(full_z_positions_cm, source_probabilities):
            source_term = openmc.IndependentSource()
            source_term.particle = 'neutron'
            source_term.angle = openmc.stats.Isotropic()
            source_term.strength = probability # scale source strength accordingly
            source_term.energy = energy_distribution

            # Set discrete axial (z) location
            z_discrete_distribution = openmc.stats.Discrete([z_pos], [1.0])
            source_term.space = openmc.stats.CylindricalIndependent(
                r_distribution, 
                phi_distribution, 
                z_discrete_distribution
            )

            sources_1d.append(source_term)

        return sources_1d
    
    def describe(self):
        print(f"1D Source Distribution from file: {self.file_name}")
        print(f" Strength (total neutrons/sec): {self.strength:.3e}")
        print(f" Radial extent: {self.radius} cm")
        print(f" Loaded source distribution from {self.file_name}")

# Subclass representing a 2D varying neutron source, based on distribution provided from a file
class Source2D(Source):
    def __init__(self, power_output, file_name):
        super().__init__(power_output)  # Initialize parent class (calculates neutron strength)
        self.file_name = file_name      # Path to file containing source distribution data

    def create_openmc_source(self):
        """
        Creates a list of OpenMC neutron sources with spatial distribution based on a 2D profile
        provided in an external CSV file.
        """
        # Extract spatial data (r, z) and source strength from the input file
        z_coords, r_coords, data_matrix = extract_coordinates_and_data(self.file_name)

        # Extend and process 2D source distribution data
        r_coords_2d, z_coords_2d, source_strength_2d = extend_and_plot_r_vs_z(r_coords, z_coords, data_matrix)

        # Monoenergetic neutron energy (14.1 MeV)
        energy_distribution = openmc.stats.Discrete([14.1e6], [1.0])

        # Uniform angular distribution around cylindrical axis
        phi_distribution = openmc.stats.Uniform(0, 2 * np.pi)

        sources_2d = []
        # Create individual discrete OpenMC sources for each spatial point
        for j in range(len(r_coords_2d)):
            for i in range(len(z_coords_2d)):
                z_pos = z_coords_2d[i]  # Axial position (z)
                r_pos = r_coords_2d[j]  # Radial position (r)
                probability = source_strength_2d[j][i]  # Source strength probability
                
                # Define an OpenMC independent neutron source
                source_term = openmc.IndependentSource()
                source_term.particle = 'neutron'
                source_term.angle = openmc.stats.Isotropic()
                source_term.strength = probability * self.strength  # Scale source strength
                source_term.energy = energy_distribution

                # Set discrete spatial distributions for (r, z)
                z_discrete_distribution = openmc.stats.Discrete([z_pos], [1.0])
                r_discrete_distribution = openmc.stats.Discrete([r_pos], [1.0])
                
                source_term.space = openmc.stats.CylindricalIndependent(
                    r_discrete_distribution, 
                    phi_distribution, 
                    z_discrete_distribution
                )
                
                sources_2d.append(source_term)

        return sources_2d
    
    def describe(self):
        """
        Prints a description of the 2D source, including its file and computed strength.
        """
        print("2D Source Distribution:")
        print(f" Strength (total neutrons/sec): {self.strength:.3e}")
        print(f" Loaded source distribution from {self.file_name}")


def load_source_from_yaml(yaml_file):
    """
    Reads a YAML file and creates the corresponding OpenMC source object using predefined subclasses.
    
    Parameters:
        yaml_file (str): Path to the YAML file containing source information.
    
    Returns:
        openmc.Source or list of openmc.Source: The generated OpenMC source(s).
    """
    # Load YAML file
    with open(yaml_file, 'r') as f:
        source_data = yaml.safe_load(f)
    
    # Extract general source information
    source_type = source_data['source']['type']
    power_output = source_data['source']['power_output']
    
    # Create source based on type using predefined classes
    if source_type == "Uniform":
        length = source_data['source']['uniform']['length']
        radius = source_data['source']['uniform']['radius']
        return UniformSource(power_output, length, radius).create_openmc_source()
    
    elif source_type == "1D_Varying":
        file_name = source_data['source']['source_1D']['file_name']
        radius = source_data['source']['source_1D']['radius']
        return Source1D(power_output, radius, file_name).create_openmc_source()
    
    elif source_type == "2D_Varying":
        file_name = source_data['source']['source_2D']['file_name']
        return Source2D(power_output, radius, file_name).create_openmc_source()
    
    else:
        raise ValueError(f"Unsupported source type: {source_type}")



