# paratan

**Paratan** is a lightweight OpenMC wrapper designed for rapid parametric modeling of tandem mirror fusion devices. It simplifies geometry construction, material assignment, and source specification, enabling easy setup and iteration of neutronics simulations. Built for early-stage design and analysis, **Paratan** helps researchers and engineers explore design spaces efficiently.

# 🐳 Setting Up the Environment with Docker

This project includes a Dockerfile to ensure a consistent and reproducible development environment across different systems.

## 🔧 Step 1: Build the Docker Image

Run the following command from the root directory of the repository (where the Dockerfile is located):

    sudo docker build -t openmc-paratan .

This builds a Docker image named openmc-paratan using the instructions provided in the Dockerfile.

## 🚀 Step 2: Run the Docker Container

After the image is built, you can launch the container with:

    sudo docker run -v /local_paratan_directory/:/local_workspace -it openmc-paratan

The -v flag mounts the project directory into the container so that changes you make on your host system are reflected inside the container.

The -it flag runs the container in interactive mode with a terminal.

# Model Setup

The input files are written in YAML format for easy readability, structured organization, and seamless parsing within the codebase.

## Source Configuration

The source input file defines how neutrons are introduced into the system. It supports multiple source types and simulation settings:

### Source Types (choose one via the type key):

**Uniform** – A constant neutron generation rate within a cylindrical region. Requires:

length (cm): axial extent of the source

radius (cm): radial extent of the source

**1D_Varying** – A spatially varying source along one dimension (typically axial). Requires:

file_name: CSV file specifying the axial distribution

radius (cm): radius of the cylindrical source

**2D_Varying** – A spatially varying source in two dimensions (radial and axial). Requires:

file_name: CSV file specifying the 2D distribution

### Common Parameter:

power_output – Total power output of the source in megawatts (MW): Used for postprocessing

### Simulation Settings (under settings):

batches: Number of OpenMC batches to run

particles_per_batch: Number of particles per batch

weight_windows: Enable/disable variance reduction via weight windows

statepoint_frequency: Frequency (in batches) to write statepoint files

photon_transport: Enable photon transport mode






    
