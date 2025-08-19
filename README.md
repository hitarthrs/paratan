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

# 🖥️ Running Simulations with CLI

Paratan provides a command-line interface (CLI) for running neutronics simulations. The CLI supports two modes: **simple** and **tandem**.

## Basic Usage

```bash
python cli.py --mode <mode> --input <input_file> --output-dir <output_directory>
```

### Parameters:
- `--mode`: Choose between "simple", "tandem", or "build-only" simulation modes
- `--input`: Path to the YAML input file
- `--output-dir`: Directory to write output files (default: "output")

## Available Modes

### Simple Mode
```bash
python cli.py --mode simple --input input_files/simple_parametric_input.yaml --output-dir simple_output
```
- **Purpose**: Simple mirror fusion device simulation with full geometry construction
- **Input**: Uses `simple_parametric_input.yaml` for configuration
- **Features**: Vacuum vessel, central cell, LF/HF coils, end cells with configurable materials and dimensions
- **Status**: Fully functional - builds complete OpenMC models and runs simulations

### Tandem Mode
```bash
python cli.py --mode tandem --input input_files/tandem_parametric_input.yaml --output-dir tandem_output
```
- **Purpose**: ATandem mirror fusion device simulation
- **Input**: Uses `tandem_parametric_input.yaml` for configuration
- **Features**: Full geometry construction, material assignment, and source specification

### Build-Only Mode
```bash
python cli.py --mode build-only --input input_files/simple_parametric_input.yaml --output-dir model_output
```
- **Purpose**: Generate OpenMC XML files and model plot without running the simulation (Works for simple mirrors right now, tandem mode to be added)
- **Use case**: Model validation, geometry inspection, or preparation for manual OpenMC runs
- **Output**: Complete set of OpenMC input files (geometry.xml, materials.xml, tallies.xml, settings.xml)

## Input Files

The simulation configuration is defined in YAML files located in the `input_files/` directory:

- `simple_parametric_input.yaml`: Configuration for simple mirror simulations
- `tandem_parametric_input.yaml`: Configuration for tandem mirror simulations  
- `source_information.yaml`: Source configuration for neutron generation

## Example Usage

1. **Run a tandem simulation:**
   ```bash
   python cli.py --mode tandem --input input_files/tandem_parametric_input.yaml --output-dir my_tandem_run
   ```

2. **Run a simple simulation:**
   ```bash
   python cli.py --mode simple --input input_files/simple_parametric_input.yaml --output-dir my_simple_run
   ```

3. **Build model files only:**
   ```bash
   python cli.py --mode build-only --input input_files/simple_parametric_input.yaml --output-dir model_files
   ```

The simulation will create output files in the specified output directory, including OpenMC statepoint files and other results. The build-only mode generates the complete set of OpenMC input files without running the simulation.

# Model Setup

The input files are written in YAML format for easy readability, structured organization, and seamless parsing within the codebase.

## Source Configuration

The source input file defines how neutrons are introduced into the system. It supports multiple source types and simulation settings:

### Source Types

**Volumetric (Default)** – Automatically approximates the vacuum vessel shape using a series of cylinders with relative strengths set according to the vessel geometry. This is the recommended source type for most fusion device simulations.

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






    
