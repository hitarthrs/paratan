# NEWTON Magnetic Mirror Source Model

This repository builds the NEWTON magnetic mirror neutron source and the corresponding ParaTAN/OpenMC model from a YAML input file.

The main run command performs the source model calculation, generates the correlated OpenMC neutron source, and exports the OpenMC model XML files. It does **not** start an OpenMC transport calculation.

## Requirements

The code requires:

- Python 3.10 or newer
- OpenMC with working Python bindings
- NumPy
- SciPy
- PyYAML
- h5py
- Matplotlib

From the repository root, install the Python package and dependencies with:

```bash
python -m pip install -e .
```

Verify that OpenMC is available in the same environment:

```bash
python -c "import openmc; print(openmc.__version__)"
```

If OpenMC transport will be run, configure a compatible nuclear data library before starting OpenMC. For example:

```bash
export OPENMC_CROSS_SECTIONS=/path/to/cross_sections.xml
```

## Running the NEWTON model

The primary NEWTON input is:

```text
input_files/NEWTON_base_source.yaml
```

The OpenMC run settings are:

```text
input_files/source_information.yaml
```

A standard model build command is:

```bash
python scripts/run_model.py     --input input_files/NEWTON_base_source.yaml     --source-info input_files/source_information.yaml     --output-dir desired/output_path/
```

The output directory can be changed to any desired run location.

During a normal run, the console reports the electron temperature trials, beam density iterations, final numerical convergence state, neutron rate, OpenMC source location, and output directories.

A completed process is not by itself sufficient to establish that the numerical solution converged. For a production result, confirm that the final console output reports:

```text
Numerical convergence: passed
```

and check the saved run summary and metadata.

## Useful run options

Generate the ParaTAN geometry cross section during the model build:

```bash
python scripts/run_model.py     --input input_files/NEWTON_base_source.yaml     --source-info input_files/source_information.yaml     --output-dir runs/NEWTON_run     --plot-geometry
```

Disable iterative source model progress messages:

```bash
python scripts/run_model.py     --input input_files/NEWTON_base_source.yaml     --source-info input_files/source_information.yaml     --output-dir runs/NEWTON_run     --no-source-model-progress
```

Show the more detailed assessment and runtime output:

```bash
python scripts/run_model.py     --input input_files/NEWTON_base_source.yaml     --source-info input_files/source_information.yaml     --output-dir runs/NEWTON_run     --verbose-source-model-output
```

## Main inputs

`input_files/NEWTON_base_source.yaml` contains the ParaTAN device geometry and the magnetic mirror source model configuration, including the magnetic field targets, plasma state, neutral beams, numerical settings, fusion calculation, neutron calculation, and OpenMC source export settings.

`input_files/source_information.yaml` contains OpenMC run settings such as particles per batch, number of batches, statepoint frequency, photon transport, and tally settings.

`input_files/README_inputs.md` gives detailed descriptions of the items within a .yaml input and their responsibilities.

## Run outputs

A successful run writes its primary results under the selected run directory. The important files are:

```text
runs/NEWTON_run/
    input_used.yaml
    resolved_input.yaml
    run_summary.json
    source_model_metadata.json
    source_model_arrays.npz
    correlated_neutron_events.npz

    <configured OpenMC source directory>/
        source.h5
        source_file_validity.json
        source_sampling_metadata.json

    model_xml_files/
        geometry.xml
        materials.xml
        settings.xml
        tallies.xml
        source.h5
```

The main outputs are:

- `input_used.yaml` is a copy of the submitted input
- `resolved_input.yaml` records the fully parsed source model configuration
- `run_summary.json` contains the compact principal run results
- `source_model_metadata.json` contains the full source model metadata and numerical qualification information
- `source_model_arrays.npz` stores larger arrays referenced by the metadata
- `correlated_neutron_events.npz` stores the correlated neutron event bank used for detailed analysis and plotting
- `source.h5` is the correlated OpenMC FileSource
- `model_xml_files/` is the exported OpenMC model package

The `source.h5` copy inside `model_xml_files/` is packaged with the XML model so that `settings.xml` can reference the source locally.

## Generating source model plots

Generate the standard diagnostic plots for an existing run with:

```bash
python scripts/plot_source_model_diagnostics.py     --run-dir runs/NEWTON_run     --groups all
```

By default, plots are written to:

```text
runs/NEWTON_run/plots/
```

Individual plot groups can be requested by replacing `all` with the desired group names.

## Running OpenMC transport

The model build exports an OpenMC-ready model under:

```text
runs/NEWTON_run/model_xml_files/
```

After configuring the OpenMC nuclear data library, transport can be started from that directory:

```bash
cd runs/NEWTON_run/model_xml_files
openmc
```

The packaged `settings.xml` references the local `source.h5` file in the same directory.

## Repository layout

```text
input_files/                  Model input files
scripts/                      Main run and plotting entry points
src/mirror_paratan/           ParaTAN geometry, materials, tallies, and OpenMC model construction
src/source_model_revamp/      Magnetic mirror source model
tools/                        Diagnostic tools and nuclear data
```
