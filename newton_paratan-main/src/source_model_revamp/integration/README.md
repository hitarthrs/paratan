# Integration

The integration package connects parsed source model inputs to geometry, beam attenuation, modal FBIS closure, expander populations, fusion, neutron events, OpenMC export, validation, assessment, and run outputs. The direct files define the shared stage result types, final orchestration, validation and assessment rules, full device population mapping, metadata handling, output helpers, and reporting utilities used around the integration subpackages.

## Files

`pipeline.py`: Runs the complete source model stage sequence, applies direct state validation, constructs the final convergence and applicability assessment, reports represented population scope, assembles metadata, and records stage runtime.

`pipeline_types.py`: Defines the immutable stage result and continuation state dataclasses shared across the integration workflow. These include geometry, beam, kinetic, operating point, expander, fusion, neutron, OpenMC export, and complete pipeline results.

`assessment.py`: Converts typed stage metadata into canonical numerical convergence and model applicability records. A check can pass, fail, or be marked not applicable, and the final result status combines the two assessment groups.

`validation.py`: Applies direct consistency checks to calculated geometry, density, kinetic, electrostatic, current, expander, fusion, neutron, correlated event, and OpenMC source states. Invalid calculated states fail before numerical qualification is summarized.

`evaluation.py`: Defines operating point evaluation modes and optional progress events used by fixed temperature and self consistent electron temperature searches.

`full_device_populations.py`: Builds the common full device population grid, conservatively remaps confined cell averaged quantities by overlap volume, retains explicit full device fast ion distributions when available, and keeps the solved electron profile unavailable outside its confined support.

`nbi_supported_balance.py`: Evaluates separate D and T stationary particle balances using beam ionization, charge exchange transfer and target sink rates, and kinetic terminal losses.

`pipeline_common.py`: Provides shared positive profile checks, volume averaged source reduction, distribution roundoff handling, and source aligned speed grid construction.

`pipeline_reporting.py`: Reports represented fast ion populations, excluded populations, and explicit limitations for the selected plasma and expander closure state.

`metadata.py`: Converts metadata into JSON compatible values, archives large arrays in a companion NPZ file, requires canonical run keys, and writes the compact run summary.

`output_writer.py`: Writes the submitted input copy, resolved input, metadata and array outputs, correlated neutron event archive, and paths to any OpenMC source files created by the export stage.

`plasma_closure.py`: Provides the canonical plasma closure name and identifies when configured D and T profiles are initialization data only.

`progress.py`: Formats optional console progress for temperature trials, beam density iterations, kinetic solves, and final qualification.

`yaml_loader.py`: Loads YAML mappings and parses the `source_model` block into the typed run configuration.

## Subpackages

`config/`: Parses and validates the typed source model configuration, including geometry, beam, plasma, modal numerical controls, operating point coupling, electron temperature closure, fusion, neutron, export, and reporting settings.

`modal_stage/`: Connects geometry and beam sources to the modal FBIS solver, including Eq 42 basis coupling, density closure, coupled D and T operation, prompt source handling, full device fast ion mapping, and end loss power availability.

`pipeline_stages/`: Implements the stage level geometry, beam, operating point, terminal current, expander, electron energy, fusion, neutron, OpenMC export, and metadata orchestration used by `pipeline.py`.

## Important conventions

Direct validation and numerical qualification are separate. Geometry, density, kinetic, electrostatic, current, expander, neutron event, and OpenMC state validators reject internally inconsistent calculated states. `RunAssessment` then reports convergence and applicability from the surviving stage evidence.

The confined geometry covers the magnetic throat to throat kinetic domain. Full device geometry extends to selected material boundaries. Fast ion populations can be represented on the full device grid, while the solved electron profile remains limited to its confined Eq 70 scope and is returned as NaN outside that support.

Beam source arrays use `(n_z, n_speed, n_lambda)`. Local fast ion arrays use `(n_z, n_speed, n_pitch)`. Neutron rate matrices use `(n_z, n_energy)`. Units are encoded in field names, including m, T, J, keV, W, s⁻¹, and m⁻³.

For the NBI supported stationary closure, the species particle identity is `ionization + charge exchange gain − charge exchange target sink = terminal loss`. Cross isotope charge exchange transfers species inventory, while same isotope charge exchange changes the velocity distribution without changing total species inventory.

Correlated neutron event probabilities sum to one independently of the retained physical neutron rate in s⁻¹. Source constructibility and OpenMC validation check both the normalized event representation and the physical rate evidence.

Large metadata arrays are moved to a companion NPZ archive and referenced from JSON by filename and dataset. Small scalar and sequence diagnostics remain directly readable in the JSON record.

## Physics basis

The direct integration files organize and validate the physics implemented in the source model modules rather than redefining those equations.

The modal kinetic state follows the supplied Egedal et al. framework, including Eq 42 density weighting, Eq 59 fast ion velocity evolution, Eqs 61 through 64 ion losses, Eq 68 and Eq 69 end current and energy relations, and Eqs 70 through 72 quasineutral electrostatic reconstruction.

Neutral beam inputs and stopping targets feed the beam attenuation model based on the supplied deuterium beam stopping reference. Fusion stage outputs use the Bosch and Hale total cross section representation. Neutron spectra and correlated events use the supplied arbitrary reactant relativistic kinematics framework together with the evaluated center of momentum angular data handled by the nuclear data and neutron event modules.
