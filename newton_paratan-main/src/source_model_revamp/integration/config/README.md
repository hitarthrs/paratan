# Integration configuration

The config folder parses the `source_model` input block into immutable typed configuration objects used by the integration pipeline. It validates units, model names, numerical controls, geometry links, beam definitions, closure choices, downstream source settings, and cross block consistency before any physics stage is assembled.

## Files

`common.py`: Provides shared mapping, scalar, integer, boolean, array, coordinate, alias, and geometry layer parsing helpers.

`constants.py`: Defines canonical model identifiers and alias tables shared by the configuration blocks.

`geometry.py`: Parses magnetic field targets, reduced plasma geometry, ParaTAN linked device geometry, effective LF and HF coil fits, and plasma facing boundary selection controls.

`plasma.py`: Parses background D and T seed or reference densities, ion temperature, electron initial guesses, density closure choice, and axial background profile specification.

`beam.py`: Parses deuterium and tritium beam power, energy components, pitch integration controls, atomic data choice, and straight beamline geometry.

`kinetic.py`: Holds the modal FBIS, Eq 42 density weighting, Eq 59 velocity solve, Eq 70 electrostatic solve, fixed magnetic loss, current balance, expander, and convergence numerical controls.

`operating_point.py`: Defines the beam attenuation and density fixed point tolerances, relaxation, floors, and metadata representation.

`power.py`: Defines fixed or self consistent electron temperature closure controls, external electron heating, root bracket, and scalar energy balance tolerances.

`downstream.py`: Defines fusion pair integration, neutron energy and angular sampling, correlated event controls, radial source profile choice, and OpenMC file source output settings.

`output.py`: Defines optional reporting controls that do not alter the solved physics state.

`root.py`: Assembles the full configuration, checks consistency between blocks, selects enabled beams, and creates the resolved mapping written with run outputs.
