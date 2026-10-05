# Pipeline stages

The pipeline_stages folder connects typed configuration to the source model physics modules. It builds geometry and beam deposition, closes the fixed temperature operating point, applies the electron temperature and terminal current closures, evaluates fusion and neutron sources, optionally writes the correlated OpenMC source, and assembles the final run metadata.

## Files

`geometry_stage.py`: Builds the analytic or ParaTAN linked magnetic field, confined and full device axial grids, flux tube volumes, magnetic throats, first material intersections, plasma boundaries, and startup or reference background geometry.

`beam_stage.py`: Resolves beam pitch and path geometry, selects the active stopping target density, evaluates atomic stopping rates, attenuates every enabled beam independently, and combines the resulting fast ion births by species.

`kinetic_stage.py`: Thin wrapper that forwards operating point continuation state into the modal FBIS integration stage.

`operating_point_stage.py`: Solves the beam attenuation and plasma density fixed point at one electron temperature. It checks pointwise, volume L2, and absolute density changes together with beam rate, power, axial source shape, and component deposition changes before rebuilding the candidate state for final consistency.

`ion_current_helpers.py`: Rebuilds the coupled kinetic state for explicit ion current components and checks whether current closure changed the electron density state used by beam attenuation.

`expander_stage.py`: Couples fixed magnetic boundary Eq 63 throat crossing rates to source connected expander populations and iterates terminal ion current with the shared Eq 68 barrier and Eq 70 potential state.

`electron_energy_balance_stage.py`: Assembles the represented electron energy terms from Eq 59 ion to electron transfer, prescribed external heating, electron wall kinetic loss, and ambipolar ion acceleration.

`temperature_closure_stage.py`: Selects fixed or self consistent electron temperature closure. The self consistent path repeatedly closes the beam, kinetic, terminal current, expander, and represented energy states before admitting a residual to the scalar temperature solve.

`fusion_stage.py`: Builds the full device fast D and fast T reactant registry and evaluates active D D and D T fusion components with the Bosch and Hale Table IV pair kernel.

`neutron_stage.py`: Builds deterministic axial neutron energy spectra from the same full device reactant populations and assembles the correlated neutron event bank when evaluated angular data are requested.

`neutron_event_stage.py`: Converts neutron producing fusion components into correlated event specifications and builds events with evaluated center of momentum angular laws and relativistic lab kinematics.

`openmc_export_stage.py`: Validates normalized event probabilities and absolute neutron rate conservation before optionally writing the correlated OpenMC HDF5 file source bundle.

`metadata_stage.py`: Filters and combines stage metadata, records missing physics availability, and defines the plotting and population grid ownership contract.

## Important conventions

The confined magnetic grid is used by the modal kinetic solve. The full device grid extends through both throats to the selected material boundaries and is used for mapped fast ion populations, fusion, and neutron source construction.

For startup seed closure, configured D and T profiles initialize the first beam attenuation target. Once a solved kinetic density state exists, the solved electron density becomes the stopping target and solved kinetic D and T distributions supply heavy particle atomic rates.

At fixed electron temperature, beam attenuation and the modal kinetic density state are iterated together. A candidate fixed point must be rebuilt using its own solved density target before final consistency is accepted.

The terminal current closure is evaluated after the throat level kinetic state. Species terminal target rates, the Eq 68 barrier, the Eq 70 profile, and the expander potential are iterated together before a final unrelaxed consistency rebuild.

Self consistent electron temperature closure admits a temperature residual only when the beam density fixed point, kinetic state, terminal current closure, expander state, and represented electron energy terms are evaluable. The accepted root is rebuilt without warm continuation for final qualification.

Fusion and neutron stages use the same full device reactant population authority. Fusion rates come directly from the Table IV pair kernel without posterior rate rescaling.

The deterministic neutron energy matrix uses the exact isotropic center of momentum energy marginal. Evaluated angular data enter the separate correlated event path, which retains position, energy, and direction correlations.

The correlated event bank stores normalized sampling probabilities separately from its physical neutron rate. OpenMC export requires both probability normalization and equality between the event bank physical rate and the neutron stage total rate.

## Physics basis

These files organize the physics modules rather than redefining their equations.

The modal kinetic and electrostatic stages follow the supplied Egedal et al. framework, including Eq 42 density weighting, Eq 59 velocity evolution, Eqs 61 through 64 ion loss relations, Eq 68 and Eq 69 end current and energy relations, and Eqs 70 through 72 quasineutral potential and accessibility mapping.

Neutral beam attenuation uses the supplied deuterium beam stopping basis through the beam module. Fusion uses the Bosch and Hale Table IV total cross section pair kernel. Neutron spectra and correlated events use the supplied arbitrary reactant relativistic kinematics framework, with evaluated center of momentum angular data applied by the event path.
