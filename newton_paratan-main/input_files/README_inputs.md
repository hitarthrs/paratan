# NEWTON Input Configuration Reference

`input_files/NEWTON_base_source.yaml` controls both the ParaTAN/OpenMC device geometry and the magnetic mirror source calculation.

Unless a field name specifies another unit, ParaTAN geometry dimensions are given in centimeters. Source model geometry uses the units written in the field names, including meters (`_m`), tesla (`_T`), watts (`_W`), kiloelectronvolts (`_keV`), and number density in m\(^{-3}\) (`_m3`).

## Vacuum vessel

| Option | Effect |
| --- | --- |
| `vacuum_vessel.geometry_style` | Selects the vessel geometry. `perpendicular` uses straight central and bottleneck cylinders with perpendicular transitions. `conical` uses the hourglass style conical transition geometry |
| `vacuum_vessel.outer_axial_length` | Axial geometry parameter used by the conical vessel model. It is ignored when `geometry_style` is `perpendicular` |
| `vacuum_vessel.central_axial_length` | Sets the total axial length of the central vacuum vessel section. |
| `vacuum_vessel.central_radius` | Sets the inner radius of the central vacuum region |
| `vacuum_vessel.bottleneck_radius` | Sets the inner radius of the narrower vacuum regions beyond the central section |
| `vacuum_vessel.left_bottleneck_length` | Sets the axial length of the left bottleneck region |
| `vacuum_vessel.right_bottleneck_length` | Sets the axial length of the right bottleneck region |
| `vacuum_vessel.axial_midplane` | Shifts the complete device geometry along the OpenMC `z` axis |
| `vacuum_vessel.structure` | Defines successive structural layers outside the vacuum region |
| `vacuum_vessel.structure.<layer>.thickness` | Sets the radial thickness of that structural layer |
| `vacuum_vessel.structure.<layer>.material` | Selects the material assigned to that structural layer |

## Central cell

| Option | Effect |
| --- | --- |
| `central_cell.axial_length` | Sets the axial length of the central blanket and shielding package. It also contributes to the axial placement of the HF coil assemblies |
| `central_cell.layers` | Defines successive radial material layers outside the vacuum vessel structure |
| `central_cell.layers[].thickness` | Sets the radial thickness of one central cell layer |
| `central_cell.layers[].material` | Selects the material for that layer |
| `central_cell.tallies` | Defines OpenMC tallies associated with central cell regions |
| `central_cell.tallies.breeder` | Associates the enclosed tally definitions with the first central cell layer, which is treated as the breeder region |
| `mesh_tallies[].scores` | Lists the OpenMC tally scores to record, such as `heating` |
| `mesh_tallies[].dimensions` | Sets cylindrical mesh resolution in `(nr, nphi, nz)` |

The supplied `[40, 1, 25]` mesh therefore uses 40 radial cells, one azimuthal cell, and 25 axial cells.

## Low field coils

| Option | Effect |
| --- | --- |
| `lf_coil.positions` | Gives the axial centers of the LF coils relative to `vacuum_vessel.axial_midplane` |
| `lf_coil.inner_dimensions.axial_length` | Sets the axial thickness of each LF magnet region |
| `lf_coil.inner_dimensions.radial_thickness` | Sets the radial thickness of each LF magnet region |
| `lf_coil.shell_thicknesses.front` | Sets the inner side radial shield thickness surrounding each LF magnet |
| `lf_coil.shell_thicknesses.back` | Sets the outer side radial shield thickness |
| `lf_coil.shell_thicknesses.axial` | Sets the axial shield thickness added at both ends |
| `lf_coil.materials.shield` | Selects the LF coil shield material |
| `lf_coil.materials.magnet` | Selects the LF magnet material |

The LF coil positions and dimensions also enter the effective coil representation used to construct the source model magnetic field.

## High field coils

| Option | Effect |
| --- | --- |
| `hf_coil.magnet.radial_thickness` | Sets the radial thickness of the HF magnet |
| `hf_coil.magnet.axial_thickness` | Sets the axial thickness of the HF magnet and contributes to the full device axial layout |
| `hf_coil.magnet.material` | Selects the HF magnet material |
| `hf_coil.casing_layers` | Defines ordered casing layers surrounding each HF magnet |
| `hf_coil.casing_layers[].thickness` | Sets the thickness of one casing layer |
| `hf_coil.casing_layers[].material` | Selects its material |
| `hf_coil.shield.shield_central_cell_gap` | Sets the axial separation between the central cell package and HF shielding. This changes HF coil position and therefore also affects the geometry linked magnetic field |
| `hf_coil.shield.radial_thickness[0]` | Sets the inner side radial HF shield thickness |
| `hf_coil.shield.radial_thickness[1]` | Sets the outer side radial HF shield thickness |
| `hf_coil.shield.axial_thickness[0]` | Sets the central facing axial shield thickness |
| `hf_coil.shield.axial_thickness[1]` | Sets the outward facing axial shield thickness |
| `hf_coil.shield.material` | Selects the HF shield material |

## End cells

| Option | Effect |
| --- | --- |
| `end_cell.axial_length` | Sets the axial length of each end cell interior |
| `end_cell.diameter` | Sets the inner diameter of the end cell volume |
| `end_cell.shell_thickness` | Sets the radial and end cap shell thickness |
| `end_cell.shell_material` | Selects the end cell shell material |
| `end_cell.inner_material` | Selects the material filling the end cell interior, normally vacuum |

The end cell positions are derived automatically from the HF coil geometry and the full ParaTAN axial layout.

## Source model selection

| Option | Effect |
| --- | --- |
| `source_model.model` | Selects the magnetic mirror source model pipeline. The production value `single_beam_egedal_modal_fbis_openmc_source` should normally be left unchanged |

## Magnetic source geometry

| Option | Effect |
| --- | --- |
| `source_model.geometry.model` | Selects how the source model magnetic field is generated. `paratan_coil_field` derives an effective on axis field from the ParaTAN LF and HF coil geometry |
| `source_model.geometry.domain_mode` | Selects how the source model axial domain is defined. `geometry_linked` derives it from the ParaTAN device geometry |
| `source_model.geometry.plasma_radius_m` | Sets the reference plasma flux tube radius at the midplane. It controls source model volume, beam plasma overlap, and radial source extent. It does not change the solid OpenMC vessel geometry |
| `source_model.geometry.axial_bins` | Sets the number of axial cells used across the confined source model domain. Increasing it improves axial resolution but increases computational cost |
| `source_model.geometry.field_strength.midplane_B_T` | Sets the requested on axis magnetic field at the midplane |
| `source_model.geometry.field_strength.mirror_throat_B_T` | Sets the requested mirror throat field. Together with the midplane field it determines the target mirror ratio |
| `source_model.geometry.coil_field_grid_points` | Sets the axial resolution used to evaluate and fit the effective ParaTAN coil field. Larger values give finer throat and field profile resolution at increased setup cost |

## Plasma facing boundaries

| Option | Effect |
| --- | --- |
| `source_model.plasma_boundaries.detection_mode` | Controls discovery of plasma facing material surfaces. Normal implementation uses `automatic_with_overrides` |
| `default_material_interface_role` | Sets the default particle interaction assigned to automatically detected material interfaces. Usual input uses `absorbing` |
| `default_electrical_model` | Sets the default electrical treatment of detected surfaces. Supported defaults are `floating` and `grounded` |
| `overrides` | Allows selected detected surfaces to receive different particle or electrical boundary conditions |

An override can select a surface or component and specify a `particle_role`, `electrical_model`, and, for a prescribed electrical bias, `prescribed_bias_V`.

Supported particle roles include `absorbing`, `collector`, `end_ring`, `reflecting`, and `excluded`.

## Plasma closure and startup state

| Option | Effect |
| --- | --- |
| `source_model.plasma_closure.model` | Selects the plasma density closure. `nbi_supported_stationary` is the default |
| `background_ions.deuterium_midplane_density_m3` | Sets the startup deuterium seed density at the midplane |
| `background_ions.tritium_midplane_density_m3` | Sets the startup tritium seed density at the midplane |
| `background_ions.ion_temperature_keV` | Sets the startup or reference thermal ion temperature used by the closure |
| `electrons.density_mode` | Selects how electron density is obtained. `quasineutral_from_total_positive_charge` derives it from quasineutrality |
| `electrons.temperature_initial_guess_keV` | Sets the initial electron temperature guess supplied to the self consistent temperature solution |

For `nbi_supported_stationary`, the configured D and T densities are initialization values. They are not imposed as the final stationary fast ion density.

## Neutral beams

Each entry under `source_model.beams` defines one neutral beam. The NEWTON_base input contains one deuterium and one tritium beam.

| Option | Effect |
| --- | --- |
| `id` | Unique identifier used to track the beam in metadata and component accounting |
| `enabled` | Enables or disables that beam |
| `model` | Selects the beam deposition model. `geometry_linked_attenuation` follows the configured physical beamline through the geometry linked plasma |
| `species` | Selects deuterium or tritium |
| `power_W` | Sets total injected neutral beam power |
| `energy_keV` | Sets the full energy component energy |
| `injection_angle_deg` | Reference injection pitch angle relative to the magnetic field direction |
| `injection_pitch_source` | Selects whether pitch is derived from the beamline geometry or deliberately taken from the configured angle. Production uses `beamline_endpoints` |
| `injection_angle_tolerance_deg` | Maximum allowed disagreement between `injection_angle_deg` and the angle derived from the beamline endpoints |
| `pitch_distribution_model` | Selects the treatment of finite beam pitch width. Production uses geometry linked finite width integration |
| `pitch_quadrature_points` | Sets the number of quadrature points used to resolve the finite beam pitch distribution |
| `components` | Defines full and fractional beam energy components |
| `components[].id` | Identifier for one energy component |
| `components[].power_fraction` | Fraction of total beam power assigned to that component. Fractions for a beam must sum to one |
| `components[].energy_multiplier` | Multiplies the parent `energy_keV` to obtain that component's energy |
| `atomic_data_model` | Selects the atomic stopping and charge exchange data model. Production uses `hydrogenic_dt_ground_state` |
| `channel` | Labels the physical source channel passed into the fast ion source construction. Production uses `fast_ion_birth` |
| `beamline.start_m` | Three dimensional beamline starting coordinate in meters |
| `beamline.end_m` | Three dimensional beamline ending coordinate in meters |
| `beamline.radius_m` | Initial physical beam radius |
| `beamline.divergence_half_angle_deg` | Half angle used to expand the beam radius along the path |
| `beamline.radial_overlap_model` | Controls how the finite beam footprint is intersected with the plasma flux tube. Production uses `hard_edge_overlap` |

Changing beam power, energy, component fractions, geometry, or species changes the deposited fast ion source and therefore propagates through the stationary plasma, fusion, and neutron calculations.

## Kinetic and electrostatic model

The `kinetic_electrostatic` block contains both physical model selections and numerical controls for the modal FBIS calculation. Model selection fields should normally remain at their production values. Resolution and tolerance fields can be changed for convergence studies.

### Physical model selections

| Option | Effect |
| --- | --- |
| `electrostatic_feedback_model` | Selects the coupling between the modal fast ion calculation and the self consistent electrostatic potential |
| `modal_eq42_density_weighting_model` | Selects the density weighting used in the orbit averaged scattering operator |
| `ion_loss_closure_model` | Selects the ion loss boundary treatment. Production uses the fixed magnetic boundary |
| `lost_ion_parallel_temperature_model` | Selects the closure used to obtain the lost ion parallel temperature at the throat |
| `modal_velocity_solution_model` | Selects the modal velocity space solve. Production uses the hot ion Rosenbluth Eq. 59 model |
| `expander_radial_profile_model` | Selects the prescribed radial weighting used when routing the source connected population through the expander |

### Main velocity space grids

| Option | Effect |
| --- | --- |
| `speed_bins` | Number of cells in the main fast ion speed grid |
| `lambda_bins` | Number of cells used for the magnetic invariant `Lambda` representation |
| `pitch_bins` | Number of pitch cells used for reconstructed local distributions |
| `speed_max_factor` | Sets the upper speed domain extent relative to the highest configured beam component speed |
| `modal_speed_core_cell_fraction` | Controls how the stretched speed grid distributes resolution between the source or core region and high speed tail |
| `modal_speed_tail_stretch_power` | Controls stretching of the high speed tail portion of the speed grid |

Increasing grid resolution generally improves numerical resolution but increases runtime and memory use.

### Modal basis construction

| Option | Effect |
| --- | --- |
| `modal_n_lambda_grid` | Resolution of the `Lambda` grid used to construct the physical modal basis |
| `modal_n_eta_grid` | Resolution of the transformed `eta` coordinate used by the modal basis |
| `modal_n_square_basis_modes` | Number of square mirror basis modes constructed before transformation to the physical basis |
| `modal_n_physical_modes` | Number of physical modes retained in the production kinetic solution |
| `modal_square_scan_points` | Resolution used when locating square mirror eigenvalues |
| `modal_basis_convergence_levels` | Number of grid refinement levels used to test basis convergence |
| `modal_basis_eigenvalue_relative_tolerance` | Maximum accepted relative eigenvalue change across basis refinements |
| `modal_basis_eigenfunction_overlap_tolerance` | Minimum accepted eigenfunction overlap across refinements |

More retained modes and finer basis grids generally improve modal resolution at increased computational cost.

### Eq. 59 Rosenbluth velocity solve

| Option | Effect |
| --- | --- |
| `modal_rosenbluth_max_iterations` | Maximum number of nonlinear Rosenbluth iterations |
| `modal_rosenbluth_min_iterations` | Minimum number performed before convergence may be accepted |
| `modal_rosenbluth_relative_tolerance` | Relative convergence criterion for the Eq. 59 iteration |
| `modal_rosenbluth_absolute_tolerance` | Absolute convergence criterion. A value of zero leaves the relative criterion dominant |
| `modal_rosenbluth_relaxation` | Under relaxation applied between Eq. 59 iterations. Smaller values damp updates more strongly |
| `modal_rosenbluth_high_speed_tail_population_tolerance` | Maximum allowed unresolved population fraction in the high speed tail |
| `modal_rosenbluth_high_speed_tail_energy_tolerance` | Maximum allowed unresolved energy fraction in the high speed tail |
| `modal_rosenbluth_speed_convergence_relative_tolerance` | Required agreement under speed grid refinement |
| `modal_retained_mode_relative_tolerance` | Required agreement when assessing retained mode convergence |
| `modal_local_speed_refinement_relative_tolerance` | Required agreement under refinement of the local reconstructed speed grid |

Reducing a tolerance makes the associated qualification stricter and can require additional resolution or iteration.

### Local reconstruction checks

| Option | Effect |
| --- | --- |
| `negative_distribution_roundoff_tolerance` | Defines the magnitude of small negative distribution values that can be treated as numerical roundoff |
| `modal_reconstruction_negative_particle_fraction_tolerance` | Maximum accepted particle fraction associated with negative values during reconstruction |
| `modal_reconstruction_energy_correction_relative_tolerance` | Maximum accepted relative energy correction introduced by reconstruction |
| `modal_local_velocity_quadrature_order` | Quadrature order used for local velocity space integrations |

### Electrostatic potential and density closure

| Option | Effect |
| --- | --- |
| `modal_phi_z_iterations` | Maximum iterations used to converge the axial electrostatic potential profile |
| `modal_phi_z_root_scan_points` | Number of scan points used to locate local quasineutrality roots. The value must be odd and at least five |
| `modal_phi_z_relative_tolerance` | Relative convergence threshold for the axial potential solution |
| `modal_phi_z_relaxation` | Under relaxation applied to potential profile updates |
| `modal_phi_z_low_energy_weight_fraction_tolerance` | Maximum allowed distribution weight associated with the low energy region where the potential reconstruction requires special treatment |
| `modal_electrostatic_feedback_iterations` | Maximum outer iterations coupling the electrostatic solution back into the modal kinetic solution |
| `modal_electrostatic_feedback_relative_tolerance` | Convergence criterion for that outer electrostatic feedback |
| `modal_density_closure_iterations` | Maximum iterations used for the self consistent density closure |
| `modal_density_closure_relative_tolerance` | Overall density closure convergence threshold |
| `modal_density_iteration_relative_tolerance` | Additional numerical consistency criterion applied to density updates |
| `modal_density_closure_relaxation` | Relaxation applied when updating the self consistent density |

### Density weighted modal basis

| Option | Effect |
| --- | --- |
| `modal_eq42_density_basis_iterations` | Maximum iterations coupling the solved density profile back into the density weighted modal basis |
| `modal_eq42_density_shape_relative_tolerance` | Relative convergence criterion for the density shape used by the orbit average |
| `modal_eq42_density_shape_relaxation` | Relaxation applied to density shape updates |
| `modal_eq42_basis_eigenvalue_relative_tolerance` | Required eigenvalue agreement during the density weighted basis iteration |
| `modal_eq42_basis_eigenfunction_overlap_tolerance` | Required eigenfunction overlap during that iteration |
| `modal_eq42_density_symmetry_tolerance` | Maximum accepted left or right density asymmetry for the symmetric model |

### Loss and source consistency

| Option | Effect |
| --- | --- |
| `modal_loss_convention_relative_tolerance` | Tolerance used when checking consistency of the represented ion loss formulations |
| `modal_source_projection_relative_tolerance` | Maximum relative error accepted when comparing the physical source rate with its modal projection |
| `fast_fusion_burnup_relative_tolerance` | Maximum fast ion fusion burnup fraction that may be neglected as an explicit sink in the Eq. 59 fast ion solve |

### Current balance and expander solution

| Option | Effect |
| --- | --- |
| `total_current_balance_iterations` | Maximum iterations used to match the ion and electron terminal currents |
| `total_current_balance_relative_tolerance` | Relative ion or electron current mismatch allowed at convergence |
| `total_current_balance_relaxation` | Relaxation applied to current balance updates |
| `total_current_balance_end_asymmetry_relative_tolerance` | Maximum accepted difference between the left and right end solutions for the symmetric end model |
| `expander_potential_iterations` | Maximum iterations used to solve the expander electrostatic potential |
| `expander_potential_relative_tolerance` | Relative convergence criterion for the expander potential |
| `expander_fast_throat_pitch_quantile_nodes` | Number of pitch quantile nodes used to map the lost fast ion distribution through the expander |
| `expander_fast_nonterminal_fraction_tolerance` | Maximum fast ion fraction allowed to remain without a terminal classification |
| `expander_unclosed_population_fraction_tolerance` | Maximum total expander population fraction allowed to remain outside the represented closure |

## Beam density fixed point

The `beam_density_coupling` block controls the iteration between beam attenuation and deposition and the plasma density supported by that deposited beam source.

| Option | Effect |
| --- | --- |
| `max_iterations` | Maximum number of beam density fixed point iterations |
| `electron_profile_relative_tolerance` | Pointwise relative convergence threshold for the electron density profile |
| `electron_profile_volume_L2_relative_tolerance` | Volume weighted L2 convergence threshold for the electron density profile |
| `electron_profile_absolute_reference_tolerance` | Additional reference normalized profile convergence criterion |
| `total_birth_rate_relative_tolerance` | Required relative agreement in total fast ion birth rate |
| `deposited_power_relative_tolerance` | Required relative agreement in total deposited beam power |
| `axial_birth_profile_relative_tolerance` | Required relative agreement in the axial fast ion birth profile |
| `energy_component_deposition_relative_tolerance` | Required relative agreement in power deposition by beam energy component |
| `target_density_relaxation` | Relaxation applied to density updates between fixed point iterations |
| `density_floor_absolute_m3` | Absolute density floor used to regularize relative convergence metrics near zero density |
| `density_floor_reference_fraction` | Additional density floor expressed as a fraction of the reference density |
| `birth_profile_floor_absolute_m3_s` | Absolute floor used when comparing low valued axial birth rate density cells |
| `birth_profile_floor_reference_fraction` | Reference scaled floor for the axial birth rate profile comparison |
| `component_power_floor_absolute_W` | Absolute power floor used when comparing very small beam components |
| `component_power_floor_reference_fraction` | Reference scaled floor for beam component power comparisons |

## Electron temperature power balance

| Option | Effect |
| --- | --- |
| `electron_temperature_mode` | Selects how electron temperature is obtained. `self_consistent_electron_energy` solves it from the represented electron energy balance. `fixed_closure` uses the initial electron temperature guess as a prescribed value and skips the self-consistent solve. This option will reduce runtime significantly but is not as physically accurate and may cause run failures.|
| `external_electron_heating_W` | Adds prescribed external electron heating to the represented energy balance, can be used to simulate electron cyclotron resonance heating or other electron heating systems. |
| `bracket_min_keV` | Lower bound of the electron temperature root search |
| `bracket_max_keV` | Upper bound of the root search |
| `bracket_scan_points` | Number of temperatures initially evaluated while locating a valid root bracket |
| `max_iterations` | Maximum iterations of the temperature root solution |
| `residual_tolerance_W` | Absolute accepted electron energy residual |
| `relative_residual_tolerance` | Relative accepted energy balance residual |
| `root_temperature_relative_tolerance` | Relative temperature tolerance used by the root solution |

The self consistent solution must remain inside the configured temperature bracket. The configured `plasma_closure.electrons.temperature_initial_guess_keV` is an initialization value rather than the final imposed temperature.

## Fusion calculation

| Option | Effect |
| --- | --- |
| `fusion.include_fast_fast` | Includes fusion reactions between represented fast ion populations |
| `fusion.num_gyroangle_points` | Sets the relative gyroangle quadrature resolution used by the fusion pair kernel |
| `fusion.max_pair_states_per_batch` | Limits the number of reactant pair states processed at once. This mainly controls memory use and batching rather than the physical model |

## Neutron calculation

| Option | Effect |
| --- | --- |
| `neutrons.model` | Selects the neutron spectrum model. Production uses `distribution_kinematics_spectrum` |
| `energy_min_MeV` | Lower edge of the saved deterministic neutron energy histogram |
| `energy_max_MeV` | Upper edge of the neutron energy histogram |
| `energy_bins` | Number of bins between the configured energy limits |
| `angular_model` | Selects the CM angular emission model. Production uses evaluated ENDF/B VIII.1 CM angular data for correlated events |
| `allow_out_of_range` | Controls whether weighted neutron events outside the configured energy histogram are permitted. When `false`, significant out of range rate causes the calculation to fail |
| `max_pair_states_per_batch` | Limits the reactant pair states processed simultaneously during neutron spectrum construction |
| `correlated_event_count` | Number of correlated neutron events generated for the event bank and OpenMC FileSource. Increasing it reduces event sampling noise and increases file size and generation time |
| `correlated_event_seed` | Random number seed used for reproducible correlated event generation |
| `radial_source_profile_model` | Selects the prescribed normalized radial source profile used to place neutron births across the local flux tube cross section |

## OpenMC source export

| Option | Effect |
| --- | --- |
| `openmc_export.source_output_directory` | Sets the directory, relative to the selected run directory unless an absolute path is supplied, where the validated OpenMC source bundle is written |
| `source_filename` | Sets the HDF5 OpenMC FileSource filename |
| `validity_filename` | Sets the filename containing source validity and qualification information |
| `sampling_metadata_filename` | Sets the filename containing source sampling metadata |
| `source_time_s` | Sets the time value written into exported OpenMC source particles |
| `overwrite_existing` | Controls whether an existing source bundle may be replaced |

The three source bundle filenames must be distinct. The HDF5 filename must end in `.h5`, and the two metadata filenames must end in `.json`.

## OpenMC run settings

`input_files/source_information.yaml` controls the OpenMC settings written into the exported model. These settings do not change the magnetic mirror source model physics.

| Option | Effect |
| --- | --- |
| `source.type` | Identifies the source integration as `source_model_revamp`. This should remain unchanged for the production workflow |
| `settings.particles_per_batch` | Sets `settings.particles` in the exported OpenMC model |
| `settings.batches` | Sets the number of OpenMC batches |
| `settings.statepoint_frequency` | Controls the batch interval at which statepoint files are requested |
| `settings.weight_windows` | Enables or disables the configured OpenMC weight window generation path |
| `settings.photon_transport` | Enables or disables coupled photon transport |
| `settings.tallies` | Controls whether the ParaTAN tally definitions are included in the exported OpenMC model |

These settings are written during model construction, but the NEWTON build command does not itself run OpenMC transport.