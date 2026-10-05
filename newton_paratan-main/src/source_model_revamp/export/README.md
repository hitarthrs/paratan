# Source export

The export folder serializes correlated neutron event banks and converts them into the OpenMC HDF5 FileSource representation used by the transport model. It preserves correlated position direction and energy information while keeping normalized sampling probability separate from the absolute physical neutron rate.

## Files

`event_bank.py`: Writes and reads the versioned compressed NPZ sidecar for correlated neutron events. The sidecar stores position direction energy normalized probability reaction labels axial cell indices center of mass audit quantities reactant velocities and the separate physical neutron rate.

`openmc_file_source.py`: Converts positive probability correlated events into an OpenMC HDF5 source bank, performs measured OpenMC read back checks, writes validity and sampling metadata JSON sidecars, and returns the corresponding `openmc.FileSource` object.

## Important conventions

The correlated event bank stores positions in m and neutron kinetic energies in J. OpenMC source particles are written with positions in cm and energies in eV. Lab directions are unit vectors in both representations.

`normalized_weights` represents event sampling probability and sums to one. Zero probability audit events remain in the NPZ event bank but are excluded from the OpenMC HDF5 source bank.

For `N` exported positive probability events with probabilities `p_i`, the OpenMC source particle weight is `w_i = N p_i`. The resulting source particle weights have mean one under uniform row selection and preserve the normalized event probability distribution.

The physical neutron rate is stored separately as `physical_total_rate_s` in units of s⁻¹. The OpenMC `FileSource` strength is one, so transport results use the separate physical rate for absolute normalization rather than interpreting source particle weights as neutrons per second.

The OpenMC bundle is staged and validated before the destination directory is replaced. Read back checks cover HDF5 structure, position, direction, energy, time, particle identity, probability weights, weighted energy direction covariance, and rowwise lab energy direction correlation.
