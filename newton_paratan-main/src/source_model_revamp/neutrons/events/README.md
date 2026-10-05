# Correlated neutron events

The events folder builds correlated D D and D T neutron source particles from the represented fusion reactant populations. It samples reacting ion pairs inside fixed component and axial rate strata, applies evaluated center of mass angular laws, transforms the two body final state to the lab frame, samples source positions, and retains the event fields needed for transport export and validation.

## Files

`angular_sampling.py`: Samples normalized evaluated center of mass angular distributions. Tabulated `p(mu)` data use piecewise linear inverse CDF sampling, while Legendre representations use coefficient interpolation and rejection sampling.

`event_kinematics.py`: Connects sampled reactant pairs to invariant `s`, the equivalent deuteron incident energy used by the evaluated angular law, center of mass emission directions, relativistic lab neutron states, and residual mass shell diagnostics.

`generator.py`: Allocates events across component and axial strata, samples pair and angular state, applies conditional `sigma * g` importance weights, samples physical positions, assembles the correlated event bank, and records sampling diagnostics.

`kinematics.py`: Computes the reactant pair invariant and the deuteron projectile direction in the pair center of mass frame. It also constructs emission directions around arbitrary projectile axes.

`pair_sampling.py`: Samples represented gyrotropic reactant velocity cells with the `f d^3v` measure, assigns gyro angles, evaluates `E_cm` and `sigma(E_cm) * g`, and handles D D projectile exchange symmetry when requested.

`spatial_sampling.py`: Samples axisymmetric positions inside selected axial cells using the configured radial probability law and source connected radial support.

`types.py`: Defines validated component specifications and correlated event bank records, including physical rate, normalized probability weights, reaction labels, sampled center of mass variables, and both reactant velocities.

## Important conventions

Each fusion component and axial cell forms a stratum with physical rate `R_s`. Positive strata receive at least one event. Remaining samples are apportioned by `R_s`, but the physical rate itself is not estimated from the event count.

Reactant velocity cells are sampled with probability proportional to `f d^3v`. The remaining conditional reaction importance is `sigma(E_cm) * g`. Event probability inside each stratum is normalized from this importance, then multiplied by the stratum fraction of the total physical neutron rate.

The event bank therefore separates probability and normalization. `normalized_weights` sum to one and describe event sampling probability. `physical_total_rate_s` stores the absolute neutron rate. Their product gives the physical rate represented by each stored event.

D T components use the deuteron as reactant a for the evaluated projectile axis. Identical D D populations randomly exchange the two sampled deuterons before angular sampling so the projectile assignment is exchange symmetric.

The evaluated angular variable is `mu = cos(theta_cm)` relative to the deuteron projectile direction in the reactant pair center of mass frame. The equivalent incident energy is derived from invariant `s` for a deuteron projectile on a stationary target. The sampled neutron four momentum is then transformed to the lab frame.

Events with negligible `sigma * g` importance outside the evaluated angular energy range can be retained with zero probability weight for audit. A material excluded importance fraction is rejected. Downstream OpenMC file source creation uses only positive probability events.

Positions are sampled after the component and axial stratum is selected. Axial position is uniform within the selected cell, azimuth is uniform, and the radial law is conditioned on the available inner and source connected outer support.

Event positions are in m, neutron kinetic energies are in J, reactant velocities are in m s⁻¹, equivalent deuteron incident energies are in eV, and directions are unit vectors.

## Physics basis

The correlated event chain follows the supplied Eriksson et al. arbitrary reactant spectrum method. Reactant velocities define the incoming four momentum and invariant `s`, the differential reaction weighting contains relative speed and angular cross section information, and the two body final state is transformed between the center of mass and lab frames.

The total D D and D T reaction weighting uses the Bosch and Hale cross section path from the fusion module. The conditional center of mass angular law comes from the packaged ENDF B VIII.1 data in `nuclear_data/endf_b_viii1/`. The angular law changes event direction probability but does not replace the total fusion rate supplied by the fusion calculation.
