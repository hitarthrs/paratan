# Neutron spectra

The neutrons folder converts represented D D and D T reactant distributions into neutron energy spectra and relativistic lab kinematics. The direct files provide the deterministic spectrum path and shared two body kinematics, while the events subpackage builds correlated neutron event banks with evaluated angular data.

## Files

`kinematics.py`: Defines D T and D D neutron branch masses and evaluates relativistic two body neutron kinematics. It provides scalar and batched Lorentz transformations, exact isotropic center of momentum lab energy endpoints, and deterministic isotropic direction sets.

`lambda_spectrum.py`: Connects invariant or local FBIS reactant distributions to axial neutron energy matrices. It can conservatively map `F(v, Λ)` to local `f(z, v, ξ)`, attach cell volumes, retain raw and out of range rate diagnostics, and form globally integrated spectra.

`spectrum.py`: Evaluates deterministic neutron spectra from local gyrotropic reactant distributions. The pair kernel uses `f_a d³v_a f_b d³v_b σ(E_cm) |v_a − v_b|`, includes the identical population `1 / 2` factor when required, and supports exact isotropic center of momentum energy marginal integration.

## Subpackages

`events/`: Builds correlated D D and D T neutron events from sampled reactant pairs, evaluated center of momentum angular laws, relativistic two body kinematics, and sampled source positions. It is annotated separately from the direct neutron files.

## Important conventions

Reactant velocity distributions are gyrotropic about the local magnetic field and use pitch cosine `ξ = v_parallel / v`. Local distribution arrays have shape `(n_z, n_speed, n_pitch)`, while invariant distributions use speed and `Λ = μ B0 / E`.

Total reaction weighting uses the same total fusion cross section path as the fusion rate calculation. For identical reactant populations the pair contribution is multiplied by `1 / 2` to avoid double counting.

The deterministic energy matrix stores neutron rate density integrated over each energy bin, not `dR / dE`. Dividing by the energy bin width gives the spectral density.

All neutron kinetic energies and energy edges in these direct helpers are in J. Velocities are in m s⁻¹ and cell volumes are in m³.

## Physics basis

The relativistic two body kinematics and arbitrary reactant spectrum structure follow the supplied Eriksson et al. arbitrary reactant methodology. Reactant velocities determine the incoming four momentum and invariant `s`, the neutron and residual are solved as a two body final state in the center of momentum frame, and the neutron four momentum is transformed to the lab frame.

Reaction weighting follows the supplied Bosch and Hale total cross section representation through the fusion module. The evaluated angular dependence used by correlated events is handled in `nuclear_data/endf_b_viii1/` and `neutrons/events/`, not in the deterministic direct spectrum files.
