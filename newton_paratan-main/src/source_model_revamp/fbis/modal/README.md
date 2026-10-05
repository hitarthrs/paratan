# Modal FBIS backend

The modal/ folder contains the Egedal modal fast beam ion solver(FBIS) used by the source model. It connects the orbit averaged magnetic basis to attenuated beam births, the cold reference and hot eq 59 velocity solves, ion loss closures, the eq 70-72 electrostatic reconstruction and species level solver results.

## Files

`basis.py`: Builds the magnetic Λ to η mapping, square mirror basis, density weighted orbit averaged Lorentz eigenbasis, Eq 61 geometry data, and basis convergence diagnostics. It also evaluates convergence with retained physical mode count.

`density_weighting.py`: Builds the even axial density ratio used in Egedal Eq 42. Profiles are normalized by the physical flux tube cell volume average and can represent the uniform reference, a prescribed profile, or the self consistent quasineutral profile.

`electrostatic_feedback.py`: Implements the published Egedal low energy first eigenvalue substitution. It also constructs the same shaped reference magnetic profile at mirror ratio 1.5 used to obtain the substituted first eigenvalue.

`full_device_lost_reconstruction.py`: Reconstructs the production fixed magnetic boundary Eq 63 lost ion branches in local speed and pitch after the electrostatic potential has been solved. Central cell pitch integrals are evaluated analytically and expander reconstruction is deferred to the expander model.

`loss_geometry.py`: Evaluates the dimensionless Eq 62 geometry factor and the magnetic flux tube volume integral from midpoint field data. The loss geometry can also include the active Eq 42 density ratio.

`lost_ion_distribution.py`: Builds Eq 61 through Eq 64 loss quantities. It calculates the production throat parallel temperature state, matches left and right Eq 63 branches independently to their one end Eq 61 speed spectra, and retains magnetic reference loss reconstructions for diagnostics.

`models.py`: Defines and validates canonical selectable velocity solution and electrostatic feedback model names.

`source_projection.py`: Partitions attenuated beam births into prompt magnetic losses and confined births, projects the confined source onto the physical eigenfunctions, and audits signed modal truncation without forcing a rate normalization.

`species_state.py`: Defines one species solver requests, solved or inactive species states, and the shared multi species fast ion system state.

`types.py`: Defines the numerical controls, modal basis state, Rosenbluth coefficient state, electrostatic profile, local reconstruction, beam component result, and completed modal FBIS result.

`utils.py`: Provides small integer validation and row interpolation helpers shared by modal modules.

`velocity_solve_cold.py`: Implements the Egedal Eq 14 and Eq 15 cold ion reference solution and the discrete monoenergetic source cell representation used by velocity space operators.

## Subpackages

`eq59/`: Hot ion Eq 59 velocity space collision operator, nonlinear iteration, Rosenbluth state, energy moments, and speed convergence.

`local/`: Mapping from the invariant modal solution to local speed, Λ, and pitch representations after an electrostatic potential is applied.

`electrostatic/`: Eq 70 electron closure, Eq 71 and Eq 72 ion accessibility mapping, exact axial node handling, and coupled quasineutral potential iteration.

`solver/`: High level one species and coupled species orchestration that assembles sources, collisions, velocity solves, loss closures, electrostatics, reconstruction, convergence, and metadata.

## Coordinate and array conventions

The invariant pitch coordinate is `Λ = μ B0 / E`. For the magnetic reference problem it spans zero through one and the magnetic trapped passing boundary is `Λ = 1 / RM`. The Eq 47 coordinate η is built from the normalized bounce time. It runs in the opposite direction to Λ, and the confined physical eigenbasis is evaluated from `η = 0` through the trapped passing value `ηTP`.

Modal velocity distributions use shape `(n_mode, n_speed)`. Physical invariant distributions generally use `(n_speed, n_eta)` or `(n_speed, n_lambda)`. Local reconstructions add axial position first, and pitch reconstructions therefore use `(n_z, n_local_speed, n_pitch)`.

The source projection uses the folded even parallel sign measure represented by η. Physical beam birth rates are first calculated from the axial birth rate density and cell volumes. Only births mapped inside the confined η interval enter the modal right hand side. Prompt magnetic loss births remain outside the confined solve.

## Eq 42 density weighting

Egedal Eq 42 weights the orbit average by `n(z) / <n>` because the scattering frequency varies with plasma density. The implementation normalizes this ratio with physical flux tube cell volumes and symmetrizes it about the midplane. Density weighting changes the orbit averaged Lorentz coefficients. It does not change the magnetic Eq 47 η mapping.

## Loss closure

The production ion loss path uses the fixed magnetic boundary `ΛM = 1 / RM`. The Eq 61 one end loss spectrum provides the target flux for each directed throat branch. The Eq 63 branch construction solves for `H(U)` independently in each invariant speed cell using the continuous loss cone flux integral, so no posterior rate rescaling is applied.

After the quasineutral potential is known, `full_device_lost_reconstruction.py` maps the invariant total energy to local kinetic coordinates and reconstructs left moving and right moving lost populations on the local speed and pitch grid. The Eq 71 and Eq 72 moving electrostatic boundary remains part of the confined population mapping and is not substituted for the production fixed magnetic loss boundary.

## Physics basis

The modal chain follows the supplied Egedal et al. 2022 formulation. The direct files use the square mirror modal expansion and cold reference solution from Eqs 9 through 15, the orbit averaged Lorentz operator and density weighting from Eqs 41 and 42, the η coordinate and physical eigenbasis from Eqs 47 and 50, the hot ion velocity equation from Eq 59, the ion loss relations from Eqs 61 through 64, and the electrostatic reconstruction from Eqs 70 through 72.
