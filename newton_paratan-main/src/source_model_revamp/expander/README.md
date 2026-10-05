# Expander mapping and electrostatic closure

The `expander/` folder maps directed fast ion loss populations from each mirror throat through the physical expander geometry. It solves the side resolved quasineutral potential, reconstructs finite transit ion distributions on the full device grid, routes source connected ions to material surfaces, and reports terminal current and power quantities used by the integration stage.

## Files

`mapping.py`: Propagates directed throat samples using conserved total energy and magnetic moment, classifies terminal and reflected trajectories, assigns finite transit inventory to local speed and pitch cells, applies radial material accessibility, and audits particle, power, and inventory conservation.

`solver.py`: Builds left and right throat source samples from the Eq 63 loss state, constructs the physical expander paths, solves the cell center Eq 70 quasineutral potential on each side, maps the resulting branches, and assembles terminal current, power, material loss, and closure diagnostics.

`types.py`: Defines and validates branch, potential, and coupled expander states together with their geometry, distribution, conservation, and status identities.

## Important conventions

Each path is ordered outward from a mirror throat to the selected terminal surface. Edge arrays therefore have length `n_cell + 1`, while cell centered arrays have length `n_cell`. The mapped velocity distribution has shape `(n_full_z, n_speed, n_pitch)` and integrates with the gyrotropic velocity cell measure to number density in m⁻³.

Fast ion inputs are built from the left and right Eq 63 one end speed spectra. Within each speed bin the pitch variable is `x = 1 − R_M Λ`, with `Λ = (1 − x) / R_M`, `U = 0.5 m v²`, and `μ B0 = U Λ`. Midpoint quantiles of the truncated Eq 63 exponential flux preserve the supplied one end speed spectrum without an additional rate rescaling.

Source connected motion follows the implemented parallel energy relation `K_parallel = U + Z ΔΦ − (μ B0) B_tilde`. A state is terminal only when it remains accessible from the throat through the full path. A first blocked segment produces a reflected branch, while positive parallel energy beyond that block marks local well topology without restoring source connectivity.

The side potential is solved at physical expander cell centers. Throat and wall potential energies are fixed, intermediate edge values are linearly reconstructed, and each physical cell is split into two synthetic mapping segments during the Eq 70 residual evaluation. Ion charge density is the charge number weighted sum of the source connected fast populations. Electron density uses the parent Maxwellian normalization from the shared central cell electrostatic state.

Radial accessibility is represented by `rho = r / a_nominal(z)` and a survival probability from the configured radial profile. Loss of radial support is assigned to the first material surface encountered along the path. The selected axial terminal surface receives the remaining source connected population that reaches the end of the path.

Prompt magnetic losses are retained as unresolved when only a total prompt rate is available and no side resolved throat distribution exists. Such a rate contributes to the throat total but cannot be mapped into an expander distribution or terminal current branch.

## Physics basis

The fast ion throat source follows the supplied Egedal et al. 2022 Eq 63 lost ion form and the Eq 61 one end loss spectrum. The expander mapping preserves the total energy `U` and magnetic moment used by that formulation. The shared wall barrier comes from the Eq 68 current balance, and the local electron density closure follows Eq 70. The additional path discretization, radial material accessibility, finite transit inventory mapping, and conservation audits are numerical and geometric structures of this implementation.
