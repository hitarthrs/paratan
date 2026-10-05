# Modal electrostatic reconstruction

The electrostatic/ folder solves the axial Egedal Eq 70 to Eq 72 closure used by the modal FBIS backend. It reconstructs confined fast ion charge from the invariant distributions, includes fixed magnetic boundary Eq 63 lost ion charge when supplied, solves the shared electron and ion quasineutrality profile, and records exact midplane and throat diagnostics.

## Files

`accessibility.py`: Checks whether the Eq 71 trapped passing boundary is actually set by the mirror throats and measures the modal inventory in the low energy interval where the Eq 71 confined domain closes.

`eq70.py`: Evaluates the Eq 70 confined Maxwellian electron density fraction, limits potential clipping to floating point roundoff, and calculates local and volume integrated quasineutrality residual metrics.

`iteration.py`: Implements the single species Eq 70 to Eq 72 profile iteration with an exact zero potential midplane reference, direct throat roots, local invariant reconstruction, and convergence diagnostics.

`nodes.py`: Inserts exact left throat, midplane, and right throat nodes into cell centered profiles and recomputes the anchor electron densities and residuals directly.

`system_iteration.py`: Solves the shared multi species electrostatic state. It couples the Eq 70 parent Maxwellian amplitude to left and right throat barriers, includes confined and Eq 63 lost ion charge, follows bounded axial roots, solves the exact midplane separately, and returns species resolved local reconstructions.

## Important conventions

The electron closure uses a nonnegative potential energy coordinate bounded by the wall barrier. The Eq 70 helper follows the sign convention `ΔΦ = −e ϕ`, while the reported electrical potential relative to the exact midplane is formed from the difference in potential energy coordinates divided by `−e`.

Eq 70 uses a parent Maxwellian density `n0`. In the shared solve, when an electron midplane density is not prescribed, `n0` floats with the maximum positive charge state available at zero potential drop. The relaxed scalar iteration is followed by an exact reference reconciliation before the axial roots are finalized.

Confined ions preserve total energy `U` and magnetic moment through the local mapping. Eq 71 sets the accessible trapped passing interval and Eq 72 maps that interval back to the magnetic reference coordinate. Fixed magnetic boundary Eq 63 lost populations enter the ion charge closure separately from the confined population.

Fast ion local reconstructions use species specific speed grids. Local `Λ` arrays have shape `(n_z, n_local_speed, n_lambda)`, local pitch arrays have shape `(n_z, n_local_speed, n_pitch)`, and density profiles have shape `(n_z,)`. Charge closure uses charge number weighted number density in m⁻³, not Coulombs per cubic metre.

Exact `ζ = −1`, `ζ = 0`, and `ζ = 1` nodes are retained for throat and midplane checks even when the axial solve is cell centered. Symmetric inputs may be reflected and averaged, but the exact anchor states are still evaluated directly.

## Physics basis

The electron density closure follows Egedal Eq 70. Confined ion accessibility uses the Eq 71 trapped passing boundary and Eq 72 invariant compression. The fixed magnetic boundary lost ion contribution follows the Eq 63 branch construction supplied by the modal loss closure. The solver combines these densities to enforce axial quasineutrality and retains the Eq 71 effective potential check to verify that the fitted mirror throats set the passing boundary.
