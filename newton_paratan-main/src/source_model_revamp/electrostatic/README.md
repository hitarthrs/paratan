# Electrostatic current balance

The electrostatic/ folder provides the ambipolar current balance closure used by the source model. It converts ion particle losses into a positive ion loss current and solves for the electron confining barrier that gives the matching electron loss current.

## Files

`current_balance.py`: Builds species resolved ion loss currents and inverts the Egedal electron tail refilling relation for the barrier energy, normalized barrier, wall potential, electron particle loss rate, and wall power loss. It supports both the actual midplane electron density normalization and a fixed Eq 70 parent Maxwellian density.

## Important conventions

Ion and electron currents are stored as positive loss magnitudes. For ion species `s`, `I_s = e Z_s Γ_s`, and species values are summed over axis zero.

The electron barrier is represented by `y = E_barrier / T_e`. `E_barrier` and `T_e` are energies in J. Before the full axial profile is reconciled, `Φ_wall − Φ_midplane = −E_barrier / e`. With a floating Eq 70 reference, the stored wall potential uses `Φ_wall − Φ_midplane = −(E_barrier − ΔE_midplane) / e`.

The actual midplane density path uses `n_e(0) = n0 P(3/2, y)` and therefore updates the parent Maxwellian density `n0` for every trial barrier. The fixed parent density path holds `n0` constant and is used when the coupled electrostatic profile has already supplied the Eq 70 normalization.

`number_of_ends` explicitly scales the summed electron loss rate. The modal solver uses two symmetric ends when matching the electron current to the total ion loss current.

The electron collision density is retained with the balance state, while the electron loss prefactor uses the supplied electron collision frequency directly.

## Physics basis

The current balance follows the supplied Egedal et al. 2022 treatment of electron end loss. The electron loss rate uses the exact tail integral that precedes the high barrier Eq 68 approximation, and the wall power uses the exact escaping energy integral whose high barrier limit gives Eq 69. The Eq 70 electron density relation supplies the connection between the confined midplane density and the parent Maxwellian density `n0`.

The full axial quasineutral potential profile is solved separately in `fbis/modal/electrostatic/`. This module supplies the scalar electron barrier and current balance state used by that coupled solve.