# Electron temperature coupling

The `coupling/` folder assembles the represented electron energy balance and solves the self consistent electron temperature used by the source model. It connects Eq 59 fast ion energy transfer, electron and ion end loss powers, current balance, expander terminal quantities, and repeated operating point solves into one scalar temperature closure.

## Files

`electron_energy_balance.py`: Collects fast ion to electron energy transfer, prescribed external electron heating, electron wall kinetic loss, and ambipolar ion acceleration power. It validates the independent current and terminal power identities and reports whether the represented terms are available, numerically evaluable, and qualified.

`electron_temperature_power_balance.py`: Solves the scalar electron temperature closure. Each admitted `Te` trial rebuilds an operating point, the bracket search uses geometrically spaced temperatures, rejected trials do not enter a sign bracket, and a valid sign change is refined with SciPy Brent root finding.

`electron_temperature_preflight.py`: Checks configuration level requirements before the temperature search and separates checks that can be decided immediately from checks that require one solved operating point.

## Important conventions

The temperature root coordinate is in keV, while `ElectronEnergyBalanceState.electron_temperature_J` and barrier energies are in J. Powers are in W and currents are in A.

The represented residual is `R = P_heat − P_loss`. Heating is the sum of direct Eq 59 fast ion to electron transfer and `external_electron_heating_W`. Loss is the electron wall kinetic power plus ambipolar ion acceleration power. The latter is evaluated as `P_amb = (I_i / e) E_barrier`, using the ion current matched by the shared current closure.

When an expander terminal state is available, terminal ion current and terminal electrostatic power replace the central cell estimates. Prompt losses that lack a resolved expander branch use the shared barrier as a fallback contribution and remain explicitly marked in integration metadata.

A temperature residual is admitted only when the operating point and total current fixed point have converged and the represented energy terms are evaluable. Static configuration failures terminate the search. Trial specific physical failures are retained as rejected temperatures and are never used to infer a root bracket.

The scan is geometric in `Te`. Brent refinement starts only between adjacent admitted trials with opposite residual signs. If the first admitted negative residual follows an invalid lower scan point, the solver geometrically refines that gap until it finds a converged trial, a valid sign bracket, or reaches the configured temperature tolerance.

## Physics basis

The temperature closure follows the energy balance logic of Egedal et al. 2022 section 3.3, where electron heating and mirror end losses determine `Te`. The implementation evaluates the represented powers directly rather than using the reduced normalized Eq 24 expression. Fast ion to electron transfer comes from the discrete Eq 59 collision energy moment. Electron wall loss and the matched electron ion current relation follow the Eq 68 and Eq 69 end loss treatment, while the ambipolar ion acceleration term uses the same wall barrier energy carried by the shared current closure.
