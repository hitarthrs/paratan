# Fusion rates and source profiles

The `fusion/` folder evaluates D D and D T reaction rates from the represented fast ion populations and converts those rates into axial neutron and power source profiles. It uses Bosch and Hale total cross sections for arbitrary local gyrotropic distributions and keeps the independent Maxwellian thermal reactivity fit available for validation.

## Files

`bosch_hale.py`: Stores the Bosch and Hale Table IV S factor coefficients and evaluates total fusion cross sections as functions of center of mass energy.

`cross_sections.py`: Provides reduced mass and center of mass energy helpers, the Table VII Maxwellian thermal reactivity fit, and the `σ(E_cm) g` pair kernel wrapper.

`pair_kernel.py`: Applies the Table IV total cross section to local gyrotropic reactant distributions. It returns the axial reaction rate, reaction weighted kinetic energy removed from each reactant population, and the contribution from pair energies inside and outside the stated Table IV fit interval. It also compares Maxwellian Table IV quadrature with the independent Table VII fit.

`populations.py`: Defines the represented fast deuterium and fast tritium populations and the allowed D D and D T population pair definitions.

`reactions.py`: Defines reaction branch energies, neutron yields, compact reaction keys, and energy, cross section, and reactivity unit conversions.

`reactivity.py`: Evaluates general gyrotropic pair integrals over speed, pitch, and relative gyrophase. It supports shared kernels, reaction weighted state moments, axial stacks, identical population counting, and density normalized `<σg>` values.

`source_profiles.py`: Converts reaction rate density profiles into neutron source and power profiles and integrates them over supplied physical cell volumes.

## Important conventions

Local reactant distributions use `f(z, v, ξ)` with shape `(n_z, n_speed, n_pitch)`, where `ξ = v_parallel / v`. The FBIS velocity cell measure represents `d^3v = 2π v^2 dv dξ`. The remaining relative azimuth is integrated by midpoint sampling of the relative gyrophase.

For a reactant pair, `g = |v_a − v_b|` and `E_cm = 0.5 * μ * g^2`. The arbitrary distribution rate kernel is `σ(E_cm) g`. A factor of `1/2` is applied only when both reactants are drawn from the same represented population so identical pairs are not counted twice.

The Table IV cross section interface uses center of mass energy and returns `σ` in `m^2`. The pair kernel evaluates the represented pair distribution without clipping it to the Table IV fit interval and separately records the reaction rate inside and outside that interval. No rate rescaling is applied by this module.

The Table VII interface accepts ion temperature energy in J and returns Maxwellian reactivity in `m^3/s`. Its D T and D D validation interval is `0.2 <= T_i <= 100 keV`.

`FusionPairKernelResult.reactant_a_kinetic_energy_removal_density_W_m3` and the corresponding reactant B field are reaction weighted kinetic energy consumption rates. They use the kinetic energy of the reacting state before the reaction and therefore have units `W/m^3`.

Reaction rate density and neutron source density use `m^−3 s^−1`. Power density uses `W/m^3`. When physical cell volumes are supplied, the source profile helpers also report integrated rates in `s^−1` and powers in W.

## Physics basis

The cross section path follows Bosch and Hale 1992 Table IV and Eqs 8 and 9, where the fitted S factor is combined with the Coulomb penetrability to obtain the total fusion cross section as a function of center of mass energy. The pair integration follows the distribution averaged reactivity in their Eq 11 and evaluates it directly for the represented non Maxwellian gyrotropic ion distributions.

The independent Maxwellian check follows Bosch and Hale Table VII and Eq 14. The Table VII fit is used only as a thermal comparison to the Table IV velocity space quadrature. Neutron angular emission and correlated two body kinematics are handled later by the `neutrons/` package.
