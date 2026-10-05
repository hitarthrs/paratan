# Orbits

This module supplies the magnetic single particle coordinates and orbit measures used by the modal FBIS calculation.
Converts local velocity geometry into the conserved pitch invariant Lambda, evaluates magnetic bounce motion and Egedal Eq 42 orbit averages, and constructs the eta coordinate used by the physical scattering basis.
Electrostatic potential effects are handled later in the source model and are not included in these magnetic orbit helpers.

## Files

invariants.py:        Converts between energy, speed, magnetic moment, pitch angle, xi, and Lambda and evaluates magnetic accessibility and trapped passing boundaries.
bounce_averages.py:   Finds magnetic orbit limits, evaluates normalized and full bounce times, handles the magnetic separatrix, and evaluates Egedal Eq 42 spatial averages.
eta_mapping.py:       Builds the Egedal eta coordinate from normalized bounce time and provides Lambda, eta, and derivative interpolators.

## How it fits into the model

fbis/beam_source_definition.py: converts the beam birth pitch angle and local normalized field into the source value of Lambda.

fbis/modal/basis.py: builds the full Lambda grid, evaluates the bounce time measure, maps Lambda to eta, locates the trapped passing boundary in eta, and supplies the mapping to the physical modal basis.

scattering/lorentz_operator.py: uses the Eq 42 orbit averages to form the density weighted coefficients of the orbit averaged Lorentz operator.
scattering/square_mirror_basis.py: uses the square mirror xi trapped passing boundary when constructing the reference basis.

## Coordinates and units

| Quantity | Convention |
| ---  | --- |
| `energy_J`               | Kinetic energy in joules |
| `mass_kg`                | Particle mass in kilograms |
| `speed_m_s`              | Particle speed in meters per second |
| `local_B_T` and `B0_T`   | Magnetic field in tesla |
| `B_tilde`                | Local field divided by the reference field B0 |
| `zeta`                   | Positive axial position normalized by the mirror half length |
| `mirror_ratio`           | Magnetic throat field divided by the reference field |
| `Lambda`                 | Magnetic pitch invariant `mu * B0 / E`, from zero through one in the central mirror convention |
| `xi`                     | Signed local parallel speed fraction `v_parallel / v` |
| `tau_tilde_b`            | Dimensionless normalized magnetic bounce time |
| `eta`                    | Egedal orbit pitch coordinate, decreasing from one at Lambda zero to zero at Lambda one |
| `density_ratio_function` | Eq 42 density ratio `n(zeta) / <n>` normalized by physical flux tube volume |

The magnetic trapped passing boundary is `Lambda = 1 / R_M`.
Trapped particles terminate at a magnetic turning point where `B_tilde = 1 / Lambda`, while passing particles extend to the throat.

## References

The orbit invariant, normalized bounce time, orbit average, phase space normalization, and eta coordinate follow J Egedal et al, *Nuclear Fusion* 62 (2022) 126053, Eqs 35 through 47.
