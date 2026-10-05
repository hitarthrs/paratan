# ENDF B VIII.1 fusion angular data

The `endf_b_viii1/` folder provides the evaluated center of mass neutron angular laws used for correlated D D and D T neutron events. It loads compact runtime arrays, validates the angular probability laws, interpolates them in incident energy, and converts arbitrary reactant pair invariant mass to the equivalent deuteron laboratory energy used by the evaluations.

## Files

`angular_distribution.py`: Loads the D D and D T runtime arrays, verifies reaction identity and angular representation, evaluates `p(mu | E)`, converts it to probability per solid angle, and measures angular normalization and nonnegativity.

`incident_energy.py`: Converts relativistic pair invariant `s` to the equivalent projectile laboratory kinetic energy for a stationary target in J or eV.

`interpolation.py`: Validates the supported ENDF incident energy interpolation metadata and provides linear interpolation in incident energy and cosine.

`angular_data_provenance.json`: Records the source evaluation identity, reaction metadata, hashes, energy coverage, angular representation, and preprocessing audit for the packaged runtime arrays.

`ddn_angular_data.npz`: Stores the `D(d,n)3He` MAT 128, MF 6, MT 50 angular data as tabulated probability density in `mu`.

`dtn_angular_data.npz`: Stores the `T(d,n)4He` MAT 131, MF 6, MT 50 angular data as Legendre coefficients.

## Important conventions

The angular cosine is `mu = cos(theta_cm)`, where `theta_cm` is measured from the incident deuteron direction to the emitted neutron direction in the center of mass frame. The stored conditional density satisfies `integral from −1 to 1 of p(mu | E) dmu = 1`. Uniform azimuth about the incident deuteron axis gives `p_Omega = p_mu / (2 pi)` in `sr^−1`.

The evaluated incident energy coordinate is deuteron laboratory kinetic energy on a stationary target. For arbitrary reactant velocities, the event model preserves the pair invariant and uses

`K_lab = [s − (m_p c^2 + m_t c^2)^2] / (2 m_t c^2)`

to obtain the equivalent incident energy before angular sampling.

The D D dataset uses `LCT = 2`, `LAW = 2`, and `LANG = 12` with a piecewise linear tabulated `p(mu)`. Its packaged angular grid covers 20 eV through 20 MeV. The D T dataset uses `LCT = 2`, `LAW = 2`, and `LANG = 0` with Legendre coefficients. Its packaged angular grid covers 100 eV through 40 MeV.

Incident energy interpolation is limited to ENDF `INT = 2`, so angular probabilities or Legendre coefficients are blended linearly between adjacent energy knots. Queries outside the packaged energy interval are rejected.

The angular distributions are normalized conditional probabilities only. The reaction rate weight comes from the Bosch and Hale Table IV total cross section in `fusion/`.

## Physics basis

The runtime arrays are derived from the supplied ENDF B VIII.1 incident deuteron evaluations for `D(d,n)3He` and `T(d,n)4He`, using MF 6, MT 50 two body center of mass angular data. The D D evaluation is represented by tabulated angular probabilities and the D T evaluation by Legendre coefficients.

The correlated neutron event method follows the supplied Eriksson et al. arbitrary reactant treatment in which reactant four vectors determine the center of mass state and evaluated angular information supplies the directional dependence. The total fusion cross section remains the Bosch and Hale Table IV model handled by `fusion/`.
