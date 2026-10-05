# Scattering module

`source_model_revamp.scattering` constructs the pitch angle scattering operators and eigenbases used by the modal FBIS calculation. 
It starts from the Lorentz operator, expresses that operator in the magnetic invariant coordinate `Lambda`, applies orbit averaging for the physical mirror field, and solves for the physical pitch eigenfunctions used to project and reconstruct the fast ion distribution.

This module is magnetic and invariant space infrastructure. It does not solve the FBIS velocity equation or the ambipolar potential. The self consistent density profile can enter through the Egedal Eq 42 scattering weight, while electrostatic corrections are handled later in the modal electrostatic and local reconstruction paths.

# Files

## `__init__.py`

Provides the package description without defining additional exports.

## `lorentz_operator.py`

Defines the `Lambda` space Lorentz coefficient container and helpers for local Eq 40 coefficients, orbit averaged Eq 41 coefficients, Eq 42 density weighted coefficients, direct construction from magnetic bounce averages, and finite difference application of a stored operator.

## `square_mirror_basis.py`

Builds the square mirror reference basis. It evaluates the even hypergeometric form of the Legendre solutions, locates the noninteger orders that satisfy the loss boundary condition, computes the corresponding eigenvalues and normalization integrals, and evaluates mode derivatives.

## `physical_eigenbasis.py`

Maps the square mirror basis to the physical `eta` interval, converts eta derivatives to `Lambda` derivatives, assembles the strong or conservative weak projected operator, solves the generalized eigenproblem, reconstructs the physical modes, and records residual, symmetry, conditioning, boundary, normalization, and orthogonality diagnostics.

# Physics coordinates

The module uses three related dimensionless pitch coordinates.

* `xi = v_parallel / v` is the square mirror pitch coordinate. Even symmetry allows the reference basis to use the positive trapped interval from `0` to `xi_TP`
* `Lambda = mu * B0 / E` is the magnetic invariant used for the nonuniform mirror. The trapped region extends from `1 / R_M` to `1`
* `eta` is the phase space coordinate constructed from the normalized bounce time in `orbits/eta_mapping.py`. On the trapped interval used here, `eta = 0` corresponds to `Lambda = 1` and `eta = eta_TP` corresponds to `Lambda = 1 / R_M`

All three coordinates, the Lorentz coefficients, and the scattering eigenvalues are dimensionless.

# Operator and basis chain

The square mirror reference problem follows the Lorentz operator

```
L = d/dxi [(1 - xi^2) d/dxi]
L M_j = -lambda_j M_j
```

with even midplane conditions and `M_j(xi_TP) = 0`. The boundary condition selects noninteger Legendre orders. These reference functions are mapped linearly from the square mirror `xi` interval onto the physical trapped `eta` interval.

For a nonuniform magnetic field, `lorentz_operator.py` implements the local `Lambda` coefficients from Egedal Eq 40 and the orbit averaged form from Eqs 41 and 42. With density weighting, the coefficients are

```
A = 4 <1 / B_tilde>_z,n - 6 Lambda <1>_z,n
D = 4 Lambda [<1 / B_tilde>_z,n - Lambda <1>_z,n]
```

where the averages are supplied by the orbit module. A uniform density profile reduces `<1>_z,n` to one.

`physical_eigenbasis.py` expands each physical mode in the mapped square mirror basis and solves the projected eigenproblem

```
K V = -lambda G V
```

The production modal basis builder selects the conservative weak projection. The strong Eq 41 projection remains available in the implementation. The resulting physical modes are normalized to `I_j(eta = 0) = 1` and vanish at the trapped passing boundary.

# Production use

The `fbis/modal/basis.py` caller builds a full `Lambda` to `eta` map, constructs the square mirror reference basis, evaluates Eq 42 weighted Lorentz coefficients on `Lambda(eta)`, and requests the conservative weak physical basis. The physical eigenfunctions are then used by modal source projection and by local reconstruction. Their boundary slopes are also converted back to `dI/dLambda` for the ion loss closure.

The basis convergence path rebuilds this chain at several `Lambda` and `eta` resolutions and compares physical eigenvalues and matched eigenfunctions. The diagnostics stored in `PhysicalEigenbasis` support those numerical qualification checks.

# Physics reference

The implemented scattering chain follows J. Egedal et al., *Fusion by beam ions in a low collisionality, high mirror ratio magnetic mirror*, Nuclear Fusion 62, 126053 (2022). The relevant development is the square mirror Lorentz eigenbasis in Eqs 4 through 8, the nonuniform mirror invariant and orbit averaged operator in Eqs 35 through 42, the `eta` coordinate in Eqs 45 through 49, and the physical eigenfunction expansion in Eqs 50 through 55.
