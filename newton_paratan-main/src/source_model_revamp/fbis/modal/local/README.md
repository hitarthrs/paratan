# Local distribution reconstruction

This subpackage reconstructs physical local fast ion distributions from the invariant modal solution. It maps between modal eta, magnetic invariant Lambda, local speed, and local pitch while accounting for the configured electrostatic potential and accessibility constraints.


## Files

`mapping.py`: Converts modal coefficients to `f(v, η)` and `f(v, Λ)`, evaluates invariant moments, builds interpolation of the magnetic reference distribution, and implements the Eq 71 and Eq 72 invariant mappings.

`pitch.py`: Maps the invariant distribution to local `f(v_local, ξ)` cell averages, evaluates local density, performs conservative speed remapping, and builds the `ϕ = 0` magnetic reference reconstruction.

`speed.py`: Extends the local physical speed domain for electrostatic acceleration, evaluates density and kinetic energy moments, and checks local speed discretization convergence.

## Coordinate conventions

The magnetic reference distribution uses invariant speed and `Λ`. Local pitch is `ξ = v_parallel / v` and local pitch arrays use shape `(n_z, n_local_speed, n_pitch)`. Modal distributions first reconstruct `f(v, η)` with shape `(n_speed, n_eta)` before sampling onto `Λ`.

## Physics basis

The magnetic reference distribution follows the modal solution in the Egedal constant of motion variables. The electrostatic reconstruction uses the trapped passing boundary and compressed invariant mapping associated with Egedal Eqs 71 and 72 while retaining the implementation's explicit accessibility checks.