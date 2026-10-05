# Losses module

The `losses/` folder contains shared end loss utilities for the magnetic mirror source model. It defines electron electrostatic barrier conventions, electron end loss rate models, and adiabatic classification of directed lost ion trajectories between a mirror throat and a selected material surface.

## Files

`electron_loss_boundary.py`: Defines the electron barrier sign convention and conversions among wall potential, barrier energy, normalized barrier, and escape speed. It also provides Maxwellian number and flux tail factors and mean escaping electron energies.

`electron_loss_rate.py`: Evaluates reference Maxwellian surface losses and the Egedal electron tail refilling model. It includes the Eq 66 collision scale, the exact `E1(y)` particle tail integral, Eq 68 and Eq 69 high barrier comparisons, and electron particle, current, and energy loss states.

`expander_trajectories.py`: Traces directed throat crossing ion bins through supplied magnetic and electrostatic exhaust paths. It classifies each invariant energy and Λ bin as deposited, reflected, locally trapped, or unresolved and calculates the particle rate and kinetic power reaching the selected surface.

## Important conventions

Electron potential is measured relative to the midplane,

`Δϕ = ϕ_wall - ϕ_midplane`

and the positive electron barrier energy is

`E_barrier = e(ϕ_midplane - ϕ_wall) = -e Δϕ`

Electron temperature is represented as an energy in joules. The normalized barrier and parallel escape speed are

`y = E_barrier / T_e`

`v_w = sqrt(2 E_barrier / m_e)`

The Maxwellian number tail `0.5 erfc(sqrt(y))` and flux tail `exp(-y)` are distinct quantities. The Egedal tail refilling path evaluates the exact particle kernel `E1(y)` while retaining the Eq 68 asymptotic form `exp(-y) / y` and the Eq 69 high barrier wall energy limit for comparison.

For directed lost ions, `invariant_total_energy_J` has shape `(N_v,)`, `global_lambda` has shape `(N_Lambda,)`, and `directed_particle_rate_v_lambda_s` has shape `(N_v, N_Lambda)`. The traced parallel kinetic energy has shape `(N_v, N_Lambda, N_z)` and is evaluated from

`K_parallel(z) = U + P(z) - Λ U B_tilde(z)`

with `P(z) = -q ϕ(z)`. The first path point is the mirror throat and the final point is the selected surface. Only bins classified as deposited contribute to the reported surface particle rate and power.

## Physics basis

The electron loss model follows the supplied Egedal et al. 2022 formulation for collisional electron loss and ambipolar confinement. The implementation uses the electron collision scale from Eq 66, the loss rate construction leading to Eq 68, and the wall energy relation leading to Eq 69. The lost ion trajectory classifier uses the same invariant total energy `U` and pitch invariant `Λ` framework used for the Egedal lost ion and expander treatment.
