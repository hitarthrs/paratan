# Eq 59 velocity solver

This subpackage implements the hot ion velocity space solve based on Egedal Eq 59. 
It assembles speed dependent pairwise collision coefficients, solves the retained modal equations, iterates the Rosenbluth state, evaluates energy moments, and checks numerical convergence in the speed coordinate.

The Python files in this folder have not been annotation reviewed in this pass. This README records their current responsibilities from the implementation and callers.

## Files

`collision_operator_state.py`: Builds the speed dependent collision operator state for each active ordered field population and combines those pair contributions for Eq 59.

`convergence.py`: Constructs nested speed grids and evaluates Eq 59 convergence using modal populations, energy moments, upper speed tails, and mapped loss terms.

`cross_species.py`: Builds reduced collision states and Rosenbluth coefficients for an external non Maxwellian fast ion field population.

`energy_moments.py`: Evaluates discrete Eq 59 particle and energy balances, pairwise energy transfer, pitch loss power, charge exchange sink power, source power, and fast ion electron heating.

`maxwellian_pairs.py`: Evaluates analytic density normalized Rosenbluth functions for Maxwellian field populations.

`operator.py`: Normalizes Rosenbluth coefficients, constructs Scharfetter Gummel face coefficients, assembles the modal tridiagonal operator, solves each mode, and evaluates nonlinear residuals.

`solve.py`: Performs the nonlinear Eq 59 iteration, updating the collision operator from the evolving fast ion distribution until the configured convergence criteria are evaluated.

`state.py`: Builds and remaps warm start states, computes modal inventory and effective energy metrics, and classifies fixed point iteration behavior.

`types.py`: Defines Eq 59 convergence, warm start, speed convergence, and internal iteration metric dataclasses.

## Numerical representation

The speed grid is a `SpeedGrid` with exact spherical shell measures. Modal distributions are stored by retained mode and speed cell. Pairwise collision coefficients are evaluated on the same speed grid so the drift, diffusion, source, and loss terms can be assembled consistently.

Warm starts use conservative shell remapping when the speed grid changes. Convergence checks distinguish nonlinear fixed point convergence from speed domain and speed resolution qualification.

## Physics basis

Egedal Eq 59 extends the modal fast ion velocity equation to include energy diffusion through density normalized Rosenbluth potentials. The implementation keeps the modal pitch scattering eigenvalues separate from the speed dependent collision coefficients and retains explicit pairwise field population provenance.
