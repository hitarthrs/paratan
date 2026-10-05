# Modal FBIS solver orchestration

The `solver/` folder assembles the modal FBIS state from the magnetic basis, attenuated beam source, Eq 14 or Eq 59 velocity solve, ion loss closure, shared current balance, electrostatic reconstruction, and numerical diagnostics. The coupled system path is also used by the one species entry point so both interfaces share the same orchestration.

## Files

`solver_balance.py`: Builds the Eq 61 ion loss state, Eq 68 electron current balance, confinement times, wall power terms, and loss consistency diagnostics for one solved modal species.

`solver_components.py`: Projects attenuated beam births onto the physical eigenbasis, separates prompt magnetic losses from confined births, solves Eq 14 or Eq 59 in speed, reconstructs `f(v, η)`, and evaluates mode and speed convergence.

`solver_convergence.py`: Provides safe ratio helpers, upper speed particle and energy tail metrics, and source rate reconstruction from modal coefficients.

`solver_core.py`: Coordinates the complete one species modal solve, including basis reuse, source and velocity solution, loss and current balance, electrostatic reconstruction, metadata, warm start packaging, and Eq 59 energy moments.

`solver_fixed_boundary.py`: Solves the single species Eq 70 local reconstruction and builds the left and right fixed magnetic boundary Eq 63 branches when that loss closure is active.

`solver_loss_rates.py`: Reconstructs Eq 61 boundary losses, compares weak form and cell average boundary slopes, audits the modal eigenvalue sink, integrates directed Eq 63 throat flux, and selects the configured device loss representation.

`solver_metadata.py`: Collects solver settings, numerical convergence, collision state, source, loss, electrostatic, and local reconstruction diagnostics into the `ModalFBISResult` metadata.

`species_solver.py`: Routes a one species request through the shared system solver and returns the requested species state.

`system_solver.py`: Couples active fast deuterium and tritium species through a shared basis, Eq 68 wall barrier, Eq 70 potential, directed Eq 63 charge, and reduced fast cross species collision fields.

## Important conventions

The modal velocity state uses shape `(n_mode, n_speed)`. Reconstruction to the physical invariant grid uses `f(v, η)` with shape `(n_speed, n_eta)` and then `f(v, Λ)` with shape `(n_speed, n_lambda)`.

Only beam births inside the confined η interval enter the modal source. Births outside that interval are retained as prompt magnetic losses and are added to the total ion loss after the confined Eq 61 closure.

Eq 61 first reconstructs a one end boundary spectrum. The symmetric device loss uses both ends. The final stored boundary flux spectrum follows the selected device loss representation, while Eq 63 left and right branches are each matched to the corresponding one end Eq 61 spectrum.

The coupled system owns the shared electrostatic state. Eq 68 first obtains an electron barrier from the summed ion loss current. Eq 70 then determines the parent Maxwellian amplitude and axial quasineutral potential. The barrier and parent normalization are reconciled until the scalar reference is consistent with the axial profile.

Directed fixed magnetic boundary Eq 63 populations are attached before the shared Eq 70 solve so their charge can enter the quasineutrality closure. Local confined populations continue to use the Eq 71 accessibility boundary and Eq 72 invariant mapping.

## Physics basis

The solver combines the supplied Egedal et al. 2022 modal framework. Beam source projection and the cold reference use Eqs 9 through 15, the orbit averaged scattering basis uses Eqs 41, 42, 47, and 50, the hot ion speed solve uses Eq 59, and the ion loss chain uses Eqs 61 through 64. Electron current balance and axial quasineutrality use Eqs 68 through 70, while local confined ion accessibility uses Eqs 71 and 72.
