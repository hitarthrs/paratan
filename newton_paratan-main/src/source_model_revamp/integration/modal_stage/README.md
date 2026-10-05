# Modal FBIS integration stage

The modal_stage folder connects pipeline geometry and attenuated beam sources to the modal FBIS solver. It owns reusable Eq 42 basis construction, density closure, single and coupled fast ion stage orchestration, prompt only handling, full device population mapping, end loss power availability checks, and the typed density and metadata state returned to the operating point workflow.

## Files

`basis_density.py`: Converts kinetic configuration into modal numerics, builds uniform, prescribed, or self consistent Eq 42 density profiles, constructs and caches density weighted modal bases, evaluates optional basis convergence, and builds collision parameter states for supplied density compositions.

`density_closure.py`: Implements the single species deuterium density closure. It runs fixed collision state modal solves, forms the quasineutral Eq 42 density profile, compares density dependent bases, and iterates the Eq 42 density shape and basis to convergence.

`system_stage.py`: Implements the coupled fast D and fast T operating point. It iterates species densities, Eq 59 collision operators, reduced fast cross species Rosenbluth fields, the shared Eq 68 wall barrier, the shared Eq 70 profile, and the Eq 42 density weighted basis. It also assembles species source to loss particle audits and final convergence metadata.

`stage.py`: Provides the top level modal kinetic entry point. It routes NBI supported stationary closure and active tritium cases to the coupled system path and retains the deuterium benchmark path for the single species closure.

`prompt_source.py`: Represents the exact zero confined state when all deposited beam births are prompt magnetic losses. It retains prompt source and collision reference bookkeeping without inventing a confined kinetic population.

`full_device.py`: Maps confined local distributions to the full device grid and adds the active fixed magnetic boundary Eq 63 lost ion branches. It also evaluates whether directed throat rates can be traced through the expander geometry to selected material surfaces.

`end_loss_power.py`: Validates whether ion and electron end loss power terms are complete and satisfy the wall power identities required by the electron temperature energy balance.

`metadata_adapter.py`: Provides common modal stage convergence fields, ordered failure reasons, and the shared metadata representation used by single and coupled paths.

`types.py`: Defines the solved confined density state, density closure result, run local modal basis cache, and reusable modal basis bundle.

## Important conventions

The magnetic invariant coordinate is `Λ = μ B0 / E` and the active fixed magnetic ion loss boundary is `Λ_M = 1 / R_M`. Eq 71 and Eq 72 remain part of the local electrostatic accessibility mapping rather than replacing that fixed magnetic boundary.

Eq 42 weights the orbit averaged scattering operator with `n(z) / <n>`. The self consistent path begins from a uniform ratio, solves the quasineutral density state, rebuilds the density weighted physical basis, relaxes the density shape, and iterates until both density shape and basis changes satisfy their tolerances.

The coupled D plus T path shares one magnetic basis, Eq 68 current balance, and Eq 70 electrostatic profile. Fast D and fast T retain separate speed grids, Eq 59 modal distributions, collision states, local reconstructions, Eq 61 losses, and Eq 63 directed lost branches.

The coupled scalar fixed point tracks fast species midpoint and volume average densities, electron midpoint and collision densities, Eq 59 drag, energy diffusion, and pitch scattering arrays, the wall barrier, and reduced fast cross species Rosenbluth fields. Warm states can continue these quantities across nearby operating point evaluations only when their species, grids, and pair identities remain compatible.

`OperatingPointDensityState` stores confined cell and exact electrostatic node profiles. Electron and total positive charge densities are checked at the exact midplane and by confined volume integrated inventory.

Prompt magnetic losses are separated from the confined modal source before the kinetic solve. A source with no confined births returns exact zero confined arrays and does not fabricate Eq 68 or Eq 70 closure states.

Full device mapping keeps throat crossing distinct from material deposition. The central Eq 63 lost population is reconstructed from the fixed magnetic boundary state, while the later expander stage supplies the source connected expander closure.

## Physics basis

The stage organizes the supplied Egedal et al. modal FBIS framework without changing its underlying equations. The magnetic basis uses the orbit averaged Lorentz operator and Eq 42 density weighting, the velocity solve uses Eq 59, the fixed magnetic ion loss chain uses Eqs 61 through 64, electron current balance uses Eq 68 and Eq 69, and the quasineutral potential and local accessibility mapping use Eqs 70 through 72.

The coupled D plus T path applies these same modal equations to separate fast species while sharing the magnetic basis and electrostatic state. Reduced fast cross species collision fields are constructed from the solved species states and iterated with the scalar density closure.

The full device lost ion mapping uses the Eq 63 invariant total energy description together with the solved electrostatic potential before the separate expander trajectory model classifies wall connected, reflected, local well, and unresolved branches.
