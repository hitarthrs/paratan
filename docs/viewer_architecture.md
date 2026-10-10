# Viewer platform

The viewer supports simple mirrors and tandem devices through one application.
Run `python -m src.paratan.viewer --input input_files/tandem_parametric_input.yaml`
from the repository root, load another YAML with the input file control, or use
the bundled models under **Example devices**.

## Boundaries

- `state.py`: canonical `SceneState`, validated `SceneActions`, and the Trame
  projection. State and actions can be exercised without a browser or VTK.
- `devices.py`: `DeviceAdapter` protocol and simple/tandem adapters. They translate
  device-specific inputs into common components and tally configurations.
- `model.py`: caches full/half/quarter geometry. The scene preloads all cutaways
  so switching cuts only changes actor transforms. Browser array requests share
  in-flight downloads, and uploads update the live scene without a page reload.
- `scene.py`: owns VTK actors, dynamic slices, overlays, and result rendering.
  It receives ordinary Python state, not Trame state. Geometry signatures separate
  structural changes from opacity, visibility, and explode updates.
- `results.py`: immutable cylindrical coordinates and tally arrays. Array order
  is always `(r, phi, z)`. Synthetic and OpenMC data use `TallyDataset`.
- `manifest.py`: versioned, JSON-safe component, cell, and tally identities.
  Importing it from a simulation builder does not import VTK or Trame.
- `presets.py`: validates an entire saved view before it affects the scene.
- `ui.py`: widgets and presentation only; callbacks dispatch application actions.
- `app.py`: coordinates model loading, the state bridge, result loading, camera,
  presets, and export. Existing `create_app`, `build_model`, and launch entry
  points remain available.

To add another device, implement `DeviceAdapter.build_components`, register it in
`ADAPTERS`, and extend input detection in `manifest.device_kind`. Components use
stable names, common groups, materials, and optional axisymmetric profiles.
Profiles provide the existing cutaway, cage, and slice capabilities automatically.
Tally bindings are adapter metadata rather than assumptions in the renderer.

## Model identities and simulation manifests

The viewer's downloadable display manifest describes analytic meshes, their
bounds, stable device-prefixed component IDs, and configured tallies. It deliberately
contains no guessed OpenMC cell IDs.

The standard simple-mirror and tandem model-building entry points also write
`simulation-manifest.json` in their output directory. This records actual cell
IDs and actual generated tally IDs, with cylindrical coordinates and origin.
Both manifests use a canonical input-data fingerprint. Formatting a YAML file
alone does not change this fingerprint. Saved camera presets separately retain
raw-file fingerprints for compatibility with the previous viewer.

Simulation geometry and analytic display geometry are different representations.
Neighboring-component boolean exclusions are not universally reproduced by the
preview. The tandem VNS input displays its central test-region envelope; individual
VNS test modules are not modeled by that adapter.

## Loading real results

1. Rebuild a simple or standard tandem model to generate its simulation manifest.
2. Run OpenMC normally to create a statepoint.
3. Open the matching device YAML in the viewer and select a component.
4. Expand **OpenMC results**. Enter the statepoint path, simulation-manifest path,
   bound tally ID, score, and nuclide (`total` when appropriate).
5. If non-mesh filters have multiple bins, supply explicit zero-based filter-index
   to bin-index selections, e.g. `{"0": 1, "2": 0}`. This refers to the filter
   order in the tally, which is listed in the simulation manifest. No bins are
   silently summed.
6. Load the result. Use mean, standard deviation, or relative uncertainty; the
   same full-volume/cutaway and azimuthal/axial slice tools apply.

The loader validates device/input identity, component/tally binding, mesh edges,
origin, score, and nuclide. OpenMC mesh-bin tuples map file ordering into the
common `(r, phi, z)` arrays. `TallyDataset.save/load` uses compressed NPZ without
pickle; this is a transport format for future result providers, not a replacement
for an OpenMC statepoint.

Native heating is **eV per source particle**. There is no implicit conversion to
watts or W/cm³. Standard deviation has the same units; relative uncertainty is a
fraction and is undefined at zero mean (NaN), never shown as zero uncertainty.
Synthetic heating stays labeled demo data in arbitrary units and provides no
invented uncertainty estimates.

The first reader supports one cylindrical MeshFilter, full azimuth, and meshes
centered on the machine axis. Nonuniform bin edges and axial origins are supported.
Off-axis origins and partial-azimuth grids are rejected instead of mispositioned.
Large real datasets are not automatically coarsened or averaged; deployment memory
and rendering budgets still need measuring with production statepoints.

## Verification

Run `python -m unittest discover -s tests -p 'test_viewer*.py'`.
Tests cover geometry, adapter schemas, actions/reset, immutable results, nonuniform
bins, OpenMC filter ordering, HDF5 loading, manifests, presets, and actor retention.
The HDF5 test is a synthetic OpenMC-format fixture, not a transport simulation;
production-statepoint validation remains necessary before scientific use.
