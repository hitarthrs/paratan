# ParaTAN Viewer

Interactive 3D geometry studio (Trame + PyVista) for ParaTAN YAML devices.
Builds analytic revolved meshes from the input, serves a browser UI for cutaways,
component picking, and optional OpenMC statepoint overlays.

## Dependencies

Use the repo conda environment ([`environment.yml`](../../../environment.yml)).
Viewer-critical packages:

| Package | Role |
|---------|------|
| `pyvista`, `vtk` | Mesh generation / rendering |
| `trame`, `trame-vtk`, `trame-vuetify` | Browser UI + WebGL bridge |
| `pyyaml`, `numpy` | Input parsing / numerics |

```bash
# from repository root
conda env create -f environment.yml          # first time
conda activate paratan
# or: conda env update -f environment.yml --prune
```

Headless / remote sessions often need OSMesa or a display. The launcher applies
GL workarounds via `gl_env.py` before VTK imports; do not import PyVista first.

OpenMC is **not** required to browse analytic geometry. It is required only when
loading a statepoint through **OpenMC results**.

## Launch

From the repository root (so `PYTHONPATH` resolves `src`):

```bash
# default: simple_parametric_input_new.yaml on http://127.0.0.1:8080/
./view_simple_mirror.py

# explicit module entry (same CLI)
PYTHONPATH=. python -m src.paratan.viewer \
  --input input_files/simple_parametric_input_new.yaml

# tandem example
PYTHONPATH=. python -m src.paratan.viewer \
  -i input_files/tandem_parametric_input.yaml
```

`view_simple_mirror.py` forwards all CLI flags to the module entry point.

| Flag | Default | Notes |
|------|---------|--------|
| `-i` / `--input` | `input_files/simple_parametric_input_new.yaml` | Device YAML |
| `--host` / `--port` | `127.0.0.1` / `8080` | Bind address |
| `--n-theta` | `96` | Angular tessellation (lower = faster load) |
| `--render` | `trame` | `trame` (local WebGL), `server`, `client` |
| `--timeout` | `2` | Exit seconds after last tab closes (`0` = keep alive) |
| `--no-browser` | off | Do not auto-open a tab |

Blank canvas: `--render server`, then hard-refresh (`Ctrl+Shift+R`).

## Usage notes

- Switch devices with the upload control or **Example devices** in the sidebar.
- Vacuum vessel styles for simple mirrors: `geometry_style: conical` or `perpendicular`.
- Toolbar: full / half / quarter cutaways and camera controls.
- **Geometry**: group visibility, Solo / Hide, material legend.
- **Inspect**: dimensions, tally config, cylindrical cages, φ/z slices.
- **OpenMC results**: statepoint + matching `simulation-manifest.json` (exported by
  the simple/tandem builders). Heating remains eV/source; no auto power normalize.
- Presets write next to the launch YAML as `<input>.viewer.json`.

## Layout

| Module | Role |
|--------|------|
| `__main__.py` / `gl_env.py` | CLI; GL env before VTK |
| `app.py` | Trame app, uploads, presets, materials catalog |
| `ui.py` | Browser chrome + client transport helpers |
| `scene.py` | VTK actors, slices, overlays |
| `model.py` | Cutaway cache, `SceneState` / actions |
| `devices.py` | Device adapters (simple / FFH / tandem) |
| `components.py` | `ComponentMesh`, labels, material colors |
| `simple_mirror_meshes.py` | YAML → axisymmetric meshes |
| `port_openings.py` | NBI / radial port meshes + surface trims |
| `revolve.py` | Analytic `(r,z)` profiles and sections |
| `results.py` | Cyl grids, tally overlays, statepoint datasets |
| `manifest.py` | Display / simulation identity JSON (no VTK) |

Deeper boundaries: [`docs/viewer_architecture.md`](../../../docs/viewer_architecture.md).
Tests: `tests/test_viewer_*.py`.

## Scope

- **Supported:** simple mirror (conical / perpendicular VV), FFH inlet/outlet pipes, tandem display geometry
- **Partial:** tandem VNS (envelope / under construction in UI)
- **Not:** full OpenMC CSG parity for every boolean exclusion
