# ParaTAN Viewer

Interactive 3D viewer for ParaTAN simple-mirror geometry (Trame + PyVista).
Loads a parametric YAML, builds revolved CSG-style meshes, and serves a browser UI
for cutaways, component picking, and optional OpenMC statepoint overlays.

## Run

From the `paratan-development` repo root:

```bash
PYTHONPATH=. python -m src.paratan.viewer \
  --input input_files/simple_parametric_input_new.yaml
```

Then open `http://127.0.0.1:8080/`.

Useful flags:

| Flag | Default | Notes |
|------|---------|--------|
| `-i / --input` | `input_files/simple_parametric_input_new.yaml` | Simple-mirror YAML |
| `--host` / `--port` | `127.0.0.1` / `8080` | Bind address |
| `--n-theta` | `96` | Angular tessellation (lower = faster) |
| `--render` | `trame` | `trame` (snappy), `server` (if blank canvas), `client` |
| `--timeout` | `2` | Exit this many seconds after the last browser tab closes (`0` = keep alive) |
| `--no-browser` | off | Do not auto-open a tab |

If the canvas is blank: try `--render server`, then hard-refresh (`Ctrl+Shift+R`).

## Layout

| Module | Role |
|--------|------|
| `__main__.py` | CLI entry; applies GL env before VTK import |
| `app.py` | Trame app, state bridge, upload / presets |
| `model.py` | Device model, groups, visibility |
| `simple_mirror_meshes.py` | YAML → component meshes |
| `scene.py` | PyVista scene / cutaways / camera |
| `ui.py` | Browser UI |
| `results.py` / `tally_overlay.py` | Statepoint + tally overlays |
| `state.py` | Scene state + actions |
| `geometry helpers` | `revolve.py`, `cyl_grid.py`, `mesh_primitives.py`, `port_openings.py` |

Heavy viz deps are imported lazily so non-viewer code can import the package lightly.

## Develop here next

1. Prefer editing `model.py` / `simple_mirror_meshes.py` for geometry, `scene.py` / `ui.py` for UX.
2. Saved camera presets land next to the input YAML as `<input>.viewer.json`.
3. Viewer-focused tests live under `tests/test_viewer_*.py` (repo root).
4. Keep GL setup in `gl_env.py` — it must run before any PyVista/VTK import.

## Scope today

- **Supported:** simple mirror from ParaTAN YAML
- **Not yet:** full tandem-mirror parity in this UI
