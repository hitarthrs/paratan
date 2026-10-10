"""Trame presentation layer. Widgets dispatch application actions."""
from __future__ import annotations
import base64
from pathlib import Path
from trame.ui.vuetify3 import SinglePageLayout
from trame.widgets import client, html, vuetify3 as vuetify
from pyvista.trame.ui import plotter_ui
from src.paratan.viewer.model import GROUPS, MODE_LABELS
from src.paratan.viewer.labels import GROUP_LABELS
from trame_client.widgets.core import AbstractElement
from src.paratan.viewer import transport

class ViewerUI:
    def _build_ui(self, render_mode: str) -> None:
        server, ctrl, pl = self.server, self.ctrl, self.pl
        server.enable_module(transport)
        with SinglePageLayout(server) as layout:
            logo_path = Path(__file__).with_name("static") / "paratan-icon.svg"
            logo_uri = "data:image/svg+xml;base64," + base64.b64encode(logo_path.read_bytes()).decode("ascii")
            server.state.trame__favicon = logo_uri
            layout.title.set_text("")
            with layout.title:
                with html.Div(style="display:flex;align-items:center;gap:10px;white-space:nowrap;"):
                    html.Img(src=logo_uri, alt="Paratan logo", width=32, height=32,
                             style="display:block;flex-shrink:0;")
                    html.Span("PARATAN · Geometry studio")
            layout.root.style = "height:100vh; overflow:hidden; background:#f1f5f9; color:#172b4d;"
            with vuetify.VSnackbar(v_model=('upload_notice', False), timeout=8000,
                                  color=("upload_phase === 'error' ? 'error' : 'success'",), location='bottom right'):
                html.Span('{{ load_status }}', role='status', aria_live='polite')
                with vuetify.Template(v_slot_actions='{}'):
                    vuetify.VBtn('Close', click='upload_notice = false', variant='text')
            with layout.toolbar as toolbar:
                toolbar.density = "compact"
                vuetify.VSpacer()
                vuetify.VChip("{{ n_components }} parts · cm", size="small", variant="tonal")
                with vuetify.VMenu():
                    with vuetify.Template(v_slot_activator="{ props }"):
                        vuetify.VBtn(
                            "{{ slice_phi_on || slice_z_on ? 'Cut: slice' : (section_mode === 'half' ? 'Cut: half' : "
                            "(section_mode === 'quarter' ? 'Cut: quarter' : 'Cut')) }}",
                            v_bind="props", variant="text", size="small", append_icon="mdi-chevron-down",
                            loading=('!cuts_ready',), disabled=('!cuts_ready',))
                    with vuetify.VList(density="compact"):
                        for value, label in MODE_LABELS:
                            vuetify.VListItem(title=label, click=(ctrl.set_section_mode, f"['{value}']"),
                                              active=(f"section_mode === {value!r}",),
                                              disabled=("slice_phi_on || slice_z_on",))
                        html.Div("A slice plane (Inspect) sets the cut while it is on.",
                                 v_if="slice_phi_on || slice_z_on", classes="text-caption px-4 pb-2",
                                 style="color:#64748b;max-width:220px;")
                for label, action in (("Side", ctrl.cam_side), ("End", ctrl.cam_end), ("Iso", ctrl.cam_iso),
                                      ("Fit selection", ctrl.cam_fit_selected), ("Reset", ctrl.reset_view)):
                    vuetify.VBtn(label, click=action, variant="text", size="small")
            with layout.content:
                layout.content.style = "height:100vh; padding:48px 0 22px; box-sizing:border-box; overflow:hidden;"
                with html.Div(style="display:flex; height:100%;"):
                    with html.Div(style="width:380px; min-width:300px; max-width:45vw; height:100%; overflow-y:auto; "
                                        "background:white; border-right:1px solid #dbe3ed; padding:16px;"):
                        self._panel_header()
                        with vuetify.VTabs(v_model=("panel", "geometry"), density="compact", color="primary"):
                            vuetify.VTab("Geometry", value="geometry")
                            vuetify.VTab("Inspect", value="inspect")
                        self._geometry_panel()
                        self._inspect_panel()
                        self._results_panel()
                        vuetify.VDivider(classes="my-4")
                        self._views_panel()
                    with html.Div(style="flex:1; min-width:0; position:relative; height:100%;"):
                        self.view = plotter_ui(
                            pl, mode=render_mode, default_server_rendering=render_mode == "server",
                            add_menu=False, server=server, still_ratio=1.25, interactive_ratio=0.4,
                            still_quality=85, interactive_quality=40,
                            **({"EndAnimation": (self._sync_camera, f"[trame.refs['view_{pl._id_name}'].getCamera()]"),
                                "interactor_events": ("['EndAnimation']",),
                                # Local WebGL never reaches the server-side picker: pick in the browser.
                                "picking_modes": ("['click']",),
                                "click": (self._on_view_click, "[$event]")} if render_mode != "server" else {}),
                            style="width:100%;height:100%;")
                        with self.view:
                            AbstractElement('paratan-scene-transport')
                        ctrl.view_update = self.view.update
                        with html.Div(v_if='upload_busy', role='status', aria_live='polite',
                                      style='position:absolute;inset:0;z-index:10;background:rgba(241,245,249,.88);display:flex;align-items:center;justify-content:center;'):
                            with html.Div(style='text-align:center;padding:24px;max-width:420px;'):
                                vuetify.VProgressCircular(indeterminate=True, color='primary', size=36)
                                html.Div('{{ load_status }}', style='margin-top:12px;font-size:14px;font-weight:600;overflow-wrap:anywhere;')
                                html.Div('Preparing the viewer. This may take a moment.', style='margin-top:6px;font-size:12px;color:#64748b;')
                        # Server-invoked page reload: lets the browser resubscribe after a model swap.
                        self._page_reload = client.JSEval(exec="utils.get('location').reload()")
                        with html.Div(v_if="tally_summary.active", style="position:absolute;top:16px;left:16px;right:16px;"
                                      "max-width:340px;background:rgba(255,255,255,.95);border:1px solid #e2e8f0;"
                                      "padding:14px 16px;border-radius:12px;box-shadow:0 4px 20px rgba(15,23,42,.08);pointer-events:none;"):
                            html.Div('{{ tally_summary.source }}', style='font-size:10px;color:#b45309;font-weight:700;letter-spacing:.6px;')
                            html.Div('{{ selected_label }}', style='font-size:15px;font-weight:600;margin:5px 0;')
                            html.Div('{{ tally_summary.bins }}', style='font-size:12px;color:#475569;')
                            html.Div('{{ tally_summary.extent }}', style='font-size:11px;color:#64748b;margin-top:4px;')
                            html.Div('{{ tally_summary.quantity }} · {{ tally_summary.range }}', v_if='tally_summary.range',
                                     style='font-size:12px;color:#0f172a;margin-top:8px;')
                        html.Div("{{ selected_label || 'Orbit to explore · select a part to inspect' }}",
                                 style="position:absolute;bottom:18px;left:18px;background:rgba(255,255,255,.92);"
                                       "padding:8px 12px;border-radius:8px;font-size:12px;pointer-events:none;")

    def _panel_header(self) -> None:
        html.Div("{{ model_name }}", classes="text-subtitle-2 mb-1")
        html.Div("{{ device_label }} · analytic geometry", classes="text-caption mb-3")
        with html.Div(v_if="device_under_construction", classes="mb-3",
                      style="padding:10px 12px;background:#fffbeb;border:1px solid #fde68a;border-radius:6px;color:#92400e;"):
            html.Div('Under construction', style='font-size:12px;font-weight:600;')
            html.Div('Tandem and VNS previews are incomplete and still being developed.',
                     style='font-size:11px;margin-top:4px;')
        upload_event = ("const input = $event.target; const f = input.files[0]; if (f && !upload_busy) { "
                        "upload_busy = true; upload_phase = 'loading'; upload_notice = false; "
                        "load_status = 'Loading ' + f.name + '…'; "
                        "if (f.size > 2000000) { upload_busy = false; upload_phase = 'error'; "
                        "load_status = f.name + ' is larger than 2 MB.'; upload_notice = true; input.value = ''; } else { "
                        "const Controller = utils.get('AbortController'); const aborter = new Controller(); "
                        "const timer = utils.get('setTimeout')(() => aborter.abort(), 30000); "
                        "utils.get('fetch')('/paratan-upload?name=' + utils.get('encodeURIComponent')(f.name), "
                        "{method: 'POST', body: f, signal: aborter.signal})"
                        ".then(response => response.json())"
                        ".then(result => { upload_phase = result.phase; load_status = result.message; upload_notice = true; })"
                        ".catch(error => { upload_phase = 'error'; upload_notice = true; "
                        "load_status = 'Could not upload ' + f.name + ': ' + error.message + '. Please try again.'; })"
                        ".finally(() => { utils.get('clearTimeout')(timer); input.value = ''; upload_busy = false; }); } }")
        with html.Div(style='padding:10px;border:1px solid #cbd5e1;border-radius:6px;'):
            html.Label('Load input file (.yaml)', for_='device-input-file',
                       style='display:block;font-size:11px;color:#64748b;margin-bottom:6px;')
            html.Input(id='device-input-file', type='file', accept='.yaml,.yml',
                       change=upload_event, disabled=('upload_busy',),
                       style='width:100%;font-size:11px;color:#334155;')
        vuetify.VProgressLinear(v_if='upload_busy', indeterminate=True, color='primary', classes='mt-2')
        vuetify.VAlert('{{ load_status }}', v_if='load_status', density='compact', variant='tonal',
                       type=("upload_phase === 'error' ? 'error' : (upload_phase === 'success' ? 'success' : 'info')",),
                       classes='mt-2', style='font-size:12px;overflow-wrap:anywhere;', role='status', aria_live='polite')
        with html.Details(classes='mt-2'):
            html.Summary('Example devices', style='font-size:11px;color:#64748b;cursor:pointer;')
            with html.Div(classes='d-flex flex-wrap mt-2', style='gap:4px;'):
                for key, label in [('simple','Simple mirror'),('tandem','Tandem · Under construction'),('tandem_vns','Tandem VNS · Under construction')]:
                    vuetify.VBtn(label, click=(self.ctrl.load_example, f"['{key}']"), disabled=('upload_busy',), size='x-small', variant='tonal')
        html.A('Download model manifest', href=('manifest_url',), download='display-manifest.json', classes='text-caption d-block mt-2')
        html.Div(classes="mb-3")

    def _geometry_panel(self) -> None:
        ctrl = self.ctrl
        with html.Div(v_show="panel === 'geometry'", classes="pt-4"):
            vuetify.VBtn('Explore HF mesh tally', click=ctrl.view_hf_tally, variant='tonal', color='primary',
                         block=True, classes='mb-3', v_if="items_HF.length")
            vuetify.VTextField(v_model=("search", ""), label="Search components or materials", density="compact",
                               variant="outlined", clearable=True, hide_details=True)
            with html.Div(style="max-height:50vh; overflow-y:auto; margin-top:8px;"):
                for g in GROUPS:
                    with html.Div(v_if=f"items_{g}.length", classes="mb-1"):
                        with html.Div(style="display:flex;align-items:center;gap:2px;"):
                            vuetify.VCheckbox(v_model=(f"vis_{g}", True), density="compact", hide_details=True,
                                              title="Show / hide group", style="flex:none;")
                            html.Div(f"{GROUP_LABELS[g].upper()} · {{{{ items_{g}.length }}}}", classes="text-overline")
                        with html.Div(v_for=f"item in items_{g}", key="item.name"):
                            with html.Div(v_show="(item.label + ' ' + item.material).toLowerCase()"
                                                 ".includes((search || '').toLowerCase())"):
                                vuetify.VBtn("{{ item.label }}", click=(ctrl.select_component, "[item.name]"),
                                             block=True, size="small",
                                             variant=("selected_name === item.name ? 'tonal' : 'text'",),
                                             color=("selected_name === item.name ? 'primary' : undefined",),
                                             classes="text-none justify-start mb-1", style="font-size:12px;")
            with html.Div(classes="d-flex flex-wrap mt-3", style="gap:6px;"):
                for label, action in (("Solo", ctrl.toggle_solo), ("Hide", ctrl.hide_selected),
                                      ("Clear", ctrl.clear_selection)):
                    vuetify.VBtn(label, click=action, size="small", variant="tonal", disabled=("!has_selection",))
                vuetify.VBtn("Show all", click=ctrl.show_all_components, size="small", variant="text")
            with html.Details(classes="mt-4"):
                html.Summary("Appearance", style="cursor:pointer;font-weight:600;")
                vuetify.VSlider(v_model=("explode", 0), min=0, max=1, step=0.05, label="Explode", classes="mt-4",
                                hide_details=True)
                vuetify.VSlider(v_model=("selected_opacity", 1), min=0, max=1, step=0.05,
                                label="Selection opacity", disabled=("!has_selection",), hide_details=True)
                html.Div("MATERIALS", classes="text-overline mt-3")
                with html.Div(v_for="item in legend_items", key="item.name", classes="d-flex align-center mb-2",
                              style="gap:8px;"):
                    html.Div(style=("`width:10px;height:10px;border-radius:50%;flex:none;background:${item.color}`",))
                    html.Div("{{ item.name }}", classes="text-caption")

    def _inspect_panel(self) -> None:
        def caption(text: str) -> None:
            html.Div(text, style="font-size:10px;font-weight:600;color:#94a3b8;letter-spacing:.7px;margin-bottom:3px;")

        with html.Div(v_show="panel === 'inspect'", classes="pt-4"):
            with html.Div(v_if="!has_selection", style="padding:24px 16px;background:#f8fafc;border:1px dashed #cbd5e1;"
                                                     "border-radius:12px;color:#64748b;font-size:13px;line-height:1.6;"):
                html.Div("Select a component", style="font-weight:600;color:#334155;margin-bottom:4px;")
                html.Div("Choose a part from Geometry to see its dimensions and configured tallies.")
            with html.Div(v_if="has_selection"):
                with html.Div(style="padding:16px;border:1px solid #e2e8f0;border-radius:12px;background:#fff;"
                                    "box-shadow:0 2px 6px rgba(15,23,42,.03);"):
                    html.Span("{{ inspector_group }}", style="display:inline-block;font-size:10px;font-weight:700;"
                              "padding:3px 7px;border-radius:5px;background:#eff6ff;color:#2563eb;"
                              "letter-spacing:.5px;margin-bottom:8px;")
                    html.Div("{{ inspector_title }}", style="font-size:17px;font-weight:650;color:#0f172a;"
                                                           "line-height:1.35;margin-bottom:14px;")
                    caption("MATERIAL")
                    html.Div("{{ inspector_material }}", style="font-size:13px;color:#334155;margin-bottom:12px;")
                    with html.Div(v_for="row in inspector_rows", key="row.label",
                                  style="display:flex;justify-content:space-between;gap:12px;padding:7px 0;"
                                        "border-top:1px solid #f1f5f9;font-size:12px;"):
                        html.Span("{{ row.label }}", style="color:#64748b;")
                        with html.Span(style="color:#0f172a;font-weight:600;font-variant-numeric:tabular-nums;text-align:right;"):
                            html.Span("{{ row.value }}")
                            html.Span(" {{ row.unit }}", style="font-weight:400;color:#94a3b8;")
                with html.Div(style="display:flex;align-items:center;justify-content:space-between;margin:22px 0 10px;"):
                    html.Div("Tallies", style="font-size:15px;font-weight:600;color:#0f172a;")
                    html.Span("{{ inspector_tallies.length }} configured",
                              style="font-size:11px;color:#64748b;background:#f1f5f9;padding:4px 8px;border-radius:20px;")
                html.Div("Tally configuration · load a statepoint below to inspect results",
                         style="font-size:11px;line-height:1.5;color:#64748b;margin-bottom:12px;")
                html.Div("{{ inspector_path }}", v_if="inspector_path",
                         style="font-size:10px;color:#94a3b8;overflow-wrap:anywhere;margin-bottom:12px;")
                html.Div("No active tallies attached to this component.", v_if="!inspector_tallies.length",
                         style="font-size:12px;color:#64748b;padding:14px;background:#f8fafc;border-radius:8px;")
                with html.Div(v_for="tally in inspector_tallies", key="tally.title",
                              style="border:1px solid #e2e8f0;border-radius:10px;padding:14px;margin-bottom:10px;background:#fff;"):
                    html.Div("{{ tally.title }}", style="font-size:12px;font-weight:650;color:#334155;margin-bottom:10px;")
                    with html.Div(style="display:flex;flex-wrap:wrap;gap:5px;margin-bottom:12px;"):
                        html.Span("{{ score }}", v_for="score in tally.scores", key="score",
                                  style="font-size:11px;font-weight:600;color:#1d4ed8;background:#eff6ff;"
                                        "border-radius:5px;padding:4px 8px;")
                    caption("FILTERS")
                    html.Div("{{ tally.filters.join(' · ') }}", style="font-size:12px;color:#475569;line-height:1.5;margin-bottom:9px;")
                    caption("NUCLIDES")
                    html.Div("{{ tally.nuclides }}", style="font-size:12px;color:#475569;line-height:1.5;")
                    with html.Div(v_if="tally.dimensions.length === 3",
                                  style="display:grid;grid-template-columns:repeat(3,1fr);gap:6px;margin-top:12px;"):
                        for index, axis in enumerate(("Radial", "Azimuthal", "Axial")):
                            with html.Div(style="text-align:center;background:#f8fafc;border-radius:6px;padding:8px 4px;"):
                                html.Div("{{ tally.dimensions[" + str(index) + "] }}",
                                         style="font-size:16px;font-weight:600;color:#0f172a;font-variant-numeric:tabular-nums;")
                                html.Div(axis, style="font-size:9px;color:#94a3b8;margin-top:2px;")
                    vuetify.VBtn("{{ (preview_tally && tally_mesh_index === tally.index ? 'Showing ' : 'View ') + tally.scores.join(' / ') + ' mesh' }}",
                                 v_if='tally.dimensions.length === 3',
                                 click=(self.ctrl.view_cylindrical_tally, '[tally.index]'),
                                 variant='tonal', color='primary', size='small', block=True, classes='mt-3')
                self._slice_card()
                with html.Div(v_if="tally_meshes.length", style="margin-top:20px;padding-top:16px;border-top:1px solid #e2e8f0;"):
                    html.Div("Mesh coverage", style="font-size:13px;font-weight:600;color:#334155;margin-bottom:14px;")
                    vuetify.VBtn('View mesh on geometry', click=self.ctrl.view_cylindrical_tally,
                                 variant='tonal', color='primary', block=True, classes='mb-3')
                    vuetify.VSelect(v_model=("tally_mesh_index", 0), items=("tally_meshes",), label="Preview mesh",
                                    density="compact", variant="outlined", hide_details=True, classes="mb-3",
                                    style="font-size:12px;")
                    vuetify.VCheckbox(v_model=("preview_tally", False), label="Show bin grid", density="compact",
                                      hide_details=True)
                    html.Div('Viewing a tally opens a quarter cut with the mesh overlaid on the geometry.',
                             style='font-size:11px;color:#64748b;margin-top:6px;')
                    with html.Div(style="padding:12px;background:#fffbeb;border:1px solid #fde68a;border-radius:8px;margin:12px 0;"):
                        html.Div("DEMO HEATING MAP", style="font-size:10px;font-weight:700;color:#92400e;letter-spacing:.6px;")
                        html.Div("Invented bin values · arbitrary units · no OpenMC results loaded",
                                 style="font-size:11px;color:#a16207;line-height:1.5;margin-top:4px;")
                        vuetify.VCheckbox(v_model=("demo_heating", False), label="Show demo heating map",
                                          density="compact", hide_details=True, disabled=("!heating_available",))
                        html.Div("Select a mesh tally with the heating score to enable the demo.",
                                 v_if="!heating_available", style="font-size:11px;color:#a16207;")
                        vuetify.VSlider(v_model=("demo_opacity", 1.0), label="Map opacity", min=.1, max=1, step=.1,
                                        hide_details=True, density="compact", v_if="result_source",
                                        classes="mt-2")
                    html.Div("Grid lines mark bin boundaries (large meshes are thinned for display; configured bin "
                             "counts are unchanged).", style="font-size:11px;line-height:1.6;color:#94a3b8;margin:8px 0 18px;")
            with html.Details():
                html.Summary("Tally setup example", style="cursor:pointer;font-weight:600;")
                html.Pre("{{ tally_example }}", style="white-space:pre-wrap;font-size:11px;background:#f1f5f9;padding:10px;")
                html.Div("Cell tallies aggregate over a component; mesh tallies resolve space. Dimensions are radial, "
                         "azimuthal, axial. Filters select particles or energy groups. Scores are per source particle; "
                         "physical rates require source normalization. Rebuild and run the model to obtain results.",
                         classes="text-caption")

    def _slice_card(self) -> None:
        with html.Div(v_if="has_slice", style="margin-top:20px;padding-top:16px;border-top:1px solid #e2e8f0;"):
            html.Div("Slice planes", style="font-size:13px;font-weight:600;color:#334155;margin-bottom:4px;")
            html.Div("Isolates this component and opens it at a chosen azimuth or height; the bins on that plane are shown.",
                     style="font-size:11px;line-height:1.5;color:#64748b;margin-bottom:8px;")
            for key, title, readout in (("phi", "Azimuthal plane (r–z map)", "Number(slice_phi).toFixed(1) + '°'"),
                                        ("z", "Axial plane (r–φ map)", "Number(slice_z).toFixed(1) + ' cm'")):
                vuetify.VSwitch(v_model=(f"slice_{key}_on", False), label=title, density="compact",
                                hide_details=True, color="primary")
                with html.Div(v_if=f"slice_{key}_on", style="padding:0 4px 8px;"):
                    limits = ({"min": 0, "max": 360, "step": 1} if key == "phi" else
                              {"min": ("slice_z_min", 0), "max": ("slice_z_max", 1), "step": 0.5})
                    vuetify.VSlider(v_model=(f"slice_{key}", 0), hide_details=True, density="compact",
                                    color="primary", **limits)
                    with html.Div(style="display:flex;align-items:center;gap:6px;"):
                        vuetify.VBtn("◀ bin", click=(self.ctrl.step_slice, f"['{key}', -1]"), size="x-small", variant="tonal")
                        vuetify.VBtn("bin ▶", click=(self.ctrl.step_slice, f"['{key}', 1]"), size="x-small", variant="tonal")
                        html.Span("{{ " + readout + " }}",
                                  style="font-size:12px;font-weight:600;font-variant-numeric:tabular-nums;")
            with html.Div(v_if="slice_phi_on || slice_z_on"):
                html.Div("{{ slice_info }}", style="white-space:pre-line;font-size:11px;color:#475569;line-height:1.6;margin:6px 0;")
                vuetify.VBtn("Face the slice", click=self.ctrl.face_slice, size="small", variant="tonal")

    def _results_panel(self):
        with html.Details(classes='mt-4'):
            html.Summary('OpenMC results', style='cursor:pointer;font-weight:600;')
            html.Div('Load a local statepoint and its simulation-manifest.json. No bins are summed automatically.', classes='text-caption my-2')
            for key, label in [('statepoint_path','Statepoint path (.h5)'), ('simulation_manifest_path','Simulation manifest path'),
                               ('statepoint_score','Score'), ('statepoint_nuclide','Nuclide'), ('statepoint_filter_bins','Extra filter bin indices (JSON)')]:
                vuetify.VTextField(v_model=(key,), label=label, density='compact', variant='outlined', hide_details=True, classes='mb-2')
            vuetify.VTextField(v_model=('statepoint_tally_id',0), label='Tally ID', type='number', density='compact', variant='outlined', hide_details=True, classes='mb-2')
            vuetify.VBtn('Load result', click=self.ctrl.load_results, size='small', variant='tonal')
            vuetify.VBtn('Clear result', click=self.ctrl.clear_results, size='small', variant='text')
            html.Div('{{ result_status }}', role='status', classes='text-caption mt-2', style='overflow-wrap:anywhere;')
        vuetify.VSelect(v_model=('result_id',''), items=('result_choices',), label='Loaded results for this component',
                           density='compact', variant='outlined', hide_details=True, v_if='result_choices.length', classes='mt-3')
        with html.Div(v_if='result_source', classes='mt-3'):
            html.Div('{{ result_summary }}', classes='text-caption mb-2')
            vuetify.VSelect(v_model=('result_quantity','mean'), label='Quantity',
                           items=([{'title':'Mean','value':'mean'}, {'title':'Standard deviation','value':'std_dev'},
                                   {'title':'Relative uncertainty','value':'relative_error'}],),
                           density='compact', variant='outlined', hide_details=True)

    def _views_panel(self) -> None:
        ctrl = self.ctrl
        with html.Details():
            html.Summary("Saved views & export", style="cursor:pointer;font-weight:600;margin-bottom:12px;")
            vuetify.VTextField(v_model=("preset_name", ""), label="View name", density="compact",
                               variant="outlined", hide_details=True)
            vuetify.VBtn("Save current view", click=ctrl.save_preset, size="small", variant="tonal", classes="my-2")
            vuetify.VSelect(v_model=("preset_selected", ""), items=("preset_names",), label="Saved views",
                            density="compact", variant="outlined", hide_details=True)
            with html.Div(classes="d-flex my-2", style="gap:6px;"):
                vuetify.VBtn("Restore", click=ctrl.restore_preset, size="small", disabled=("!preset_selected",))
                vuetify.VBtn("Delete", click=ctrl.delete_preset, size="small", variant="text", disabled=("!preset_selected",))
            html.A("Download views JSON", href=("presets_url",), download="paratan-views.json", v_if="presets_url",
                   classes="text-caption d-block mb-2")
            with html.Details():
                html.Summary("Import views JSON", classes="text-caption")
                vuetify.VTextarea(v_model=("preset_import", ""), label="Paste exported JSON", rows=3,
                                  variant="outlined", density="compact")
                vuetify.VBtn("Import", click=ctrl.import_presets, size="small")
            vuetify.VBtn("Prepare PNG", click=ctrl.export_png, size="small", variant="tonal", classes="mt-3")
            html.A("Download PNG", href=("png_url",), download="paratan-view.png", v_if="png_url",
                   classes="text-caption d-block mt-2")
            html.Div("{{ action_status }}", role="status", classes="text-caption mt-2", style="overflow-wrap:anywhere;")
