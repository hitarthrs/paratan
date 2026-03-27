"""
Paratan 2 — Simple Mirror Machine Builder (Desktop GUI)
Run: python gui.py
"""

import os
import sys
import tempfile
import threading
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

plt.rcParams.update({
    "font.family":    "serif",
    "font.serif":     ["Liberation Serif", "DejaVu Serif"],
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

from paratan2.config.models import (
    SimpleMachineConfig, VacuumVesselConfig, VacuumVesselStructureLayer,
    CentralCellConfig, MaterialLayer,
    LFCoilConfig, LFCoilShellThicknesses, LFCoilInnerDimensions, LFCoilMaterials,
    HFCoilConfig, HFCoilMagnet, HFCoilShield,
    EndCellConfig, RoomConfig,
)

# ---------------------------------------------------------------------------
# Palette
# ---------------------------------------------------------------------------

BG      = "#1e1e1e"
BG2     = "#252526"
BG3     = "#2d2d2d"
ACCENT  = "#2d7dd2"
FG      = "#e0e0e0"
FG_DIM  = "#888"
BORDER  = "#444"
SUCCESS = "#4caf50"
ERROR   = "#f44336"
WARN    = "#ff9800"

FONT_BODY   = ("Liberation Serif", 9)
FONT_LABEL  = ("Liberation Serif", 8)
FONT_BOLD   = ("Liberation Serif", 10, "bold")
FONT_BTN    = ("Liberation Serif", 12, "bold")
FONT_TBAR   = ("Liberation Serif", 9)

# ---------------------------------------------------------------------------
# Dynamic material list from materials module
# ---------------------------------------------------------------------------

def _load_all_mats():
    import openmc
    from paratan2.materials import material as _m
    return sorted(name for name, obj in vars(_m).items() if isinstance(obj, openmc.Material))

ALL_MATS = _load_all_mats()

# ---------------------------------------------------------------------------
# Shared widget helpers
# ---------------------------------------------------------------------------

def make_spin(parent, value, lo=0.1, hi=9999.0, step=1.0):
    v = tk.DoubleVar(value=value)
    w = tk.Spinbox(parent, from_=lo, to=hi, textvariable=v, increment=step,
                   width=12, bg=BG3, fg=FG, insertbackground=FG,
                   buttonbackground=BG3, relief="flat",
                   highlightthickness=1, highlightcolor=ACCENT,
                   highlightbackground=BORDER, font=FONT_BODY)
    return w, v


def make_int_spin(parent, value, lo=1, hi=1_000_000, step=1):
    v = tk.IntVar(value=value)
    w = tk.Spinbox(parent, from_=lo, to=hi, textvariable=v, increment=step,
                   width=12, bg=BG3, fg=FG, insertbackground=FG,
                   buttonbackground=BG3, relief="flat",
                   highlightthickness=1, highlightcolor=ACCENT,
                   highlightbackground=BORDER, font=FONT_BODY)
    return w, v


def make_combo(parent, options, default=None):
    v = tk.StringVar(value=default or (options[0] if options else ""))
    w = ttk.Combobox(parent, textvariable=v, values=options,
                     style="Dark.TCombobox", width=26, state="readonly",
                     font=FONT_BODY)
    return w, v


def make_entry(parent, value="", width=28):
    v = tk.StringVar(value=value)
    w = tk.Entry(parent, textvariable=v,
                 bg=BG3, fg=FG, insertbackground=FG,
                 relief="flat", highlightthickness=1,
                 highlightbackground=BORDER, font=FONT_BODY, width=width)
    return w, v


def section_label(parent, text):
    tk.Frame(parent, bg=ACCENT, height=1).pack(fill="x", pady=(12, 0))
    tk.Label(parent, text=text, bg=BG2, fg=FG,
             font=FONT_BOLD).pack(anchor="w", padx=6, pady=(2, 4))


def param_row(parent, label_text, widget, unit=""):
    """Vertical: label above, widget below."""
    row = tk.Frame(parent, bg=BG2)
    row.pack(fill="x", padx=8, pady=(0, 5))
    lbl = f"{label_text}  ({unit})" if unit else label_text
    tk.Label(row, text=lbl, bg=BG2, fg=FG_DIM,
             font=FONT_LABEL, anchor="w").pack(anchor="w")
    widget.pack(fill="x", anchor="w")


def scrollable_frame(parent):
    """Returns (outer_frame, inner_frame, canvas) — pack outer_frame."""
    outer = tk.Frame(parent, bg=BG2)
    canvas = tk.Canvas(outer, bg=BG2, highlightthickness=0)
    sb = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=sb.set)
    sb.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)
    inner = tk.Frame(canvas, bg=BG2)
    _win = canvas.create_window((0, 0), window=inner, anchor="nw")

    def _resize(e):
        canvas.configure(scrollregion=canvas.bbox("all"))
        canvas.itemconfig(_win, width=e.width)
    canvas.bind("<Configure>", _resize)
    inner.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
    canvas.bind_all("<Button-4>", lambda e: canvas.yview_scroll(-1, "units"))
    canvas.bind_all("<Button-5>", lambda e: canvas.yview_scroll( 1, "units"))
    return outer, inner, canvas


# ---------------------------------------------------------------------------
# Main App
# ---------------------------------------------------------------------------

class ParatanApp:
    def __init__(self, root):
        self.root = root
        root.title("Paratan 2 — Simple Mirror Machine Builder")
        root.configure(bg=BG)
        root.geometry("1500x920")

        self._builder  = None
        self._geometry = None

        self._style_ttk()
        self._build_layout()

    # ------------------------------------------------------------------
    # TTK styles
    # ------------------------------------------------------------------

    def _style_ttk(self):
        s = ttk.Style()
        s.theme_use("clam")
        s.configure(".", background=BG, foreground=FG, fieldbackground=BG3,
                     troughcolor=BG3, bordercolor=BORDER, font=FONT_BODY)
        s.configure("TNotebook", background=BG2, bordercolor=BORDER, tabmargins=[0, 0, 0, 0])
        s.configure("TNotebook.Tab", background=BG3, foreground=FG_DIM,
                     padding=(14, 6), font=FONT_BODY)
        s.map("TNotebook.Tab",
              background=[("selected", ACCENT)],
              foreground=[("selected", "white")])
        s.configure("TScrollbar", background=BG3, troughcolor=BG2,
                     arrowcolor=FG_DIM, bordercolor=BG2)
        s.configure("Dark.TCombobox",
                    fieldbackground=BG3, background=BG3,
                    foreground=FG, arrowcolor=FG,
                    selectbackground=ACCENT, selectforeground=FG)

    # ------------------------------------------------------------------
    # Top-level layout
    # ------------------------------------------------------------------

    def _build_layout(self):
        pane = tk.PanedWindow(self.root, orient="horizontal",
                              bg=BG, sashwidth=4, sashrelief="flat")
        pane.pack(fill="both", expand=True, padx=6, pady=(6, 0))

        # ---- Left side: tabbed panels ----
        left_outer = tk.Frame(pane, bg=BG, width=380)
        pane.add(left_outer, minsize=340)

        # Tab notebook on the left
        self.left_nb = ttk.Notebook(left_outer)
        self.left_nb.pack(fill="both", expand=True)

        # Build button below the notebook
        btn_row = tk.Frame(left_outer, bg=BG, pady=6)
        btn_row.pack(fill="x", padx=4, side="bottom")
        self.build_btn = tk.Button(
            btn_row, text="  Build Geometry",
            command=self._on_build,
            bg=ACCENT, fg="white", activebackground="#3a8ee0",
            activeforeground="white", relief="flat",
            font=FONT_BTN, cursor="hand2", padx=10, pady=8)
        self.build_btn.pack(fill="x")

        # ---- Geometry tab ----
        geo_tab = tk.Frame(self.left_nb, bg=BG2)
        self.left_nb.add(geo_tab, text="  Geometry  ")
        geo_scroll_outer, self.form, _ = scrollable_frame(geo_tab)
        geo_scroll_outer.pack(fill="both", expand=True)
        self._build_geometry_form()

        # ---- Settings tab ----
        set_tab = tk.Frame(self.left_nb, bg=BG2)
        self.left_nb.add(set_tab, text="  Settings  ")
        set_scroll_outer, self.set_form, _ = scrollable_frame(set_tab)
        set_scroll_outer.pack(fill="both", expand=True)
        self._build_settings_form()

        # ---- Right side: toolbar + plot notebook ----
        self.right = tk.Frame(pane, bg=BG)
        pane.add(self.right, minsize=700)
        self._build_right_panel()

        # ---- Status bar ----
        self.status_var = tk.StringVar(value="Configure the machine and click Build Geometry.")
        self.status_bar = tk.Label(
            self.root, textvariable=self.status_var,
            bg=BG3, fg=FG_DIM, anchor="w", padx=8, pady=3,
            font=FONT_LABEL)
        self.status_bar.pack(fill="x", side="bottom")

    # ------------------------------------------------------------------
    # Geometry form
    # ------------------------------------------------------------------

    def _build_geometry_form(self):
        f = self.form

        # ── Vacuum Vessel ─────────────────────────────────────────────
        section_label(f, "Vacuum Vessel")
        self.vv_outer,  self.vv_outer_v  = make_spin(f, 450, 10, 5000)
        self.vv_cent,   self.vv_cent_v   = make_spin(f, 600, 10, 5000)
        self.vv_cr,     self.vv_cr_v     = make_spin(f, 75,  1,  500)
        self.vv_br,     self.vv_br_v     = make_spin(f, 25,  1,  500)
        self.vv_lbn,    self.vv_lbn_v    = make_spin(f, 300, 1,  2000)
        self.vv_rbn,    self.vv_rbn_v    = make_spin(f, 300, 1,  2000)
        for lbl, w, u in [
            ("Outer axial length",   self.vv_outer, "cm"),
            ("Central axial length", self.vv_cent,  "cm"),
            ("Central radius",       self.vv_cr,    "cm"),
            ("Bottleneck radius",    self.vv_br,    "cm"),
            ("Left BN length",       self.vv_lbn,   "cm"),
            ("Right BN length",      self.vv_rbn,   "cm"),
        ]:
            param_row(f, lbl, w, u)

        section_label(f, "  VV Structure")
        self.fw_t,   self.fw_t_v   = make_spin(f, 1.0, 0.1, 50)
        self.fw_mat, self.fw_mat_v = make_combo(f, ALL_MATS, "tungsten")
        self.vs_t,   self.vs_t_v   = make_spin(f, 2.0, 0.1, 50)
        self.vs_mat, self.vs_mat_v = make_combo(f, ALL_MATS, "FW_FNSF")
        param_row(f, "First wall thickness", self.fw_t,   "cm")
        param_row(f, "First wall material",  self.fw_mat, "")
        param_row(f, "VV shell thickness",   self.vs_t,   "cm")
        param_row(f, "VV shell material",    self.vs_mat, "")

        # ── Central Cell ──────────────────────────────────────────────
        section_label(f, "Central Cell")
        self.cc_axial, self.cc_axial_v = make_spin(f, 550, 10, 5000)
        param_row(f, "Axial length", self.cc_axial, "cm")

        layers = [
            ("Breeding blanket", "bl",  77,  2, "eutectic_breeding_material"),
            ("Back wall",        "bw",  2,   2, "BW_FNSF"),
            ("HT Shield",        "sh",  30,  2, "HT_Shield_filler"),
            ("Shield back plate","sbp", 2.5, 2, "water_cooled_wc"),
        ]
        for title, attr, t_def, t_lo, mat_def in layers:
            section_label(f, f"  {title}")
            t_w, t_v = make_spin(f, t_def, t_lo, 500)
            m_w, m_v = make_combo(f, ALL_MATS, mat_def)
            setattr(self, f"{attr}_t",     t_w);  setattr(self, f"{attr}_t_v",   t_v)
            setattr(self, f"{attr}_mat",   m_w);  setattr(self, f"{attr}_mat_v", m_v)
            param_row(f, "Thickness", t_w, "cm")
            param_row(f, "Material",  m_w, "")

        # ── LF Coils ──────────────────────────────────────────────────
        section_label(f, "LF (Solenoid) Coils")
        self.lf_pos_entry, self.lf_pos_v = make_entry(f, "-75, 75")
        param_row(f, "Positions (comma-separated)", self.lf_pos_entry, "cm")

        self.lf_al,  self.lf_al_v  = make_spin(f, 20, 1, 200)
        self.lf_rt,  self.lf_rt_v  = make_spin(f, 10, 1, 200)
        self.lf_fs,  self.lf_fs_v  = make_spin(f, 10, 0.1, 100)
        self.lf_bs,  self.lf_bs_v  = make_spin(f, 5,  0.1, 100)
        self.lf_as,  self.lf_as_v  = make_spin(f, 10, 0.1, 100)
        self.lf_sm,  self.lf_sm_v  = make_combo(f, ALL_MATS, "stainless")
        self.lf_mm,  self.lf_mm_v  = make_combo(f, ALL_MATS, "Magnet_Winding_Pack_2")
        for lbl, w, u in [
            ("Winding axial length",     self.lf_al, "cm"),
            ("Winding radial thickness", self.lf_rt, "cm"),
            ("Shell front thickness",    self.lf_fs, "cm"),
            ("Shell back thickness",     self.lf_bs, "cm"),
            ("Shell axial thickness",    self.lf_as, "cm"),
            ("Shield material",          self.lf_sm, ""),
            ("Magnet material",          self.lf_mm, ""),
        ]:
            param_row(f, lbl, w, u)

        # ── HF Coils ──────────────────────────────────────────────────
        section_label(f, "HF (Mirror) Coils — Magnet")
        self.hf_br,  self.hf_br_v  = make_spin(f, 60,  1,  500)
        self.hf_rt,  self.hf_rt_v  = make_spin(f, 150, 10, 1000)
        self.hf_at,  self.hf_at_v  = make_spin(f, 200, 10, 1000)
        self.hf_mm,  self.hf_mm_v  = make_combo(f, ALL_MATS, "Magnet_Winding_Pack_2")
        for lbl, w, u in [
            ("Bore radius",      self.hf_br, "cm"),
            ("Radial thickness", self.hf_rt, "cm"),
            ("Axial thickness",  self.hf_at, "cm"),
            ("Magnet material",  self.hf_mm, ""),
        ]:
            param_row(f, lbl, w, u)

        section_label(f, "  HF Casing layers (inner → outer)")
        # Dynamic list — each entry is (thickness_DoubleVar, material_StringVar, row_Frame)
        self._hf_casing_rows: list[tuple[tk.DoubleVar, tk.StringVar, tk.Frame]] = []
        self._hf_casing_container = tk.Frame(f, bg=BG2)
        self._hf_casing_container.pack(fill="x")
        # "Add layer" button
        tk.Button(f, text="+ Add casing layer",
                  command=self._hf_casing_add,
                  bg=BG3, fg=FG, activebackground=BG2, activeforeground=FG,
                  relief="flat", font=FONT_LABEL, cursor="hand2",
                  padx=6, pady=3).pack(anchor="w", padx=8, pady=(2, 6))
        # Seed with the default two layers
        self._hf_casing_add(thickness=1.0, material="vacuum")
        self._hf_casing_add(thickness=1.0, material="ss316ln")

        section_label(f, "  HF Shield")
        self.hf_sg,  self.hf_sg_v  = make_spin(f, 5,  0, 200)
        self.hf_rg,  self.hf_rg_v  = make_spin(f, 1,  0, 50)
        self.hf_sri, self.hf_sri_v = make_spin(f, 20, 1, 300)
        self.hf_sro, self.hf_sro_v = make_spin(f, 30, 1, 300)
        self.hf_sai, self.hf_sai_v = make_spin(f, 20, 1, 300)
        self.hf_sao, self.hf_sao_v = make_spin(f, 10, 1, 300)
        self.hf_sm,  self.hf_sm_v  = make_combo(f, ALL_MATS, "cooled_w2b5")
        for lbl, w, u in [
            ("Shield–CC axial gap",      self.hf_sg,  "cm"),
            ("Radial gap before casing", self.hf_rg,  "cm"),
            ("Shield radial inner",      self.hf_sri, "cm"),
            ("Shield radial outer",      self.hf_sro, "cm"),
            ("Shield axial (midplane)",  self.hf_sai, "cm"),
            ("Shield axial (away)",      self.hf_sao, "cm"),
            ("Shield material",          self.hf_sm,  ""),
        ]:
            param_row(f, lbl, w, u)

        # ── End Cells ─────────────────────────────────────────────────
        section_label(f, "End Cells")
        self.ec_ax, self.ec_ax_v = make_spin(f, 500, 10, 2000)
        self.ec_st, self.ec_st_v = make_spin(f, 3,  0.1, 50)
        self.ec_d,  self.ec_d_v  = make_spin(f, 350, 10, 1000)
        self.ec_sm, self.ec_sm_v = make_combo(f, ALL_MATS, "stainless")
        for lbl, w, u in [
            ("Axial length",    self.ec_ax, "cm"),
            ("Shell thickness", self.ec_st, "cm"),
            ("Diameter",        self.ec_d,  "cm"),
            ("Shell material",  self.ec_sm, ""),
        ]:
            param_row(f, lbl, w, u)

        # bottom padding
        tk.Frame(f, bg=BG2, height=12).pack()

    # ------------------------------------------------------------------
    # Settings form
    # ------------------------------------------------------------------

    def _build_settings_form(self):
        f = self.set_form

        # ── Simulation ────────────────────────────────────────────────
        section_label(f, "Simulation")
        self.sim_particles, self.sim_particles_v = make_int_spin(f, 1000, 100, 10_000_000, 100)
        self.sim_batches,   self.sim_batches_v   = make_int_spin(f, 10,   1,   10_000,     1)
        self.sim_sp_freq,   self.sim_sp_freq_v   = make_int_spin(f, 5,    1,   10_000,     1)
        param_row(f, "Particles per batch",    self.sim_particles, "")
        param_row(f, "Number of batches",      self.sim_batches,   "")
        param_row(f, "Statepoint frequency",   self.sim_sp_freq,   "batches")

        # Photon transport toggle
        tk.Frame(f, bg=BG2, height=6).pack()
        self.photon_var = tk.BooleanVar(value=False)
        photon_row = tk.Frame(f, bg=BG2)
        photon_row.pack(fill="x", padx=8, pady=(0, 4))
        tk.Checkbutton(
            photon_row, text="  Photon transport",
            variable=self.photon_var,
            bg=BG2, fg=FG, selectcolor=BG3,
            activebackground=BG2, activeforeground=FG,
            font=FONT_BODY
        ).pack(anchor="w")
        tk.Label(photon_row, text="Enables coupled neutron–photon heating tallies",
                 bg=BG2, fg=FG_DIM, font=FONT_LABEL).pack(anchor="w", padx=20)

        # Weight windows toggle
        tk.Frame(f, bg=BG2, height=2).pack()
        self.ww_var = tk.BooleanVar(value=False)
        ww_row = tk.Frame(f, bg=BG2)
        ww_row.pack(fill="x", padx=8, pady=(0, 4))
        tk.Checkbutton(
            ww_row, text="  Weight windows (MAGIC)",
            variable=self.ww_var,
            bg=BG2, fg=FG, selectcolor=BG3,
            activebackground=BG2, activeforeground=FG,
            font=FONT_BODY
        ).pack(anchor="w")
        tk.Label(ww_row, text="Variance reduction for HF coil region",
                 bg=BG2, fg=FG_DIM, font=FONT_LABEL).pack(anchor="w", padx=20)

        # ── Neutron Source ────────────────────────────────────────────
        section_label(f, "Neutron Source")

        # Source type selector
        self.src_type_v = tk.StringVar(value="Volumetric (plasma)")
        src_types = ["Volumetric (plasma)", "Uniform cylinder", "From YAML file"]
        src_combo, _ = make_combo(f, src_types, "Volumetric (plasma)")
        src_combo.configure(textvariable=self.src_type_v)
        param_row(f, "Source type", src_combo, "")

        # Volumetric params
        self.src_frame_vol = tk.Frame(f, bg=BG2)
        self.src_frame_vol.pack(fill="x")

        self.src_r_v   = tk.DoubleVar(value=50.0)
        self.src_zhl_v = tk.DoubleVar(value=275.0)
        r_spin  = tk.Spinbox(self.src_frame_vol, from_=1, to=500, textvariable=self.src_r_v,
                             increment=1, width=12, bg=BG3, fg=FG, insertbackground=FG,
                             buttonbackground=BG3, relief="flat",
                             highlightthickness=1, highlightcolor=ACCENT,
                             highlightbackground=BORDER, font=FONT_BODY)
        z_spin  = tk.Spinbox(self.src_frame_vol, from_=1, to=2000, textvariable=self.src_zhl_v,
                             increment=5, width=12, bg=BG3, fg=FG, insertbackground=FG,
                             buttonbackground=BG3, relief="flat",
                             highlightthickness=1, highlightcolor=ACCENT,
                             highlightbackground=BORDER, font=FONT_BODY)
        param_row(self.src_frame_vol, "Plasma radius",       r_spin,  "cm")
        param_row(self.src_frame_vol, "Plasma half-length",  z_spin,  "cm")

        # YAML source path (hidden until "From YAML file" selected)
        self.src_frame_yaml = tk.Frame(f, bg=BG2)
        yaml_entry_row = tk.Frame(self.src_frame_yaml, bg=BG2)
        yaml_entry_row.pack(fill="x", padx=8, pady=(0, 4))
        tk.Label(yaml_entry_row, text="Source YAML path", bg=BG2, fg=FG_DIM,
                 font=FONT_LABEL, anchor="w").pack(anchor="w")
        path_row = tk.Frame(yaml_entry_row, bg=BG2)
        path_row.pack(fill="x")
        self.src_yaml_v = tk.StringVar(value="")
        yaml_e = tk.Entry(path_row, textvariable=self.src_yaml_v,
                          bg=BG3, fg=FG, insertbackground=FG,
                          relief="flat", highlightthickness=1,
                          highlightbackground=BORDER, font=FONT_BODY)
        yaml_e.pack(side="left", fill="x", expand=True)
        tk.Button(path_row, text="Browse", command=self._browse_yaml,
                  bg=BG3, fg=FG, relief="flat", font=FONT_LABEL,
                  cursor="hand2", padx=6, pady=2).pack(side="left", padx=(4, 0))

        def _update_src_frames(*_):
            sel = self.src_type_v.get()
            if sel == "From YAML file":
                self.src_frame_vol.pack_forget()
                self.src_frame_yaml.pack(fill="x")
            else:
                self.src_frame_yaml.pack_forget()
                self.src_frame_vol.pack(fill="x")
        self.src_type_v.trace_add("write", _update_src_frames)

        # ── Room (bounding box) ───────────────────────────────────────
        section_label(f, "Room (bounding box)")
        self.room_xy, self.room_xy_v = make_spin(f, 400,  50, 5000)
        self.room_z,  self.room_z_v  = make_spin(f, 1350, 50, 10000)
        param_row(f, "XY half-extent", self.room_xy, "cm")
        param_row(f, "Z half-extent",  self.room_z,  "cm")

        tk.Frame(f, bg=BG2, height=12).pack()

    def _browse_yaml(self):
        path = filedialog.askopenfilename(
            title="Select source YAML",
            filetypes=[("YAML files", "*.yaml *.yml"), ("All files", "*.*")],
        )
        if path:
            self.src_yaml_v.set(path)

    # ------------------------------------------------------------------
    # Right panel — toolbar + plot notebook
    # ------------------------------------------------------------------

    def _hf_casing_add(self, thickness=1.0, material="vacuum"):
        c = self._hf_casing_container
        t_v = tk.DoubleVar(value=thickness)
        m_v = tk.StringVar(value=material)

        row = tk.Frame(c, bg=BG2, pady=2)
        row.pack(fill="x", padx=8)

        idx_label = tk.Label(row, bg=BG2, fg=FG_DIM, font=FONT_LABEL, width=3, anchor="e")
        idx_label.pack(side="left")

        t_spin = tk.Spinbox(row, from_=0.1, to=200, textvariable=t_v, increment=0.5,
                            width=7, bg=BG3, fg=FG, insertbackground=FG,
                            buttonbackground=BG3, relief="flat",
                            highlightthickness=1, highlightcolor=ACCENT,
                            highlightbackground=BORDER, font=FONT_BODY)
        t_spin.pack(side="left", padx=(2, 4))

        tk.Label(row, text="cm", bg=BG2, fg=FG_DIM, font=FONT_LABEL).pack(side="left")

        m_combo = ttk.Combobox(row, textvariable=m_v, values=ALL_MATS,
                               style="Dark.TCombobox", width=20, state="readonly",
                               font=FONT_BODY)
        m_combo.pack(side="left", padx=(6, 4))

        entry = (t_v, m_v, row)
        self._hf_casing_rows.append(entry)

        def _delete(e=entry):
            self._hf_casing_delete(e)

        tk.Button(row, text="×", command=_delete,
                  bg=BG2, fg=ERROR, activebackground=BG3, activeforeground=ERROR,
                  relief="flat", font=("Liberation Serif", 11, "bold"),
                  cursor="hand2", padx=4, pady=0).pack(side="left", padx=2)

        self._hf_casing_renumber()

    def _hf_casing_delete(self, entry):
        if len(self._hf_casing_rows) <= 1:
            return   # keep at least one layer
        t_v, m_v, row = entry
        self._hf_casing_rows.remove(entry)
        row.destroy()
        self._hf_casing_renumber()

    def _hf_casing_renumber(self):
        for i, (t_v, m_v, row) in enumerate(self._hf_casing_rows):
            # First child of row is the index label
            row.winfo_children()[0].config(text=f"{i+1}.")

    def _build_right_panel(self):
        toolbar = tk.Frame(self.right, bg=BG3, pady=4)
        toolbar.pack(fill="x")

        def _tbtn(text, cmd, state="disabled"):
            b = tk.Button(toolbar, text=text, command=cmd,
                          bg=BG3, fg=FG, activebackground=BG2,
                          activeforeground=FG, relief="flat",
                          font=FONT_TBAR, cursor="hand2",
                          padx=10, pady=4, disabledforeground=BORDER)
            b.pack(side="left", padx=4)
            b.config(state=state)
            return b

        self.btn_save_img   = _tbtn("Save Images…",      self._on_save_image)
        self.btn_export_xml = _tbtn("Export Model XML…", self._on_export_xml)
        self.btn_check_geom = _tbtn("Check Geometry",    self._on_check_geometry)

        self.plot_nb = ttk.Notebook(self.right)
        self.plot_nb.pack(fill="both", expand=True)

        # ── Legend panel below the notebook ───────────────────────────
        legend_outer = tk.Frame(self.right, bg=BG3)
        legend_outer.pack(fill="x", side="bottom")
        tk.Frame(legend_outer, bg=BORDER, height=1).pack(fill="x")
        hdr = tk.Frame(legend_outer, bg=BG3, pady=4)
        hdr.pack(fill="x", padx=8)
        tk.Label(hdr, text="Materials", bg=BG3, fg=FG,
                 font=FONT_BOLD).pack(side="left")
        # Scrollable inner row for swatches
        leg_canvas = tk.Canvas(legend_outer, bg=BG3, height=54, highlightthickness=0)
        leg_hbar = ttk.Scrollbar(legend_outer, orient="horizontal",
                                  command=leg_canvas.xview)
        leg_canvas.configure(xscrollcommand=leg_hbar.set)
        leg_hbar.pack(fill="x", side="bottom")
        leg_canvas.pack(fill="x")
        self._legend_inner = tk.Frame(leg_canvas, bg=BG3)
        _lw = leg_canvas.create_window((0, 0), window=self._legend_inner, anchor="nw")
        self._legend_inner.bind(
            "<Configure>",
            lambda e: leg_canvas.configure(scrollregion=leg_canvas.bbox("all"))
        )
        # placeholder text
        self._legend_placeholder = tk.Label(
            self._legend_inner, text="Build geometry to populate the legend.",
            bg=BG3, fg=FG_DIM, font=FONT_LABEL)
        self._legend_placeholder.pack(padx=8, pady=6)

        # ── XZ tab: fixed-size figure in a scrollable canvas ──────────
        # width=(1000,2800) → portrait ratio ≈ 1:2.8
        # figsize=(8, 22.4) at 100 dpi → 800 × 2240 px on screen
        tab_xz = tk.Frame(self.plot_nb, bg=BG)
        self.plot_nb.add(tab_xz, text="  Side view (XZ)  ")

        xz_scroll = tk.Canvas(tab_xz, bg=BG, highlightthickness=0)
        xz_vbar = ttk.Scrollbar(tab_xz, orient="vertical",   command=xz_scroll.yview)
        xz_hbar = ttk.Scrollbar(tab_xz, orient="horizontal", command=xz_scroll.xview)
        xz_scroll.configure(yscrollcommand=xz_vbar.set, xscrollcommand=xz_hbar.set)
        xz_vbar.pack(side="right",  fill="y")
        xz_hbar.pack(side="bottom", fill="x")
        xz_scroll.pack(fill="both", expand=True)

        self.fig_xz, self.ax_xz = plt.subplots(figsize=(8, 22.4), facecolor=BG, dpi=100)
        self.ax_xz.set_facecolor(BG)
        self.ax_xz.set_title("Build geometry to see the machine", color=FG_DIM)
        self.canvas_xz = FigureCanvasTkAgg(self.fig_xz, master=xz_scroll)
        xz_widget = self.canvas_xz.get_tk_widget()
        xz_scroll.create_window((0, 0), window=xz_widget, anchor="nw")
        xz_widget.bind("<Configure>",
                       lambda e: xz_scroll.configure(scrollregion=xz_scroll.bbox("all")))
        # Mouse-wheel scrolling on the XZ tab
        xz_scroll.bind("<Button-4>", lambda e: xz_scroll.yview_scroll(-1, "units"))
        xz_scroll.bind("<Button-5>", lambda e: xz_scroll.yview_scroll( 1, "units"))

        # ── XY tab: fills the panel (square, no scroll needed) ────────
        tab_xy = tk.Frame(self.plot_nb, bg=BG)
        self.plot_nb.add(tab_xy, text="  End-on view (XY)  ")
        self.fig_xy, self.ax_xy = plt.subplots(figsize=(7, 7), facecolor=BG)
        self.ax_xy.set_facecolor(BG)
        self.ax_xy.set_title("Build geometry to see the machine", color=FG_DIM)
        self.canvas_xy = FigureCanvasTkAgg(self.fig_xy, master=tab_xy)
        self.canvas_xy.get_tk_widget().pack(fill="both", expand=True)

    # ------------------------------------------------------------------
    # Config assembly
    # ------------------------------------------------------------------

    def _get_config(self):
        lf_positions = [float(x.strip()) for x in self.lf_pos_v.get().split(",") if x.strip()]
        r = self.room_xy_v.get()
        z = self.room_z_v.get()
        return SimpleMachineConfig(
            vacuum_vessel=VacuumVesselConfig(
                outer_axial_length=self.vv_outer_v.get(),
                central_axial_length=self.vv_cent_v.get(),
                central_radius=self.vv_cr_v.get(),
                bottleneck_radius=self.vv_br_v.get(),
                left_bottleneck_length=self.vv_lbn_v.get(),
                right_bottleneck_length=self.vv_rbn_v.get(),
                structure={
                    "first_wall": VacuumVesselStructureLayer(thickness=self.fw_t_v.get(), material=self.fw_mat_v.get()),
                    "vv_shell":   VacuumVesselStructureLayer(thickness=self.vs_t_v.get(), material=self.vs_mat_v.get()),
                },
            ),
            central_cell=CentralCellConfig(
                axial_length=self.cc_axial_v.get(),
                layers=[
                    MaterialLayer(thickness=self.bl_t_v.get(),  material=self.bl_mat_v.get()),
                    MaterialLayer(thickness=self.bw_t_v.get(),  material=self.bw_mat_v.get()),
                    MaterialLayer(thickness=self.sh_t_v.get(),  material=self.sh_mat_v.get()),
                    MaterialLayer(thickness=self.sbp_t_v.get(), material=self.sbp_mat_v.get()),
                ],
            ),
            lf_coil=LFCoilConfig(
                shell_thicknesses=LFCoilShellThicknesses(front=self.lf_fs_v.get(), back=self.lf_bs_v.get(), axial=self.lf_as_v.get()),
                inner_dimensions=LFCoilInnerDimensions(axial_length=self.lf_al_v.get(), radial_thickness=self.lf_rt_v.get()),
                positions=lf_positions,
                materials=LFCoilMaterials(shield=self.lf_sm_v.get(), magnet=self.lf_mm_v.get()),
            ),
            hf_coil=HFCoilConfig(
                magnet=HFCoilMagnet(
                    bore_radius=self.hf_br_v.get(), radial_thickness=self.hf_rt_v.get(),
                    axial_thickness=self.hf_at_v.get(), material=self.hf_mm_v.get(),
                ),
                casing_layers=[
                    MaterialLayer(thickness=t_v.get(), material=m_v.get())
                    for t_v, m_v, _ in self._hf_casing_rows
                ],
                shield=HFCoilShield(
                    radial_gap_before_casing=self.hf_rg_v.get(),
                    shield_central_cell_gap=self.hf_sg_v.get(),
                    radial_thickness=[self.hf_sri_v.get(), self.hf_sro_v.get()],
                    axial_thickness=[self.hf_sai_v.get(), self.hf_sao_v.get()],
                    material=self.hf_sm_v.get(),
                ),
            ),
            end_cell=EndCellConfig(
                axial_length=self.ec_ax_v.get(), shell_thickness=self.ec_st_v.get(),
                diameter=self.ec_d_v.get(), shell_material=self.ec_sm_v.get(),
                inner_material="vacuum",
            ),
            room=RoomConfig(
                xmin=-r, xmax=r, ymin=-r, ymax=r,
                zmin=-z, zmax=z,
            ),
        )

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------

    def _set_status(self, text, color=FG_DIM):
        self.status_var.set(text)
        self.status_bar.config(fg=color)

    def _on_build(self):
        try:
            config = self._get_config()
        except Exception as e:
            self._set_status(f"Config error: {e}", ERROR)
            return

        self.build_btn.config(state="disabled", text="  Building…")
        self._set_status("Building geometry…", FG_DIM)

        def _worker():
            try:
                from paratan2.models.simple_machine_builder import SimpleMachineBuilder, _load_materials
                m  = _load_materials()
                builder = SimpleMachineBuilder(config, m)
                orig = os.getcwd()
                with tempfile.TemporaryDirectory() as tmp:
                    os.chdir(tmp)
                    try:
                        builder.build()
                    finally:
                        os.chdir(orig)
                geom = builder.get_geometry()
                self.root.after(0, lambda: self._on_build_done(geom, builder))
            except Exception as e:
                import traceback
                tb = traceback.format_exc()
                self.root.after(0, lambda: self._on_build_error(tb))

        threading.Thread(target=_worker, daemon=True).start()

    def _on_build_error(self, msg):
        self._set_status("Build failed", ERROR)
        self.build_btn.config(state="normal", text="  Build Geometry")
        messagebox.showerror("Build failed", msg)

    def _on_build_done(self, geom, builder):
        n = len(builder.all_cells)
        self._set_status(f"Built successfully — {n} cells", SUCCESS)
        self.build_btn.config(state="normal", text="  Build Geometry")
        self._builder  = builder
        self._geometry = geom
        self.btn_save_img.config(state="normal")
        self.btn_export_xml.config(state="normal")
        self.btn_check_geom.config(state="normal")
        self._render_plots(geom, builder)

    # ------------------------------------------------------------------
    # Plot rendering — use full room extent so nothing is cut off
    # ------------------------------------------------------------------

    def _plot_extents(self, builder):
        """Returns (w_x, w_z) — the full physical width in each direction."""
        rc = builder.config.room
        w_x = rc.xmax - rc.xmin   # radial extent  (e.g. 800 cm)
        w_z = rc.zmax - rc.zmin   # axial extent   (e.g. 2700 cm)
        return w_x, w_z

    @staticmethod
    def _build_color_map(builder):
        """Assign a distinct color to every unique material in the built cells.
        Returns {openmc.Material: (R, G, B)} with integer values in [0, 255]
        as required by openmc's colors= parameter."""
        import openmc
        import matplotlib

        seen = {}
        for cell in builder.all_cells:
            mat = cell.fill
            if isinstance(mat, openmc.Material) and mat not in seen:
                seen[mat] = None

        mats = list(seen.keys())
        palette = matplotlib.colormaps.get_cmap("tab20")
        for i, mat in enumerate(mats):
            rgba = palette(i % 20)   # floats [0,1]
            # OpenMC expects integer RGB in [0, 255]
            seen[mat] = tuple(int(v * 255) for v in rgba[:3])

        return seen   # {Material: (r, g, b)  integers 0-255}

    def _populate_legend(self, color_map):
        """Render the color_map as tkinter swatches in the legend strip."""
        for w in self._legend_inner.winfo_children():
            w.destroy()

        def _to_hex(rgb):
            # rgb values are already integers 0-255
            return "#{:02x}{:02x}{:02x}".format(
                max(0, min(255, rgb[0])),
                max(0, min(255, rgb[1])),
                max(0, min(255, rgb[2])),
            )

        for mat, rgb in color_map.items():
            hex_col = _to_hex(rgb)
            item = tk.Frame(self._legend_inner, bg=BG3)
            item.pack(side="left", padx=(6, 2), pady=4)
            tk.Label(item, bg=hex_col, width=2, height=1,
                     relief="flat").pack(side="left")
            tk.Label(item, text=mat.name, bg=BG3, fg=FG,
                     font=FONT_LABEL, anchor="w").pack(side="left", padx=(3, 0))

    def _render_plots(self, geom, builder):
        w_x, w_z = self._plot_extents(builder)

        # Build a consistent color map for all materials so both plots and
        # the legend all use the same colours.
        color_map = self._build_color_map(builder)
        # openmc expects {Material: (r,g,b)} — same format we built
        omc_colors = color_map

        # ── XZ side view — same dimensions as paratan v1 ──────────────
        self.ax_xz.clear()
        geom.root_universe.plot(
            basis="xz", origin=(0, 0, 0),
            width=(1000, 2800),
            pixels=(800, 2200),
            color_by="material",
            colors=omc_colors,
            axes=self.ax_xz,
        )
        self.ax_xz.set_title("Side view (XZ)", color=FG)
        self.ax_xz.set_xlabel("X (cm)", color=FG)
        self.ax_xz.set_ylabel("Z (cm)", color=FG)
        self.ax_xz.tick_params(colors=FG)
        for sp in self.ax_xz.spines.values():
            sp.set_edgecolor(BORDER)
        self.fig_xz.set_facecolor(BG)
        self.ax_xz.set_facecolor(BG)
        self.fig_xz.tight_layout()
        self.canvas_xz.draw()

        # ── XY end-on view ────────────────────────────────────────────
        self.ax_xy.clear()
        geom.root_universe.plot(
            basis="xy", origin=(0, 0, 0),
            width=(w_x, w_x),
            pixels=(800, 800),
            color_by="material",
            colors=omc_colors,
            axes=self.ax_xy,
        )
        self.ax_xy.set_aspect("auto")
        self.ax_xy.set_title("End-on view (XY) at midplane", color=FG)
        self.ax_xy.set_xlabel("X (cm)", color=FG)
        self.ax_xy.set_ylabel("Y (cm)", color=FG)
        self.ax_xy.tick_params(colors=FG)
        for sp in self.ax_xy.spines.values():
            sp.set_edgecolor(BORDER)
        self.fig_xy.set_facecolor(BG)
        self.ax_xy.set_facecolor(BG)
        self.fig_xy.tight_layout()
        self.canvas_xy.draw()

        # ── Legend strip ──────────────────────────────────────────────
        self._populate_legend(color_map)

        self.plot_nb.select(0)

    # ------------------------------------------------------------------
    # Check geometry — runs openmc -g in a temp dir, parses overlap output
    # ------------------------------------------------------------------

    def _on_check_geometry(self):
        self.btn_check_geom.config(state="disabled", text="Checking…")
        self._set_status("Running geometry debug check (openmc -g, 5000 particles)…", FG_DIM)

        geom = self._geometry

        def _worker():
            import subprocess
            import openmc

            try:
                with tempfile.TemporaryDirectory() as tmpdir:
                    orig = os.getcwd()
                    os.chdir(tmpdir)
                    try:
                        # Use a single dummy air material for all cells —
                        # material accuracy doesn't matter for overlap checking.
                        air = openmc.Material(name="air")
                        air.add_nuclide("N14", 0.78, "ao")
                        air.add_nuclide("O16", 0.22, "ao")
                        air.set_density("g/cm3", 0.001205)

                        all_cells   = geom.root_universe.get_all_cells().values()
                        orig_fills  = {}
                        for cell in all_cells:
                            if isinstance(cell.fill, openmc.Material):
                                orig_fills[cell] = cell.fill
                                cell.fill = air

                        try:
                            geom.export_to_xml("geometry.xml")
                            openmc.Materials([air]).export_to_xml("materials.xml")
                        finally:
                            for cell, mat in orig_fills.items():
                                cell.fill = mat

                        settings = openmc.Settings()
                        settings.run_mode  = "fixed source"
                        settings.particles = 5000
                        settings.batches   = 1
                        settings.source    = openmc.IndependentSource(
                            space=openmc.stats.Point((0, 0, 0)),
                            energy=openmc.stats.Discrete([14.1e6], [1.0]),
                            particle="neutron",
                        )
                        settings.export_to_xml("settings.xml")

                        proc = subprocess.run(
                            ["openmc", "-g"],
                            capture_output=True, text=True, timeout=300,
                        )
                        output = proc.stdout + proc.stderr
                    finally:
                        os.chdir(orig)

                self.root.after(0, lambda: self._on_check_done(output))

            except FileNotFoundError:
                self.root.after(0, lambda: self._on_check_error(
                    "openmc executable not found.\n\nMake sure OpenMC is on your PATH."))
            except Exception:
                import traceback
                tb = traceback.format_exc()
                self.root.after(0, lambda: self._on_check_error(tb))

        threading.Thread(target=_worker, daemon=True).start()

    def _on_check_done(self, output):
        self.btn_check_geom.config(state="normal", text="Check Geometry")

        lines = output.splitlines()

        # Real overlap detection — OpenMC prints this if two cells share a point
        overlap_errors = [l for l in lines if "overlapping cells detected" in l.lower()]

        # Coverage line — "There were N cells with less than 10 overlap checks"
        coverage_line = next((l.strip() for l in lines
                              if "cells with less than 10 overlap checks" in l), None)
        try:
            under_checked = int(coverage_line.split()[2]) if coverage_line else 0
        except (IndexError, ValueError):
            under_checked = 0

        if overlap_errors:
            status = f"Overlaps detected — {len(overlap_errors)} issue(s)"
            self._set_status(status, ERROR)
            title = "Geometry Check — OVERLAPS FOUND"
            ok    = False
        elif under_checked > 0:
            status = f"Geometry OK but {under_checked} cell(s) had < 10 checks — try more particles"
            self._set_status(status, WARN)
            title = "Geometry Check — Low Coverage Warning"
            ok    = True   # no overlaps, just a coverage advisory
        else:
            self._set_status("Geometry clean — 0 overlapping cells, all cells well-checked", SUCCESS)
            title = "Geometry Check Passed"
            ok    = True

        self._show_output_dialog(title=title, output=output, ok=ok)

    def _show_output_dialog(self, title, output, ok):
        """Scrollable text window showing full openmc -g output."""
        win = tk.Toplevel(self.root)
        win.title(title)
        win.configure(bg=BG)
        win.geometry("900x600")

        hdr_color = SUCCESS if ok else ERROR
        hdr_text  = "No overlapping cells detected." if ok else "Overlapping cells detected — see output below."
        tk.Label(win, text=hdr_text, bg=BG, fg=hdr_color,
                 font=FONT_BOLD, pady=6).pack(fill="x", padx=12)

        frame = tk.Frame(win, bg=BG)
        frame.pack(fill="both", expand=True, padx=12, pady=(0, 12))

        txt = tk.Text(frame, bg=BG3, fg=FG, font=("Courier New", 9),
                      wrap="none", relief="flat", bd=0)
        vsb = ttk.Scrollbar(frame, orient="vertical",   command=txt.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=txt.xview)
        txt.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        vsb.pack(side="right",  fill="y")
        hsb.pack(side="bottom", fill="x")
        txt.pack(side="left",   fill="both", expand=True)

        txt.insert("1.0", output)
        for i, line in enumerate(output.splitlines(), start=1):
            ll = line.lower()
            if "overlapping cells detected" in ll:
                txt.tag_add("overlap_err", f"{i}.0", f"{i}.end")
            elif "overlap check summary" in ll or "overlap checks" in ll:
                txt.tag_add("summary", f"{i}.0", f"{i}.end")
            elif "there were 0 cells" in ll:
                txt.tag_add("good", f"{i}.0", f"{i}.end")
        txt.tag_config("overlap_err", foreground=ERROR,   font=("Courier New", 9, "bold"))
        txt.tag_config("summary",     foreground=ACCENT)
        txt.tag_config("good",        foreground=SUCCESS, font=("Courier New", 9, "bold"))
        txt.config(state="disabled")

        tk.Button(win, text="Close", command=win.destroy,
                  bg=ACCENT, fg=FG, font=FONT_BTN,
                  relief="flat", padx=16, pady=4).pack(pady=(0, 12))

    def _on_check_error(self, msg):
        self.btn_check_geom.config(state="normal", text="Check Geometry")
        self._set_status("Geometry check failed", ERROR)
        messagebox.showerror("Geometry Check Failed", msg)

    # ------------------------------------------------------------------
    # Save images — two files, correct aspect ratios
    # ------------------------------------------------------------------

    def _on_save_image(self):
        path = filedialog.asksaveasfilename(
            title="Save images — base name (_xz and _xy will be appended)",
            defaultextension=".png",
            filetypes=[("PNG image", "*.png"), ("PDF", "*.pdf"), ("All files", "*.*")],
            initialfile="paratan2_geometry",
        )
        if not path:
            return

        base, ext = os.path.splitext(path)
        ext = ext or ".png"

        w_x, w_z = self._plot_extents(self._builder)
        omc_colors = self._build_color_map(self._builder)

        # ── XZ: same dimensions as paratan v1 ────────────────────────
        fig_xz = plt.figure(figsize=(10, 28), facecolor=BG)
        ax_xz  = fig_xz.add_subplot(111, facecolor=BG)
        self._geometry.root_universe.plot(
            basis="xz", origin=(0, 0, 0),
            width=(1000, 2800), pixels=(2400, 4000),
            color_by="material", colors=omc_colors, axes=ax_xz,
        )
        ax_xz.set_aspect("auto")
        ax_xz.set_title("Side view (XZ)", color=FG, fontsize=13)
        ax_xz.set_xlabel("X (cm)", color=FG)
        ax_xz.set_ylabel("Z (cm)", color=FG)
        ax_xz.tick_params(colors=FG)
        for sp in ax_xz.spines.values():
            sp.set_edgecolor(BORDER)
        path_xz = base + "_xz" + ext
        fig_xz.savefig(path_xz, dpi=150, bbox_inches="tight", facecolor=BG)
        plt.close(fig_xz)

        # ── XY: square ────────────────────────────────────────────────
        fig_xy = plt.figure(figsize=(9, 9), facecolor=BG)
        ax_xy  = fig_xy.add_subplot(111, facecolor=BG)
        self._geometry.root_universe.plot(
            basis="xy", origin=(0, 0, 0),
            width=(w_x, w_x), pixels=(1400, 1400),
            color_by="material", colors=omc_colors, axes=ax_xy,
        )
        ax_xy.set_aspect("auto")
        ax_xy.set_title("End-on view (XY) at midplane", color=FG, fontsize=13)
        ax_xy.set_xlabel("X (cm)", color=FG)
        ax_xy.set_ylabel("Y (cm)", color=FG)
        ax_xy.tick_params(colors=FG)
        for sp in ax_xy.spines.values():
            sp.set_edgecolor(BORDER)
        path_xy = base + "_xy" + ext
        fig_xy.savefig(path_xy, dpi=150, bbox_inches="tight", facecolor=BG)
        plt.close(fig_xy)

        self._set_status(f"Saved: …{os.path.basename(path_xz)}  +  …{os.path.basename(path_xy)}", SUCCESS)
        messagebox.showinfo("Images saved",
                            f"Side view:   {path_xz}\nEnd-on view: {path_xy}")

    # ------------------------------------------------------------------
    # Export model XML
    # ------------------------------------------------------------------

    def _on_export_xml(self):
        out_dir = filedialog.askdirectory(title="Select output directory for model XML files")
        if not out_dir:
            return

        self.btn_export_xml.config(state="disabled", text="Exporting…")
        self._set_status("Exporting model XML files…", FG_DIM)

        builder         = self._builder
        geometry        = self._geometry
        photon          = self.photon_var.get()
        ww_on           = self.ww_var.get()
        particles       = self.sim_particles_v.get()
        batches         = self.sim_batches_v.get()
        sp_freq         = self.sim_sp_freq_v.get()
        src_type        = self.src_type_v.get()
        src_yaml_path   = self.src_yaml_v.get().strip()

        def _worker():
            try:
                import openmc

                os.makedirs(out_dir, exist_ok=True)
                orig = os.getcwd()
                os.chdir(out_dir)
                try:
                    geometry.export_to_xml("geometry.xml")

                    mats = openmc.Materials(builder.material_ns.materials)
                    mats.export_to_xml("materials.xml")

                    tallies = openmc.Tallies(builder.get_tallies())
                    tallies.export_to_xml("tallies.xml")

                    settings = openmc.Settings()
                    settings.run_mode = "fixed source"
                    settings.particles = particles
                    settings.batches   = batches
                    settings.output    = {"tallies": False}
                    settings.photon_transport = photon
                    settings.statepoint = {
                        "batches": [1]
                        + list(range(sp_freq, batches, sp_freq))
                        + [batches]
                    }

                    # Source
                    if src_type == "From YAML file" and src_yaml_path:
                        from paratan2.source.volumetric import load_source_from_yaml
                        settings.source = load_source_from_yaml(src_yaml_path)
                    elif src_type == "Volumetric (plasma)":
                        from paratan2.source.volumetric import VolumetricSource
                        src_r   = self.src_r_v.get()
                        src_zhl = self.src_zhl_v.get()
                        settings.source = VolumetricSource(
                            radius=src_r, half_length=src_zhl
                        ).to_openmc_source()
                    elif src_type == "Uniform cylinder":
                        from paratan2.source.volumetric import UniformSource
                        src_r   = self.src_r_v.get()
                        src_zhl = self.src_zhl_v.get()
                        settings.source = UniformSource(
                            radius=src_r, half_length=src_zhl
                        ).to_openmc_source()

                    if ww_on:
                        wwgs = builder.get_weight_window_generators()
                        settings.weight_window_generators = wwgs
                        settings.weight_windows_on = True
                        settings.weight_window_checkpoints = {"collision": True, "surface": True}

                    settings.export_to_xml("settings.xml")
                finally:
                    os.chdir(orig)

                self.root.after(0, lambda: self._on_export_done(out_dir))
            except Exception as e:
                import traceback
                tb = traceback.format_exc()
                self.root.after(0, lambda: self._on_export_error(tb))

        threading.Thread(target=_worker, daemon=True).start()

    def _on_export_done(self, out_dir):
        self.btn_export_xml.config(state="normal", text="Export Model XML…")
        self._set_status(f"XML files written to: {out_dir}", SUCCESS)
        messagebox.showinfo("Export complete",
                            f"geometry.xml  ·  materials.xml  ·  tallies.xml  ·  settings.xml\n\n{out_dir}")

    def _on_export_error(self, msg):
        self.btn_export_xml.config(state="normal", text="Export Model XML…")
        self._set_status("Export failed", ERROR)
        messagebox.showerror("Export failed", msg)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    try:
        root = tk.Tk()
    except tk.TclError as e:
        err = str(e).lower()
        if "display" in err:
            print(
                "This GUI needs a display. Examples:\n"
                "  • From a desktop session, run: python gui.py\n"
                "  • Over SSH with forwarding: ssh -Y user@host  then  python gui.py\n"
                "  • Headless / no GPU session: xvfb-run -a python gui.py\n",
                file=sys.stderr,
            )
            sys.exit(1)
        raise
    app = ParatanApp(root)
    root.mainloop()
