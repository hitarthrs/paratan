"""
Paratan 2 — Simple Mirror Machine Builder (Streamlit UI).

  pip install streamlit                    # or: pip install -e ".[gui]" here
  PYTHONPATH=src streamlit run app.py      # omit PYTHONPATH if paratan2 is installed editable
"""

import os
import tempfile
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

from paratan2.config.models import (
    SimpleMachineConfig,
    VacuumVesselConfig, VacuumVesselStructureLayer,
    CentralCellConfig, MaterialLayer,
    LFCoilConfig, LFCoilShellThicknesses, LFCoilInnerDimensions, LFCoilMaterials,
    HFCoilConfig, HFCoilMagnet, HFCoilShield,
    EndCellConfig,
    RoomConfig,
)

# ---------------------------------------------------------------------------
# Page setup
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="Paratan 2 — Simple Mirror Builder",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("Paratan 2 — Simple Mirror Machine Builder")
st.caption("Configure each component in the sidebar, then hit **Build Geometry**.")

# ---------------------------------------------------------------------------
# Material options (common ones from the library)
# ---------------------------------------------------------------------------

STRUCTURAL_MATERIALS = [
    "stainless", "ss316ln", "aluminum_6061", "tungsten", "rafm_steel",
    "cooled_rafm_steel", "cooled_tungsten", "cooled_tungsten_carbide",
]
BREEDER_MATERIALS = [
    "eutectic_breeding_material", "LiPb_breeder", "flibe_breeding_material",
    "hcpb_breeding_material", "dcll_enriched_lithium_breeding_material",
    "lithium_metatitanate", "lithium_orthosilicate",
]
SHIELD_MATERIALS = [
    "HT_Shield_filler", "LT_Shield_filler", "cooled_w2b5", "cooled_wb_shield",
    "cooled_wb2_shield", "cooled_wc_shield", "water_cooled_wc", "cooled_B4C",
    "B4C", "tungsten_carbide", "cooled_tungsten_boride",
]
MAGNET_MATERIALS = ["Magnet_Winding_Pack_2", "Magnet_Winding_Pack", "rebco", "copper"]
WALL_MATERIALS = ["tungsten", "FW_FNSF", "BW_FNSF", "stainless"]
ALL_MATERIALS = sorted(set(
    STRUCTURAL_MATERIALS + BREEDER_MATERIALS + SHIELD_MATERIALS
    + MAGNET_MATERIALS + WALL_MATERIALS
    + ["vacuum", "air", "deuterium"]
))


def mat_select(label, default, options=None, key=None):
    opts = options or ALL_MATERIALS
    idx = opts.index(default) if default in opts else 0
    return st.selectbox(label, opts, index=idx, key=key)


# ---------------------------------------------------------------------------
# Sidebar — component inputs
# ---------------------------------------------------------------------------

sb = st.sidebar

# ---- Vacuum Vessel ----
sb.header("Vacuum Vessel")
vv_outer_axial   = sb.number_input("Outer axial length (cm)",   value=450.0,  min_value=10.0,  key="vv_outer")
vv_central_axial = sb.number_input("Central axial length (cm)", value=600.0,  min_value=10.0,  key="vv_central")
vv_central_r     = sb.number_input("Central radius (cm)",       value=75.0,   min_value=1.0,   key="vv_cr")
vv_bottle_r      = sb.number_input("Bottleneck radius (cm)",    value=25.0,   min_value=1.0,   key="vv_br")
vv_left_bn       = sb.number_input("Left bottleneck length (cm)",  value=300.0, min_value=1.0, key="vv_lbn")
vv_right_bn      = sb.number_input("Right bottleneck length (cm)", value=300.0, min_value=1.0, key="vv_rbn")

sb.markdown("**VV Structure layers** (inner → outer)")
with sb.expander("First wall"):
    fw_on  = st.checkbox("Include", value=True, key="fw_on")
    fw_t   = st.number_input("Thickness (cm)", value=1.0, min_value=0.1, key="fw_t")
    fw_mat = mat_select("Material", "tungsten", WALL_MATERIALS, key="fw_mat")
with sb.expander("VV shell"):
    vs_on  = st.checkbox("Include", value=True, key="vs_on")
    vs_t   = st.number_input("Thickness (cm)", value=2.0, min_value=0.1, key="vs_t")
    vs_mat = mat_select("Material", "FW_FNSF", WALL_MATERIALS, key="vs_mat")

# ---- Central Cell ----
sb.header("Central Cell")
cc_axial = sb.number_input("Axial length (cm)", value=550.0, min_value=10.0, key="cc_axial")
sb.markdown("**Radial layers** (inner → outer)")
with sb.expander("Breeding blanket"):
    bl_t   = sb.number_input("Thickness (cm)", value=77.0, min_value=1.0, key="bl_t")
    bl_mat = mat_select("Material", "eutectic_breeding_material", BREEDER_MATERIALS, key="bl_mat")
with sb.expander("Back wall"):
    bw_t   = sb.number_input("Thickness (cm)", value=2.0, min_value=0.1, key="bw_t")
    bw_mat = mat_select("Material", "BW_FNSF", WALL_MATERIALS, key="bw_mat")
with sb.expander("HT Shield"):
    sh_t   = sb.number_input("Thickness (cm)", value=30.0, min_value=1.0, key="sh_t")
    sh_mat = mat_select("Material", "HT_Shield_filler", SHIELD_MATERIALS, key="sh_mat")
with sb.expander("Shield back plate"):
    sbp_t   = sb.number_input("Thickness (cm)", value=2.5, min_value=0.1, key="sbp_t")
    sbp_mat = mat_select("Material", "water_cooled_wc", SHIELD_MATERIALS, key="sbp_mat")

# ---- LF Coils ----
sb.header("LF (Solenoid) Coils")
lf_positions_str = sb.text_input("Coil positions (cm, comma-separated)", value="-75, 75", key="lf_pos")
lf_axial_len   = sb.number_input("Winding pack axial length (cm)",   value=20.0,  min_value=1.0, key="lf_al")
lf_radial_t    = sb.number_input("Winding pack radial thickness (cm)", value=10.0, min_value=1.0, key="lf_rt")
lf_front_sh    = sb.number_input("Shell front thickness (cm)", value=10.0, min_value=0.1, key="lf_fs")
lf_back_sh     = sb.number_input("Shell back thickness (cm)",  value=5.0,  min_value=0.1, key="lf_bs")
lf_axial_sh    = sb.number_input("Shell axial thickness (cm)", value=10.0, min_value=0.1, key="lf_as")
lf_shield_mat  = mat_select("Shield material", "stainless", STRUCTURAL_MATERIALS, key="lf_sm")
lf_magnet_mat  = mat_select("Magnet material", "Magnet_Winding_Pack_2", MAGNET_MATERIALS, key="lf_mm")

# ---- HF Coils ----
sb.header("HF (Mirror) Coils")
hf_bore_r   = sb.number_input("Bore radius (cm)",           value=60.0,  min_value=1.0,  key="hf_br")
hf_rad_t    = sb.number_input("Radial thickness (cm)",       value=150.0, min_value=10.0, key="hf_rt")
hf_ax_t     = sb.number_input("Axial thickness (cm)",        value=200.0, min_value=10.0, key="hf_at")
hf_mag_mat  = mat_select("Magnet material", "Magnet_Winding_Pack_2", MAGNET_MATERIALS, key="hf_mm")
hf_sh_gap   = sb.number_input("Shield–CC axial gap (cm)",    value=5.0,   min_value=0.0,  key="hf_sg")
hf_rad_gap  = sb.number_input("Radial gap before casing (cm)", value=1.0, min_value=0.0,  key="hf_rg")
hf_sh_ri    = sb.number_input("Shield radial thickness inner (cm)", value=20.0, min_value=1.0, key="hf_sri")
hf_sh_ro    = sb.number_input("Shield radial thickness outer (cm)", value=30.0, min_value=1.0, key="hf_sro")
hf_sh_ai    = sb.number_input("Shield axial thickness toward midplane (cm)", value=20.0, min_value=1.0, key="hf_sai")
hf_sh_ao    = sb.number_input("Shield axial thickness away (cm)", value=10.0, min_value=1.0, key="hf_sao")
hf_sh_mat   = mat_select("Shield material", "cooled_w2b5", SHIELD_MATERIALS, key="hf_sm")

# ---- End Cells ----
sb.header("End Cells")
ec_axial   = sb.number_input("Axial length (cm)",    value=500.0, min_value=10.0,  key="ec_ax")
ec_shell_t = sb.number_input("Shell thickness (cm)", value=3.0,   min_value=0.1,   key="ec_st")
ec_diam    = sb.number_input("Diameter (cm)",         value=350.0, min_value=10.0,  key="ec_d")
ec_sh_mat  = mat_select("Shell material", "stainless", STRUCTURAL_MATERIALS, key="ec_sm")

# ---- Room ----
with sb.expander("Room (bounding box)"):
    room_xy  = sb.number_input("XY half-extent (cm)", value=400.0, min_value=100.0, key="room_xy")
    room_z   = sb.number_input("Z half-extent (cm)",  value=1350.0, min_value=100.0, key="room_z")

# ---------------------------------------------------------------------------
# Assemble config from sidebar values
# ---------------------------------------------------------------------------

def build_config() -> SimpleMachineConfig:
    vv_structure = {}
    if fw_on:
        vv_structure["first_wall"] = VacuumVesselStructureLayer(thickness=fw_t, material=fw_mat)
    if vs_on:
        vv_structure["vv_shell"] = VacuumVesselStructureLayer(thickness=vs_t, material=vs_mat)

    lf_positions = [float(x.strip()) for x in lf_positions_str.split(",") if x.strip()]

    return SimpleMachineConfig(
        vacuum_vessel=VacuumVesselConfig(
            outer_axial_length=vv_outer_axial,
            central_axial_length=vv_central_axial,
            central_radius=vv_central_r,
            bottleneck_radius=vv_bottle_r,
            left_bottleneck_length=vv_left_bn,
            right_bottleneck_length=vv_right_bn,
            structure=vv_structure,
        ),
        central_cell=CentralCellConfig(
            axial_length=cc_axial,
            layers=[
                MaterialLayer(thickness=bl_t,  material=bl_mat),
                MaterialLayer(thickness=bw_t,  material=bw_mat),
                MaterialLayer(thickness=sh_t,  material=sh_mat),
                MaterialLayer(thickness=sbp_t, material=sbp_mat),
            ],
        ),
        lf_coil=LFCoilConfig(
            shell_thicknesses=LFCoilShellThicknesses(front=lf_front_sh, back=lf_back_sh, axial=lf_axial_sh),
            inner_dimensions=LFCoilInnerDimensions(axial_length=lf_axial_len, radial_thickness=lf_radial_t),
            positions=lf_positions,
            materials=LFCoilMaterials(shield=lf_shield_mat, magnet=lf_magnet_mat),
        ),
        hf_coil=HFCoilConfig(
            magnet=HFCoilMagnet(
                bore_radius=hf_bore_r,
                radial_thickness=hf_rad_t,
                axial_thickness=hf_ax_t,
                material=hf_mag_mat,
            ),
            casing_layers=[
                MaterialLayer(thickness=1.0, material="vacuum"),
                MaterialLayer(thickness=1.0, material="ss316ln"),
            ],
            shield=HFCoilShield(
                radial_gap_before_casing=hf_rad_gap,
                shield_central_cell_gap=hf_sh_gap,
                radial_thickness=[hf_sh_ri, hf_sh_ro],
                axial_thickness=[hf_sh_ai, hf_sh_ao],
                material=hf_sh_mat,
            ),
        ),
        end_cell=EndCellConfig(
            axial_length=ec_axial,
            shell_thickness=ec_shell_t,
            diameter=ec_diam,
            shell_material=ec_sh_mat,
            inner_material="vacuum",
        ),
        room=RoomConfig(
            xmin=-room_xy, xmax=room_xy,
            ymin=-room_xy, ymax=room_xy,
            zmin=-room_z,  zmax=room_z,
        ),
    )


# ---------------------------------------------------------------------------
# Main area — validation banner + build button
# ---------------------------------------------------------------------------

# Live config validation
try:
    config = build_config()
    st.success("Config valid")
except Exception as e:
    st.error(f"Config error: {e}")
    st.stop()

col1, col2 = st.columns([1, 3])
with col1:
    build_clicked = st.button("Build Geometry", type="primary", use_container_width=True)

# Config summary
with st.expander("Config summary", expanded=False):
    st.json(config.model_dump())

# ---------------------------------------------------------------------------
# Build + plot
# ---------------------------------------------------------------------------

if build_clicked:
    from paratan2.models.simple_machine_builder import SimpleMachineBuilder, _load_materials

    with st.spinner("Building geometry..."):
        try:
            m = _load_materials()
            builder = SimpleMachineBuilder(config, m)

            with tempfile.TemporaryDirectory() as tmpdir:
                orig = os.getcwd()
                os.chdir(tmpdir)
                try:
                    builder.build()
                finally:
                    os.chdir(orig)

            geometry = builder.get_geometry()
            n_cells = len(builder.all_cells)
            st.success(f"Built successfully — {n_cells} cells")

        except Exception as e:
            st.error(f"Build failed: {e}")
            st.exception(e)
            st.stop()

    with st.spinner("Rendering cross-sections..."):
        fig, axes = plt.subplots(1, 2, figsize=(20, 9))

        total_z = vv_central_axial / 2 + max(vv_left_bn, vv_right_bn)
        total_x = (vv_central_r
                   + sum([fw_t if fw_on else 0, vs_t if vs_on else 0])
                   + bl_t + bw_t + sh_t + sbp_t
                   + hf_rad_t + hf_sh_ri + hf_sh_ro + 50)

        geometry.root_universe.plot(
            basis="xz",
            origin=(0, 0, 0),
            width=(total_x * 2.2, total_z * 2.2),
            pixels=(1600, 2400),
            color_by="material",
            axes=axes[0],
        )
        axes[0].set_title("Side view (XZ)", fontsize=13)
        axes[0].set_xlabel("Z (cm)")
        axes[0].set_ylabel("X (cm)")

        geometry.root_universe.plot(
            basis="xy",
            origin=(0, 0, 0),
            width=(total_x * 2.2, total_x * 2.2),
            pixels=(800, 800),
            color_by="material",
            axes=axes[1],
        )
        axes[1].set_title("End-on view (XY) at midplane", fontsize=13)
        axes[1].set_xlabel("X (cm)")
        axes[1].set_ylabel("Y (cm)")

        plt.suptitle("Simple Mirror — Geometry", fontsize=15, fontweight="bold")
        plt.tight_layout()
        st.pyplot(fig)
        plt.close(fig)

    # Metrics
    st.divider()
    m1, m2, m3, m4 = st.columns(4)
    total_vv = vv_left_bn + vv_central_axial + vv_right_bn
    cc_r = (vv_central_r
            + sum([fw_t if fw_on else 0, vs_t if vs_on else 0])
            + bl_t + bw_t + sh_t + sbp_t)
    m1.metric("Total machine length", f"{total_vv:.0f} cm")
    m2.metric("CC blanket outer radius", f"{cc_r:.1f} cm")
    m3.metric("LF coil count", len([x.strip() for x in lf_positions_str.split(",") if x.strip()]))
    m4.metric("Total cells", n_cells)
