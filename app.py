"""MD trajectory viewer: Streamlit + MDAnalysis + py3Dmol."""
from pathlib import Path

import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

from mdviewer import analysis, io

DEMO = Path(__file__).parent / "demo"
DEMO_TOP, DEMO_TRAJ = DEMO / "egfr_T790M_osimertinib.pdb", DEMO / "egfr_T790M_osimertinib.xtc"

st.set_page_config(page_title="MD trajectory viewer", layout="wide")
st.title("MD trajectory viewer")


@st.cache_resource(show_spinner="Loading trajectory…")
def get_universe(top: str, traj: str):
    return io.load_universe(top, traj)


@st.cache_data(show_spinner="Computing…")
def cached_series(top, traj, kind, select, stride):
    u = get_universe(top, traj)
    if kind == "rmsd":
        return analysis.rmsd_series(u, select, stride=stride)
    if kind == "rg":
        return analysis.rg_series(u, select, stride=stride)
    raise ValueError(kind)


@st.cache_data(show_spinner="Computing RMSF…")
def cached_rmsf(top, traj, select, stride):
    # RMSF fits the trajectory in memory, so it works on a private copy.
    return analysis.rmsf_per_residue(io.load_universe(top, traj), select, stride)


def viewer_html(pdb: str, style: str, animate: bool, highlight: str | None) -> str:
    spec = {"cartoon": "{cartoon:{color:'spectrum'}}", "stick": "{stick:{}}", "line": "{line:{}}"}[style]
    anim = "v.animate({loop:'forward',interval:120});" if animate else ""
    hl = f"v.addStyle({{resn:'{highlight}'}},{{stick:{{colorscheme:'orangeCarbon',radius:0.25}}}});" if highlight else ""
    return f"""
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"></script>
<div id="v" style="width:100%;height:520px;position:relative"></div>
<script>
const v = $3Dmol.createViewer('v', {{backgroundColor:'white'}});
v.addModelsAsFrames(`{pdb}`, 'pdb');
v.setStyle({{}}, {spec});
{hl}
v.zoomTo(); v.render(); {anim}
</script>"""


with st.sidebar:
    st.header("Data")
    source = st.radio("Source", ["Demo (EGFR T790M + osimertinib)", "Upload"])
    if source.startswith("Demo"):
        top, traj = str(DEMO_TOP), str(DEMO_TRAJ)
    else:
        up_top = st.file_uploader("Topology", type=io.TOPOLOGY_TYPES)
        up_traj = st.file_uploader("Trajectory", type=io.TRAJECTORY_TYPES)
        if not (up_top and up_traj):
            st.info("Upload a topology and a trajectory (200 MB limit each).")
            st.stop()
        wd = st.session_state.setdefault("workdir", io.new_workdir())
        top, traj = str(io.save_upload(up_top, wd)), str(io.save_upload(up_traj, wd))

    st.header("Analysis")
    stride = st.number_input("Stride (every Nth frame)", 1, 1000, 1)
    sel_view = st.text_input("View selection", "protein or resname MOL")
    sel_ana = st.text_input("Analysis selection", "protein and name CA")
    style = st.selectbox("Style", ["cartoon", "stick", "line"])
    hl = st.text_input("Highlight residue name (stick)", "MOL")

try:
    u = get_universe(top, traj)
    u.select_atoms(sel_view)
except Exception as e:  # invalid files or selection syntax
    st.error(f"Could not load or parse selection: {e}")
    st.stop()

n = len(u.trajectory)
dt = getattr(u.trajectory, "dt", 0.0)
st.caption(f"{u.atoms.n_atoms} atoms · {n} frames · {dt:.1f} ps/frame")

col3d, colplot = st.columns([3, 2])
with col3d:
    frame = st.slider("Frame", 0, n - 1, 0)
    animate = st.checkbox(f"Play in browser (≤{io.MAX_ANIMATION_FRAMES} frames, subsampled)")
    try:
        pdb = io.multi_model_pdb(u, sel_view, stride) if animate else io.frame_pdb(u, frame, sel_view)
        components.html(viewer_html(pdb, style, animate, hl or None), height=540)
    except Exception as e:
        st.error(f"Selection error: {e}")

with colplot:
    try:
        rmsd = cached_series(top, traj, "rmsd", sel_ana, stride)
        rg = cached_series(top, traj, "rg", sel_view, stride)
        t = [i * stride * dt / 1000 for i in range(len(rmsd))]
        for name, y, unit in [("RMSD", rmsd, "Å"), ("Radius of gyration", rg, "Å")]:
            fig = go.Figure(go.Scatter(x=t, y=y, mode="lines"))
            fig.add_vline(x=frame * dt / 1000, line_dash="dot")
            fig.update_layout(title=name, xaxis_title="time (ns)", yaxis_title=unit,
                              height=250, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig, use_container_width=True)
    except Exception as e:
        st.error(f"Analysis error: {e}")

st.subheader("Per-residue RMSF")
try:
    resid, f = cached_rmsf(top, traj, sel_ana, stride)
    fig = go.Figure(go.Scatter(x=resid, y=f, mode="lines"))
    fig.update_layout(xaxis_title="residue", yaxis_title="RMSF (Å)", height=260,
                      margin=dict(l=10, r=10, t=10, b=10))
    st.plotly_chart(fig, use_container_width=True)
except Exception as e:
    st.error(f"RMSF error: {e}")
