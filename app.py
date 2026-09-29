"""Streamlit frontend. All molecular work is delegated to martini_solv.builder."""
from __future__ import annotations

import io
import json
from pathlib import Path
import tempfile
import zipfile

import streamlit as st

from martini_solv.builder import build
from martini_solv.models import BuildConfig, Solute, reline_counts
from martini_solv.theme import STYLE


st.set_page_config(page_title="MartiniSolv", page_icon="🧪", layout="wide")
st.markdown(STYLE, unsafe_allow_html=True)
st.markdown('''<div class="hero"><div class="eyebrow">MARTINI 3 · PROTEIN SYSTEM BUILDER</div>
<h1>MartiniSolv</h1><p>From a protein structure to a solvated coarse-grained GROMACS system.</p>
<span class="pill">Protein</span><span class="pill">Water / Reline</span>
<span class="pill">Additional molecules</span></div>''', unsafe_allow_html=True)

with st.sidebar:
    st.markdown("### Build a system")
    step = st.radio("Workflow", ["Protein", "Martini model", "Solvent & molecules", "Review & build", "Results"])
    st.caption("Martini 3 · A small, independent project inspired by MartiniSurf")

if step == "Protein":
    st.header("01 / Protein")
    uploaded = st.file_uploader("Upload protein PDB", type=["pdb"])
    if uploaded:
        st.session_state["pdb_bytes"] = uploaded.getvalue()
        st.session_state["pdb_name"] = uploaded.name
    if st.session_state.get("pdb_bytes"):
        lines = st.session_state["pdb_bytes"].decode("utf-8", "replace").splitlines()
        atom_count = sum(x.startswith("ATOM  ") for x in lines)
        st.metric("ATOM records", atom_count)
        if atom_count == 0:
            st.error("This file contains no protein ATOM records.")
        st.caption("The original PDB is preserved. Resolve missing backbone atoms and chain breaks before building.")

elif step == "Martini model":
    st.header("02 / Martini 3 model")
    st.markdown('<div class="panel"><h3>Protein coarse graining</h3><p class="subtle">martinize2 · martini3001 · DSSP · elastic network</p></div>', unsafe_allow_html=True)
    st.info("The first version uses one conservative protein model. Its generated topology and coordinates are included in the download.")
    st.number_input("Box margin (nm)", min_value=0.5, max_value=5.0, value=1.2, step=0.1, key="box_distance")

elif step == "Solvent & molecules":
    st.header("03 / Solvent & molecules")
    st.selectbox("Solvent", ["Water", "Reline · ChCl:urea 1:2"], key="solvent_ui")
    if st.session_state.solvent_ui.startswith("Water"):
        st.number_input("Added NaCl (M)", min_value=0.0, max_value=2.0, value=0.15, step=0.05, key="salt")
    else:
        st.slider("Water mole fraction in reline", min_value=0.0, max_value=0.4,
                  value=0.0, step=0.01, key="water_fraction")
        st.caption("0 = dry. 0.10 is an example of hydrated reline; other fractions require independent validation.")
        st.caption("One Martini W bead represents four water molecules; the actual composition is reported in the manifest.")
    st.subheader("Additional free molecules")
    st.caption("One species per row. Models generated from SMILES require scientific validation.")
    default = [{"name": "", "smiles": "", "count": 1}]
    st.data_editor(default, num_rows="dynamic", key="solute_rows",
                   column_config={"name": "Name", "smiles": "SMILES", "count": st.column_config.NumberColumn("Copies", min_value=1, step=1)})

elif step == "Review & build":
    st.header("04 / Review & build")
    solvent = "reline" if st.session_state.get("solvent_ui", "Water").startswith("Reline") else "water"
    rows = st.session_state.get("solute_rows", {})
    solute_rows = rows.get("edited_rows", {})
    original = [{"name": "", "smiles": "", "count": 1}]
    for index, patch in solute_rows.items():
        original[int(index)].update(patch)
    original.extend(rows.get("added_rows", []))
    for index in sorted(rows.get("deleted_rows", []), reverse=True):
        original.pop(index)
    try:
        solutes = [Solute(str(r.get("name", "")).strip(), str(r.get("smiles", "")).strip(),
                          int(r.get("count", 1))) for r in original if r.get("name") or r.get("smiles")]
        config = BuildConfig(solvent=solvent, water_fraction=float(st.session_state.get("water_fraction", 0)) if solvent == "reline" else 0,
                             salt_m=float(st.session_state.get("salt", 0.15)) if solvent == "water" else 0,
                             box_distance_nm=float(st.session_state.get("box_distance", 1.2)), solutes=solutes)
        config.validate()
        st.json(json.loads(config.to_json()))
        if not st.session_state.get("pdb_bytes"):
            st.warning("Upload a PDB in the Protein step.")
        elif st.button("Build system", type="primary"):
            with st.status("Building system", expanded=True):
                root = Path(tempfile.mkdtemp(prefix="martinisolv_"))
                input_path = root / "input.pdb"
                input_path.write_bytes(st.session_state.pdb_bytes)
                try:
                    built = build(input_path, root / "system", config)
                    data = io.BytesIO()
                    with zipfile.ZipFile(data, "w", zipfile.ZIP_DEFLATED) as archive:
                        for file in sorted(built.rglob("*")):
                            if file.is_file():
                                archive.write(file, file.relative_to(built))
                    st.session_state["result_zip"] = data.getvalue()
                    st.session_state["result_manifest"] = (built / "manifest.json").read_text()
                    st.success("GROMACS grompp passed. The prepared system is ready for minimization.")
                except Exception as exc:
                    st.error(str(exc))
                    st.caption(f"Diagnostic files: {root / 'system.incomplete'}")
    except (ValueError, RuntimeError) as exc:
        st.error(str(exc))

else:
    st.header("05 / Results")
    if st.session_state.get("result_zip"):
        st.json(json.loads(st.session_state.result_manifest))
        st.download_button("Download GROMACS package", st.session_state.result_zip,
                           file_name="martini-solv-system.zip", mime="application/zip", type="primary")
        st.caption("Includes system.gro, system.top, protein CG files, force-field includes, minimization.mdp, build.log and provenance.")
    else:
        st.info("Build a system to access the download.")
