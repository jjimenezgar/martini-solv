"""Protein-in-solvent Streamlit workflow; molecular work lives in the builder."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import urllib.request

import streamlit as st

from martini_solv.builder import _map_solute, build
from martini_solv.models import BuildConfig, Solute
from martini_solv.short_md import archive, run_short_md
from martini_solv.theme import STYLE

STEPS = ["Structure", "Model", "Environment", "Review & Build", "Short MD"]
st.set_page_config(page_title="MartiniSolv", page_icon="🧪", layout="wide")
st.markdown(STYLE, unsafe_allow_html=True)
st.markdown('''<div class="hero"><div class="eyebrow">MARTINI 3 · PROTEIN IN SOLUTION</div>
<h1>MartiniSolv</h1><p>Prepare a coarse-grained protein in water or hydrated reline, with optional freely dissolved molecules.</p>
<span class="pill">Protein</span><span class="pill">Water / Reline</span><span class="pill">Free molecules</span></div>''', unsafe_allow_html=True)
st.session_state.setdefault("solutes", [])
st.session_state.setdefault("active_step", "Structure")

with st.sidebar:
    st.markdown("### Build a system")
    step = st.radio("Workflow step", STEPS, key="active_step", label_visibility="collapsed")
    st.caption("Structure → Model → Environment → Review & Build → Short MD")
    st.caption("Martini 3 · Inspired by MartiniSurf, without surfaces or linkers")


def _protein_summary(pdb: bytes) -> tuple[int, list[dict[str, object]]]:
    residues: dict[str, set[tuple[str, str]]] = {}
    atoms = 0
    for line in pdb.decode("utf-8", "replace").splitlines():
        if line.startswith("ATOM  "):
            atoms += 1
            residues.setdefault(line[21:22].strip() or "A", set()).add((line[22:27], line[17:20]))
    return atoms, [{"Chain": chain, "Residues": len(items)} for chain, items in sorted(residues.items())]


def _show_mapping(gro: Path) -> None:
    lines = gro.read_text().splitlines()
    bead_count = int(lines[1])
    st.dataframe([{"Bead": i + 1, "Residue": row[5:10].strip(), "Name": row[10:15].strip(),
                   "x (nm)": row[20:28].strip(), "y (nm)": row[28:36].strip(),
                   "z (nm)": row[36:44].strip()}
                  for i, row in enumerate(lines[2:2 + bead_count])], hide_index=True, use_container_width=True)


def _show_molecule(smiles: str) -> None:
    from rdkit import Chem
    from rdkit.Chem import Draw
    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        st.image(Draw.MolToImage(molecule, size=(460, 280)), caption="2D structure from SMILES")


def _solute_generator() -> None:
    with st.expander("Free molecule generator · SMILES", expanded=True):
        st.caption("Generate a Martini 3 model for a freely dissolved molecule. It is never bonded or oriented toward the protein.")
        left, right = st.columns([2.2, 0.8])
        smiles = left.text_input("SMILES", key="new_smiles", placeholder="CCO")
        name = right.text_input("Molecule name", key="new_name", placeholder="ETOH")
        copies = st.number_input("Number of copies", min_value=1, max_value=1000, value=1, key="new_copies")
        if smiles:
            try:
                _show_molecule(smiles)
            except ImportError:
                st.warning("RDKit is unavailable; SMILES preview cannot be shown.")
        if st.button("Generate molecule", type="primary", use_container_width=True):
            try:
                species = Solute(name.strip(), smiles.strip(), int(copies))
                species.validate()
                if any(row["name"].upper() == species.name.upper() for row in st.session_state.solutes):
                    raise ValueError("Choose a unique molecule name")
                with st.spinner("Running Martini Mapper"):
                    root = Path(tempfile.mkdtemp(prefix="martinisolv_mapper_"))
                    gro, itp = _map_solute(root, species.name, species.smiles, root / "mapper.log")
                st.session_state.solutes.append({"name": species.name, "smiles": species.smiles,
                                                 "count": species.count, "gro": str(gro), "itp": str(itp)})
                st.success(f"{species.name}: Martini Mapper generated .gro and .itp files.")
            except (ValueError, RuntimeError, OSError) as exc:
                st.error(str(exc))
        for index, row in enumerate(st.session_state.solutes):
            with st.container(border=True):
                title, action = st.columns([3, 1])
                title.markdown(f"**{row['name']}** · {row['count']} copies · `{row['smiles']}`")
                if action.button("Remove", key=f"remove_{index}"):
                    st.session_state.solutes.pop(index)
                    st.rerun()
                with st.expander(f"Preview {row['name']} coarse-grained beads"):
                    gro = Path(row["gro"])
                    if gro.is_file():
                        _show_mapping(gro)
                    else:
                        st.info("Preview expired; regenerate the molecule. It will also be generated during build.")
        st.caption("Review the generated mapping and topology before scientific use. Only neutral free molecules are accepted.")


def _config() -> BuildConfig:
    solvent = "reline" if st.session_state.get("solvent_ui", "Water") == "Reline · ChCl:urea 1:2" else "water"
    return BuildConfig(
        solvent=solvent,
        salt_m=float(st.session_state.get("salt", 0.15)) if solvent == "water" else 0,
        water_fraction=float(st.session_state.get("water_fraction", 0)) if solvent == "reline" else 0,
        box_distance_nm=float(st.session_state.get("box_distance", 1.2)),
        elastic=bool(st.session_state.get("elastic", True)),
        elastic_force=int(st.session_state.get("elastic_force", 700)),
        solutes=[Solute(row["name"], row["smiles"], int(row["count"])) for row in st.session_state.solutes],
    )


if step == "Structure":
    st.header("01 / Structure")
    left, right = st.columns([0.9, 1.1], gap="large")
    with left:
        st.subheader("Protein structure")
        st.text_input("Project name", key="project_name", placeholder="My protein")
        uploaded = st.file_uploader("Protein PDB", type=["pdb"], key="structure_upload")
        if uploaded:
            st.session_state.pdb_bytes = uploaded.getvalue()
            st.session_state.pdb_name = uploaded.name
        pdb_id = st.text_input("Or fetch a PDB ID", placeholder="1UBQ", key="pdb_id").strip().upper()
        if st.button("Fetch PDB", disabled=not bool(pdb_id)):
            if len(pdb_id) != 4 or not pdb_id.isalnum():
                st.error("Enter a four-character PDB ID")
            else:
                try:
                    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=20) as response:
                        st.session_state.pdb_bytes = response.read()
                    st.session_state.pdb_name = f"{pdb_id}.pdb"
                except Exception as exc:
                    st.error(f"Could not fetch {pdb_id}: {exc}")
    with right:
        st.subheader("Structure preview")
        if st.session_state.get("pdb_bytes"):
            atoms, chains = _protein_summary(st.session_state.pdb_bytes)
            a, b = st.columns(2)
            a.metric("ATOM records", atoms)
            b.metric("Chains", len(chains))
            st.dataframe(chains, hide_index=True, use_container_width=True)
            st.caption(f"Input: {st.session_state.get('pdb_name', 'protein.pdb')}. Water and ligands are excluded from the protein sent to martinize2.")
            if not atoms:
                st.error("No protein ATOM records found")
        else:
            st.info("Upload a PDB or fetch a PDB ID. Try 1UBQ as a small example.")

elif step == "Model":
    st.header("02 / Model")
    st.subheader("Martini 3 protein model")
    st.caption("martinize2 · martini3001 · DSSP secondary structure")
    st.toggle("Elastic network", value=True, key="elastic", help="Optional martinize2 elastic network for protein structure preservation.")
    if st.session_state.elastic:
        st.number_input("Elastic force constant (kJ/mol/nm²)", min_value=100, max_value=1500,
                        value=700, step=50, key="elastic_force")
    st.info("DSSP is enabled. The generated topology and CG coordinates can be inspected in the download. GōMartini is not configured in this compact workflow.")

elif step == "Environment":
    st.header("03 / Environment")
    st.selectbox("Solvent", ["Water", "Reline · ChCl:urea 1:2"], key="solvent_ui")
    st.number_input("Protein-to-box distance (nm)", min_value=0.5, max_value=5.0,
                    value=1.2, step=0.1, key="box_distance")
    if st.session_state.solvent_ui == "Water":
        st.number_input("NaCl concentration (M)", min_value=0.0, max_value=2.0,
                        value=0.15, step=0.05, key="salt")
        st.caption("INSANE places Martini 3 water and ions around the protein and neutralizes its charge.")
    else:
        st.slider("Water mole fraction in reline", min_value=0.0, max_value=0.4,
                  value=0.0, step=0.01, key="water_fraction")
        st.caption("0 = dry. Hydration is adjustable; user-defined compositions require validation. Added NaCl is unavailable in reline.")
    _solute_generator()

elif step == "Review & Build":
    st.header("04 / Review & Build")
    try:
        config = _config()
        config.validate()
        st.json(json.loads(config.to_json()))
        pdb = st.session_state.get("pdb_bytes")
        if not pdb or not _protein_summary(pdb)[0]:
            st.warning("Provide a protein PDB in Structure first.")
        if st.button("Build system", type="primary", disabled=not bool(pdb and _protein_summary(pdb)[0])):
            root = Path(tempfile.mkdtemp(prefix="martinisolv_"))
            input_path = root / "input.pdb"
            input_path.write_bytes(pdb)
            with st.status("Building Martini 3 system", expanded=True) as status:
                try:
                    built = build(input_path, root / "system", config)
                    st.session_state.result_dir = str(built)
                    st.session_state.result_manifest = (built / "manifest.json").read_text()
                    status.update(label="GROMACS-ready system generated", state="complete")
                    st.success("gmx grompp passed with zero warnings. The system still needs minimization and equilibration.")
                except Exception as exc:
                    status.update(label="Build failed", state="error")
                    st.error(str(exc))
                    log = root / "system.incomplete" / "build.log"
                    if log.is_file():
                        st.code(log.read_text(errors="replace")[-8000:], language="text")
        built = Path(st.session_state.get("result_dir", ""))
        if built.is_dir() and (built / "manifest.json").is_file():
            st.download_button("Download GROMACS package", archive(built), file_name="martini-solv-system.zip",
                               mime="application/zip", type="primary")
            st.caption("Includes system.gro, system.top, Martini force field, molecular ITPs, minimization.mdp and provenance log.")
    except (ValueError, RuntimeError) as exc:
        st.error(str(exc))

elif step == "Short MD":
    st.header("05 / Short MD")
    st.caption("Optional local sanity check: energy minimization followed by a short NVT run at 300 K. No production MD.")
    built = Path(st.session_state.get("result_dir", ""))
    if not (built.is_dir() and (built / "system.gro").is_file()):
        st.info("Build a system in Review & Build first.")
    else:
        a, b = st.columns(2)
        duration = a.number_input("NVT duration (ps)", min_value=1.0, max_value=100.0, value=20.0, step=1.0)
        threads = b.number_input("CPU threads", min_value=1, max_value=4, value=2)
        if st.button("Run Short MD", type="primary"):
            with st.status("Minimization → NVT", expanded=True) as status:
                try:
                    run_short_md(built, duration, threads)
                    status.update(label="Short MD completed", state="complete")
                    st.success("Sanity check completed. Inspect the energies, log and structure before further simulation.")
                except (ValueError, RuntimeError, FileExistsError) as exc:
                    status.update(label="Short MD failed", state="error")
                    st.error(str(exc))
        log = built / "short_md" / "short_md.log"
        if log.is_file():
            with st.expander("Short MD execution log"):
                st.code(log.read_text(errors="replace")[-10000:], language="text")
        if (built / "short_md" / "summary.json").is_file():
            st.json(json.loads((built / "short_md" / "summary.json").read_text()))
        st.download_button("Download system and Short MD files", archive(built),
                           file_name="martini-solv-system.zip", mime="application/zip")
