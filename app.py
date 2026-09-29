"""Friendly MartiniSolv Streamlit workflow modelled on MartiniSurf."""
from __future__ import annotations

import json
from pathlib import Path
import shlex
import shutil
import tempfile
import urllib.request

import streamlit as st

from martini_solv.builder import _map_solute, build
from martini_solv.models import BuildConfig, Solute
from martini_solv.molecular_viewer import render_cg_structure, render_structure_preview, render_trajectory
from martini_solv.short_md import (
    DEFAULT_GROMPP_MAXWARN,
    DEFAULT_STAGE_SETTINGS,
    DEFAULT_XTC_WRITE_EVERY_PS,
    STAGE_ORDER,
    ShortMDConfig,
    StageSettings,
    archive,
    run_short_md,
    selected_stages,
    validate_stage_order,
)
from martini_solv.theme import STYLE


STEPS = ["Structure", "Model", "Environment", "Review & Build", "Short MD"]
CHAIN_LABELS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

st.set_page_config(page_title="MartiniSolv", page_icon="MS", layout="wide", initial_sidebar_state="expanded")
st.markdown(STYLE, unsafe_allow_html=True)

DEFAULTS = {
    "active_step": "Structure",
    "project_name": "MartiniSolv Protein - 1UBQ",
    "pdb_id": "1UBQ",
    "molecule_name": "Protein",
    "merge_chain_count": 1,
    "dssp": True,
    "go": True,
    "go_eps": 9.414,
    "position_restraints": "backbone",
    "elastic": False,
    "elastic_force": 700,
    "maxwarn": 1,
    "martinize_extra": "",
    "solvent_ui": "Water",
    "box_distance": 1.0,
    "salt": 0.15,
    "water_fraction": 0.0,
    "solutes": [],
    "short_md_output_tag": "Protein MD",
    "short_md_xtc_write_every_ps": DEFAULT_XTC_WRITE_EVERY_PS,
    "short_md_grompp_maxwarn": DEFAULT_GROMPP_MAXWARN,
    "short_md_threads": 2,
    "short_md_view_stage": "production",
    "short_md_view_stride": 1,
}
for key, value in DEFAULTS.items():
    st.session_state.setdefault(key, value)
for stage, defaults in DEFAULT_STAGE_SETTINGS.items():
    st.session_state.setdefault(f"short_md_run_{stage}", bool(defaults["enabled"]))
    st.session_state.setdefault(f"short_md_{stage}_dt", float(defaults["dt_ps"]))
    st.session_state.setdefault(f"short_md_{stage}_time", float(defaults["time_ns"]))


def _sync_pdb_id_from_widget() -> None:
    """Keep the chosen PDB ID even when the Structure widget is not rendered."""
    st.session_state["pdb_id"] = str(st.session_state.get("_pdb_id_input", "")).strip().upper()


@st.cache_data(show_spinner=False, ttl=3600)
def _fetch_pdb(pdb_id: str) -> bytes:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=20) as response:
        return response.read()


def _protein_summary(pdb: bytes) -> tuple[int, list[dict[str, object]]]:
    residues: dict[str, set[tuple[str, str]]] = {}
    atoms = 0
    for line in pdb.decode("utf-8", "replace").splitlines():
        if line.startswith("ATOM  "):
            atoms += 1
            chain = line[21:22].strip() or "A"
            residues.setdefault(chain, set()).add((line[22:27], line[17:20]))
    rows = [{"Chain": chain, "Amino acids": len(items)} for chain, items in sorted(residues.items())]
    return atoms, rows


def _ensure_remote_pdb() -> None:
    pdb_id = str(st.session_state.get("pdb_id", "")).strip().upper()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        return
    if st.session_state.get("pdb_source_id") == pdb_id and st.session_state.get("pdb_bytes"):
        return
    try:
        data = _fetch_pdb(pdb_id)
    except Exception as exc:
        st.session_state["pdb_fetch_error"] = str(exc)
        return
    st.session_state["pdb_bytes"] = data
    st.session_state["pdb_name"] = f"{pdb_id}.pdb"
    st.session_state["pdb_source_id"] = pdb_id
    st.session_state["pdb_fetch_error"] = ""


def _show_mapping(gro: Path) -> None:
    lines = gro.read_text().splitlines()
    bead_count = int(lines[1])
    rows = [
        {
            "Bead": i + 1,
            "Residue": row[5:10].strip(),
            "Name": row[10:15].strip(),
            "x (nm)": row[20:28].strip(),
            "y (nm)": row[28:36].strip(),
            "z (nm)": row[36:44].strip(),
        }
        for i, row in enumerate(lines[2:2 + bead_count])
    ]
    st.dataframe(rows, hide_index=True, use_container_width=True)


def _show_molecule(smiles: str) -> None:
    from rdkit import Chem
    from rdkit.Chem import Draw

    molecule = Chem.MolFromSmiles(smiles)
    if molecule:
        st.image(Draw.MolToImage(molecule, size=(460, 280)), caption="2D structure from SMILES")


def _solute_generator() -> None:
    with st.expander("Free molecule generator · SMILES", expanded=True):
        st.caption("Add freely dissolved Martini 3 molecules. No linker, anchor or orientation logic is used.")
        left, middle, right = st.columns([2.0, 0.8, 0.65])
        smiles = left.text_input("SMILES", key="new_smiles", placeholder="CCO")
        name = middle.text_input("Molecule name", key="new_name", placeholder="ETOH")
        copies = right.number_input("Copies", min_value=1, max_value=1000, value=1, key="new_copies")
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
                st.session_state.solutes.append(
                    {"name": species.name, "smiles": species.smiles, "count": species.count,
                     "gro": str(gro), "itp": str(itp)}
                )
                st.success(f"{species.name}: Martini Mapper generated the coarse-grained model.")
            except (ValueError, RuntimeError, OSError) as exc:
                st.error(str(exc))

        for index, row in enumerate(list(st.session_state.solutes)):
            with st.container(border=True):
                title, action = st.columns([3, 0.6])
                title.markdown(f"**{row['name']}** · {row['count']} copies · `{row['smiles']}`")
                if action.button("Remove", key=f"remove_{index}", use_container_width=True):
                    st.session_state.solutes.pop(index)
                    st.rerun()
                gro = Path(row["gro"])
                if gro.is_file():
                    with st.expander(f"Preview {row['name']} CG beads"):
                        _show_mapping(gro)


def _chain_labels() -> list[str]:
    pdb = st.session_state.get("pdb_bytes")
    if not pdb:
        return ["A"]
    _, rows = _protein_summary(pdb)
    return [str(row["Chain"]) for row in rows] or ["A"]


def _config() -> BuildConfig:
    solvent = "reline" if st.session_state.solvent_ui == "Reline · ChCl:urea 1:2" else "water"
    labels = _chain_labels()
    count = max(1, min(int(st.session_state.merge_chain_count), len(labels)))
    merge_chains = ",".join(labels[:count])
    extras = shlex.split(str(st.session_state.get("martinize_extra", "")))
    return BuildConfig(
        solvent=solvent,
        salt_m=float(st.session_state.salt) if solvent == "water" else 0.0,
        water_fraction=float(st.session_state.water_fraction) if solvent == "reline" else 0.0,
        box_distance_nm=float(st.session_state.box_distance),
        molecule_name=str(st.session_state.molecule_name).strip() or "Protein",
        merge_chains=merge_chains,
        dssp=bool(st.session_state.dssp),
        go=bool(st.session_state.go),
        go_eps=float(st.session_state.go_eps),
        elastic=bool(st.session_state.elastic),
        elastic_force=int(st.session_state.elastic_force),
        position_restraints=str(st.session_state.position_restraints),
        maxwarn=int(st.session_state.maxwarn),
        martinize_extra_args=extras,
        solutes=[Solute(row["name"], row["smiles"], int(row["count"])) for row in st.session_state.solutes],
    )


def _short_md_config() -> ShortMDConfig:
    return ShortMDConfig(
        output_tag=str(st.session_state.short_md_output_tag),
        xtc_write_every_ps=float(st.session_state.short_md_xtc_write_every_ps),
        grompp_maxwarn=int(st.session_state.short_md_grompp_maxwarn),
        threads=int(st.session_state.short_md_threads),
        stages=tuple(
            StageSettings(
                stage,
                bool(st.session_state[f"short_md_run_{stage}"]),
                float(st.session_state[f"short_md_{stage}_dt"]),
                float(st.session_state[f"short_md_{stage}_time"]),
            )
            for stage in STAGE_ORDER
        ),
    )


def _store_short_md_result(result) -> None:
    st.session_state["short_md_work_dir"] = str(result.work_dir)
    st.session_state["short_md_log"] = str(result.log)
    st.session_state["short_md_stage_results"] = [
        {
            "name": stage.name,
            "gro": str(stage.gro),
            "tpr": str(stage.tpr),
            "xtc": str(stage.xtc or ""),
            "elapsed_s": stage.elapsed_s,
            "ns_day": stage.ns_day,
        }
        for stage in result.stages
    ]
    available = [row["name"] for row in st.session_state.short_md_stage_results if row["xtc"]]
    if available:
        st.session_state.short_md_view_stage = "production" if "production" in available else available[-1]


def _reset_build() -> None:
    result = st.session_state.get("result_dir")
    if result:
        path = Path(str(result))
        if path.exists():
            shutil.rmtree(path, ignore_errors=True)
    for key in [
        "result_dir", "result_manifest", "short_md_work_dir", "short_md_log",
        "short_md_stage_results",
    ]:
        st.session_state.pop(key, None)


with st.sidebar:
    st.markdown('<div class="ms-side-title">MartiniSolv</div>', unsafe_allow_html=True)
    st.caption("Protein in solution · Martini 3")
    st.radio("Workflow step", STEPS, key="active_step", label_visibility="collapsed")
    st.markdown('<div class="ms-side-flow">Structure → Model → Environment → Review & Build → Short MD</div>', unsafe_allow_html=True)

step = st.session_state.active_step

if step == "Structure":
    left, right = st.columns([0.9, 1.1], gap="large")
    with left:
        st.markdown('<div class="ms-panel-title">Structure input</div>', unsafe_allow_html=True)
        st.text_input("Project name", key="project_name")
        uploaded = st.file_uploader("Protein structure", type=["pdb"], key="structure_upload")
        if uploaded is not None:
            st.session_state.pdb_bytes = uploaded.getvalue()
            st.session_state.pdb_name = uploaded.name
            st.session_state.pdb_source_id = ""
        if "_pdb_id_input" not in st.session_state:
            st.session_state["_pdb_id_input"] = str(st.session_state.get("pdb_id", "1UBQ"))
        st.text_input(
            "PDB ID",
            key="_pdb_id_input",
            on_change=_sync_pdb_id_from_widget,
            help="Enter a four-character RCSB PDB ID, for example 1UBQ.",
        )
        # Keep the persistent value synchronized on the first render too.
        _sync_pdb_id_from_widget()
        with st.expander("Optional auxiliary files"):
            st.caption("No surface/linker files are needed in MartiniSolv.")

    if uploaded is None:
        _ensure_remote_pdb()

    with right:
        st.markdown('<div class="ms-panel-title">Structure preview</div>', unsafe_allow_html=True)
        pdb = st.session_state.get("pdb_bytes")
        if pdb:
            atoms, chains = _protein_summary(pdb)
            source = st.session_state.get("pdb_source_id") or st.session_state.get("pdb_name", "protein.pdb")
            a, b = st.columns(2)
            a.metric("Source", source)
            b.metric("Chains", len(chains))
            render_structure_preview(pdb.decode("utf-8", "replace"), height=430)
            if chains:
                st.dataframe(chains, hide_index=True, use_container_width=True)
            if not atoms:
                st.error("No protein ATOM records were found.")
        else:
            st.markdown('<div class="ms-empty-preview"><strong>Preview pending</strong><span>Upload a PDB or enter an RCSB PDB ID.</span></div>', unsafe_allow_html=True)
        if st.session_state.get("pdb_fetch_error"):
            st.warning(f"Remote structure could not be loaded: {st.session_state.pdb_fetch_error}")

elif step == "Model":
    st.markdown('<div class="ms-panel-title">Model decisions</div>', unsafe_allow_html=True)
    labels = _chain_labels()
    st.text_input("Molecule name", key="molecule_name")
    max_count = max(1, len(labels))
    if int(st.session_state.merge_chain_count) > max_count:
        st.session_state.merge_chain_count = max_count
    options = list(range(1, max_count + 1))
    st.selectbox(
        "Merge chains",
        options,
        key="merge_chain_count",
        format_func=lambda value: f"{value} chain" if value == 1 else f"{value} chains",
    )
    c1, c2 = st.columns(2)
    c1.toggle("DSSP", key="dssp", help="Use DSSP secondary-structure assignment in martinize2.")
    c2.toggle("GōMartini", key="go", help="Enable the martinize2 Gō model.")
    st.number_input("Go_Epsilon", min_value=0.0, step=0.001, format="%.3f", key="go_eps", disabled=not st.session_state.go)
    st.segmented_control(
        "Position restraints",
        ["backbone", "all", "none"],
        key="position_restraints",
        default=st.session_state.position_restraints,
    )
    with st.expander("Advanced model options"):
        a, b = st.columns(2)
        a.toggle("Elastic network", key="elastic")
        b.number_input(
            "Elastic force constant (kJ/mol/nm²)",
            min_value=100,
            max_value=1500,
            step=50,
            key="elastic_force",
            disabled=not st.session_state.elastic,
        )
        st.number_input("martinize2 max warnings", min_value=0, max_value=20, step=1, key="maxwarn")
        st.text_area(
            "Extra martinize2 args",
            key="martinize_extra",
            height=80,
            help="Advanced users only. Arguments are parsed safely without a shell.",
        )

elif step == "Environment":
    st.markdown('<div class="ms-panel-title">Environment</div>', unsafe_allow_html=True)
    a, b = st.columns(2)
    a.selectbox("Solvent", ["Water", "Reline · ChCl:urea 1:2"], key="solvent_ui")
    b.number_input("Protein-to-box distance (nm)", min_value=0.5, max_value=5.0, step=0.1, key="box_distance")
    if st.session_state.solvent_ui == "Water":
        st.number_input("NaCl concentration (M)", min_value=0.0, max_value=2.0, step=0.05, key="salt")
        st.caption("INSANE solvates the coarse-grained protein with Martini 3 water and ions.")
    else:
        st.slider("Water mole fraction in reline", min_value=0.0, max_value=0.4, step=0.01, key="water_fraction")
        st.caption("Reline uses ChCl:urea 1:2. Added NaCl is not enabled for this path.")
    _solute_generator()

elif step == "Review & Build":
    st.markdown('<div class="ms-panel-title">Review & Build</div>', unsafe_allow_html=True)
    try:
        config = _config()
        config.validate()
        pdb = st.session_state.get("pdb_bytes")
        valid_pdb = bool(pdb and _protein_summary(pdb)[0])

        summary_a, summary_b, summary_c, summary_d = st.columns(4)
        summary_a.metric("Protein", config.molecule_name)
        summary_b.metric("Model", "GōMartini" if config.go else ("Elastic" if config.elastic else "Martini 3"))
        summary_c.metric("Solvent", "Water" if config.solvent == "water" else "Reline")
        summary_d.metric("Free species", len(config.solutes))

        left, right = st.columns([1.15, 0.85], gap="large")
        with right:
            build_clicked = st.button("Build system", type="primary", use_container_width=True, disabled=not valid_pdb)
            if st.button("Reset working folder", use_container_width=True):
                _reset_build()
                st.rerun()
            if not valid_pdb:
                st.warning("Provide a valid protein structure in Structure first.")

        if build_clicked:
            root = Path(tempfile.mkdtemp(prefix="martinisolv_"))
            input_path = root / "input.pdb"
            input_path.write_bytes(pdb)
            with st.status("Building Martini 3 system", expanded=True) as status:
                try:
                    built = build(input_path, root / "system", config)
                    st.session_state.result_dir = str(built)
                    st.session_state.result_manifest = (built / "manifest.json").read_text()
                    status.update(label="GROMACS-ready system generated", state="complete")
                except Exception as exc:
                    status.update(label="Build failed", state="error")
                    st.error(str(exc))
                    log = root / "system.incomplete" / "build.log"
                    if log.is_file():
                        with st.expander("Build execution log", expanded=True):
                            st.code(log.read_text(errors="replace")[-10000:], language="text")

        built = Path(str(st.session_state.get("result_dir", "")))
        with left:
            st.markdown('<div class="ms-panel-title">Generated structure</div>', unsafe_allow_html=True)
            if built.is_dir():
                protein_cg = built / "protein_cg.pdb"
                system_gro = built / "system.gro"
                view_options = []
                if protein_cg.is_file():
                    view_options.append("Martini protein")
                if system_gro.is_file():
                    view_options.append("Full solvated system")
                if view_options:
                    selected = st.selectbox("Structure to view", view_options, label_visibility="collapsed")
                    render_cg_structure(protein_cg if selected == "Martini protein" else system_gro, height=620)
                else:
                    st.info("No generated structure is available for preview.")
            else:
                st.markdown('<div class="ms-empty-preview"><strong>No generated files yet.</strong><span>Build the system to inspect the Martini structure here.</span></div>', unsafe_allow_html=True)

        if built.is_dir() and (built / "manifest.json").is_file():
            st.success("System generated and checked with GROMACS.")
            with st.expander("Build details"):
                manifest = json.loads((built / "manifest.json").read_text())
                rows = [
                    {"Setting": "Protein model", "Value": "GōMartini" if config.go else "Martini 3"},
                    {"Setting": "DSSP", "Value": "On" if config.dssp else "Off"},
                    {"Setting": "Position restraints", "Value": config.position_restraints},
                    {"Setting": "Solvent", "Value": config.solvent},
                    {"Setting": "Box distance", "Value": f"{config.box_distance_nm:.2f} nm"},
                    {"Setting": "Salt", "Value": f"{config.salt_m:.2f} M" if config.solvent == "water" else "—"},
                    {"Setting": "Build check", "Value": manifest.get("status", "")},
                ]
                st.dataframe(rows, hide_index=True, use_container_width=True)
            st.download_button(
                "Download GROMACS package",
                archive(built),
                file_name="martini-solv-system.zip",
                mime="application/zip",
                type="primary",
                use_container_width=True,
            )
    except (ValueError, RuntimeError) as exc:
        st.error(str(exc))

elif step == "Short MD":
    st.markdown('<div class="ms-panel-title">Short MDs</div>', unsafe_allow_html=True)
    st.caption("Optional GROMACS fast check: minimization → NVT → NPT → Production.")

    built = Path(str(st.session_state.get("result_dir", "")))
    if not (built.is_dir() and (built / "system.gro").is_file()):
        st.info("Build the system in Review & Build first.")
    else:
        top_a, top_b, top_c = st.columns(3)
        top_a.text_input("Simulation name", key="short_md_output_tag")
        top_b.number_input(
            "XTC write every (ps)", min_value=0.001, step=1.0, format="%.4f",
            key="short_md_xtc_write_every_ps",
        )
        top_c.number_input("GROMPP max warnings", min_value=0, step=1, key="short_md_grompp_maxwarn")

        with st.expander("Short MD protocol", expanded=True):
            for stage in STAGE_ORDER:
                defaults = DEFAULT_STAGE_SETTINGS[stage]
                cols = st.columns([0.85, 1, 1])
                cols[0].toggle(stage.upper() if stage != "production" else "Production", key=f"short_md_run_{stage}")
                cols[1].number_input(
                    f"{stage.upper() if stage != 'production' else 'Production'} timestep (ps)",
                    min_value=0.000001,
                    step=float(defaults["dt_ps"]),
                    format="%.4f",
                    key=f"short_md_{stage}_dt",
                )
                cols[2].number_input(
                    f"{stage.upper() if stage != 'production' else 'Production'} time (ns)",
                    min_value=0.000001,
                    step=float(defaults["time_ns"]),
                    format="%.4f",
                    key=f"short_md_{stage}_time",
                )
            st.number_input("CPU threads", min_value=1, max_value=4, step=1, key="short_md_threads")

        md_config = _short_md_config()
        errors = validate_stage_order(selected_stages(md_config))
        for error in errors:
            st.error(error)

        if st.button("Run Short MD", type="primary", use_container_width=True, disabled=bool(errors)):
            with st.status("Running Short MD", expanded=True) as status:
                try:
                    result = run_short_md(built, md_config)
                    _store_short_md_result(result)
                    status.update(label="Short MD completed", state="complete")
                except Exception as exc:
                    status.update(label="Short MD failed", state="error")
                    st.error(str(exc))

        rows = st.session_state.get("short_md_stage_results") or []
        if rows:
            st.success("Short MD finished successfully.")
            table_rows = []
            for row in rows:
                table_rows.append({
                    "Stage": str(row["name"]).upper(),
                    "Wall time (s)": round(float(row["elapsed_s"]), 2),
                    "Performance (ns/day)": round(float(row["ns_day"]), 3) if row["ns_day"] is not None else "—",
                })
            st.dataframe(table_rows, hide_index=True, use_container_width=True)

            trajectory_rows = [row for row in rows if row.get("xtc")]
            if trajectory_rows:
                stage_names = [str(row["name"]) for row in trajectory_rows]
                if st.session_state.short_md_view_stage not in stage_names:
                    st.session_state.short_md_view_stage = stage_names[-1]
                controls = st.columns([1.2, 0.8])
                controls[0].selectbox(
                    "Stage to view",
                    stage_names,
                    key="short_md_view_stage",
                    format_func=lambda value: value.upper() if value != "production" else "Production",
                )
                controls[1].number_input("Frame stride", min_value=1, step=1, key="short_md_view_stride")

                selected = next(row for row in trajectory_rows if row["name"] == st.session_state.short_md_view_stage)
                gro = Path(selected["gro"])
                xtc = Path(selected["xtc"])
                if gro.is_file() and xtc.is_file():
                    try:
                        frames = render_trajectory(gro, xtc, stride=int(st.session_state.short_md_view_stride), height=680)
                        st.caption(f"Animated trajectory preview · {frames} displayed frames.")
                    except Exception as exc:
                        st.warning(f"Trajectory preview is unavailable: {exc}")

            st.download_button(
                "Download system and Short MD files",
                archive(built),
                file_name="martini-solv-system.zip",
                mime="application/zip",
                use_container_width=True,
            )
            log_path = Path(str(st.session_state.get("short_md_log", "")))
            if log_path.is_file():
                with st.expander("Short MD execution log"):
                    st.code(log_path.read_text(errors="replace")[-12000:], language="text")
