"""Friendly MartiniSolv Streamlit workflow modelled on MartiniSurf."""
from __future__ import annotations

import json
from pathlib import Path
import shlex
import shutil
import tempfile
import urllib.request

import streamlit as st

from martini_solv.builder import _itp_net_charge, _map_solute, _prepare_uploaded_solute, build
from martini_solv.models import BuildConfig, Solute
from martini_solv.molecular_viewer import render_build_viewer, render_free_molecule_mapping, render_structure_preview, render_trajectory
from martini_solv.short_md import (
    DEFAULT_GROMPP_MAXWARN,
    DEFAULT_STAGE_SETTINGS,
    DEFAULT_XTC_WRITE_EVERY_PS,
    STAGE_ORDER,
    ShortMDConfig,
    StageSettings,
    archive,
    run_short_md,
    run_short_md_analysis,
    selected_stages,
    validate_stage_order,
)
from martini_solv.theme import STYLE


STEPS = ["Home", "Structure", "Model", "Environment", "Review & Build", "Short MD"]
CHAIN_LABELS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
APP_STATE_VERSION = 3

st.set_page_config(page_title="MartiniSolv", page_icon="MS", layout="wide", initial_sidebar_state="expanded")
st.markdown(STYLE, unsafe_allow_html=True)

DEFAULTS = {
    "active_step": "Home",
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
    "new_smiles": "",
    "new_name": "",
    "new_copies": 1,
    "upload_name": "",
    "upload_copies": 1,
    "build_view": "Full solvated system",
    "viewer_show_connectivity": True,
    "viewer_bead_radius": 0.85,
    "viewer_bond_radius": 0.20,
    "viewer_topology_bond_max_nm": 0.75,
    "viewer_show_protein": True,
    "viewer_show_solute": True,
    "viewer_show_solvent": True,
    "viewer_show_ions": True,
    "short_md_output_tag": "Protein MD",
    "short_md_xtc_write_every_ps": DEFAULT_XTC_WRITE_EVERY_PS,
    "short_md_grompp_maxwarn": DEFAULT_GROMPP_MAXWARN,
    "short_md_threads": 2,
    "short_md_view_stage": "production",
    "short_md_view_stride": 1,
    "short_md_view_protein": True,
    "short_md_view_solute": False,
    "short_md_view_solvent": False,
    "short_md_view_ions": False,
}
if st.session_state.get("_martinisolv_state_version") != APP_STATE_VERSION:
    # One-time migration: make the requested defaults visible even in browser
    # sessions created before this release.
    for key, value in DEFAULTS.items():
        if key != "active_step":
            st.session_state[key] = value
    st.session_state["_martinisolv_state_version"] = APP_STATE_VERSION
else:
    for key, value in DEFAULTS.items():
        st.session_state.setdefault(key, value)

for stage, defaults in DEFAULT_STAGE_SETTINGS.items():
    st.session_state.setdefault(f"short_md_run_{stage}", bool(defaults["enabled"]))
    st.session_state.setdefault(f"short_md_{stage}_dt", float(defaults["dt_ps"]))
    st.session_state.setdefault(f"short_md_{stage}_time", float(defaults["time_ns"]))


BUILD_INPUT_KEYS = {
    "project_name", "pdb_id", "molecule_name", "merge_chain_count", "dssp", "go",
    "go_eps", "position_restraints", "elastic", "elastic_force", "maxwarn",
    "martinize_extra", "solvent_ui", "box_distance", "salt", "water_fraction",
    "new_smiles", "new_name", "new_copies",
}


def _invalidate_build_outputs() -> None:
    """Prevent stale water/Reline systems from being shown after configuration changes."""
    for key in [
        "result_dir", "result_manifest", "short_md_work_dir", "short_md_log",
        "short_md_stage_results",
    ]:
        st.session_state.pop(key, None)


def _widget_key(name: str) -> str:
    return f"_widget_{name}"


def _prime_widget(name: str) -> str:
    """Restore a widget from persistent workflow state when its page reappears."""
    key = _widget_key(name)
    if key not in st.session_state:
        st.session_state[key] = st.session_state.get(name)
    return key


def _store_widget(name: str) -> None:
    """Persist widget state and invalidate any build made from older settings."""
    previous = st.session_state.get(name)
    current = st.session_state.get(_widget_key(name))
    st.session_state[name] = current
    if name in BUILD_INPUT_KEYS and previous != current:
        _invalidate_build_outputs()


def _store_pdb_id() -> None:
    previous = st.session_state.get("pdb_id")
    current = str(st.session_state.get(_widget_key("pdb_id"), "")).strip().upper()
    st.session_state["pdb_id"] = current
    if previous != current:
        _invalidate_build_outputs()


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


def _itp_bead_types(path: Path) -> dict[int, str]:
    """Read Martini bead types from an ITP [ atoms ] section."""
    rows: dict[int, str] = {}
    if not path.is_file():
        return rows
    section = ""
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.startswith("[") and "]" in line:
            section = line.strip("[]").strip().lower()
            continue
        if section != "atoms" or not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            rows[int(parts[0])] = parts[1]
        except ValueError:
            continue
    return rows


def _solute_charge_state(row: dict[str, object]) -> tuple[int, int | None]:
    charge = int(row.get("net_charge", 0) or 0)
    bead = row.get("charged_bead")
    return charge, (int(bead) if bead is not None else None)


def _itp_bead_charges(path: Path) -> dict[int, float]:
    """Read per-bead charges from an ITP [ atoms ] section."""
    rows: dict[int, float] = {}
    if not path.is_file():
        return rows
    section = ""
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.startswith("[") and "]" in line:
            section = line.strip("[]").strip().lower()
            continue
        if section != "atoms" or not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        try:
            rows[int(parts[0])] = float(parts[6])
        except ValueError:
            continue
    return rows


def _solute_generator() -> None:
    with st.expander("Free molecule generator · SMILES", expanded=True):
        st.caption("Generate a freely dissolved Martini 3 molecule from a neutral SMILES.")
        left, middle, right = st.columns([2.0, 0.8, 0.65])
        smiles = left.text_input(
            "SMILES", key=_prime_widget("new_smiles"), placeholder="CCO",
            on_change=_store_widget, args=("new_smiles",),
        )
        name = middle.text_input(
            "Molecule name", key=_prime_widget("new_name"), placeholder="ETOH",
            on_change=_store_widget, args=("new_name",),
        )
        copies = right.number_input(
            "Copies", min_value=1, max_value=1000,
            key=_prime_widget("new_copies"),
            on_change=_store_widget, args=("new_copies",),
        )
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
                    {
                        "name": species.name,
                        "smiles": species.smiles,
                        "count": species.count,
                        "net_charge": 0,
                        "charged_bead": None,
                        "source": "smiles",
                        "template_itp": "",
                        "gro": str(gro),
                        "itp": str(itp),
                    }
                )
                _invalidate_build_outputs()
                st.success(f"{species.name}: Martini Mapper generated the coarse-grained model.")
            except (ValueError, RuntimeError, OSError) as exc:
                st.error(str(exc))

    with st.expander("Free molecule · Upload Martini topology", expanded=False):
        st.caption(
            "Already have a Martini model? Upload only its .itp topology. "
            "MartiniSolv reconstructs a starting .gro geometry from [ atoms ] plus "
            "[ bonds ]/[ constraints ], then inserts the requested number of copies."
        )
        up_left, up_mid, up_right = st.columns([1.35, 1.0, 0.65])
        uploaded_itp = up_left.file_uploader(
            "Martini topology (.itp)",
            type=["itp"],
            key="free_molecule_itp_upload",
        )
        upload_name = up_mid.text_input(
            "Molecule name",
            key=_prime_widget("upload_name"),
            placeholder="LIG",
            on_change=_store_widget,
            args=("upload_name",),
        )
        upload_copies = up_right.number_input(
            "Copies",
            min_value=1,
            max_value=1000,
            key=_prime_widget("upload_copies"),
            on_change=_store_widget,
            args=("upload_copies",),
        )

        source_itp = None
        if uploaded_itp is not None:
            try:
                root = Path(tempfile.mkdtemp(prefix="martinisolv_upload_"))
                source_itp = root / Path(uploaded_itp.name).name
                source_itp.write_bytes(uploaded_itp.getvalue())
                bead_types = _itp_bead_types(source_itp)
                topology_charge = _itp_net_charge(source_itp)
                info_a, info_b = st.columns(2)
                info_a.metric("Detected beads", len(bead_types))
                info_b.metric("Topology charge", f"{topology_charge:+.2f} e")
                if abs(topology_charge - round(topology_charge)) > 0.01:
                    st.warning(
                        "The uploaded topology has a non-integer net charge. "
                        "Automatic counterion neutralization requires an approximately integer molecular charge."
                    )
            except Exception as exc:
                st.warning(f"Could not inspect uploaded ITP: {exc}")
                source_itp = None

        if st.button("Add uploaded molecule", use_container_width=True):
            try:
                if source_itp is None:
                    raise ValueError("Upload a Martini .itp file first")
                species = Solute(
                    name=str(upload_name).strip(),
                    smiles="",
                    count=int(upload_copies),
                    source="upload",
                    template_itp=str(source_itp),
                )
                species.validate()
                if any(row["name"].upper() == species.name.upper() for row in st.session_state.solutes):
                    raise ValueError("Choose a unique molecule name")

                preview_root = Path(tempfile.mkdtemp(prefix="martinisolv_uploaded_preview_"))
                gro, itp = _prepare_uploaded_solute(
                    preview_root, species.name, source_itp
                )
                st.session_state.solutes.append(
                    {
                        "name": species.name,
                        "smiles": "",
                        "count": species.count,
                        "net_charge": 0,
                        "charged_bead": None,
                        "source": "upload",
                        "template_itp": str(source_itp),
                        "gro": str(gro),
                        "itp": str(itp),
                    }
                )
                _invalidate_build_outputs()
                st.success(
                    f"{species.name}: topology added and a starting GRO template was generated automatically."
                )
            except (ValueError, RuntimeError, OSError) as exc:
                st.error(str(exc))

        st.caption(
            "ITP-only reconstruction is intended for small connected Martini molecules. "
            "If the topology does not contain enough connectivity to place the beads safely, "
            "MartiniSolv stops and asks for a more complete topology instead of guessing."
        )

    for index, row in enumerate(list(st.session_state.solutes)):
        with st.container(border=True):
            title, action = st.columns([3, 0.6])
            source = str(row.get("source", "smiles"))
            preview_itp = Path(row["itp"])
            if source == "upload":
                actual_charge = _itp_net_charge(preview_itp) if preview_itp.is_file() else 0.0
                source_label = f"Uploaded ITP · topology charge {actual_charge:+.2f} e"
                charge_text = f"{actual_charge:+.2f}"
            else:
                row_charge, _ = _solute_charge_state(row)
                charge_text = f"{row_charge:+d}" if row_charge else "0"
                source_label = f"`{row['smiles']}`"

            title.markdown(
                f"**{row['name']}** · {row['count']} copies · charge {charge_text} · {source_label}"
            )
            if action.button("Remove", key=f"remove_{index}", use_container_width=True):
                st.session_state.solutes.pop(index)
                _invalidate_build_outputs()
                st.rerun()

            gro = Path(row["gro"])
            itp = Path(row["itp"])
            if not gro.is_file():
                continue

            with st.expander(f"Martini topology preview · {row['name']}", expanded=True):
                preview_col, beads_col = st.columns([1.35, 1], gap="large")
                with preview_col:
                    bead_rows = render_free_molecule_mapping(gro, height=360)
                    if source == "upload":
                        st.caption(
                            "Starting geometry reconstructed automatically from the uploaded ITP. "
                            "The build minimization will relax this template."
                        )
                    else:
                        st.caption(
                            "Generated coarse-grained topology. Bead labels start at 1, "
                            "matching the topology table."
                        )

                with beads_col:
                    bead_types = _itp_bead_types(itp)
                    bead_charges = _itp_bead_charges(itp)
                    bead_options = sorted(bead_types)

                    if source == "upload":
                        chosen_charge = None
                        chosen_bead = None
                        st.markdown("##### Topology charge")
                        st.metric("Net charge from uploaded ITP", f"{_itp_net_charge(itp):+.2f} e")
                        st.caption(
                            "Charges already present in the uploaded topology are preserved and "
                            "used automatically when MartiniSolv calculates Na⁺/Cl⁻ counterions."
                        )
                    else:
                        st.markdown("##### Charge assignment")
                        charge_key = f"solute_charge_{index}"
                        if charge_key not in st.session_state:
                            st.session_state[charge_key] = int(row.get("net_charge", 0) or 0)
                        chosen_charge = st.select_slider(
                            "Net charge",
                            options=[-2, -1, 0, 1, 2],
                            key=charge_key,
                            help="Applied after neutral-SMILES mapping. Default is 0.",
                        )
                        chosen_bead = None
                        if chosen_charge != 0 and bead_options:
                            bead_key = f"solute_charged_bead_{index}"
                            existing_bead = row.get("charged_bead")
                            if bead_key not in st.session_state:
                                st.session_state[bead_key] = (
                                    int(existing_bead) if existing_bead in bead_options else bead_options[0]
                                )
                            chosen_bead = st.selectbox(
                                "Charged bead",
                                bead_options,
                                key=bead_key,
                                format_func=lambda value: f"{value}: {bead_types.get(value, 'bead')}",
                                help="The selected bead receives the integer molecular charge in the generated ITP.",
                            )
                            st.caption(
                                "Manual charge assignment changes the ITP charge column only; "
                                "the Martini bead type assigned by Martini Mapper is preserved."
                            )
                        if (
                            int(row.get("net_charge", 0) or 0) != int(chosen_charge)
                            or row.get("charged_bead") != chosen_bead
                        ):
                            st.session_state.solutes[index]["net_charge"] = int(chosen_charge)
                            st.session_state.solutes[index]["charged_bead"] = chosen_bead
                            _invalidate_build_outputs()

                    table = []
                    for bead_index, bead in enumerate(bead_rows, start=1):
                        if source == "upload":
                            displayed_charge = bead_charges.get(bead_index, 0.0)
                        else:
                            displayed_charge = (
                                int(chosen_charge)
                                if chosen_charge != 0 and bead_index == chosen_bead
                                else 0
                            )
                        table.append({
                            "Bead": bead["Bead"],
                            "Martini type": bead_types.get(bead_index, "—"),
                            "Charge": displayed_charge,
                            "Residue": bead["Residue"],
                            "x (nm)": bead["x (nm)"],
                            "y (nm)": bead["y (nm)"],
                            "z (nm)": bead["z (nm)"],
                        })

                    if table:
                        st.markdown("##### Bead mapping")
                        st.table(table)
                    elif bead_types:
                        st.markdown("##### Bead mapping")
                        st.table([
                            {
                                "Bead": bead_index,
                                "Martini type": bead_type,
                                "Charge": bead_charges.get(bead_index, 0.0)
                                if source == "upload"
                                else (
                                    int(chosen_charge)
                                    if chosen_charge != 0 and bead_index == chosen_bead
                                    else 0
                                ),
                            }
                            for bead_index, bead_type in sorted(bead_types.items())
                        ])
                    st.caption(
                        "This Martini topology will be inserted as a free molecule "
                        "in the simulation box."
                    )

def _chain_labels() -> list[str]:
    pdb = st.session_state.get("pdb_bytes")
    if not pdb:
        return ["A"]
    _, rows = _protein_summary(pdb)
    return [str(row["Chain"]) for row in rows] or ["A"]


def _config() -> BuildConfig:
    solvent_label = str(st.session_state.solvent_ui)
    if solvent_label == "Reline · ChCl:urea 1:2":
        solvent = "reline"
    elif solvent_label == "ChCl:sorbitol 1:1":
        solvent = "chcl_sorbitol"
    else:
        solvent = "water"
    labels = _chain_labels()
    count = max(1, min(int(st.session_state.merge_chain_count), len(labels)))
    merge_chains = ",".join(labels[:count])
    extras = shlex.split(str(st.session_state.get("martinize_extra", "")))
    return BuildConfig(
        solvent=solvent,
        salt_m=float(st.session_state.salt) if solvent == "water" else 0.0,
        water_fraction=float(st.session_state.water_fraction) if solvent != "water" else 0.0,
        chcl_sorbitol_density_g_cm3=1.20,
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
        solutes=[
            Solute(
                row["name"],
                row.get("smiles", ""),
                int(row["count"]),
                int(row.get("net_charge", 0) or 0),
                (int(row["charged_bead"]) if row.get("charged_bead") is not None else None),
                str(row.get("source", "smiles")),
                str(row.get("template_itp", "")),
            )
            for row in st.session_state.solutes
        ],
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
            "edr": str(stage.edr or ""),
            "elapsed_s": stage.elapsed_s,
            "ns_day": stage.ns_day,
        }
        for stage in result.stages
    ]
    available = [row["name"] for row in st.session_state.short_md_stage_results if row["xtc"]]
    if available:
        st.session_state.short_md_view_stage = "production" if "production" in available else available[-1]


def _render_stage_analysis(selected: dict[str, object]) -> None:
    """Render analyses for exactly the stage selected in 'Stage to view'."""
    stage = str(selected["name"])
    signature = (
        stage,
        str(selected.get("tpr") or ""),
        str(selected.get("xtc") or ""),
        str(selected.get("edr") or ""),
    )
    st.markdown("#### Analysis")
    st.caption(
        f"Analyses use the complete **{stage.upper()}** trajectory selected above. "
        "RMSD/RMSF use protein BB beads; Density is the total simulation-box mass density."
    )
    a, b, d = st.columns(3)
    requested = None
    if a.button("Analyse RMSD", key=f"analysis_rmsd_{stage}", use_container_width=True):
        requested = "rmsd"
    if b.button("Analyse RMSF", key=f"analysis_rmsf_{stage}", use_container_width=True):
        requested = "rmsf"
    if d.button("Analyse Density", key=f"analysis_density_{stage}", use_container_width=True):
        requested = "density"

    if requested:
        tpr = Path(str(selected.get("tpr") or ""))
        xtc_text = str(selected.get("xtc") or "")
        edr_text = str(selected.get("edr") or "")
        xtc = Path(xtc_text) if xtc_text else None
        edr = Path(edr_text) if edr_text else None
        work = Path(str(st.session_state.get("short_md_work_dir", "."))) / "analysis" / stage
        with st.spinner(f"Calculating {requested.upper()} for {stage.upper()}…"):
            try:
                result = run_short_md_analysis(requested, tpr, xtc, edr, work)
            except Exception as exc:
                st.session_state["short_md_analysis_signature"] = signature
                st.session_state["short_md_analysis_error"] = str(exc)
                st.session_state["short_md_analysis_result"] = {}
            else:
                values = [value for _, value in result.points]
                st.session_state["short_md_analysis_signature"] = signature
                st.session_state["short_md_analysis_error"] = ""
                st.session_state["short_md_analysis_result"] = {
                    "kind": result.kind,
                    "x_label": result.x_label,
                    "y_label": result.y_label,
                    "points": result.points,
                    "mean": sum(values) / len(values),
                    "maximum": max(values),
                    "minimum": min(values),
                    "output_path": str(result.output_path),
                }

    if st.session_state.get("short_md_analysis_signature") != signature:
        return
    error = str(st.session_state.get("short_md_analysis_error") or "")
    if error:
        st.error(error)
        return
    result = st.session_state.get("short_md_analysis_result") or {}
    if not result:
        return

    kind = str(result["kind"]).upper()
    metrics = st.columns(3)
    metrics[0].metric("Analysis", kind)
    if kind == "DENSITY":
        metrics[1].metric("Mean density", f"{float(result['mean']):.1f} kg/m³")
        metrics[2].metric("Range", f"{float(result['minimum']):.1f}–{float(result['maximum']):.1f}")
    else:
        metrics[1].metric("Mean", f"{float(result['mean']):.3f} nm")
        metrics[2].metric("Maximum", f"{float(result['maximum']):.3f} nm")

    rows = [
        {str(result["x_label"]): x, str(result["y_label"]): y}
        for x, y in result["points"]
    ]
    st.line_chart(rows, x=str(result["x_label"]), y=str(result["y_label"]), height=320)
    output = Path(str(result["output_path"]))
    if output.is_file():
        st.download_button(
            "Download XVG data",
            output.read_bytes(),
            file_name=output.name,
            mime="text/plain",
            key=f"download_analysis_{stage}_{result['kind']}",
            use_container_width=True,
        )


def _go_to_structure() -> None:
    """Navigate from Home before Streamlit instantiates the sidebar radio."""
    st.session_state["active_step"] = "Structure"


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
    st.markdown('<div class="ms-side-flow">Home · Structure → Model → Environment → Review & Build → Short MD</div>', unsafe_allow_html=True)

step = st.session_state.active_step

if step == "Home":
    st.markdown(
        """
        <section class="ms-home-hero">
          <div class="ms-home-kicker">MARTINI 3 · PROTEIN IN SOLUTION</div>
          <h1>MartiniSolv</h1>
          <p class="ms-home-lead">
            A streamlined workflow to prepare, solvate, inspect and validate
            coarse-grained protein systems with Martini 3.
          </p>
          <div class="ms-home-chips">
            <span>Structure preparation</span>
            <span>Water & DES solvents</span>
            <span>Free molecules</span>
            <span>Short MD validation</span>
            <span>Simulation-ready files</span>
          </div>
          <div class="ms-home-author">
            <div>
              <div class="ms-home-author-label">Developed by</div>
              <div class="ms-home-author-name">Juan Carlos Jiménez-García</div>
            </div>
            <a class="ms-home-github" href="https://github.com/jjimenezgar" target="_blank" rel="noopener noreferrer">
              github.com/jjimenezgar ↗
            </a>
          </div>
        </section>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("### Workflow")
    flow_a, flow_b, flow_c = st.columns(3)
    with flow_a:
        st.markdown(
            """
            <div class="ms-home-card">
              <div class="ms-home-card-number">01</div>
              <strong>Prepare the protein</strong>
              <span>Load a PDB, generate the Martini 3 model and configure GōMartini or elastic-network options.</span>
            </div>
            """,
            unsafe_allow_html=True,
        )
    with flow_b:
        st.markdown(
            """
            <div class="ms-home-card">
              <div class="ms-home-card-number">02</div>
              <strong>Build the environment</strong>
              <span>Choose water, Reline or ChCl:sorbitol and optionally add mapped free molecules.</span>
            </div>
            """,
            unsafe_allow_html=True,
        )
    with flow_c:
        st.markdown(
            """
            <div class="ms-home-card">
              <div class="ms-home-card-number">03</div>
              <strong>Validate & export</strong>
              <span>Inspect the full system, run a short MD check and download organised Simulation_Files.</span>
            </div>
            """,
            unsafe_allow_html=True,
        )

    st.markdown("")
    st.button(
        "Start a MartiniSolv project",
        type="primary",
        use_container_width=True,
        on_click=_go_to_structure,
    )

elif step == "Structure":
    left, right = st.columns([0.9, 1.1], gap="large")
    with left:
        st.markdown('<div class="ms-panel-title">Structure input</div>', unsafe_allow_html=True)
        st.text_input(
            "Project name", key=_prime_widget("project_name"),
            on_change=_store_widget, args=("project_name",),
        )
        uploaded = st.file_uploader("Protein structure", type=["pdb"], key="structure_upload")
        if uploaded is not None:
            uploaded_bytes = uploaded.getvalue()
            if st.session_state.get("pdb_bytes") != uploaded_bytes:
                _invalidate_build_outputs()
            st.session_state.pdb_bytes = uploaded_bytes
            st.session_state.pdb_name = uploaded.name
            st.session_state.pdb_source_id = ""
        st.text_input(
            "PDB ID",
            key=_prime_widget("pdb_id"),
            on_change=_store_pdb_id,
            help="Enter a four-character RCSB PDB ID, for example 1UBQ.",
        )
        _store_pdb_id()
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
    st.text_input(
        "Molecule name", key=_prime_widget("molecule_name"),
        on_change=_store_widget, args=("molecule_name",),
    )
    max_count = max(1, len(labels))
    if int(st.session_state.merge_chain_count) > max_count:
        st.session_state.merge_chain_count = max_count
    options = list(range(1, max_count + 1))
    merge_key = _prime_widget("merge_chain_count")
    if int(st.session_state.get(merge_key, 1)) not in options:
        st.session_state[merge_key] = int(st.session_state.merge_chain_count)
    st.selectbox(
        "Merge chains",
        options,
        key=merge_key,
        on_change=_store_widget,
        args=("merge_chain_count",),
        format_func=lambda value: f"{value} chain" if value == 1 else f"{value} chains",
    )
    c1, c2 = st.columns(2)
    c1.toggle(
        "DSSP", key=_prime_widget("dssp"),
        on_change=_store_widget, args=("dssp",),
        help="Use DSSP secondary-structure assignment in martinize2.",
    )
    c2.toggle(
        "GōMartini", key=_prime_widget("go"),
        on_change=_store_widget, args=("go",),
        help="Enable the martinize2 Gō model.",
    )
    st.number_input(
        "Go_Epsilon", min_value=0.0, step=0.001, format="%.3f",
        key=_prime_widget("go_eps"),
        on_change=_store_widget, args=("go_eps",),
        disabled=not bool(st.session_state.get(_widget_key("go"), st.session_state.go)),
    )
    st.segmented_control(
        "Position restraints",
        ["backbone", "all", "none"],
        key=_prime_widget("position_restraints"),
        on_change=_store_widget,
        args=("position_restraints",),
    )
    with st.expander("Advanced model options"):
        a, b = st.columns(2)
        a.toggle(
            "Elastic network", key=_prime_widget("elastic"),
            on_change=_store_widget, args=("elastic",),
        )
        b.number_input(
            "Elastic force constant (kJ/mol/nm²)",
            min_value=100,
            max_value=1500,
            step=50,
            key=_prime_widget("elastic_force"),
            on_change=_store_widget,
            args=("elastic_force",),
            disabled=not bool(st.session_state.get(_widget_key("elastic"), st.session_state.elastic)),
        )
        st.number_input(
            "martinize2 max warnings", min_value=0, max_value=20, step=1,
            key=_prime_widget("maxwarn"),
            on_change=_store_widget, args=("maxwarn",),
        )
        st.text_area(
            "Extra martinize2 args",
            key=_prime_widget("martinize_extra"),
            on_change=_store_widget, args=("martinize_extra",),
            height=80,
            help="Advanced users only. Arguments are parsed safely without a shell.",
        )

elif step == "Environment":
    st.markdown('<div class="ms-panel-title">Environment</div>', unsafe_allow_html=True)
    a, b = st.columns(2)
    a.selectbox(
        "Solvent", ["Water", "Reline · ChCl:urea 1:2", "ChCl:sorbitol 1:1"],
        key=_prime_widget("solvent_ui"),
        on_change=_store_widget, args=("solvent_ui",),
    )
    b.number_input(
        "Protein-to-box distance (nm)", min_value=0.5, max_value=5.0, step=0.1,
        key=_prime_widget("box_distance"),
        on_change=_store_widget, args=("box_distance",),
    )
    active_solvent = str(st.session_state.get(_widget_key("solvent_ui"), st.session_state.solvent_ui))
    if active_solvent == "Water":
        st.number_input(
            "NaCl concentration (M)", min_value=0.0, max_value=2.0, step=0.05,
            key=_prime_widget("salt"),
            on_change=_store_widget, args=("salt",),
        )
        st.caption("INSANE solvates the coarse-grained protein with Martini 3 water and ions.")
    else:
        wet_label = "Water mole fraction in reline" if active_solvent.startswith("Reline") else "Water mole fraction in ChCl:sorbitol"
        st.slider(
            wet_label, min_value=0.0, max_value=0.4, step=0.01,
            key=_prime_widget("water_fraction"),
            on_change=_store_widget, args=("water_fraction",),
        )
        if active_solvent.startswith("Reline"):
            st.caption(
                "Reline uses ChCl:urea 1:2. Initial packing is estimated from the experimental "
                "density 1.20 g/cm³. Added NaCl is not enabled for this path."
            )
        else:
            st.caption(
                "ChCl:sorbitol uses a 1:1 molar ratio. Sorbitol is represented by the supplied "
                "three-P4-bead topology. Initial packing uses a fixed density of 1.20 g/cm³. "
                "Added NaCl is not enabled for this path."
            )
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
        solvent_name = {
            "water": "Water",
            "reline": "Reline",
            "chcl_sorbitol": "ChCl:sorbitol 1:1",
        }[config.solvent]
        summary_c.metric("Solvent", solvent_name)
        summary_d.metric("Free molecules", sum(spec.count for spec in config.solutes))

        if config.solvent in {"reline", "chcl_sorbitol"}:
            des_name = "Reline · ChCl:urea = 1:2" if config.solvent == "reline" else "ChCl:sorbitol = 1:1"
            wet_prefix = "Wet" if config.water_fraction > 0 else "Dry"
            wet_text = (
                f"target water mole fraction = {config.water_fraction:.2f}"
                if config.water_fraction > 0 else "no water"
            )
            st.info(
                f"Environment selected: {wet_prefix} {des_name} · {wet_text}. "
                "A new build is required after changing the environment."
            )
        else:
            st.info(
                f"Environment selected: Martini water · NaCl = {config.salt_m:.2f} M. "
                "A new build is required after changing the environment."
            )

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
                    if st.session_state.build_view not in view_options:
                        st.session_state.build_view = view_options[0]
                    view_key = _prime_widget("build_view")
                    if st.session_state.get(view_key) not in view_options:
                        st.session_state[view_key] = st.session_state.build_view
                    selected = st.selectbox(
                        "Structure to view", view_options, label_visibility="collapsed",
                        key=view_key, on_change=_store_widget, args=("build_view",),
                    )
                    with st.expander("Viewer Options", expanded=False):
                        vis_a, vis_b, vis_c, vis_d = st.columns(4)
                        show_protein = vis_a.toggle(
                            "Protein",
                            key=_prime_widget("viewer_show_protein"),
                            on_change=_store_widget,
                            args=("viewer_show_protein",),
                        )
                        show_solute = vis_b.toggle(
                            "Free molecules",
                            key=_prime_widget("viewer_show_solute"),
                            on_change=_store_widget,
                            args=("viewer_show_solute",),
                            help="Turn off the other components to inspect every requested copy clearly.",
                        )
                        show_solvent = vis_c.toggle(
                            "Solvent",
                            key=_prime_widget("viewer_show_solvent"),
                            on_change=_store_widget,
                            args=("viewer_show_solvent",),
                        )
                        show_ions = vis_d.toggle(
                            "Ions",
                            key=_prime_widget("viewer_show_ions"),
                            on_change=_store_widget,
                            args=("viewer_show_ions",),
                        )
                        show_connectivity = st.toggle(
                            "Show protein topology connectivity",
                            key=_prime_widget("viewer_show_connectivity"),
                            on_change=_store_widget,
                            args=("viewer_show_connectivity",),
                            help="Draws short protein bonds as thin white cylinders, matching MartiniSurf.",
                        )
                        va, vb = st.columns(2)
                        bead_radius = va.number_input(
                            "Bead radius", min_value=0.05, max_value=3.0, step=0.05,
                            key=_prime_widget("viewer_bead_radius"),
                            on_change=_store_widget, args=("viewer_bead_radius",),
                        )
                        bond_radius = vb.number_input(
                            "Bond radius", min_value=0.01, max_value=1.0, step=0.01,
                            key=_prime_widget("viewer_bond_radius"),
                            on_change=_store_widget, args=("viewer_bond_radius",),
                        )
                        topology_bond_max_nm = st.number_input(
                            "Topology bond max (nm)", min_value=0.05, max_value=2.0, step=0.05,
                            key=_prime_widget("viewer_topology_bond_max_nm"),
                            on_change=_store_widget, args=("viewer_topology_bond_max_nm",),
                        )
                    view_path = system_gro if selected == "Full solvated system" else protein_cg
                    stats = render_build_viewer(
                        view_path,
                        built,
                        height=800,
                        show_connectivity=bool(show_connectivity),
                        bead_radius=float(bead_radius),
                        bond_radius=float(bond_radius),
                        topology_bond_max_nm=float(topology_bond_max_nm),
                        show_protein=bool(show_protein),
                        show_solute=bool(show_solute),
                        show_solvent=bool(show_solvent),
                        show_ions=bool(show_ions),
                    )
                    if view_path.suffix.lower() == ".gro":
                        st.caption(
                            f"Viewer: {stats['bonds']} protein bonds drawn, "
                            f"{stats['skipped_long']} long topology contacts skipped."
                        )
                else:
                    st.info("No generated structure is available for preview.")
            else:
                st.markdown('<div class="ms-empty-preview"><strong>No generated files yet.</strong><span>Build the system to inspect the Martini structure here.</span></div>', unsafe_allow_html=True)

        if built.is_dir() and (built / "manifest.json").is_file():
            st.success("System generated and checked with GROMACS.")
            manifest = json.loads((built / "manifest.json").read_text())
            composition = manifest.get("composition") or {}
            if config.solvent in {"reline", "chcl_sorbitol"} and composition:
                comp_cols = st.columns(5)
                comp_cols[0].metric("Choline", int(composition.get("CHOL", 0)))
                hbd_name = "Urea" if config.solvent == "reline" else "Sorbitol"
                hbd_key = "UREA" if config.solvent == "reline" else "SOR"
                comp_cols[1].metric(hbd_name, int(composition.get(hbd_key, 0)))
                comp_cols[2].metric("Chloride", int(composition.get("CL", 0)))
                comp_cols[3].metric("Sodium", int(composition.get("NA", 0)))
                comp_cols[4].metric("Water beads", int(composition.get("W", 0)))
                neutralizing_cl = int(composition.get("neutralizing_CL", 0))
                neutralizing_na = int(composition.get("neutralizing_NA", 0))
                protein_charge = float(composition.get("protein_net_charge", 0.0))
                if neutralizing_cl or neutralizing_na:
                    st.caption(
                        f"Protein net charge: {protein_charge:+.0f} e · automatic neutralization: "
                        f"{neutralizing_cl} extra CL⁻ and {neutralizing_na} extra NA⁺."
                    )
                actual_water = composition.get("x_water_actual")
                if actual_water is not None:
                    ratio = "ChCl:urea = 1:2" if config.solvent == "reline" else "ChCl:sorbitol = 1:1"
                    st.caption(
                        f"Built DES composition · {ratio} · "
                        f"target dry density = {float(composition.get('target_density_g_cm3', 1.20)):.2f} g/cm³ · "
                        f"actual water mole fraction ≈ {float(actual_water):.3f}."
                    )
            with st.expander("Build details"):
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
                free_rows = manifest.get("free_molecules") or []
                if free_rows:
                    st.markdown("##### Free molecules included")
                    st.dataframe(
                        [
                            {
                                "Molecule": row.get("name", ""),
                                "Source": "Uploaded ITP" if row.get("source") == "upload" else "SMILES",
                                "Copies requested": int(row.get("copies_requested", row.get("count", 0))),
                                "Copies included": int(row.get("copies_included", row.get("count", 0))),
                                "Beads / molecule": int(row.get("beads_per_molecule", 0)),
                                "Charge / molecule": int(row.get("net_charge", 0)),
                                "Charged bead": row.get("charged_bead") or "—",
                                "SMILES": row.get("smiles", ""),
                            }
                            for row in free_rows
                        ],
                        hide_index=True,
                        use_container_width=True,
                    )
                    if all(
                        int(row.get("copies_included", row.get("count", 0)))
                        == int(row.get("copies_requested", row.get("count", 0)))
                        for row in free_rows
                    ):
                        st.caption("✓ Requested free-molecule copy counts are present in the final topology.")
            st.download_button(
                "Download Simulation_Files",
                archive(built),
                file_name="Simulation_Files.zip",
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
        top_a.text_input(
            "Simulation name", key=_prime_widget("short_md_output_tag"),
            on_change=_store_widget, args=("short_md_output_tag",),
        )
        top_b.number_input(
            "XTC write every (ps)", min_value=0.001, step=1.0, format="%.4f",
            key=_prime_widget("short_md_xtc_write_every_ps"),
            on_change=_store_widget, args=("short_md_xtc_write_every_ps",),
        )
        top_c.number_input(
            "GROMPP max warnings", min_value=0, step=1,
            key=_prime_widget("short_md_grompp_maxwarn"),
            on_change=_store_widget, args=("short_md_grompp_maxwarn",),
        )

        with st.expander("Short MD protocol", expanded=True):
            for stage in STAGE_ORDER:
                defaults = DEFAULT_STAGE_SETTINGS[stage]
                cols = st.columns([0.85, 1, 1])
                run_name = f"short_md_run_{stage}"
                cols[0].toggle(
                    stage.upper() if stage != "production" else "Production",
                    key=_prime_widget(run_name),
                    on_change=_store_widget, args=(run_name,),
                )
                cols[1].number_input(
                    f"{stage.upper() if stage != 'production' else 'Production'} timestep (ps)",
                    min_value=0.000001,
                    step=float(defaults["dt_ps"]),
                    format="%.4f",
                    key=_prime_widget(f"short_md_{stage}_dt"),
                    on_change=_store_widget, args=(f"short_md_{stage}_dt",),
                )
                cols[2].number_input(
                    f"{stage.upper() if stage != 'production' else 'Production'} time (ns)",
                    min_value=0.000001,
                    step=float(defaults["time_ns"]),
                    format="%.4f",
                    key=_prime_widget(f"short_md_{stage}_time"),
                    on_change=_store_widget, args=(f"short_md_{stage}_time",),
                )
            st.number_input(
                "CPU threads", min_value=1, max_value=4, step=1,
                key=_prime_widget("short_md_threads"),
                on_change=_store_widget, args=("short_md_threads",),
            )

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
                    key=_prime_widget("short_md_view_stage"),
                    on_change=_store_widget,
                    args=("short_md_view_stage",),
                    format_func=lambda value: value.upper() if value != "production" else "Production",
                )
                controls[1].number_input(
                    "Frame stride", min_value=1, step=1,
                    key=_prime_widget("short_md_view_stride"),
                    on_change=_store_widget, args=("short_md_view_stride",),
                )

                toggle_a, toggle_b, toggle_c, toggle_d = st.columns(4)
                toggle_a.toggle(
                    "Protein", key=_prime_widget("short_md_view_protein"),
                    on_change=_store_widget, args=("short_md_view_protein",),
                )
                toggle_b.toggle(
                    "Free molecules", key=_prime_widget("short_md_view_solute"),
                    on_change=_store_widget, args=("short_md_view_solute",),
                )
                toggle_c.toggle(
                    "Solvent", key=_prime_widget("short_md_view_solvent"),
                    on_change=_store_widget, args=("short_md_view_solvent",),
                    help="Shows the active solvent: water, Reline (choline + urea + intrinsic chloride), or ChCl:sorbitol (choline + sorbitol + intrinsic chloride). Wet DES also includes water.",
                )
                toggle_d.toggle(
                    "Counterions", key=_prime_widget("short_md_view_ions"),
                    on_change=_store_widget, args=("short_md_view_ions",),
                    help="Shows only neutralizing/free ions. Chloride belonging to the selected ChCl-based solvent is shown with Solvent.",
                )

                selected = next(row for row in trajectory_rows if row["name"] == st.session_state.short_md_view_stage)
                gro = Path(selected["gro"])
                xtc = Path(selected["xtc"])
                if gro.is_file() and xtc.is_file():
                    try:
                        manifest_path = built / "manifest.json"
                        reline_chloride_count = 0
                        if manifest_path.is_file():
                            manifest = json.loads(manifest_path.read_text())
                            composition = manifest.get("composition") or {}
                            reline_chloride_count = int(
                                composition.get("intrinsic_chloride", composition.get("reline_chloride", 0))
                            )
                        frames = render_trajectory(
                            gro,
                            xtc,
                            stride=int(st.session_state.short_md_view_stride),
                            height=700,
                            show_protein=bool(st.session_state.short_md_view_protein),
                            show_solute=bool(st.session_state.short_md_view_solute),
                            show_solvent=bool(st.session_state.short_md_view_solvent),
                            show_ions=bool(st.session_state.short_md_view_ions),
                            reline_chloride_count=reline_chloride_count,
                        )
                        st.caption(
                            f"{str(st.session_state.short_md_view_stage).upper()} trajectory preview: "
                            f"{frames} displayed frames."
                        )
                    except Exception as exc:
                        st.warning(f"Trajectory preview is unavailable: {exc}")

                _render_stage_analysis(selected)

            st.download_button(
                "Download Simulation_Files + Short MD",
                archive(built),
                file_name="Simulation_Files.zip",
                mime="application/zip",
                use_container_width=True,
            )
            log_path = Path(str(st.session_state.get("short_md_log", "")))
            if log_path.is_file():
                with st.expander("Short MD execution log"):
                    st.code(log_path.read_text(errors="replace")[-12000:], language="text")
