"""3D structure and trajectory viewers used by the Streamlit interface."""
from __future__ import annotations

import json
import re
from pathlib import Path

import streamlit.components.v1 as components


BLUE = "#42C7D5"
BLUE_LIGHT = "#8FEAF2"
BG = "#07131C"
TRAJ_BG = "#0E0D11"

WATER_RESN = {"W", "WF", "SW", "TW", "SOL"}
CHOLINE_RESN = {"CHOL", "CHO"}
UREA_RESN = {"UREA", "URE"}
SORBITOL_RESN = {"SOR"}
ION_RESN = {"NA", "CL", "ION", "K", "CA", "MG", "ZN", "LI", "RB", "CS", "BA", "SR", "F", "BR", "I"}


def render_structure_preview(pdb_text: str, height: int = 430) -> None:
    """Render an atomistic PDB in the MartiniSolv blue accent."""
    script = f"""
    <div class="viewer-shell"><div id="viewer" class="viewer"></div></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: "{BG}"}});
      viewer.addModel({json.dumps(pdb_text)}, "pdb");
      viewer.setStyle({{}}, {{cartoon: {{color: "{BLUE}"}}}});
      viewer.addStyle({{hetflag: true}}, {{stick: {{radius: 0.16, colorscheme: "Jmol"}}}});
      viewer.zoomTo();
      viewer.render();
    </script>
    <style>
      html, body {{ margin:0; padding:0; overflow:hidden; background:{BG}; }}
      .viewer-shell {{
        width:100%; height:{height}px; box-sizing:border-box;
        border:1px solid rgba(116,152,170,.28); border-radius:16px;
        overflow:hidden; background:{BG};
        box-shadow: inset 0 0 42px rgba(53,201,211,.06);
      }}
      .viewer {{ width:100%; height:{height}px; overflow:hidden; }}
    </style>
    """
    components.html(script, height=height + 2)


def _parse_gro_atoms(path: Path) -> list[dict[str, float | int | str]]:
    """Parse standard GRO plus the more compact GRO formatting emitted by some mappers."""
    lines = path.read_text(errors="replace").splitlines()
    atoms: list[dict[str, float | int | str]] = []
    if len(lines) < 3:
        return atoms
    try:
        count = int(lines[1].strip())
    except ValueError:
        return atoms

    for index, raw in enumerate(lines[2:2 + count], start=1):
        # Standard fixed-width GRO formatting.
        if len(raw) >= 44:
            try:
                atoms.append({
                    "serial": int(raw[15:20]),
                    "resid": int(raw[0:5]),
                    "resn": raw[5:10].strip(),
                    "name": raw[10:15].strip(),
                    "x": float(raw[20:28]),
                    "y": float(raw[28:36]),
                    "z": float(raw[36:44]),
                })
                continue
            except ValueError:
                pass

        # Martini Mapper can emit compact whitespace-delimited GRO records.
        # Parse coordinates from the right and recover residue/name/serial
        # without requiring exact fixed-width alignment.
        parts = raw.split()
        if len(parts) < 6:
            continue
        try:
            x, y, z = map(float, parts[-3:])
            serial = int(parts[-4])
            name = parts[-5]
            residue_token = "".join(parts[:-5]) or "1MOL"
            match = re.match(r"^(\d+)(.*)$", residue_token)
            resid = int(match.group(1)) if match else 1
            resn = (match.group(2) if match and match.group(2) else "MOL").strip()
            atoms.append({
                "serial": serial or index,
                "resid": resid,
                "resn": resn,
                "name": name,
                "x": x,
                "y": y,
                "z": z,
            })
        except (TypeError, ValueError):
            continue
    return atoms



def render_free_molecule_mapping(path: Path, height: int = 360) -> list[dict[str, object]]:
    """Render a generated Martini molecule with bead labels, as in MartiniSurf."""
    atoms = _parse_gro_atoms(path)
    labels = [
        {
            "label": f"{index}: {str(atom['name']).strip()}",
            "x": 10.0 * float(atom["x"]),
            "y": 10.0 * float(atom["y"]),
            "z": 10.0 * float(atom["z"]),
        }
        for index, atom in enumerate(atoms, start=1)
    ]
    data = path.read_text(errors="replace")
    script = f"""
    <div class="viewer-shell mapping">
      <div id="viewer_free_molecule" class="viewer"></div>
      <div class="viewer-badge">Martini bead topology</div>
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer(
        document.getElementById("viewer_free_molecule"),
        {{backgroundColor: "{BG}"}}
      );
      viewer.addModel({json.dumps(data)}, "gro");
      viewer.setStyle({{}}, {{sphere: {{radius: 0.78, color: "{BLUE}"}}}});
      const labels = {json.dumps(labels)};
      for (const bead of labels) {{
        viewer.addLabel(bead.label, {{
          position: {{x: bead.x, y: bead.y, z: bead.z}},
          fontColor: "#07131C",
          backgroundColor: "#F7FBFF",
          borderColor: "{BLUE}",
          borderThickness: 1,
          fontSize: 13,
          inFront: true
        }});
      }}
      viewer.zoomTo();
      viewer.render();
    </script>
    <style>
      html, body {{margin:0;padding:0;overflow:hidden;background:{BG};}}
      .viewer-shell {{
        position:relative;width:100%;height:{height}px;box-sizing:border-box;
        border:1px solid rgba(66,199,213,.28);border-radius:16px;
        overflow:hidden;background:{BG};box-shadow:inset 0 0 42px rgba(53,201,211,.06);
      }}
      .viewer {{width:100%;height:{height}px;overflow:hidden;}}
      .viewer-badge {{
        position:absolute;left:14px;bottom:14px;padding:8px 10px;
        border:1px solid rgba(66,199,213,.30);border-radius:999px;
        background:rgba(7,19,28,.78);color:#F3F7FA;
        font:700 12px/1.2 sans-serif;pointer-events:none;
      }}
    </style>
    """
    components.html(script, height=height + 2)
    return [
        {
            "Bead": f"{index}: {str(atom['name']).strip()}",
            "Residue": str(atom["resn"]).strip(),
            "x (nm)": round(float(atom["x"]), 4),
            "y (nm)": round(float(atom["y"]), 4),
            "z (nm)": round(float(atom["z"]), 4),
        }
        for index, atom in enumerate(atoms, start=1)
    ]

def _component_resnames(gro_path: Path | None) -> dict[str, list[str]]:
    groups = {
        "protein": set(),
        "water": set(),
        "choline": set(),
        "urea": set(),
        "sorbitol": set(),
        "reline_chloride": set(),
        "ions": set(),
        "solute": set(),
    }
    if not gro_path or not gro_path.exists():
        return {key: [] for key in groups}

    by_residue: dict[tuple[int, str], set[str]] = {}
    for atom in _parse_gro_atoms(gro_path):
        resid = int(atom["resid"])
        resn = str(atom["resn"]).strip()
        name = str(atom["name"]).strip().upper()
        by_residue.setdefault((resid, resn), set()).add(name)

    for (_resid, resn), names in by_residue.items():
        upper = resn.upper()
        if upper in WATER_RESN:
            groups["water"].add(resn)
        elif upper in CHOLINE_RESN or {"N1", "OH"}.issubset(names):
            groups["choline"].add(resn)
        elif upper in UREA_RESN or {"N1", "UP", "UN"}.issubset(names):
            groups["urea"].add(resn)
        elif upper in SORBITOL_RESN or {"S1", "S2", "S3"}.issubset(names):
            groups["sorbitol"].add(resn)
        elif upper in ION_RESN:
            groups["ions"].add(resn)
        elif "BB" in names or any(name.startswith("SC") for name in names) or any(
            name.startswith("BB") and name[2:].isdigit() for name in names
        ):
            groups["protein"].add(resn)
        else:
            groups["solute"].add(resn)
    return {key: sorted(value) for key, value in groups.items()}


def _parse_itp_bonds(itp_path: Path) -> list[tuple[int, int]]:
    bonds: list[tuple[int, int]] = []
    section = None
    skip_block = False
    for raw in itp_path.read_text(errors="ignore").splitlines():
        stripped = raw.strip()
        directive = stripped.upper()
        if directive.startswith("#IFDEF GO_VIRT") or directive.startswith("#IFDEF RUBBER_BANDS"):
            skip_block = True
            continue
        if skip_block:
            if directive.startswith("#ENDIF"):
                skip_block = False
            continue

        if stripped.startswith("[") and "]" in stripped:
            section = stripped.strip("[]").strip().lower()
            continue

        line = raw.split(";", 1)[0].strip()
        if not line or line.startswith("#") or section not in {"bonds", "constraints"}:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            left, right = int(parts[0]), int(parts[1])
        except ValueError:
            continue
        if left != right:
            bonds.append(tuple(sorted((left, right))))
    return sorted(set(bonds))


def _protein_bonds(system_dir: Path) -> list[tuple[int, int]]:
    preferred = system_dir / "Protein.itp"
    if preferred.is_file():
        return _parse_itp_bonds(preferred)

    skip = ("martini_", "go_", "posre", "choline", "urea")
    for path in sorted(system_dir.glob("*.itp")):
        if any(token in path.name.lower() for token in skip):
            continue
        bonds = _parse_itp_bonds(path)
        if bonds:
            return bonds
    return []


def _short_bond_cylinders(
    gro_path: Path,
    bonds: list[tuple[int, int]],
    max_distance_nm: float = 0.75,
) -> tuple[list[dict[str, dict[str, float]]], int]:
    atoms = _parse_gro_atoms(gro_path)
    by_serial = {int(atom["serial"]): atom for atom in atoms}
    cylinders: list[dict[str, dict[str, float]]] = []
    skipped = 0
    max_d2 = max_distance_nm ** 2
    for left, right in bonds:
        a = by_serial.get(left)
        b = by_serial.get(right)
        if a is None or b is None:
            continue
        d2 = sum((float(a[k]) - float(b[k])) ** 2 for k in ("x", "y", "z"))
        if d2 > max_d2:
            skipped += 1
            continue
        cylinders.append({
            "start": {k: 10.0 * float(a[k]) for k in ("x", "y", "z")},
            "end": {k: 10.0 * float(b[k]) for k in ("x", "y", "z")},
        })
    return cylinders, skipped


def render_build_viewer(
    path: Path,
    system_dir: Path,
    *,
    height: int = 800,
    show_connectivity: bool = True,
    bead_radius: float = 0.85,
    bond_radius: float = 0.20,
    topology_bond_max_nm: float = 0.75,
    show_protein: bool = True,
    show_solute: bool = True,
    show_solvent: bool = True,
    show_ions: bool = True,
) -> dict[str, int]:
    """MartiniSurf-style build viewer: large beads plus protein connectivity."""
    data = path.read_text(errors="replace")
    fmt = "gro" if path.suffix.lower() == ".gro" else "pdb"
    components_map = _component_resnames(path if path.suffix.lower() == ".gro" else None)
    cylinders: list[dict[str, dict[str, float]]] = []
    skipped = 0
    if show_connectivity and show_protein and path.suffix.lower() == ".gro":
        cylinders, skipped = _short_bond_cylinders(
            path, _protein_bonds(system_dir), topology_bond_max_nm
        )

    script = f"""
    <div class="viewer-shell step4">
      <div id="viewer" class="viewer"></div>
      <div class="viewer-badge">Visual quality check</div>
      <div class="viewer-legend">
        <span><i style="background:{BLUE}"></i>Protein</span>
        <span><i style="background:#B0BEC5"></i>Water</span>
        <span><i style="background:#7E57C2"></i>Choline</span>
        <span><i style="background:#FFB74D"></i>Urea</span>
        <span><i style="background:#FF7043"></i>Sorbitol</span>
        <span><i style="background:limegreen"></i>Ions</span>
        <span><i style="background:yellow"></i>Free molecule</span>
      </div>
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: "{BG}"}});
      viewer.addModel({json.dumps(data)}, {json.dumps(fmt)});
      viewer.setStyle({{}}, {{sphere: {{hidden: true}}}});
      const components = {json.dumps(components_map)};
      if ({json.dumps(bool(show_solvent))} && components.water.length) {{
        viewer.setStyle({{resn: components.water}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "#B0BEC5", opacity: 0.62}}}});
      }}
      if ({json.dumps(bool(show_solvent))} && components.choline.length) {{
        viewer.setStyle({{resn: components.choline}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "#7E57C2"}}}});
      }}
      if ({json.dumps(bool(show_solvent))} && components.urea.length) {{
        viewer.setStyle({{resn: components.urea}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "#FFB74D"}}}});
      }}
      if ({json.dumps(bool(show_solvent))} && components.sorbitol.length) {{
        viewer.setStyle({{resn: components.sorbitol}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "#FF7043"}}}});
      }}
      if ({json.dumps(bool(show_ions))} && components.ions.length) {{
        viewer.setStyle({{resn: components.ions}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "limegreen"}}}});
      }}
      if ({json.dumps(bool(show_solute))} && components.solute.length) {{
        viewer.setStyle({{resn: components.solute}}, {{sphere: {{radius: {float(bead_radius) * 1.12:.4f}, color: "yellow"}}}});
      }}
      if ({json.dumps(bool(show_protein))} && components.protein.length) {{
        viewer.setStyle({{resn: components.protein}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "{BLUE_LIGHT}"}}}});
      }}
      if ({json.dumps(bool(show_protein))}) {{
        viewer.setStyle({{atom: "BB"}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "{BLUE}"}}}});
        viewer.setStyle({{atom: /^BB\\d+$/}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "{BLUE}"}}}});
        viewer.setStyle({{atom: /^SC/}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "{BLUE_LIGHT}"}}}});
      }}
      for (const cylinder of {json.dumps(cylinders)}) {{
        viewer.addCylinder({{
          start: cylinder.start,
          end: cylinder.end,
          radius: {float(bond_radius):.4f},
          color: "#F4F5F7",
          fromCap: 1,
          toCap: 1
        }});
      }}
      viewer.zoomTo();
      viewer.render();
    </script>
    <style>
      html, body {{margin:0;padding:0;overflow:hidden;background:{BG};}}
      .viewer-shell {{
        position:relative;width:100%;height:{height}px;box-sizing:border-box;
        border:1px solid rgba(116,152,170,.28);border-radius:16px;
        overflow:hidden;background:{BG};box-shadow:inset 0 0 46px rgba(53,201,211,.06);
      }}
      .viewer {{width:100%;height:{height}px;overflow:hidden;}}
      .viewer-badge {{
        position:absolute;left:14px;bottom:14px;padding:8px 10px;
        border:1px solid rgba(53,201,211,.30);border-radius:999px;
        background:rgba(7,19,28,.78);color:#F3F7FA;
        font:700 12px/1.2 sans-serif;pointer-events:none;
      }}
      .viewer-legend {{
        position:absolute;right:14px;bottom:14px;display:flex;gap:10px;flex-wrap:wrap;
        max-width:70%;padding:7px 10px;border:1px solid rgba(116,152,170,.25);
        border-radius:12px;background:rgba(7,19,28,.78);color:#F3F7FA;
        font:600 11px/1.2 sans-serif;pointer-events:none;
      }}
      .viewer-legend span {{display:flex;align-items:center;gap:5px;}}
      .viewer-legend i {{width:9px;height:9px;border-radius:50%;display:inline-block;}}
    </style>
    """
    components.html(script, height=height + 2)
    return {"bonds": len(cylinders), "skipped_long": skipped}


def _trajectory_as_multimodel_pdb(
    gro_path: Path,
    xtc_path: Path,
    stride: int,
    max_frames: int = 80,
    reline_chloride_count: int = 0,
) -> tuple[str, int]:
    try:
        import mdtraj as md
    except ImportError as exc:
        raise RuntimeError("mdtraj is required for trajectory visualization") from exc

    traj = md.load(str(xtc_path), top=str(gro_path), stride=max(1, int(stride)))
    if traj.n_frames == 0:
        raise RuntimeError("The trajectory contains no frames")
    if traj.n_frames > max_frames:
        step = max(1, traj.n_frames // max_frames)
        traj = traj[::step][:max_frames]

    lines = gro_path.read_text(errors="replace").splitlines()
    atom_count = int(lines[1])
    atoms = []
    intrinsic_cl_seen = 0
    for index, raw in enumerate(lines[2:2 + atom_count], start=1):
        original_resn = raw[5:10].strip() or "MOL"
        atom_name = raw[10:15].strip() or "B"
        upper = original_resn.upper()
        namesafe = original_resn[:3]

        # PDB residue names are limited to three characters. Use explicit
        # viewer-only labels so CHOL/UREA remain selectable after GRO->PDB
        # conversion, and separate Reline chloride from neutralizing ions.
        if upper in CHOLINE_RESN:
            namesafe = "CHO"
        elif upper in UREA_RESN:
            namesafe = "URE"
        elif upper in SORBITOL_RESN:
            namesafe = "SOR"
        elif upper in WATER_RESN:
            namesafe = "WAT"
        elif upper in {"CL", "CL-"}:
            if intrinsic_cl_seen < max(0, int(reline_chloride_count)):
                namesafe = "RCL"
                intrinsic_cl_seen += 1
            else:
                namesafe = "ICL"
        elif upper in {"NA", "NA+"}:
            namesafe = "INA"

        atoms.append({
            "serial": index,
            "resid": int(raw[0:5]),
            "resn": namesafe,
            "name": atom_name,
        })

    blocks: list[str] = []
    for frame_index, xyz in enumerate(traj.xyz, start=1):
        blocks.append(f"MODEL     {frame_index:4d}")
        for atom, coord in zip(atoms, traj.xyz[frame_index - 1]):
            x, y, z = (float(value) * 10.0 for value in coord)
            blocks.append(
                f"ATOM  {atom['serial']:5d} {atom['name'][:4]:>4} {atom['resn'][:3]:>3} A"
                f"{atom['resid'] % 10000:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          C"
            )
        blocks.append("ENDMDL")
    blocks.append("END")
    return "\n".join(blocks) + "\n", traj.n_frames

def render_trajectory(
    gro_path: Path,
    xtc_path: Path,
    *,
    stride: int = 1,
    height: int = 700,
    show_protein: bool = True,
    show_solute: bool = False,
    show_solvent: bool = False,
    show_ions: bool = False,
    reline_chloride_count: int = 0,
) -> int:
    """MartiniSurf-style trajectory viewer with large, component-aware beads."""
    pdb, frames = _trajectory_as_multimodel_pdb(
        gro_path, xtc_path, stride, reline_chloride_count=reline_chloride_count
    )
    # The GRO->PDB trajectory conversion uses stable three-character viewer
    # labels so solvent components remain selectable in 3Dmol.
    gro_components = _component_resnames(gro_path)
    components_map = {
        "water": ["WAT"],
        "choline": ["CHO"],
        "urea": ["URE"],
        "sorbitol": ["SOR"],
        "reline_chloride": ["RCL"],
        "ions": ["ICL", "INA", *[name[:3] for name in gro_components.get("ions", []) if name.upper() not in {"CL", "CL-", "NA", "NA+"}]],
        "solute": sorted({name[:3] for name in gro_components.get("solute", [])}),
    }
    script = f"""
    <div class="viewer-shell short-md">
      <div id="viewer" class="viewer"></div>
      <div class="viewer-badge">Production trajectory</div>
      <div class="viewer-legend">
        <span><i style="background:{BLUE}"></i>Protein</span>
        <span><i style="background:#B0BEC5"></i>Water</span>
        <span><i style="background:#7E57C2"></i>Choline</span>
        <span><i style="background:#FFB74D"></i>Urea</span>
        <span><i style="background:#FF7043"></i>Sorbitol</span>
        <span><i style="background:limegreen"></i>Ions</span>
        <span><i style="background:yellow"></i>Free molecule</span>
      </div>
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: "{TRAJ_BG}"}});
      viewer.addModelsAsFrames({json.dumps(pdb)}, "pdb");
      viewer.setStyle({{}}, {{sphere: {{hidden: true}}}});
      const components = {json.dumps(components_map)};
      if ({json.dumps(bool(show_solvent))} && components.water.length) {{
        viewer.setStyle({{resn: components.water}}, {{sphere: {{radius: 0.462, color: "#B0BEC5", opacity: 0.62}}}});
      }}
      if ({json.dumps(bool(show_solvent))} && components.choline.length) {{
        viewer.setStyle({{resn: components.choline}}, {{sphere: {{radius: 0.74, color: "#7E57C2"}}}});
      }}
      if ({json.dumps(bool(show_solvent))} && components.urea.length) {{
        viewer.setStyle({{resn: components.urea}}, {{sphere: {{radius: 0.70, color: "#FFB74D"}}}});
      }}
      if ({json.dumps(bool(show_solvent))} && components.sorbitol.length) {{
        viewer.setStyle({{resn: components.sorbitol}}, {{sphere: {{radius: 0.72, color: "#FF7043"}}}});
      }}
      if ({json.dumps(bool(show_solvent))} && components.reline_chloride.length) {{
        viewer.setStyle({{resn: components.reline_chloride}}, {{sphere: {{radius: 0.56, color: "#EC407A"}}}});
      }}
      if ({json.dumps(bool(show_ions))} && components.ions.length) {{
        viewer.setStyle({{resn: components.ions}}, {{sphere: {{radius: 0.605, color: "limegreen"}}}});
      }}
      if ({json.dumps(bool(show_solute))} && components.solute.length) {{
        viewer.setStyle({{resn: components.solute}}, {{sphere: {{radius: 0.836, color: "yellow"}}}});
      }}
      if ({json.dumps(bool(show_protein))}) {{
        viewer.setStyle({{atom: "BB"}}, {{sphere: {{radius: 0.902, color: "{BLUE}"}}}});
        viewer.setStyle({{atom: /^BB\\d+$/}}, {{sphere: {{radius: 0.902, color: "{BLUE}"}}}});
        viewer.setStyle({{atom: /^SC/}}, {{sphere: {{radius: 0.902, color: "{BLUE_LIGHT}"}}}});
      }}
      viewer.zoomTo();
      viewer.animate({{loop: "forward", reps: 0}});
      viewer.render();
    </script>
    <style>
      html, body {{margin:0;padding:0;overflow:hidden;background:{TRAJ_BG};}}
      .viewer-shell {{
        position:relative;width:100%;height:{height}px;box-sizing:border-box;
        border:1px solid rgba(66,199,213,.28);border-radius:16px;
        overflow:hidden;background:{TRAJ_BG};box-shadow:inset 0 0 46px rgba(66,199,213,.08);
      }}
      .viewer {{width:100%;height:{height}px;overflow:hidden;}}
      .viewer-badge {{
        position:absolute;left:14px;bottom:14px;padding:8px 10px;
        border:1px solid rgba(66,199,213,.35);border-radius:999px;
        background:rgba(14,13,17,.78);color:#F5FAFA;
        font:700 12px/1.2 sans-serif;pointer-events:none;
      }}
      .viewer-legend {{
        position:absolute;right:14px;bottom:14px;display:flex;gap:10px;flex-wrap:wrap;
        max-width:70%;padding:7px 10px;border:1px solid rgba(66,199,213,.22);
        border-radius:12px;background:rgba(14,13,17,.78);color:#F5FAFA;
        font:600 11px/1.2 sans-serif;pointer-events:none;
      }}
      .viewer-legend span {{display:flex;align-items:center;gap:5px;}}
      .viewer-legend i {{width:9px;height:9px;border-radius:50%;display:inline-block;}}
    </style>
    """
    components.html(script, height=height + 2)
    return frames
