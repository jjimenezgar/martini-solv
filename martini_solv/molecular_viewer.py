"""3D structure and trajectory viewers used by the Streamlit interface."""
from __future__ import annotations

import json
from pathlib import Path

import streamlit.components.v1 as components


PINK = "#FF4FA3"
PINK_LIGHT = "#FF9DCC"
BG = "#07131C"


def _viewer_html(data: str, fmt: str, *, height: int, trajectory: bool = False) -> str:
    loader = (
        f"viewer.addModelsAsFrames({json.dumps(data)}, {json.dumps(fmt)});"
        if trajectory
        else f"viewer.addModel({json.dumps(data)}, {json.dumps(fmt)});"
    )
    animation = """
      viewer.animate({loop: "forward", interval: 160});
    """ if trajectory else ""
    return f"""
    <div class="viewer-shell"><div id="viewer" class="viewer"></div></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: "{BG}"}});
      {loader}
      viewer.setStyle({{}}, {{sphere: {{radius: 0.18, color: "{PINK_LIGHT}"}}}});
      viewer.setStyle({{atom: "BB"}}, {{sphere: {{radius: 0.22, color: "{PINK}"}}}});
      viewer.setStyle({{atom: /^BB\\d+$/}}, {{sphere: {{radius: 0.22, color: "{PINK}"}}}});
      viewer.setStyle({{atom: /^SC/}}, {{sphere: {{radius: 0.18, color: "{PINK_LIGHT}"}}}});
      viewer.setStyle({{resn: ["W","SW","TW","SOL"]}}, {{sphere: {{radius: 0.08, color: "lightgray", opacity: 0.38}}}});
      viewer.setStyle({{resn: ["NA","CL"]}}, {{sphere: {{radius: 0.12, colorscheme: "Jmol"}}}});
      viewer.zoomTo();
      viewer.render();
      {animation}
    </script>
    <style>
      html, body {{ margin:0; padding:0; overflow:hidden; background:{BG}; }}
      .viewer-shell {{
        width:100%; height:{height}px; box-sizing:border-box;
        border:1px solid rgba(116,152,170,.28); border-radius:14px;
        overflow:hidden; background:{BG};
        box-shadow: inset 0 0 42px rgba(53,201,211,.06);
      }}
      .viewer {{ width:100%; height:{height}px; }}
    </style>
    """


def render_structure_preview(pdb_text: str, height: int = 430) -> None:
    """Render an atomistic PDB as the pink cartoon used by MartiniSurf."""
    script = f"""
    <div class="viewer-shell"><div id="viewer" class="viewer"></div></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: "{BG}"}});
      viewer.addModel({json.dumps(pdb_text)}, "pdb");
      viewer.setStyle({{}}, {{cartoon: {{color: "{PINK}"}}}});
      viewer.addStyle({{hetflag: true}}, {{stick: {{radius: 0.16, colorscheme: "Jmol"}}}});
      viewer.zoomTo();
      viewer.render();
    </script>
    <style>
      html, body {{ margin:0; padding:0; overflow:hidden; background:{BG}; }}
      .viewer-shell {{
        width:100%; height:{height}px; box-sizing:border-box;
        border:1px solid rgba(116,152,170,.28); border-radius:14px;
        overflow:hidden; background:{BG};
        box-shadow: inset 0 0 42px rgba(53,201,211,.06);
      }}
      .viewer {{ width:100%; height:{height}px; }}
    </style>
    """
    components.html(script, height=height + 2)


def render_cg_structure(path: Path, height: int = 650) -> None:
    text = path.read_text(errors="replace")
    fmt = "gro" if path.suffix.lower() == ".gro" else "pdb"
    components.html(_viewer_html(text, fmt, height=height), height=height + 2)


def _trajectory_as_multimodel_pdb(gro_path: Path, xtc_path: Path, stride: int, max_frames: int = 80) -> tuple[str, int]:
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
    for index, raw in enumerate(lines[2:2 + atom_count], start=1):
        atoms.append({
            "serial": index,
            "resid": int(raw[0:5]),
            "resn": raw[5:10].strip() or "MOL",
            "name": raw[10:15].strip() or "B",
        })

    blocks: list[str] = []
    for frame_index, xyz in enumerate(traj.xyz, start=1):
        blocks.append(f"MODEL     {frame_index:4d}")
        for atom, coord in zip(atoms, xyz):
            x, y, z = (float(value) * 10.0 for value in coord)
            blocks.append(
                f"ATOM  {atom['serial']:5d} {atom['name'][:4]:>4} {atom['resn'][:3]:>3} A"
                f"{atom['resid'] % 10000:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          C"
            )
        blocks.append("ENDMDL")
    blocks.append("END")
    return "\n".join(blocks) + "\n", traj.n_frames


def render_trajectory(gro_path: Path, xtc_path: Path, *, stride: int = 1, height: int = 650) -> int:
    pdb, frames = _trajectory_as_multimodel_pdb(gro_path, xtc_path, stride)
    components.html(_viewer_html(pdb, "pdb", height=height, trajectory=True), height=height + 2)
    return frames
