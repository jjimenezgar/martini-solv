"""Targeted visualization improvements for MartiniSolv.

This module patches the public viewer functions at package import time without
changing the simulation workflow.  It keeps the existing GRO/ITP parsing
helpers in :mod:`martini_solv.molecular_viewer` as the single source of truth.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from . import molecular_viewer as mv


LEGEND_SPEC = (
    ("protein", "Protein", mv.BLUE),
    ("water", "Water", "#B0BEC5"),
    ("choline", "Choline", "#7E57C2"),
    ("urea", "Urea", "#FFB74D"),
    ("sorbitol", "Sorbitol", "#FF7043"),
    ("ions", "Ions", "limegreen"),
    ("solute", "Free molecule", "yellow"),
)


def _dynamic_legend_html(
    components: dict[str, list[str]],
    *,
    show_protein: bool,
    show_solute: bool,
    show_solvent: bool,
    show_ions: bool,
) -> str:
    """Return legend items only for components present *and* currently shown."""
    enabled = {
        "protein": show_protein,
        "water": show_solvent,
        "choline": show_solvent,
        "urea": show_solvent,
        "sorbitol": show_solvent,
        "ions": show_ions,
        "solute": show_solute,
    }
    items = []
    for key, label, color in LEGEND_SPEC:
        if enabled[key] and components.get(key):
            items.append(f'<span><i style="background:{color}"></i>{label}</span>')
    return "\n        ".join(items)


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
    """Build viewer with a legend derived from the actual GRO composition."""
    data = path.read_text(errors="replace")
    fmt = "gro" if path.suffix.lower() == ".gro" else "pdb"
    components_map = mv._component_resnames(path if path.suffix.lower() == ".gro" else None)
    cylinders: list[dict[str, dict[str, float]]] = []
    skipped = 0
    if show_connectivity and show_protein and path.suffix.lower() == ".gro":
        cylinders, skipped = mv._short_bond_cylinders(
            path, mv._protein_bonds(system_dir), topology_bond_max_nm
        )

    legend = _dynamic_legend_html(
        components_map,
        show_protein=show_protein,
        show_solute=show_solute,
        show_solvent=show_solvent,
        show_ions=show_ions,
    )
    legend_block = f'<div class="viewer-legend">{legend}</div>' if legend else ""

    script = f"""
    <div class="viewer-shell step4">
      <div id="viewer" class="viewer"></div>
      <div class="viewer-badge">Visual quality check</div>
      {legend_block}
    </div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: "{mv.BG}"}});
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
        viewer.setStyle({{resn: components.protein}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "{mv.BLUE_LIGHT}"}}}});
      }}
      if ({json.dumps(bool(show_protein))}) {{
        viewer.setStyle({{atom: "BB"}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "{mv.BLUE}"}}}});
        viewer.setStyle({{atom: /^BB\\d+$/}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "{mv.BLUE}"}}}});
        viewer.setStyle({{atom: /^SC/}}, {{sphere: {{radius: {float(bead_radius):.4f}, color: "{mv.BLUE_LIGHT}"}}}});
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
      html, body {{margin:0;padding:0;overflow:hidden;background:{mv.BG};}}
      .viewer-shell {{
        position:relative;width:100%;height:{height}px;box-sizing:border-box;
        border:1px solid rgba(116,152,170,.28);border-radius:16px;
        overflow:hidden;background:{mv.BG};box-shadow:inset 0 0 46px rgba(53,201,211,.06);
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
    mv._components_html(script, height=height + 2)
    return {"bonds": len(cylinders), "skipped_long": skipped}


def _detailed_component(atom: dict[str, float | int | str], broad: str) -> str:
    resn = str(atom["resn"]).strip().upper()
    if broad == "protein":
        return "protein"
    if broad == "ions":
        return "ions"
    if broad == "solute":
        return "solute"
    if resn in mv.WATER_RESN:
        return "water"
    if resn in mv.CHOLINE_RESN:
        return "choline"
    if resn in mv.UREA_RESN:
        return "urea"
    if resn in mv.SORBITOL_RESN:
        return "sorbitol"
    return "solvent"


def _stable_subsample(indices: list[int], limit: int) -> list[int]:
    """Deterministically thin dense solvent while keeping the same atoms every frame."""
    if len(indices) <= limit:
        return indices
    if limit <= 1:
        return indices[:1]
    return sorted({indices[round(i * (len(indices) - 1) / (limit - 1))] for i in range(limit)})


def _kabsch_align(xyz, atom_indices: list[int]):
    """Rigidly align every frame to frame 0 using selected protein atoms."""
    import numpy as np

    if not atom_indices:
        return xyz.copy()
    aligned = xyz.copy()
    reference = xyz[0, atom_indices, :]
    ref_center = reference.mean(axis=0)
    ref0 = reference - ref_center
    for frame_index in range(aligned.shape[0]):
        mobile = xyz[frame_index, atom_indices, :]
        mobile_center = mobile.mean(axis=0)
        mob0 = mobile - mobile_center
        covariance = mob0.T @ ref0
        u, _s, vt = np.linalg.svd(covariance)
        rotation = u @ vt
        if np.linalg.det(rotation) < 0:
            u[:, -1] *= -1
            rotation = u @ vt
        aligned[frame_index] = (xyz[frame_index] - mobile_center) @ rotation + ref_center
    return aligned


def _backbone_edges(atoms, xyz0, bb_indices: list[int]) -> list[tuple[int, int]]:
    """Connect consecutive Martini backbone beads only when physically adjacent."""
    edges: list[tuple[int, int]] = []
    for left, right in zip(bb_indices, bb_indices[1:]):
        a = atoms[left]
        b = atoms[right]
        # Prevent drawing a line between separate chains / distant fragments.
        dx = xyz0[left] - xyz0[right]
        distance_nm = float((dx * dx).sum() ** 0.5)
        resid_gap = abs(int(a["resid"]) - int(b["resid"]))
        if distance_nm <= 0.75 and resid_gap <= 2:
            edges.append((left, right))
    return edges


def generate_trajectory_gif(
    gro_path: Path,
    xtc_path: Path,
    output_path: Path,
    *,
    show_protein: bool = True,
    show_solute: bool = False,
    show_solvent: bool = False,
    show_ions: bool = False,
    reline_chloride_count: int = 0,
    width: int = 720,
    height: int = 540,
    fps: int = 10,
    max_frames: int = mv.MAX_TRAJECTORY_GIF_FRAMES,
) -> tuple[Path, int]:
    """Render a protein-first trajectory GIF using only rigid alignment of real MD frames."""
    try:
        import mdtraj as md
        from PIL import Image, ImageDraw
    except ImportError as exc:
        raise RuntimeError("Trajectory GIF generation requires mdtraj and Pillow") from exc

    if not gro_path.is_file() or not xtc_path.is_file():
        raise FileNotFoundError("The selected GRO/XTC trajectory files are unavailable")
    if not any((show_protein, show_solute, show_solvent, show_ions)):
        raise ValueError("Select at least one trajectory component before generating the GIF")
    if fps < 1:
        raise ValueError("GIF FPS must be positive")

    traj = md.load(str(xtc_path), top=str(gro_path))
    if traj.n_frames == 0:
        raise RuntimeError("The selected trajectory contains no frames")
    frame_indices = mv._gif_frame_indices(traj.n_frames, max_frames=max_frames)
    traj = traj[frame_indices]

    atoms = mv._parse_gro_atoms(gro_path)
    broad_components = mv._gif_component_map(gro_path, reline_chloride_count)
    if len(atoms) != traj.n_atoms or len(broad_components) != traj.n_atoms:
        raise RuntimeError("Trajectory atom count does not match the GRO topology")
    detailed = [_detailed_component(atom, broad) for atom, broad in zip(atoms, broad_components)]

    visible_broad = {
        "protein": bool(show_protein),
        "solute": bool(show_solute),
        "solvent": bool(show_solvent),
        "ions": bool(show_ions),
    }
    visible_indices = [
        i for i, broad in enumerate(broad_components) if visible_broad.get(broad, False)
    ]
    if not visible_indices:
        raise ValueError("The selected trajectory components contain no atoms")

    # Dense water is deliberately thinned for legibility; coordinates remain real
    # and the exact same molecules are kept throughout the animation.
    if show_solvent:
        water = [i for i in visible_indices if detailed[i] == "water"]
        other_solvent = [i for i in visible_indices if broad_components[i] == "solvent" and detailed[i] != "water"]
        keep_water = set(_stable_subsample(water, 500))
        keep_other = set(_stable_subsample(other_solvent, 900))
        visible_indices = [
            i for i in visible_indices
            if broad_components[i] != "solvent" or i in keep_water or i in keep_other
        ]

    protein_indices = [i for i, broad in enumerate(broad_components) if broad == "protein"]
    bb_indices = [
        i for i in protein_indices
        if str(atoms[i]["name"]).strip().upper().startswith("BB")
    ]
    align_indices = bb_indices or protein_indices
    aligned_xyz = _kabsch_align(traj.xyz, align_indices if show_protein else visible_indices)
    backbone_edges = _backbone_edges(atoms, aligned_xyz[0], bb_indices) if show_protein else []

    center_indices = align_indices if show_protein and align_indices else visible_indices
    center = aligned_xyz[:, center_indices, :].mean(axis=(0, 1))
    coords = aligned_xyz - center

    ay = math.radians(24.0)
    ax = math.radians(-18.0)
    cy, sy = math.cos(ay), math.sin(ay)
    cx, sx = math.cos(ax), math.sin(ax)

    projected_frames: list[dict[int, tuple[float, float, float]]] = []
    all_xy: list[tuple[float, float]] = []
    for frame in coords:
        projected: dict[int, tuple[float, float, float]] = {}
        for atom_index in visible_indices:
            x, y, z = map(float, frame[atom_index])
            x1 = cy * x + sy * z
            z1 = -sy * x + cy * z
            y2 = cx * y - sx * z1
            z2 = sx * y + cx * z1
            projected[atom_index] = (x1, y2, z2)
            all_xy.append((x1, y2))
        projected_frames.append(projected)

    min_x = min(x for x, _ in all_xy)
    max_x = max(x for x, _ in all_xy)
    min_y = min(y for _, y in all_xy)
    max_y = max(y for _, y in all_xy)
    span_x = max(max_x - min_x, 0.1)
    span_y = max(max_y - min_y, 0.1)
    margin = 44.0
    scale = min((width - 2 * margin) / span_x, (height - 2 * margin) / span_y)

    palette = {
        "protein_bb": "#168D9B",
        "protein_sc": "#65C9D4",
        "water": "#D5DCE1",
        "choline": "#7E57C2",
        "urea": "#FFB74D",
        "sorbitol": "#FF7043",
        "solvent": "#B8C0C6",
        "solute": "#E2B600",
        "ions": "#48A868",
    }
    radii = {
        "protein_bb": 5,
        "protein_sc": 3,
        "water": 2,
        "choline": 3,
        "urea": 3,
        "sorbitol": 3,
        "solvent": 3,
        "solute": 5,
        "ions": 3,
    }

    def screen(point: tuple[float, float, float]) -> tuple[float, float]:
        px, py, _pz = point
        return margin + (px - min_x) * scale, height - (margin + (py - min_y) * scale)

    images = []
    for projected in projected_frames:
        image = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(image)

        # A continuous backbone makes the coarse-grained protein read as a
        # molecular structure rather than as an unrelated collection of beads.
        for left, right in backbone_edges:
            if left in projected and right in projected:
                draw.line([screen(projected[left]), screen(projected[right])], fill="#168D9B", width=8)

        # Draw beads back-to-front; side chains remain subtle around the tube.
        for atom_index in sorted(projected, key=lambda idx: projected[idx][2]):
            broad = broad_components[atom_index]
            atom_name = str(atoms[atom_index]["name"]).strip().upper()
            if broad == "protein":
                kind = "protein_bb" if atom_name.startswith("BB") else "protein_sc"
            else:
                kind = detailed[atom_index]
            color = palette.get(kind, palette.get(broad, "#B8C0C6"))
            radius = radii.get(kind, radii.get(broad, 3))
            sxp, syp = screen(projected[atom_index])
            draw.ellipse(
                [sxp - radius, syp - radius, sxp + radius, syp + radius],
                fill=color,
                outline=None if kind == "water" else "#667078",
                width=1,
            )
        images.append(image)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    duration_ms = max(40, round(1000 / fps))
    images[0].save(
        output_path,
        format="GIF",
        save_all=True,
        append_images=images[1:],
        duration=duration_ms,
        loop=0,
        optimize=True,
    )
    return output_path, len(images)


def apply() -> None:
    """Install the improved public functions on the existing viewer module."""
    mv.render_build_viewer = render_build_viewer
    mv.generate_trajectory_gif = generate_trajectory_gif
