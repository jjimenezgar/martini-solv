"""Render MartiniSolv trajectory GIFs with a clean backbone-only protein view."""
from __future__ import annotations

import math
from pathlib import Path

from . import molecular_viewer as mv
from . import trajectory_pbc as tp
from . import visualization_overrides as vo


def _backbone_only_visible_indices(atoms, components, visible_indices):
    """Keep selected non-protein atoms, but only BB beads for protein."""
    kept = []
    for index in visible_indices:
        if components[index] != "protein":
            kept.append(index)
            continue
        name = str(atoms[index]["name"]).strip().upper()
        if name.startswith("BB"):
            kept.append(index)
    return kept


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
    """Render the PBC-safe GIF while drawing protein backbone only."""
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
    traj = traj[mv._gif_frame_indices(traj.n_frames, max_frames=max_frames)]

    atoms = mv._parse_gro_atoms(gro_path)
    components = mv._gif_component_map(gro_path, reline_chloride_count)
    if len(atoms) != traj.n_atoms or len(components) != traj.n_atoms:
        raise RuntimeError("Trajectory atom count does not match the GRO topology")
    detailed = [vo._detailed_component(atom, broad) for atom, broad in zip(atoms, components)]

    enabled = {
        "protein": bool(show_protein),
        "solute": bool(show_solute),
        "solvent": bool(show_solvent),
        "ions": bool(show_ions),
    }
    visible = [i for i, broad in enumerate(components) if enabled.get(broad, False)]
    visible = _backbone_only_visible_indices(atoms, components, visible)
    if not visible:
        raise ValueError("The selected trajectory components contain no atoms")

    if show_solvent:
        water = [i for i in visible if detailed[i] == "water"]
        other_solvent = [i for i in visible if components[i] == "solvent" and detailed[i] != "water"]
        keep_water = set(vo._stable_subsample(water, 500))
        keep_other = set(vo._stable_subsample(other_solvent, 900))
        visible = [
            i for i in visible
            if components[i] != "solvent" or i in keep_water or i in keep_other
        ]

    protein_indices = [i for i, broad in enumerate(components) if broad == "protein"]
    bb_indices = [
        i for i in protein_indices
        if str(atoms[i]["name"]).strip().upper().startswith("BB")
    ]
    align_indices = bb_indices if show_protein else []
    xyz = tp._prepared_xyz(traj, atoms, components, align_indices=align_indices)
    backbone_edges = vo._backbone_edges(atoms, xyz[0], bb_indices) if show_protein else []

    center_indices = bb_indices if show_protein and bb_indices else visible
    center = xyz[:, center_indices, :].mean(axis=(0, 1))
    coords = xyz - center

    ay = math.radians(24.0)
    ax = math.radians(-18.0)
    cy, sy = math.cos(ay), math.sin(ay)
    cx, sx = math.cos(ax), math.sin(ax)
    projected_frames = []
    all_xy = []
    for frame in coords:
        projected = {}
        for atom_index in visible:
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
        "protein": "#168D9B",
        "water": "#D5DCE1",
        "choline": "#7E57C2",
        "urea": "#FFB74D",
        "sorbitol": "#FF7043",
        "solvent": "#B8C0C6",
        "solute": "#E2B600",
        "ions": "#48A868",
    }
    radii = {
        "protein": 5,
        "water": 2,
        "choline": 3,
        "urea": 3,
        "sorbitol": 3,
        "solvent": 3,
        "solute": 5,
        "ions": 3,
    }

    def screen(point):
        px, py, _pz = point
        return margin + (px - min_x) * scale, height - (margin + (py - min_y) * scale)

    images = []
    for projected in projected_frames:
        image = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(image)
        for left, right in backbone_edges:
            if left in projected and right in projected:
                draw.line([screen(projected[left]), screen(projected[right])], fill="#168D9B", width=8)

        for atom_index in sorted(projected, key=lambda idx: projected[idx][2]):
            broad = components[atom_index]
            kind = "protein" if broad == "protein" else detailed[atom_index]
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
    mv.generate_trajectory_gif = generate_trajectory_gif
