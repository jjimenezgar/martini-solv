"""PBC-safe trajectory visualization for MartiniSolv.

The simulation coordinates stored in XTC are wrapped by GROMACS.  A protein that
crosses a periodic boundary can therefore look torn apart even though bonded
beads remain physically adjacent.  This module reconstructs the protein with a
minimum-image convention before any visual alignment, then removes global
translation/rotation only for display.  No MD coordinates are deformed beyond
periodic imaging and a rigid-body transform.
"""
from __future__ import annotations

import math
from pathlib import Path

from . import molecular_viewer as mv
from . import visualization_overrides as vo


def _minimum_image(delta, box):
    """Return the minimum-image displacement for an arbitrary triclinic box."""
    import numpy as np

    matrix = np.asarray(box, dtype=float)
    if matrix.shape != (3, 3) or not np.isfinite(matrix).all():
        return delta
    try:
        inv = np.linalg.inv(matrix)
    except np.linalg.LinAlgError:
        return delta
    fractional = np.asarray(delta, dtype=float) @ inv
    fractional -= np.round(fractional)
    return fractional @ matrix


def _protein_residue_blocks(atoms, components):
    """Return consecutive protein residue blocks in GRO atom order."""
    blocks: list[list[int]] = []
    current: list[int] = []
    current_key = None
    for index, (atom, component) in enumerate(zip(atoms, components)):
        if component != "protein":
            continue
        key = (int(atom["resid"]), str(atom["resn"]).strip())
        if current and key != current_key:
            blocks.append(current)
            current = []
        current.append(index)
        current_key = key
    if current:
        blocks.append(current)
    return blocks


def _bb_index(block, atoms):
    for index in block:
        if str(atoms[index]["name"]).strip().upper().startswith("BB"):
            return index
    return block[0] if block else None


def _make_protein_whole(xyz, unitcell_vectors, atoms, components):
    """Reconstruct protein residues and backbone across periodic boundaries.

    The first backbone bead of each continuous chain fragment is kept in its
    wrapped image.  Consecutive residues are then placed using minimum-image BB
    displacements, and all beads in each residue are imaged around that BB.
    """
    import numpy as np

    whole = np.array(xyz, dtype=float, copy=True)
    if unitcell_vectors is None:
        return whole
    boxes = np.asarray(unitcell_vectors, dtype=float)
    if boxes.ndim != 3 or boxes.shape[0] != whole.shape[0]:
        return whole

    blocks = _protein_residue_blocks(atoms, components)
    if not blocks:
        return whole

    for frame_index in range(whole.shape[0]):
        box = boxes[frame_index]
        if not np.isfinite(box).all() or abs(float(np.linalg.det(box))) < 1e-8:
            continue

        previous_raw_bb = None
        previous_whole_bb = None
        previous_resid = None

        for block in blocks:
            bb = _bb_index(block, atoms)
            if bb is None:
                continue
            raw_bb = np.asarray(xyz[frame_index, bb], dtype=float)
            resid = int(atoms[bb]["resid"])

            continuous = (
                previous_raw_bb is not None
                and previous_whole_bb is not None
                and previous_resid is not None
                and 0 < resid - previous_resid <= 2
            )
            if continuous:
                bb_delta = _minimum_image(raw_bb - previous_raw_bb, box)
                # Martini peptide BB neighbours are ~0.35 nm apart.  A generous
                # cutoff prevents accidentally joining distinct chains/fragments.
                if float(np.linalg.norm(bb_delta)) <= 0.80:
                    whole_bb = previous_whole_bb + bb_delta
                else:
                    whole_bb = raw_bb
            else:
                whole_bb = raw_bb

            for atom_index in block:
                local = _minimum_image(
                    np.asarray(xyz[frame_index, atom_index], dtype=float) - raw_bb,
                    box,
                )
                whole[frame_index, atom_index] = whole_bb + local

            previous_raw_bb = raw_bb
            previous_whole_bb = whole_bb
            previous_resid = resid

    return whole


def _rigid_align(xyz, atom_indices):
    """Align all coordinates to frame 0 using a rigid Kabsch transform."""
    import numpy as np

    if not atom_indices:
        return np.array(xyz, dtype=float, copy=True)
    aligned = np.array(xyz, dtype=float, copy=True)
    reference = aligned[0, atom_indices, :]
    ref_center = reference.mean(axis=0)
    ref0 = reference - ref_center

    for frame_index in range(aligned.shape[0]):
        mobile = aligned[frame_index, atom_indices, :]
        mobile_center = mobile.mean(axis=0)
        mob0 = mobile - mobile_center
        covariance = mob0.T @ ref0
        u, _s, vt = np.linalg.svd(covariance)
        rotation = u @ vt
        if np.linalg.det(rotation) < 0:
            u[:, -1] *= -1
            rotation = u @ vt
        aligned[frame_index] = (aligned[frame_index] - mobile_center) @ rotation + ref_center
    return aligned


def _prepared_xyz(traj, atoms, components, *, align_indices=None):
    """Apply PBC reconstruction, then optional rigid display alignment."""
    xyz = _make_protein_whole(traj.xyz, traj.unitcell_vectors, atoms, components)
    if align_indices:
        xyz = _rigid_align(xyz, list(align_indices))
    return xyz


def _trajectory_as_multimodel_pdb(
    gro_path: Path,
    xtc_path: Path,
    stride: int,
    max_frames: int = 80,
    reline_chloride_count: int = 0,
) -> tuple[str, int]:
    """Convert XTC to a PBC-safe, protein-aligned multi-model PDB for 3Dmol."""
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

    parsed = mv._parse_gro_atoms(gro_path)
    components = mv._gif_component_map(gro_path, reline_chloride_count)
    if len(parsed) != traj.n_atoms or len(components) != traj.n_atoms:
        raise RuntimeError("Trajectory atom count does not match the GRO topology")

    protein_indices = [i for i, component in enumerate(components) if component == "protein"]
    bb_indices = [
        i for i in protein_indices
        if str(parsed[i]["name"]).strip().upper().startswith("BB")
    ]
    xyz = _prepared_xyz(traj, parsed, components, align_indices=bb_indices or protein_indices)

    atoms = []
    intrinsic_cl_seen = 0
    for index, atom in enumerate(parsed, start=1):
        original_resn = str(atom["resn"]).strip() or "MOL"
        atom_name = str(atom["name"]).strip() or "B"
        upper = original_resn.upper()
        namesafe = original_resn[:3]
        if upper in mv.CHOLINE_RESN:
            namesafe = "CHO"
        elif upper in mv.UREA_RESN:
            namesafe = "URE"
        elif upper in mv.SORBITOL_RESN:
            namesafe = "SOR"
        elif upper in mv.WATER_RESN:
            namesafe = "WAT"
        elif upper in {"CL", "CL-"}:
            if intrinsic_cl_seen < max(0, int(reline_chloride_count)):
                namesafe = "RCL"
                intrinsic_cl_seen += 1
            else:
                namesafe = "ICL"
        elif upper in {"NA", "NA+"}:
            namesafe = "INA"
        atoms.append(
            {
                "serial": index,
                "resid": int(atom["resid"]),
                "resn": namesafe,
                "name": atom_name,
            }
        )

    blocks: list[str] = []
    for frame_index, frame in enumerate(xyz, start=1):
        blocks.append(f"MODEL     {frame_index:4d}")
        for atom, coord in zip(atoms, frame):
            x, y, z = (float(value) * 10.0 for value in coord)
            blocks.append(
                f"ATOM  {atom['serial']:5d} {atom['name'][:4]:>4} {atom['resn'][:3]:>3} A"
                f"{atom['resid'] % 10000:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          C"
            )
        blocks.append("ENDMDL")
    blocks.append("END")
    return "\n".join(blocks) + "\n", traj.n_frames


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
    """Render a PBC-safe protein-first GIF from real MD coordinates."""
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
    detailed = [vo._detailed_component(atom, broad) for atom, broad in zip(atoms, broad_components)]

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

    if show_solvent:
        water = [i for i in visible_indices if detailed[i] == "water"]
        other_solvent = [
            i for i in visible_indices
            if broad_components[i] == "solvent" and detailed[i] != "water"
        ]
        keep_water = set(vo._stable_subsample(water, 500))
        keep_other = set(vo._stable_subsample(other_solvent, 900))
        visible_indices = [
            i for i in visible_indices
            if broad_components[i] != "solvent" or i in keep_water or i in keep_other
        ]

    protein_indices = [i for i, broad in enumerate(broad_components) if broad == "protein"]
    bb_indices = [
        i for i in protein_indices
        if str(atoms[i]["name"]).strip().upper().startswith("BB")
    ]
    align_indices = (bb_indices or protein_indices) if show_protein else []
    prepared_xyz = _prepared_xyz(traj, atoms, broad_components, align_indices=align_indices)
    backbone_edges = vo._backbone_edges(atoms, prepared_xyz[0], bb_indices) if show_protein else []

    center_indices = (bb_indices or protein_indices) if show_protein and protein_indices else visible_indices
    center = prepared_xyz[:, center_indices, :].mean(axis=(0, 1))
    coords = prepared_xyz - center

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
        for left, right in backbone_edges:
            if left in projected and right in projected:
                draw.line(
                    [screen(projected[left]), screen(projected[right])],
                    fill="#168D9B",
                    width=8,
                )

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
    """Install PBC-safe trajectory conversion and GIF rendering."""
    mv._trajectory_as_multimodel_pdb = _trajectory_as_multimodel_pdb
    mv.generate_trajectory_gif = generate_trajectory_gif
