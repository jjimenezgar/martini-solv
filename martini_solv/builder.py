"""Headless Martini 3 builder. External tools are invoked with argv, never a shell."""
from __future__ import annotations

import importlib.util
import json
import math
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

from .models import BuildConfig, chcl_sorbitol_counts, reline_counts


FF_COMMIT = "784591ebdc91d762ed4df986c4650546c938f776"
DES_COMMIT = "c648037a5517d3a6e2ff5a83620cb11a98395030"
FF_NAMES = (
    "martini_v3.0.0.itp",
    "martini_v3.0.0_ions_v1.itp",
    "martini_v3.0.0_solvents_v1.itp",
)
DES_FILES = ("itps/choline.itp", "itps/urea.itp", "gro-files/choline.gro", "gro-files/urea.gro")


def run(args: list[str], cwd: Path, log: Path, stdin: str | None = None) -> None:
    with log.open("a") as fh:
        fh.write("$ " + " ".join(args) + "\n")
        fh.flush()
        # Packmol rewinds its input; a Python subprocess pipe is not seekable.
        input_path = cwd / "packmol.inp"
        if stdin is not None:
            input_path.write_text(stdin)
        with input_path.open() if stdin is not None else open("/dev/null") as input_file:
            done = subprocess.run(args, cwd=cwd, stdin=input_file, text=True,
                                  stdout=fh, stderr=subprocess.STDOUT, check=False)
        if done.returncode:
            raise RuntimeError(f"Command failed (exit {done.returncode}): {args[0]}. See {log}")


def download_models(work: Path, include_des: bool) -> dict[str, str]:
    """Fetch published files at immutable GitHub commits; fail if unreachable."""
    urls = {}
    for name in FF_NAMES:
        url = ("https://raw.githubusercontent.com/marrink-lab/martini-forcefields/"
               f"{FF_COMMIT}/martini_forcefields/regular/v3.0.0/gmx_files/{name}")
        (work / name).write_bytes(urllib.request.urlopen(url, timeout=30).read())
        urls[name] = url
    if include_des:
        for name in DES_FILES:
            url = f"https://raw.githubusercontent.com/vainikanpete/martini3-DES-models/{DES_COMMIT}/{name}"
            (work / Path(name).name).write_bytes(urllib.request.urlopen(url, timeout=30).read())
            urls[Path(name).name] = url
    return urls


def _molecules(top: Path) -> list[tuple[str, int]]:
    inside, entries = False, []
    for line in top.read_text().splitlines():
        line = line.split(";", 1)[0].strip()
        if line.startswith("["):
            inside = bool(re.match(r"\[\s*molecules\s*\]", line, re.I))
        elif inside and line and not line.startswith("#"):
            name, count, *_ = line.split()
            entries.append((name, int(count)))
    if not entries:
        raise ValueError(f"No [ molecules ] entries in {top}")
    return entries


def _gro_count(gro: Path) -> int:
    return int(gro.read_text().splitlines()[1].strip())


def _gro_atom_identities(gro: Path, limit: int | None = None) -> list[tuple[str, str]]:
    """Return (residue, atom) identities in coordinate order."""
    lines = gro.read_text(errors="replace").splitlines()
    count = int(lines[1].strip())
    if limit is not None:
        count = min(count, int(limit))
    return [(row[5:10].strip(), row[10:15].strip()) for row in lines[2:2 + count]]


def _verify_existing_coordinates_preserved(before: Path, after: Path) -> None:
    """Ensure a solvation step did not silently discard inserted solute molecules."""
    before_ids = _gro_atom_identities(before)
    after_ids = _gro_atom_identities(after, len(before_ids))
    if len(after_ids) < len(before_ids) or after_ids != before_ids:
        raise RuntimeError(
            "Solvation did not preserve the complete protein/free-molecule coordinate block"
        )


def _reline_coordinate_composition(gro: Path) -> dict[str, int]:
    """Count DES components and neutralizing ions from GRO coordinates."""
    counts = {"CHOL": 0, "UREA": 0, "SOR": 0, "CL": 0, "NA": 0, "W": 0}
    for residue, atom in _gro_atom_identities(gro):
        res = residue.upper()
        name = atom.upper()
        if res in {"CHOL", "CHO"} and name == "N1":
            counts["CHOL"] += 1
        elif res in {"UREA", "URE"} and name == "N1":
            counts["UREA"] += 1
        elif res == "SOR" and name == "S1":
            counts["SOR"] += 1
        elif res in {"CL", "CL-"} or name in {"CL", "CL-"}:
            counts["CL"] += 1
        elif res in {"NA", "NA+"} or name in {"NA", "NA+"}:
            counts["NA"] += 1
        elif res in {"W", "WF", "SW", "TW", "SOL"} and name in {"W", "WF", "SW", "TW", "OW"}:
            counts["W"] += 1
    return counts


def _verify_reline_composition(gro: Path, expected: dict[str, float | int]) -> dict[str, int]:
    actual = _reline_coordinate_composition(gro)
    for name in ("CHOL", "UREA", "SOR", "CL", "NA", "W"):
        if actual[name] != int(expected.get(name, 0)):
            raise RuntimeError(
                f"Reline composition mismatch for {name}: requested {int(expected.get(name, 0))}, "
                f"found {actual[name]} in system.gro"
            )
    return actual


def _normalize_insane_ions(gro: Path) -> None:
    """INSANE still writes Martini 2 ion labels; Martini 3 uses NA and CL."""
    lines = gro.read_text().splitlines()
    count = int(lines[1])
    for index in range(2, 2 + count):
        row = lines[index]
        residue, atom = row[5:10].strip(), row[10:15].strip()
        if residue in {"NA+", "CL-"} or atom in {"NA+", "CL-"}:
            residue = {"NA+": "NA", "CL-": "CL"}.get(residue, residue)
            atom = {"NA+": "NA", "CL-": "CL"}.get(atom, atom)
            lines[index] = row[:5] + f"{residue:>5}{atom:>5}" + row[15:]
    gro.write_text("\n".join(lines) + "\n")


def _clean_protein_pdb(source: Path, target: Path) -> None:
    """Keep the first model's protein atoms; exclude waters, ligands and alternate B sites."""
    lines = []
    saw_model = False
    for raw in source.read_text(errors="replace").splitlines():
        if raw.startswith("MODEL "):
            if saw_model:
                break
            saw_model = True
        elif raw.startswith("ENDMDL"):
            break
        elif raw.startswith("ATOM  ") and raw[16:17] in {" ", "A"}:
            lines.append(raw[:16] + " " + raw[17:])
        elif raw.startswith("TER"):
            lines.append(raw)
    if not any(row.startswith("ATOM  ") for row in lines):
        raise ValueError("No protein ATOM records after PDB cleanup")
    target.write_text("\n".join(lines + ["END"]) + "\n")


def _molecule_type(itp: Path) -> str:
    inside = False
    for line in itp.read_text().splitlines():
        line = line.split(";", 1)[0].strip()
        if line.startswith("["):
            inside = bool(re.match(r"\[\s*moleculetype\s*\]", line, re.I))
        elif inside and line:
            return line.split()[0]
    raise ValueError(f"No [ moleculetype ] in {itp}")


def _name_molecule_type(itp: Path, destination: Path, name: str) -> None:
    """Martini Mapper calls every generated type 'res'; make it unique per species."""
    inside, replaced, lines = False, False, []
    for raw in itp.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.startswith("["):
            inside = bool(re.match(r"\[\s*moleculetype\s*\]", line, re.I))
        elif inside and line and not replaced:
            parts = raw.split(maxsplit=1)
            raw = name + (" " + parts[1] if len(parts) > 1 else "")
            replaced = True
            inside = False
        lines.append(raw)
    if not replaced:
        raise ValueError(f"No molecule type found in {itp}")
    destination.write_text("\n".join(lines) + "\n")


def _itp_net_charge(name: Path) -> float:
    total = 0.0
    inside = False
    for raw in name.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.startswith("["):
            inside = bool(re.match(r"\[\s*atoms\s*\]", line, re.I))
        elif inside and line:
            fields = line.split()
            if len(fields) >= 7:
                total += float(fields[6])
    return total


def _assign_itp_net_charge(itp: Path, net_charge: int, bead_index: int | None) -> None:
    """Adjust one Martini bead so the molecule has the requested integer net charge."""
    target = int(net_charge)
    if target == 0:
        return
    if bead_index is None or int(bead_index) < 1:
        raise ValueError("A charged free molecule needs a charged bead selection")

    current_total = _itp_net_charge(itp)
    delta = float(target) - current_total
    lines = itp.read_text().splitlines()
    inside = False
    updated = False
    output: list[str] = []
    for raw in lines:
        bare = raw.split(";", 1)[0].strip()
        if bare.startswith("["):
            inside = bool(re.match(r"\[\s*atoms\s*\]", bare, re.I))
            output.append(raw)
            continue
        if inside and bare and not bare.startswith("#"):
            fields = bare.split()
            try:
                atom_id = int(fields[0])
            except (ValueError, IndexError):
                atom_id = -1
            if atom_id == int(bead_index):
                if len(fields) < 7:
                    raise ValueError("Generated ITP atom line does not contain a charge column")
                fields[6] = f"{float(fields[6]) + delta:.6f}"
                comment = ""
                if ";" in raw:
                    comment = " ;" + raw.split(";", 1)[1]
                raw = "  " + "  ".join(fields) + comment
                updated = True
        output.append(raw)

    if not updated:
        raise ValueError(f"Charged bead {bead_index} was not found in {itp.name}")
    itp.write_text("\n".join(output) + "\n")
    if abs(_itp_net_charge(itp) - target) > 1e-6:
        raise RuntimeError(f"Could not assign net charge {target:+d} to {itp.name}")


def _protein_net_charge(work: Path) -> float:
    return sum(_itp_net_charge(path) for path in work.glob("*.itp") if path.name not in FF_NAMES
               and path.name not in {"choline.itp", "urea.itp"})


def _verify_free_molecule_topology_counts(
    system_top: Path,
    templates: list[tuple[object, Path, Path]],
) -> list[dict[str, object]]:
    """Verify requested free-molecule copy counts are present in system.top."""
    top_counts: dict[str, int] = {}
    for molecule_name, count in _molecules(system_top):
        top_counts[molecule_name] = top_counts.get(molecule_name, 0) + int(count)

    verified: list[dict[str, object]] = []
    for spec, _gro, itp in templates:
        molecule_type = _molecule_type(itp)
        requested = int(spec.count)
        actual = int(top_counts.get(molecule_type, 0))
        if actual != requested:
            raise RuntimeError(
                f"Free molecule {spec.name}: requested {requested} copies but system.top contains {actual}"
            )
        verified.append({
            "name": spec.name,
            "molecule_type": molecule_type,
            "requested": requested,
            "included": actual,
            "beads_per_molecule": _gro_count(_gro),
        })
    return verified


def _insert(work: Path, current: Path, template: Path, count: int, name: str, seed: int, log: Path) -> Path:
    if count == 0:
        return current
    out = work / f"with_{name}.gro"
    old_atoms = _gro_count(current)
    run(["gmx", "insert-molecules", "-f", str(current), "-ci", str(template), "-nmol", str(count),
         "-o", str(out), "-try", "1000", "-radius", "0.21", "-seed", str(seed)], work, log)
    expected = old_atoms + count * _gro_count(template)
    if _gro_count(out) != expected:
        raise RuntimeError(f"Only some {name} molecules fit; increase box size or reduce density")
    return out


def _pack_reline(work: Path, boxed: Path, templates: list[tuple[Path, int]], seed: int, log: Path) -> None:
    """Pack all DES components at once, keeping the protein fixed in the box."""
    box = [float(v) for v in boxed.read_text().splitlines()[-1].split()[:3]]
    if len(box) != 3 or min(box) <= 0:
        raise ValueError("Invalid GRO box")
    blocks = ["tolerance 2.5", "filetype pdb", "output packed.pdb", f"seed {seed}", ""]
    run(["gmx", "editconf", "-f", str(boxed), "-o", "boxed.pdb"], work, log)
    blocks += ["structure boxed.pdb", "  number 1", "  resnumbers 1",
               "  fixed 0. 0. 0. 0. 0. 0.", "end structure", ""]
    for index, (template, count) in enumerate(templates):
        if not count:
            continue
        pdb_name = f"pack_template_{index}.pdb"
        run(["gmx", "editconf", "-f", str(template), "-o", pdb_name], work, log)
        bounds = " ".join(f"{length * 10:.3f}" for length in box)
        blocks += [f"structure {pdb_name}", f"  number {count}", "  resnumbers 3",
                   f"  inside box 0. 0. 0. {bounds}", "end structure", ""]
    run(["packmol"], work, log, stdin="\n".join(blocks) + "\n")
    run(["gmx", "editconf", "-f", "packed.pdb", "-o", "system.gro", "-box", *map(str, box)], work, log)
    expected = _gro_count(boxed) + sum(_gro_count(template) * count for template, count in templates)
    if _gro_count(work / "system.gro") != expected:
        raise RuntimeError("Packmol output atom count does not match the requested composition")


def _topology_prefix_from_martinize(work: Path) -> list[str]:
    """Preserve martinize2 preprocessor/include ordering up to [ system ].

    This is important for GōMartini: martinize2 can emit #define directives and
    auxiliary atom-type/nonbonded includes that must appear before Protein.itp.
    Reconstructing the topology from a sorted list of ITP files loses that order
    and leads to errors such as "Atomtype Protein_1 not found".
    """
    top = work / "protein.top"
    if not top.is_file():
        raise RuntimeError("martinize2 did not produce protein.top")

    prefix: list[str] = []
    for raw in top.read_text(errors="replace").splitlines():
        stripped = raw.strip()
        if stripped.startswith("[") and re.match(r"\[\s*(system|molecules)\s*\]", stripped, re.I):
            break
        prefix.append(raw)
    while prefix and not prefix[-1].strip():
        prefix.pop()
    return prefix


def _normalize_martinize_includes(lines: list[str]) -> list[str]:
    """Map martinize2's generic Martini include to the pinned local FF file.

    Some martinize2 releases emit '#include "martini.itp"' even when the
    selected force field is martini3001. MartiniSolv downloads the pinned
    Martini 3 force field as martini_v3.0.0.itp, so leaving the generic include
    untouched makes grompp fail with 'Topology include file martini.itp not found'.
    """
    normalized: list[str] = []
    seen_includes: set[str] = set()
    for raw in lines:
        match = re.match(r'(\s*#include\s+["<])([^">]+)([">].*)', raw)
        if match:
            target = match.group(2)
            name = Path(target).name
            if name == "martini.itp":
                target = FF_NAMES[0]
                name = FF_NAMES[0]
                raw = f'{match.group(1)}{target}{match.group(3)}'
            if name in seen_includes:
                continue
            seen_includes.add(name)
        normalized.append(raw)
    return normalized


def _included_itps(lines: list[str]) -> set[str]:
    names: set[str] = set()
    for raw in lines:
        match = re.match(r'\s*#include\s+["<]([^">]+)[">]', raw)
        if match:
            names.add(Path(match.group(1)).name)
    return names


def _insert_after_forcefield(lines: list[str], additions: list[str]) -> list[str]:
    if not additions:
        return lines
    insert_at = 0
    for index, raw in enumerate(lines):
        match = re.match(r'\s*#include\s+["<]([^">]+)[">]', raw)
        if match and Path(match.group(1)).name == FF_NAMES[0]:
            insert_at = index + 1
            break
    return lines[:insert_at] + additions + lines[insert_at:]


def _relocate_go_includes(lines: list[str], work: Path, go_enabled: bool) -> list[str]:
    """Place Gō global atom-type/nonbonded directives immediately after the FF.

    GROMACS requires [atomtypes]/[nonbond_params] to be read before any
    [moleculetype]. Martini 3 ion/solvent ITPs already contain molecule types,
    so the Gō files must precede those includes as well as Protein.itp.
    """
    go_names = ("go_atomtypes.itp", "go_nbparams.itp")
    stripped: list[str] = []
    for raw in lines:
        match = re.match(r'\s*#include\s+["<]([^">]+)[">]', raw)
        if match and Path(match.group(1)).name in go_names:
            continue
        stripped.append(raw)

    if not go_enabled:
        return stripped

    missing = [name for name in go_names if not (work / name).is_file()]
    if missing:
        raise RuntimeError(
            "GōMartini was requested but martinize2 did not generate: "
            + ", ".join(missing)
        )

    go_lines = [f'#include "{name}"' for name in go_names]
    return _insert_after_forcefield(stripped, go_lines)


def _topology(
    work: Path,
    protein: list[tuple[str, int]],
    species: list[tuple[str, int]],
    *,
    go_enabled: bool = False,
) -> None:
    prefix = _normalize_martinize_includes(_topology_prefix_from_martinize(work))

    # Some martinize2 releases add this define themselves; keep it exactly where
    # they placed it. If a release generated Gō files but omitted the define,
    # add it before all includes.
    if go_enabled and not any(re.match(r"\s*#define\s+GO_VIRT\b", row) for row in prefix):
        prefix.insert(0, "#define GO_VIRT")

    included = _included_itps(prefix)

    # Ensure the Martini 3 force field is present first.
    if FF_NAMES[0] not in included:
        prefix = [f'#include "{FF_NAMES[0]}"', *prefix]

    # Gō atomtypes/nonbonded parameters are global directives and MUST be read
    # before ions, solvents, Protein.itp or any other moleculetype definition.
    prefix = _relocate_go_includes(prefix, work, go_enabled)
    included = _included_itps(prefix)

    auxiliary = []
    for name in FF_NAMES[1:]:
        if name not in included:
            auxiliary.append(f'#include "{name}"')
            included.add(name)
    for name in ("choline.itp", "urea.itp", "sorbitol.itp"):
        if (work / name).exists() and name not in included:
            auxiliary.append(f'#include "{name}"')
            included.add(name)
    if auxiliary:
        last_global = -1
        global_names = {FF_NAMES[0], "go_atomtypes.itp", "go_nbparams.itp"}
        for index, raw in enumerate(prefix):
            match = re.match(r'\s*#include\s+["<]([^">]+)[">]', raw)
            if match and Path(match.group(1)).name in global_names:
                last_global = index
        insert_at = last_global + 1 if last_global >= 0 else 0
        prefix = prefix[:insert_at] + auxiliary + prefix[insert_at:]

    # protein.top normally includes the protein ITP. Add only the actual
    # molecule ITP as a fallback; do not blindly include Gō auxiliary files,
    # because their order/ifdef placement is controlled by martinize2.
    for molecule_name, _ in protein:
        candidate = f"{molecule_name}.itp"
        if (work / candidate).is_file() and candidate not in included:
            prefix.append(f'#include "{candidate}"')
            included.add(candidate)

    # Free SMILES-derived species are independent molecule types and can be
    # included after the martinize2-controlled protein block.
    for molecule_name, _ in species:
        candidate = f"{molecule_name}.itp"
        if (work / candidate).is_file() and candidate not in included:
            prefix.append(f'#include "{candidate}"')
            included.add(candidate)

    rows = [*prefix, "", "[ system ]", "Protein in Martini 3 solvent", "", "[ molecules ]"]
    rows += [f"{name:<16} {count}" for name, count in protein + species if count]
    (work / "system.top").write_text("\n".join(rows) + "\n")

def _parse_itp_atoms(itp: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    section = ""
    for raw in itp.read_text(errors="replace").splitlines():
        bare = raw.split(";", 1)[0].strip()
        if bare.startswith("[") and "]" in bare:
            section = bare.strip("[]").strip().lower()
            continue
        if section != "atoms" or not bare or bare.startswith("#"):
            continue
        fields = bare.split()
        if len(fields) < 7:
            continue
        try:
            atom_id = int(fields[0])
            resid = int(fields[2])
        except ValueError:
            continue
        atoms.append({
            "id": atom_id,
            "type": fields[1],
            "resid": resid,
            "resname": fields[3],
            "name": fields[4],
        })
    if not atoms:
        raise ValueError(f"No [ atoms ] entries found in {itp.name}")
    return atoms


def _parse_itp_edges(itp: Path) -> list[tuple[int, int, float | None]]:
    edges: list[tuple[int, int, float | None]] = []
    section = ""
    for raw in itp.read_text(errors="replace").splitlines():
        bare = raw.split(";", 1)[0].strip()
        if bare.startswith("[") and "]" in bare:
            section = bare.strip("[]").strip().lower()
            continue
        if section not in {"bonds", "constraints"} or not bare or bare.startswith("#"):
            continue
        fields = bare.split()
        if len(fields) < 2:
            continue
        try:
            left, right = int(fields[0]), int(fields[1])
        except ValueError:
            continue
        length = None
        # Common GROMACS bond/constraint forms place equilibrium length after funct.
        for candidate in fields[3:]:
            try:
                value = float(candidate)
            except ValueError:
                continue
            if 0.15 <= value <= 1.5:
                length = value
                break
        edges.append((left, right, length))
    return edges


def _gro_from_itp(itp: Path, destination: Path, name: str) -> Path:
    """Generate a deterministic non-overlapping CG template from an ITP connectivity graph."""
    atoms = _parse_itp_atoms(itp)
    ids = {int(atom["id"]) for atom in atoms}
    edges = _parse_itp_edges(itp)
    adjacency: dict[int, list[tuple[int, float | None]]] = {atom_id: [] for atom_id in ids}
    for left, right, length in edges:
        if left in adjacency and right in adjacency:
            adjacency[left].append((right, length))
            adjacency[right].append((left, length))

    if len(atoms) > 1 and not edges:
        raise ValueError(
            f"{itp.name} has multiple beads but no [ bonds ] or [ constraints ]; "
            "a starting geometry cannot be generated safely from this ITP alone"
        )

    coords: dict[int, tuple[float, float, float]] = {}
    root = int(atoms[0]["id"])
    coords[root] = (0.0, 0.0, 0.0)
    queue = [root]
    directions = [
        (1.0, 0.0, 0.0), (-0.5, 0.866, 0.0), (-0.5, -0.866, 0.0),
        (0.0, 0.5, 0.866), (0.0, -0.5, 0.866), (0.0, 0.5, -0.866),
    ]
    placed_edges = 0
    while queue:
        parent = queue.pop(0)
        px, py, pz = coords[parent]
        children = [(child, length) for child, length in adjacency[parent] if child not in coords]
        for offset, (child, length) in enumerate(children):
            distance = float(length) if length is not None else 0.47
            dx, dy, dz = directions[(placed_edges + offset) % len(directions)]
            trial = (px + distance * dx, py + distance * dy, pz + distance * dz)
            # Nudge if a branched/cyclic topology would place two beads too close.
            shift = 0
            while any(
                sum((trial[i] - other[i]) ** 2 for i in range(3)) < 0.16 ** 2
                for other in coords.values()
            ):
                shift += 1
                trial = (trial[0], trial[1] + 0.12 * shift, trial[2] + 0.07 * shift)
            coords[child] = trial
            queue.append(child)
        placed_edges += max(1, len(children))

    missing = ids.difference(coords)
    if missing:
        raise ValueError(
            f"{itp.name} contains disconnected beads ({', '.join(map(str, sorted(missing)))}); "
            "a single-molecule GRO cannot be generated safely"
        )

    mins = [min(point[i] for point in coords.values()) for i in range(3)]
    padding = 0.5
    shifted = {
        atom_id: tuple(value - mins[i] + padding for i, value in enumerate(point))
        for atom_id, point in coords.items()
    }
    maxs = [max(point[i] for point in shifted.values()) + padding for i in range(3)]
    atom_by_id = {int(atom["id"]): atom for atom in atoms}
    lines = [f"{name} generated from ITP", str(len(atoms))]
    for serial, atom_id in enumerate(sorted(atom_by_id), start=1):
        atom = atom_by_id[atom_id]
        x, y, z = shifted[atom_id]
        resname = str(atom["resname"])[:5] or name[:5]
        atom_name = str(atom["name"])[:5] or f"B{serial}"
        resid = int(atom["resid"]) % 100000
        lines.append(
            f"{resid:5d}{resname:<5}{atom_name:>5}{serial:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    lines.append("".join(f"{max(1.0, length):10.5f}" for length in maxs))
    destination.write_text("\n".join(lines) + "\n")
    return destination


def _prepare_uploaded_solute(work: Path, name: str, source_itp: Path) -> tuple[Path, Path]:
    destination_itp = work / f"{name}.itp"
    _name_molecule_type(source_itp, destination_itp, name)
    destination_gro = work / f"{name}.gro"
    _gro_from_itp(destination_itp, destination_gro, name)
    return destination_gro, destination_itp


def _map_solute(work: Path, name: str, smiles: str, log: Path) -> tuple[Path, Path]:
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {name}")
    # Martini Mapper supports this command form in MartiniSurf; inspect outputs explicitly.
    if importlib.util.find_spec("martini_mapper") is None:
        raise RuntimeError("SMILES molecules require the Martini Mapper module; recreate the Conda environment")
    out = work / f"mapping_{name}"
    out.mkdir()
    run([sys.executable, "-m", "martini_mapper", name, smiles, "--no-xtb", "--out-dir", str(out)], work, log)
    gro = out / f"{name}.gro"
    itp = out / f"{name}.itp"
    if not gro.is_file() or not itp.is_file():
        raise RuntimeError(f"Martini Mapper did not produce {gro.name} and {itp.name}")
    _name_molecule_type(itp, work / itp.name, name)
    return gro, work / itp.name


def build(pdb: Path, output: Path, config: BuildConfig) -> Path:
    config.validate()
    pdb = pdb.resolve()
    if not pdb.is_file() or pdb.suffix.lower() != ".pdb":
        raise ValueError("Provide an existing protein PDB file")
    if not any(row.startswith("ATOM  ") for row in pdb.read_text(errors="replace").splitlines()):
        raise ValueError("PDB does not contain protein ATOM records")
    for program in ("martinize2", "gmx", "insane" if config.solvent == "water" else "packmol"):
        if not shutil.which(program):
            raise RuntimeError(f"Missing executable: {program}")
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"Choose a new output directory: {output}")
    work = output.with_name(output.name + ".incomplete")
    if work.exists():
        raise FileExistsError(f"Previous incomplete build exists: {work}")
    work.mkdir(parents=True)
    log = work / "build.log"
    try:
        shutil.copy2(pdb, work / "input.pdb")
        _clean_protein_pdb(work / "input.pdb", work / "protein_clean.pdb")
        sources = download_models(work, config.solvent in {"reline", "chcl_sorbitol"})
        if config.solvent == "chcl_sorbitol":
            (work / "sorbitol.itp").write_text(
                "; Sorbitol Martini 3 model\n"
                "; Mapping: 3 P4 beads\n\n"
                "[ moleculetype ]\n; name   nrexcl\n  SOR      1\n\n"
                "[ atoms ]\n"
                "; id   type   resnr  residu  atom   cgnr   charge   mass\n"
                "  1    P4     1      SOR     S1       1      0.0  72.000\n"
                "  2    P4     1      SOR     S2       2      0.0  72.000\n"
                "  3    P4     1      SOR     S3       3      0.0  72.000\n\n"
                "[ constraints ]\n"
                "; i   j   funct   length\n"
                "  1   2   1       0.33412\n"
                "  2   3   1       0.32440\n"
            )
            (work / "sorbitol.gro").write_text(
                "Sorbitol Martini 3\n3\n"
                "    1SOR    S1    1   0.100   0.100   0.100\n"
                "    1SOR    S2    2   0.434   0.100   0.100\n"
                "    1SOR    S3    3   0.758   0.100   0.100\n"
                "   1.00000   1.00000   1.00000\n"
            )
            sources["sorbitol.itp"] = "user-provided Martini 3 sorbitol topology"
            sources["sorbitol.gro"] = "generated starting geometry from supplied constraints"
        if importlib.util.find_spec("mdtraj") is None:
            raise RuntimeError("mdtraj is required for protein secondary structure")
        martinize = [
            "martinize2", "-f", "protein_clean.pdb", "-x", "protein_cg.pdb", "-o", "protein.top",
            "-ff", "martini3001", "-name", config.molecule_name, "-maxwarn", str(config.maxwarn), "-ignh",
        ]
        merge = config.merge_chains.strip()
        if merge and "," in merge:
            martinize += ["-merge", merge]
        if config.position_restraints != "none":
            martinize += ["-p", config.position_restraints, "-pf", str(config.position_restraint_force)]
        if config.elastic:
            martinize += ["-elastic", "-ef", str(config.elastic_force), "-el", "0.5", "-eu", "0.9"]
        if config.go:
            martinize += ["-go", "-go-eps", str(config.go_eps)]
        if config.dssp:
            martinize += ["-dssp"]
        martinize += [str(token) for token in config.martinize_extra_args if str(token).strip()]

        # Match MartiniSurf's resilient DSSP behavior: retry once without DSSP
        # when the runtime DSSP setup is incompatible.
        try:
            run(martinize, work, log)
        except RuntimeError:
            if config.dssp and "-dssp" in martinize:
                retry = [token for token in martinize if token != "-dssp"]
                with log.open("a") as handle:
                    handle.write("\nDSSP-enabled martinize2 failed; retrying without -dssp.\n")
                run(retry, work, log)
                martinize = retry
            else:
                raise
        protein = _molecules(work / "protein.top")
        protein_charge = _protein_net_charge(work)
        rounded_protein_charge = int(round(protein_charge))
        if abs(protein_charge - rounded_protein_charge) > 0.01:
            raise ValueError(
                f"Protein net charge {protein_charge:.3f} is not close to an integer; "
                "automatic counterion neutralization would be ambiguous"
            )
        run(["gmx", "editconf", "-f", "protein_cg.pdb", "-o", "boxed.gro", "-c", "-d",
             str(config.box_distance_nm), "-bt", "cubic"], work, log)
        templates = []
        for spec in config.solutes:
            if spec.source == "upload":
                gro, itp = _prepare_uploaded_solute(
                    work, spec.name, Path(spec.template_itp)
                )
            else:
                gro, itp = _map_solute(work, spec.name, spec.smiles, log)
            _assign_itp_net_charge(itp, spec.net_charge, spec.charged_bead)
            templates.append((spec, gro, itp))
        current = work / "boxed.gro"
        species = [(_molecule_type(itp), spec.count) for spec, _, itp in templates]
        free_molecule_charge = sum(int(spec.net_charge) * int(spec.count) for spec, _, _ in templates)
        total_non_solvent_charge = rounded_protein_charge + free_molecule_charge
        composition = {
            "protein_net_charge": protein_charge,
            "free_molecule_net_charge": free_molecule_charge,
            "total_non_solvent_charge": total_non_solvent_charge,
        }
        if config.solvent == "water":
            for index, (spec, gro, _) in enumerate(templates):
                current = _insert(work, current, gro, spec.count, spec.name, config.seed + index, log)
            run(["insane", "-f", str(current), "-o", "system.gro", "-p", "insane.top", "-pbc", "cubic",
                 "-d", "0", "-sol", "W", "-salt", str(config.salt_m), "-charge",
                 str(total_non_solvent_charge)], work, log)
            _normalize_insane_ions(work / "system.gro")
            _verify_existing_coordinates_preserved(current, work / "system.gro")
            entries = _molecules(work / "insane.top")
            # INSANE's output may list protein and solutes; remove only the already recorded copies.
            for name, count in entries:
                if name in {"W", "NA", "CL", "NA+", "CL-"}:
                    canonical = {"NA+": "NA", "CL-": "CL"}.get(name, name)
                    species.append((canonical, count))
                    composition[canonical] = composition.get(canonical, 0) + count
        else:
            lines = current.read_text().splitlines()
            box = [float(v) for v in lines[-1].split()[:3]]
            if len(box) != 3 or min(box) <= 0:
                raise ValueError("Invalid GRO box")
            if config.solvent == "reline":
                counts = reline_counts(
                    math.prod(box) ** (1 / 3),
                    config.water_fraction,
                    config.reline_density_g_cm3,
                )
            else:
                counts = chcl_sorbitol_counts(
                    math.prod(box) ** (1 / 3),
                    config.water_fraction,
                    config.chcl_sorbitol_density_g_cm3,
                )
            # The dry DES is neutral before adding the protein. Add only the
            # extra ions required to neutralize the protein.
            neutralizing_na = max(0, -total_non_solvent_charge)
            neutralizing_cl = max(0, total_non_solvent_charge)
            counts["NA"] = neutralizing_na
            counts["CL"] = int(counts["CL"]) + neutralizing_cl
            counts["intrinsic_chloride"] = int(counts["CHOL"])
            counts["reline_chloride"] = int(counts["CHOL"])
            counts["neutralizing_NA"] = neutralizing_na
            counts["neutralizing_CL"] = neutralizing_cl
            counts["protein_net_charge"] = protein_charge
            counts["free_molecule_net_charge"] = free_molecule_charge
            counts["total_non_solvent_charge"] = total_non_solvent_charge

            (work / "chloride.gro").write_text("Chloride\n1\n    1CL     CL    1   0.000   0.000   0.000\n   1.00000   1.00000   1.00000\n")
            (work / "sodium.gro").write_text("Sodium\n1\n    1NA     NA    1   0.000   0.000   0.000\n   1.00000   1.00000   1.00000\n")
            (work / "water.gro").write_text("Water\n1\n    1W       W    1   0.000   0.000   0.000\n   1.00000   1.00000   1.00000\n")
            to_pack = [(gro, spec.count) for spec, gro, _ in templates]
            solvent_templates = [("CHOL", "choline.gro")]
            solvent_templates.append(
                ("UREA", "urea.gro") if config.solvent == "reline"
                else ("SOR", "sorbitol.gro")
            )
            solvent_templates.extend([
                ("CL", "chloride.gro"),
                ("NA", "sodium.gro"),
                ("W", "water.gro"),
            ])
            for name, template in solvent_templates:
                count = int(counts.get(name, 0))
                if count:
                    to_pack.append((work / template, count))
                    species.append((name, count))
            _pack_reline(work, current, to_pack, config.seed, log)
            actual_counts = _verify_reline_composition(work / "system.gro", counts)
            composition = {**counts, "coordinate_counts": actual_counts}
        _topology(work, protein, species, go_enabled=config.go)
        verified_free_molecules = _verify_free_molecule_topology_counts(
            work / "system.top", templates
        )
        (work / "minimization.mdp").write_text(
            "integrator = steep\nnsteps = 5000\nemtol = 100\nemstep = 0.01\n"
            "cutoff-scheme = Verlet\nnstlist = 20\nrlist = 1.1\n"
            "vdwtype = cut-off\nvdw-modifier = Potential-shift-verlet\nrvdw = 1.1\n"
            "coulombtype = reaction-field\nrcoulomb = 1.1\nepsilon-r = 15\nepsilon-rf = 0\n"
        )
        report = {"config": json.loads(config.to_json()), "sources": sources,
                  "composition": composition, "protein_molecules": protein,
                  "free_molecules": [
                      {
                          "name": spec.name,
                          "count": spec.count,
                          "copies_requested": spec.count,
                          "copies_included": next(
                              (
                                  row["included"]
                                  for row in verified_free_molecules
                                  if row["name"] == spec.name
                              ),
                              0,
                          ),
                          "beads_per_molecule": next(
                              (
                                  row["beads_per_molecule"]
                                  for row in verified_free_molecules
                                  if row["name"] == spec.name
                              ),
                              0,
                          ),
                          "source": spec.source,
                          "smiles": spec.smiles,
                          "net_charge": spec.net_charge,
                          "charged_bead": spec.charged_bead,
                      }
                      for spec in config.solutes
                  ],
                  "martinize_command": martinize, "status": "grompp pending"}
        (work / "manifest.json").write_text(json.dumps(report, indent=2))
        run(["gmx", "grompp", "-f", "minimization.mdp", "-c", "system.gro", "-p", "system.top",
             "-o", "minimization.tpr", "-maxwarn", "0"], work, log)
        report["status"] = "grompp passed; minimization not run"
        (work / "manifest.json").write_text(json.dumps(report, indent=2))
        work.rename(output)
        return output
    except Exception:
        # Keep failed work and logs for diagnosis; never expose it as a successful system.
        raise
