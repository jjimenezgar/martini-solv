"""Conservative atomistic protein preparation before martinize2.

The preparation mirrors MartiniSurf's robust protein-only path:
- keep only the first model and protein residues;
- choose a single alternate conformation by occupancy;
- remove hydrogens and non-protein heterogens;
- convert MSE to MET;
- repair missing side-chain heavy atoms with PDBFixer;
- never invent missing residues or bridge unresolved peptide breaks;
- leave terminal patches to martinize2;
- keep a JSON audit report.
"""
from __future__ import annotations

import io
import json
import math
from collections import defaultdict
from pathlib import Path

STANDARD = set(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
)
NUCLEIC = set("A C G U I DA DC DG DT DI DU".split())


def validate_pdb_coordinates(path: Path) -> None:
    count = 0
    for line in Path(path).read_text(errors="replace").splitlines():
        if line[:6].strip() not in {"ATOM", "HETATM"}:
            continue
        count += 1
        try:
            xyz = [float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]
        except ValueError as exc:
            raise ValueError(f"Malformed coordinates in {path}: {line[:27]}") from exc
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"Non-finite coordinates in {path}: {line[:27]}")
    if not count:
        raise ValueError(f"No atoms found in {path}")


def _label(key: tuple[int, str, str, str]) -> str:
    segment, chain, resid, icode = key
    return f"{chain or '_'}:{resid.strip()}{icode.strip()} (segment {segment})"


def _select_protein(infile: Path, report: dict[str, object]) -> str:
    groups: dict[tuple[int, str, str, str], list[str]] = defaultdict(list)
    segment = 0
    models = 0

    for raw in Path(infile).read_text(errors="replace").splitlines():
        record = raw[:6].strip()
        if record == "MODEL":
            models += 1
            if models > 1:
                break
        if record == "ENDMDL":
            break
        if record == "TER":
            segment += 1
        if record not in {"ATOM", "HETATM"}:
            continue
        if len(raw) < 54:
            raise ValueError(f"Truncated PDB atom record: {raw}")
        line = raw.ljust(80)
        try:
            xyz = [float(line[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))]
        except ValueError as exc:
            raise ValueError(f"Malformed coordinates: {line[:27]}") from exc
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"Non-finite input coordinates: {line[:27]}")
        key = (segment, line[21].strip(), line[22:26], line[26])
        groups[key].append(line)

    output: list[str] = []
    previous_segment: tuple[int, str] | None = None

    for key, lines in groups.items():
        names = {line[17:20].strip() for line in lines}
        atom_names = {line[12:16].strip() for line in lines}

        if not names <= STANDARD | {"MSE"}:
            if names <= NUCLEIC or (
                all(line.startswith("HETATM") for line in lines)
                and not {"N", "CA", "C"} <= atom_names
            ):
                report["removed_residues"].append(
                    {"residue": _label(key), "names": sorted(names)}
                )
                continue
            raise ValueError(
                f"Unsupported protein residue {sorted(names)} at {_label(key)}; "
                "provide a reviewed standard-residue model"
            )

        conformers: dict[str, list[float]] = defaultdict(list)
        for line in lines:
            altloc = line[16]
            if altloc != " ":
                try:
                    occupancy = float(line[54:60].strip() or "0")
                except ValueError as exc:
                    raise ValueError(f"Invalid occupancy at {_label(key)}") from exc
                if not math.isfinite(occupancy):
                    raise ValueError(f"Invalid occupancy at {_label(key)}")
                conformers[altloc].append(occupancy)

        chosen = (
            min(
                conformers,
                key=lambda code: (
                    -sum(conformers[code]) / len(conformers[code]),
                    code != "A",
                    code,
                ),
            )
            if conformers
            else " "
        )
        selected = [line for line in lines if line[16] in {" ", chosen}]
        if conformers:
            report["alternate_conformations"].append(
                {"residue": _label(key), "selected": chosen}
            )

        selected_names = {line[17:20].strip() for line in selected}
        if len(selected_names) != 1:
            raise ValueError(f"Ambiguous residue identity at {_label(key)}")
        resname = next(iter(selected_names))
        if resname == "MSE":
            report["replacements"].append(
                {"residue": _label(key), "from": "MSE", "to": "MET"}
            )

        seen: set[str] = set()
        prepared: list[str] = []
        for line in selected:
            atom = line[12:16].strip()
            element = line[76:78].strip().upper()
            if element in {"H", "D"} or (
                not element and atom.lstrip("0123456789").startswith(("H", "D"))
            ):
                continue
            if atom in seen:
                raise ValueError(f"Duplicate atom {atom} at {_label(key)}")
            seen.add(atom)

            line = "ATOM  " + line[6:16] + " " + line[17:]
            if resname == "MSE":
                line = line[:17] + "MET" + line[20:]
                if atom == "SE":
                    line = line[:12] + " SD " + line[16:76] + " S" + line[78:]
            prepared.append(line + "\n")

        missing_backbone = {"N", "CA", "C", "O"} - seen
        if missing_backbone:
            raise ValueError(
                f"Incomplete backbone at {_label(key)}: missing {sorted(missing_backbone)}. "
                "Review/model this region before martinization."
            )

        segment_key = key[:2]
        if previous_segment is not None and previous_segment != segment_key:
            output.append("TER\n")
        output.extend(prepared)
        previous_segment = segment_key

    if not output:
        raise ValueError("No supported protein residues found in the selected input")
    report["model"] = "first"
    return "".join(output) + "END\n"


def _check_backbone(fixer) -> None:
    from openmm import unit

    xyz = fixer.positions.value_in_unit(unit.angstrom)
    for chain in fixer.topology.chains():
        previous = None
        for residue in chain.residues():
            atoms = {atom.name: atom.index for atom in residue.atoms()}
            if not {"N", "CA", "C", "O"} <= atoms.keys():
                raise ValueError(
                    f"Incomplete backbone in chain {chain.id or '_'} residue {residue.id}"
                )
            if previous is not None:
                prev, prev_atoms = previous
                distance = math.dist(xyz[prev_atoms["C"]], xyz[atoms["N"]])
                try:
                    gap = int(residue.id) - int(prev.id) > 1
                except ValueError:
                    gap = False
                if distance > 2.0 or distance < 0.9 or gap:
                    raise ValueError(
                        "Unresolved peptide break in chain "
                        f"{chain.id or '_'}: {prev.id}{prev.insertionCode} -> "
                        f"{residue.id}{residue.insertionCode} (C-N {distance:.2f} A). "
                        "Review missing residues/chain boundaries before martinization."
                    )
            previous = (residue, atoms)


def prepare_protein_pdb(infile: Path, outfile: Path, report_path: Path) -> Path:
    report: dict[str, object] = {
        "status": "preparing",
        "input": str(infile),
        "protein_only": True,
        "removed_residues": [],
        "alternate_conformations": [],
        "replacements": [],
        "added_atoms": [],
        "missing_residues_added": False,
        "repair_seed": 42,
        "terminal_atoms_deferred": [],
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        from pdbfixer import PDBFixer
        from openmm import Platform
        from openmm.app import PDBFile

        selected = _select_protein(infile, report)
        fixer = PDBFixer(
            pdbfile=io.StringIO(selected),
            platform=Platform.getPlatformByName("CPU"),
        )
        identities = [
            (res.chain.id, res.id, res.insertionCode, res.name)
            for res in fixer.topology.residues()
        ]

        _check_backbone(fixer)

        # Deliberately do not model missing residues/loops.
        fixer.missingResidues = {}
        fixer.findMissingAtoms()
        for residue, atoms in fixer.missingAtoms.items():
            report["added_atoms"].append(
                {
                    "chain": residue.chain.id,
                    "resid": residue.id,
                    "icode": residue.insertionCode,
                    "resname": residue.name,
                    "atoms": [atom.name for atom in atoms],
                }
            )
        for residue, atoms in fixer.missingTerminals.items():
            report["terminal_atoms_deferred"].append(
                {
                    "chain": residue.chain.id,
                    "resid": residue.id,
                    "icode": residue.insertionCode,
                    "resname": residue.name,
                    "atoms": list(atoms),
                }
            )

        # martinize2 owns N/C-terminal patching; avoid adding OXT here.
        fixer.missingTerminals = {}
        fixer.addMissingAtoms(seed=42)
        fixer.findMissingAtoms()
        if any(fixer.missingAtoms.values()):
            raise ValueError("Protein repair left missing heavy atoms")

        repaired_identities = [
            (res.chain.id, res.id, res.insertionCode, res.name)
            for res in fixer.topology.residues()
        ]
        if identities != repaired_identities:
            raise ValueError("Protein repair changed residue identities or numbering")

        with outfile.open("w") as handle:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

        validate_pdb_coordinates(outfile)
        _check_backbone(PDBFixer(filename=str(outfile), platform=Platform.getPlatformByName("CPU")))
        report["status"] = "prepared"
        return outfile

    except ImportError as exc:
        report["status"] = "failed"
        report["error"] = (
            "Protein preparation requires pdbfixer and OpenMM. "
            "Reinstall MartiniSolv dependencies."
        )
        raise RuntimeError(report["error"]) from exc
    except Exception as exc:
        report["status"] = "failed"
        report["error"] = str(exc)
        raise
    finally:
        report_path.write_text(json.dumps(report, indent=2) + "\n")
