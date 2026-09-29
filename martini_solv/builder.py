"""Headless Martini 3 builder. External tools are invoked with argv, never a shell."""
from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
import urllib.request
from pathlib import Path

from .models import BuildConfig, reline_counts


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
        done = subprocess.run(args, cwd=cwd, input=stdin, text=True, stdout=fh, stderr=subprocess.STDOUT, check=False)
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


def _molecule_type(itp: Path) -> str:
    inside = False
    for line in itp.read_text().splitlines():
        line = line.split(";", 1)[0].strip()
        if line.startswith("["):
            inside = bool(re.match(r"\[\s*moleculetype\s*\]", line, re.I))
        elif inside and line:
            return line.split()[0]
    raise ValueError(f"No [ moleculetype ] in {itp}")


def _protein_net_charge(work: Path) -> float:
    total = 0.0
    for name in (line for line in work.glob("*.itp") if line.name.startswith("molecule_")):
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


def _topology(work: Path, protein: list[tuple[str, int]], species: list[tuple[str, int]]) -> None:
    generated = sorted(work.glob("*.itp"))
    protein_itps = [p for p in generated if p.name not in FF_NAMES and p.name not in {"choline.itp", "urea.itp"}]
    if not protein_itps:
        raise RuntimeError("Martinize2 produced no protein .itp")
    includes = [*FF_NAMES]
    if (work / "choline.itp").exists():
        includes += ["choline.itp", "urea.itp"]
    includes += [p.name for p in protein_itps]
    rows = [f'#include "{name}"' for name in includes]
    rows += ["", "[ system ]", "Protein in Martini 3 solvent", "", "[ molecules ]"]
    rows += [f"{name:<16} {count}" for name, count in protein + species if count]
    (work / "system.top").write_text("\n".join(rows) + "\n")


def _map_solute(work: Path, name: str, smiles: str, log: Path) -> tuple[Path, Path]:
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for {name}")
    # Martini Mapper supports this command form in MartiniSurf; inspect outputs explicitly.
    if not shutil.which("martini_mapper"):
        raise RuntimeError("SMILES molecules require the martini_mapper executable; see README")
    out = work / f"mapping_{name}"
    out.mkdir()
    run(["martini_mapper", name, smiles, "--no-xtb", "--out-dir", str(out)], work, log)
    gro = out / f"{name}.gro"
    itp = out / f"{name}.itp"
    if not gro.is_file() or not itp.is_file():
        raise RuntimeError(f"Martini Mapper did not produce {gro.name} and {itp.name}")
    shutil.copy2(itp, work / itp.name)
    return gro, work / itp.name


def build(pdb: Path, output: Path, config: BuildConfig) -> Path:
    config.validate()
    pdb = pdb.resolve()
    if not pdb.is_file() or pdb.suffix.lower() != ".pdb":
        raise ValueError("Provide an existing protein PDB file")
    if not any(row.startswith("ATOM  ") for row in pdb.read_text(errors="replace").splitlines()):
        raise ValueError("PDB does not contain protein ATOM records")
    for program in ("martinize2", "gmx", "insane" if config.solvent == "water" else "gmx"):
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
        sources = download_models(work, config.solvent == "reline")
        dssp = shutil.which("mkdssp") or shutil.which("dssp")
        if not dssp:
            raise RuntimeError("mkdssp/dssp is required for protein secondary structure")
        run(["martinize2", "-f", "input.pdb", "-x", "protein_cg.pdb", "-o", "protein.top",
             "-ff", "martini3001", "-dssp", dssp, "-ignh", "-elastic", "-ef", "700", "-el", "0.5", "-eu", "0.9"], work, log)
        protein = _molecules(work / "protein.top")
        if config.solvent == "reline" and abs(_protein_net_charge(work)) > 0.001:
            raise ValueError("DES mode currently requires an electrically neutral protein; no counterion correction is implemented")
        run(["gmx", "editconf", "-f", "protein_cg.pdb", "-o", "boxed.gro", "-c", "-d",
             str(config.box_distance_nm), "-bt", "cubic"], work, log)
        templates = []
        for spec in config.solutes:
            gro, itp = _map_solute(work, spec.name, spec.smiles, log)
            templates.append((spec, gro, itp))
        current = work / "boxed.gro"
        for index, (spec, gro, _) in enumerate(templates):
            current = _insert(work, current, gro, spec.count, spec.name, config.seed + index, log)
        species = [(_molecule_type(itp), spec.count) for spec, _, itp in templates]
        composition = {}
        if config.solvent == "water":
            run(["insane", "-f", str(current), "-o", "system.gro", "-p", "insane.top", "-pbc", "cubic",
                 "-d", "0", "-sol", "W", "-salt", str(config.salt_m), "-charge", "auto"], work, log)
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
            counts = reline_counts(math.prod(box) ** (1 / 3), config.water_fraction, config.des_pairs_per_nm3)
            # The published choline/urea coordinates are packaged in the cited DES model repository.
            (work / "chloride.gro").write_text("Chloride\n1\n    1CL     CL    1   0.000   0.000   0.000\n   1.00000   1.00000   1.00000\n")
            (work / "water.gro").write_text("Water\n1\n    1W       W    1   0.000   0.000   0.000\n   1.00000   1.00000   1.00000\n")
            for index, (name, template) in enumerate((
                ("CHOL", "choline.gro"), ("UREA", "urea.gro"), ("CL", "chloride.gro"), ("W", "water.gro")
            )):
                count = int(counts[name])
                current = _insert(work, current, work / template, count, name, config.seed + 100 + index, log)
                species.append((name, count))
            shutil.copy2(current, work / "system.gro")
            composition = counts
        _topology(work, protein, species)
        (work / "minimization.mdp").write_text(
            "integrator = steep\nnsteps = 5000\nemtol = 100\nemstep = 0.01\n"
            "cutoff-scheme = Verlet\nnstlist = 20\nrlist = 1.1\n"
            "vdwtype = cut-off\nvdw-modifier = Potential-shift-verlet\nrvdw = 1.1\n"
            "coulombtype = reaction-field\nrcoulomb = 1.1\nepsilon-r = 15\nepsilon-rf = inf\n"
        )
        report = {"config": json.loads(config.to_json()), "sources": sources,
                  "composition": composition, "protein_molecules": protein, "status": "grompp pending"}
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
