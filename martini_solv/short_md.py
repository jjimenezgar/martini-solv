"""Optional, bounded GROMACS minimization and NVT sanity check."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import zipfile


def _command(args: list[str], cwd: Path, log: Path, timeout: int) -> None:
    with log.open("a") as handle:
        handle.write("$ " + " ".join(args) + "\n")
        handle.flush()
        try:
            result = subprocess.run(args, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT,
                                    timeout=timeout, check=False)
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(f"Short MD exceeded {timeout} seconds; see {log.name}") from exc
        if result.returncode:
            handle.flush()
            detail = log.read_text(errors="replace")[-3500:]
            raise RuntimeError(f"{args[0]} {args[1]} failed; see {log.name}\n{detail}")


def run_short_md(system: Path, time_ps: float = 20.0, threads: int = 2,
                 timeout_s: int = 300) -> Path:
    """Run a short check after build; never overwrite the input system."""
    if not 1 <= time_ps <= 100 or not 1 <= threads <= 4:
        raise ValueError("NVT must be 1–100 ps and use 1–4 threads")
    if not shutil.which("gmx"):
        raise RuntimeError("GROMACS is not available")
    system = Path(system).resolve()
    if not all((system / name).is_file() for name in ("system.gro", "system.top", "minimization.mdp")):
        raise ValueError("Build a GROMACS system first")
    work = system / "short_md"
    if work.exists():
        raise FileExistsError("Short MD already exists for this build")
    work.mkdir()
    log = work / "short_md.log"
    _command(["gmx", "grompp", "-f", "../minimization.mdp", "-c", "../system.gro",
              "-p", "../system.top", "-o", "em.tpr", "-maxwarn", "0"], work, log, timeout_s)
    _command(["gmx", "mdrun", "-deffnm", "em", "-nt", str(threads)], work, log, timeout_s)
    nsteps = round(time_ps / 0.001)
    (work / "nvt.mdp").write_text(
        "integrator = md\n" f"dt = 0.001\nnsteps = {nsteps}\n"
        "cutoff-scheme = Verlet\nnstlist = 20\nrlist = 1.1\n"
        "vdwtype = cut-off\nvdw-modifier = Potential-shift-verlet\nrvdw = 1.1\n"
        "coulombtype = reaction-field\nrcoulomb = 1.1\nepsilon-r = 15\nepsilon-rf = 0\n"
        "tcoupl = v-rescale\ntc-grps = System\ntau-t = 1.0\nref-t = 300\n"
        "pcoupl = no\npbc = xyz\nconstraints = none\n"
        "gen-vel = yes\ngen-temp = 300\ngen-seed = 2026\n"
        "nstenergy = 1000\nnstlog = 1000\nnstxout-compressed = 1000\n"
    )
    _command(["gmx", "grompp", "-f", "nvt.mdp", "-c", "em.gro", "-p", "../system.top",
              "-o", "nvt.tpr", "-maxwarn", "0"], work, log, timeout_s)
    _command(["gmx", "mdrun", "-deffnm", "nvt", "-nt", str(threads)], work, log, timeout_s)
    (work / "summary.json").write_text(json.dumps({"status": "completed", "nvt_ps": time_ps,
                                                  "threads": threads, "scope": "sanity check"}, indent=2))
    return work


def archive(system: Path) -> bytes:
    """Bundle the prepared system plus any completed short check."""
    import io
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w", zipfile.ZIP_DEFLATED) as zipped:
        for file in sorted(Path(system).rglob("*")):
            if file.is_file():
                zipped.write(file, file.relative_to(system))
    return data.getvalue()
