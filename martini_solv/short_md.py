"""Optional short GROMACS protocol used by the Streamlit interface."""
from __future__ import annotations

from dataclasses import dataclass, field
import io
import re
from pathlib import Path
import shutil
import subprocess
import time
import zipfile


DEFAULT_XTC_WRITE_EVERY_PS = 10.0
DEFAULT_GROMPP_MAXWARN = 3
DEFAULT_STAGE_SETTINGS = {
    "nvt": {"enabled": True, "dt_ps": 0.001, "time_ns": 0.020},
    "npt": {"enabled": True, "dt_ps": 0.005, "time_ns": 0.020},
    "production": {"enabled": True, "dt_ps": 0.010, "time_ns": 0.100},
}
STAGE_ORDER = ["nvt", "npt", "production"]


@dataclass(frozen=True)
class StageSettings:
    name: str
    enabled: bool
    dt_ps: float
    time_ns: float


@dataclass(frozen=True)
class ShortMDConfig:
    output_tag: str = "Protein MD"
    xtc_write_every_ps: float = DEFAULT_XTC_WRITE_EVERY_PS
    grompp_maxwarn: int = DEFAULT_GROMPP_MAXWARN
    threads: int = 2
    stages: tuple[StageSettings, ...] = field(default_factory=tuple)


@dataclass
class StageResult:
    name: str
    gro: Path
    tpr: Path
    xtc: Path | None
    elapsed_s: float
    ns_day: float | None


@dataclass
class ShortMDResult:
    work_dir: Path
    stages: list[StageResult]
    log: Path

    @property
    def last_gro(self) -> Path:
        return self.stages[-1].gro

    @property
    def last_tpr(self) -> Path:
        return self.stages[-1].tpr

    @property
    def last_xtc(self) -> Path | None:
        return self.stages[-1].xtc


def selected_stages(config: ShortMDConfig) -> list[StageSettings]:
    return [stage for stage in config.stages if stage.enabled]


def validate_stage_order(stages: list[StageSettings]) -> list[str]:
    names = [stage.name for stage in stages]
    errors: list[str] = []
    if "npt" in names and "nvt" not in names:
        errors.append("NPT requires NVT to be enabled.")
    if "production" in names and "npt" not in names:
        errors.append("Production requires NPT to be enabled.")
    return errors


def _safe_tag(value: str) -> str:
    clean = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value.strip())
    return clean or "Protein_MD"


def _command(args: list[str], cwd: Path, log: Path, timeout: int) -> tuple[float, str]:
    started = time.monotonic()
    with log.open("a") as handle:
        handle.write("$ " + " ".join(args) + "\n")
        handle.flush()
        try:
            result = subprocess.run(
                args,
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(f"Short MD exceeded {timeout} seconds; see {log.name}") from exc
        handle.write(result.stdout or "")
        handle.write("\n")
    if result.returncode:
        detail = (result.stdout or "")[-5000:]
        raise RuntimeError(f"{args[0]} {args[1]} failed (exit {result.returncode}).\n{detail}")
    return time.monotonic() - started, result.stdout or ""


def _ns_day(text: str) -> float | None:
    match = re.search(r"Performance:\s+([0-9.]+)", text)
    return float(match.group(1)) if match else None


def _common_mdp(dt_ps: float, nsteps: int, xtc_stride: int) -> list[str]:
    return [
        "integrator = md",
        f"dt = {dt_ps:.8f}",
        f"nsteps = {nsteps}",
        "cutoff-scheme = Verlet",
        "nstlist = 20",
        "rlist = 1.1",
        "vdwtype = cut-off",
        "vdw-modifier = Potential-shift-verlet",
        "rvdw = 1.1",
        "coulombtype = reaction-field",
        "rcoulomb = 1.1",
        "epsilon-r = 15",
        "epsilon-rf = 0",
        "tcoupl = v-rescale",
        "tc-grps = System",
        "tau-t = 1.0",
        "ref-t = 300",
        "pbc = xyz",
        "constraints = none",
        "nstenergy = 1000",
        "nstlog = 1000",
        f"nstxout-compressed = {xtc_stride}",
        "compressed-x-precision = 100",
    ]


def _stage_mdp(stage: StageSettings, xtc_every_ps: float) -> str:
    if stage.dt_ps <= 0 or stage.time_ns <= 0:
        raise ValueError(f"{stage.name}: timestep and duration must be positive")
    nsteps = max(1, round(stage.time_ns * 1000.0 / stage.dt_ps))
    xtc_stride = max(1, round(xtc_every_ps / stage.dt_ps))
    rows = _common_mdp(stage.dt_ps, nsteps, xtc_stride)

    if stage.name == "nvt":
        rows += [
            "define = -DPOSRES",
            "pcoupl = no",
            "gen-vel = yes",
            "gen-temp = 300",
            "gen-seed = 2026",
            "continuation = no",
        ]
    else:
        rows += [
            "pcoupl = C-rescale",
            "pcoupltype = isotropic",
            "tau-p = 12.0",
            "ref-p = 1.0",
            "compressibility = 3e-4",
            "gen-vel = no",
            "continuation = yes",
        ]
        if stage.name == "npt":
            rows.append("define = -DPOSRES")
    return "\n".join(rows) + "\n"


def run_short_md(system: Path, config: ShortMDConfig, timeout_s: int = 900) -> ShortMDResult:
    """Run minimization followed by selected NVT/NPT/Production stages."""
    if not shutil.which("gmx"):
        raise RuntimeError("GROMACS is not available")
    if not 1 <= int(config.threads) <= 4:
        raise ValueError("CPU threads must be between 1 and 4")
    if config.xtc_write_every_ps <= 0:
        raise ValueError("XTC write interval must be positive")

    stages = selected_stages(config)
    errors = validate_stage_order(stages)
    if errors:
        raise ValueError(" ".join(errors))
    if not stages:
        raise ValueError("Enable at least one Short MD stage")

    system = Path(system).resolve()
    if not all((system / name).is_file() for name in ("system.gro", "system.top", "minimization.mdp")):
        raise ValueError("Build a GROMACS system first")

    work = system / "short_md"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()
    log = work / "short_md.log"
    tag = _safe_tag(config.output_tag)

    # Energy minimization is always run before the selectable short protocol.
    em_tpr = work / f"{tag}_minimization.tpr"
    _command([
        "gmx", "grompp", "-f", "../minimization.mdp", "-c", "../system.gro", "-r", "../system.gro",
        "-p", "../system.top", "-o", em_tpr.name, "-maxwarn", str(config.grompp_maxwarn)
    ], work, log, timeout_s)
    em_base = work / f"{tag}_minimization"
    elapsed, stdout = _command(["gmx", "mdrun", "-deffnm", em_base.name, "-nt", str(config.threads)], work, log, timeout_s)
    results = [StageResult("minimization", em_base.with_suffix(".gro"), em_tpr, None, elapsed, _ns_day(stdout))]

    prev_gro = em_base.with_suffix(".gro")
    prev_cpt = em_base.with_suffix(".cpt")
    for stage in stages:
        mdp = work / f"{tag}_{stage.name}.mdp"
        mdp.write_text(_stage_mdp(stage, config.xtc_write_every_ps))
        base = work / f"{tag}_{stage.name}"
        command = [
            "gmx", "grompp", "-f", mdp.name, "-c", prev_gro.name, "-r", prev_gro.name,
            "-p", "../system.top", "-o", base.with_suffix(".tpr").name,
            "-maxwarn", str(config.grompp_maxwarn),
        ]
        if prev_cpt.is_file():
            command += ["-t", prev_cpt.name]
        _command(command, work, log, timeout_s)
        elapsed, stdout = _command(["gmx", "mdrun", "-deffnm", base.name, "-nt", str(config.threads)], work, log, timeout_s)
        xtc = base.with_suffix(".xtc")
        results.append(StageResult(
            stage.name,
            base.with_suffix(".gro"),
            base.with_suffix(".tpr"),
            xtc if xtc.is_file() else None,
            elapsed,
            _ns_day(stdout),
        ))
        prev_gro = base.with_suffix(".gro")
        prev_cpt = base.with_suffix(".cpt")

    return ShortMDResult(work, results, log)


def archive(system: Path) -> bytes:
    """Bundle the prepared system plus any completed short protocol."""
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w", zipfile.ZIP_DEFLATED) as zipped:
        for file in sorted(Path(system).rglob("*")):
            if file.is_file():
                zipped.write(file, file.relative_to(system))
    return data.getvalue()
