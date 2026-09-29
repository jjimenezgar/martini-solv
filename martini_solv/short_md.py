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
    edr: Path | None
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
    results = [StageResult("minimization", em_base.with_suffix(".gro"), em_tpr, None, em_base.with_suffix(".edr"), elapsed, _ns_day(stdout))]

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
        edr = base.with_suffix(".edr")
        results.append(StageResult(
            stage.name,
            base.with_suffix(".gro"),
            base.with_suffix(".tpr"),
            xtc if xtc.is_file() else None,
            edr if edr.is_file() else None,
            elapsed,
            _ns_day(stdout),
        ))
        prev_gro = base.with_suffix(".gro")
        prev_cpt = base.with_suffix(".cpt")

    return ShortMDResult(work, results, log)


@dataclass
class AnalysisResult:
    kind: str
    x_label: str
    y_label: str
    points: list[tuple[float, float]]
    output_path: Path


def parse_xvg(path: Path) -> list[tuple[float, float]]:
    points: list[tuple[float, float]] = []
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith(("#", "@")):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            points.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    return points


def _gmx_energy_term_index(gmx: str, edr: Path, term: str, cwd: Path) -> int:
    """Discover an energy-term index from the GROMACS menu instead of relying on name input."""
    probe = subprocess.run(
        [gmx, "energy", "-f", str(edr), "-o", str(cwd / "_probe.xvg")],
        cwd=cwd,
        input="\n",
        text=True,
        capture_output=True,
        check=False,
    )
    menu = (probe.stdout or "") + "\n" + (probe.stderr or "")
    # GROMACS prints numbered energy terms in columns; match Density exactly.
    match = re.search(rf"(?m)(?:^|\s)(\d+)\s+{re.escape(term)}(?:\s|$)", menu)
    if not match:
        available = " ".join(line.strip() for line in menu.splitlines()[-40:])
        raise RuntimeError(
            f"GROMACS energy term {term!r} was not found in the selected EDR. "
            f"Available menu output: {available[-1800:]}"
        )
    return int(match.group(1))


def run_short_md_analysis(
    kind: str,
    tpr: Path,
    xtc: Path | None,
    edr: Path | None,
    output_dir: Path,
) -> AnalysisResult:
    """Run RMSD/RMSF or box-density analysis for one selected Short MD stage."""
    normalized = kind.strip().lower()
    if normalized not in {"rmsd", "rmsf", "density"}:
        raise ValueError("Analysis must be rmsd, rmsf or density")
    gmx = shutil.which("gmx")
    if not gmx:
        raise RuntimeError("GROMACS is not available")
    output_dir.mkdir(parents=True, exist_ok=True)

    if normalized in {"rmsd", "rmsf"}:
        if not tpr.is_file() or xtc is None or not xtc.is_file():
            raise FileNotFoundError("The selected stage needs both TPR and XTC files")
        index_path = output_dir / "protein_backbone.ndx"
        select = subprocess.run(
            [gmx, "select", "-s", str(tpr), "-on", str(index_path), "-select", "name BB"],
            cwd=output_dir,
            text=True,
            capture_output=True,
            check=False,
        )
        if select.returncode != 0 or not index_path.is_file():
            detail = (select.stderr or select.stdout or "Could not select BB beads.").strip()
            raise RuntimeError(detail[-2000:])

        output = output_dir / f"protein_{normalized}.xvg"
        if normalized == "rmsd":
            command = [
                gmx, "rms", "-s", str(tpr), "-f", str(xtc), "-n", str(index_path),
                "-o", str(output), "-tu", "ns",
            ]
            selection = "0\n0\n"
            x_label, y_label = "Time (ns)", "RMSD (nm)"
        else:
            command = [
                gmx, "rmsf", "-s", str(tpr), "-f", str(xtc), "-n", str(index_path),
                "-o", str(output), "-res",
            ]
            selection = "0\n"
            x_label, y_label = "Residue", "RMSF (nm)"
        result = subprocess.run(
            command, cwd=output_dir, input=selection, text=True,
            capture_output=True, check=False,
        )
    else:
        if edr is None or not edr.is_file():
            raise FileNotFoundError("The selected stage needs an EDR file for density analysis")
        output = output_dir / "system_density.xvg"
        density_index = _gmx_energy_term_index(gmx, edr, "Density", output_dir)
        # gmx energy in GROMACS 2025.x does not support -tu. Energy XVG time
        # is emitted in ps, so convert to ns after parsing.
        command = [gmx, "energy", "-f", str(edr), "-o", str(output)]
        result = subprocess.run(
            command,
            cwd=output_dir,
            input=f"{density_index}\n",
            text=True,
            capture_output=True,
            check=False,
        )
        x_label, y_label = "Time (ns)", "Density (kg/m³)"

    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "GROMACS analysis failed").strip()
        raise RuntimeError(detail[-2500:])
    points = parse_xvg(output)
    if normalized == "density":
        points = [(x / 1000.0, y) for x, y in points]
    if not points:
        raise RuntimeError(f"{output.name} contains no numeric data")
    return AnalysisResult(normalized, x_label, y_label, points, output)

def _packaged_topology_text(system: Path) -> str:
    """Rewrite local ITP includes for the MartiniSurf-style package layout."""
    top = system / "system.top"
    if not top.is_file():
        raise FileNotFoundError("system.top is missing from the built system")
    local_itps = {path.name for path in system.glob("*.itp")}
    lines: list[str] = []
    pattern = re.compile(r'^(\s*#include\s+["<])([^">]+)([">].*)
    for raw in top.read_text(errors="replace").splitlines():
        match = pattern.match(raw)
        if match:
            name = Path(match.group(2)).name
            if name in local_itps:
                raw = f'{match.group(1)}system_itp/{name}{match.group(3)}'
        lines.append(raw)
    return "\n".join(lines) + "\n"


def _package_readme(has_short_md: bool) -> str:
    short_md_text = (
        "3_short_md/\n"
        "  Results of the optional MartiniSolv Short MD workflow, including TPR, GRO, "
        "XTC, EDR, CPT, logs and any generated analyses.\n\n"
        if has_short_md
        else ""
    )
    return (
        "MartiniSolv Simulation_Files\n"
        "=============================\n\n"
        "This package is organised like MartiniSurf so the generated system can be "
        "continued directly with GROMACS.\n\n"
        "0_topology/\n"
        "  system.top                Main GROMACS topology.\n"
        "  system_itp/               Martini force-field, protein, Go-model, solvent "
        "and free-molecule ITP files.\n\n"
        "1_mdp/\n"
        "  minimization.mdp          Build validation/minimisation input.\n"
        "  nvt.mdp                   Short NVT template.\n"
        "  npt.mdp                   Short NPT template.\n"
        "  production.mdp            Short production template; extend nsteps for a "
        "scientific production run.\n\n"
        "2_system/\n"
        "  system.gro                Final solvated/packed starting coordinates.\n"
        "  protein_cg.pdb            Coarse-grained protein, when available.\n"
        "  protein_clean.pdb         Cleaned atomistic input, when available.\n\n"
        + short_md_text +
        "metadata/\n"
        "  manifest.json             Exact MartiniSolv configuration and composition.\n"
        "  build.log                 Build commands and tool output.\n\n"
        "Example workflow from inside Simulation_Files:\n"
        "  gmx grompp -f 1_mdp/minimization.mdp -c 2_system/system.gro "
        "-r 2_system/system.gro -p 0_topology/system.top -o em.tpr\n"
        "  gmx mdrun -deffnm em\n"
        "  gmx grompp -f 1_mdp/nvt.mdp -c em.gro -r em.gro "
        "-p 0_topology/system.top -o nvt.tpr\n"
        "  gmx mdrun -deffnm nvt\n\n"
        "The supplied NVT/NPT/Production MDP files are short validation templates. "
        "Review timestep, duration, coupling groups and scientific protocol before "
        "using them for production research.\n"
    )


def archive(system: Path) -> bytes:
    """Create a MartiniSurf-style, simulation-ready Simulation_Files package."""
    system = Path(system).resolve()
    required = ("system.gro", "system.top", "minimization.mdp")
    missing = [name for name in required if not (system / name).is_file()]
    if missing:
        raise FileNotFoundError("Cannot package system; missing: " + ", ".join(missing))

    root = Path("Simulation_Files")
    data = io.BytesIO()
    short_md = system / "short_md"
    has_short_md = short_md.is_dir()

    with zipfile.ZipFile(data, "w", zipfile.ZIP_DEFLATED) as zipped:
        # 0_topology: main TOP plus every local include required by that topology.
        zipped.writestr(str(root / "0_topology" / "system.top"), _packaged_topology_text(system))
        for path in sorted(system.glob("*.itp")):
            zipped.write(path, str(root / "0_topology" / "system_itp" / path.name))

        # 1_mdp: always provide a complete minimal continuation protocol.
        zipped.write(
            system / "minimization.mdp",
            str(root / "1_mdp" / "minimization.mdp"),
        )
        for stage_name in STAGE_ORDER:
            defaults = DEFAULT_STAGE_SETTINGS[stage_name]
            stage = StageSettings(
                stage_name,
                True,
                float(defaults["dt_ps"]),
                float(defaults["time_ns"]),
            )
            zipped.writestr(
                str(root / "1_mdp" / f"{stage_name}.mdp"),
                _stage_mdp(stage, DEFAULT_XTC_WRITE_EVERY_PS),
            )

        # 2_system: canonical coordinates plus useful structure references.
        zipped.write(system / "system.gro", str(root / "2_system" / "system.gro"))
        for name in ("protein_cg.pdb", "protein_clean.pdb", "input.pdb"):
            path = system / name
            if path.is_file():
                target = "input_atomistic.pdb" if name == "input.pdb" else name
                zipped.write(path, str(root / "2_system" / target))

        # Optional Short MD results remain separate from reusable setup files.
        if has_short_md:
            for path in sorted(short_md.rglob("*")):
                if path.is_file():
                    zipped.write(
                        path,
                        str(root / "3_short_md" / path.relative_to(short_md)),
                    )

        # Reproducibility metadata and logs.
        for name in ("manifest.json", "build.log"):
            path = system / name
            if path.is_file():
                zipped.write(path, str(root / "metadata" / name))

        zipped.writestr(str(root / "README.txt"), _package_readme(has_short_md))

    return data.getvalue()
)
    for raw in top.read_text(errors="replace").splitlines():
        match = pattern.match(raw)
        if match:
            name = Path(match.group(2)).name
            if name in local_itps:
                raw = f'{match.group(1)}system_itp/{name}{match.group(3)}'
        lines.append(raw)
    return "\n".join(lines) + "\n"


def _package_readme(has_short_md: bool) -> str:
    short_md_text = (
        "3_short_md/\n"
        "  Results of the optional MartiniSolv Short MD workflow, including TPR, GRO, "
        "XTC, EDR, CPT, logs and any generated analyses.\n\n"
        if has_short_md
        else ""
    )
    return (
        "MartiniSolv Simulation_Files\n"
        "=============================\n\n"
        "This package is organised like MartiniSurf so the generated system can be "
        "continued directly with GROMACS.\n\n"
        "0_topology/\n"
        "  system.top                Main GROMACS topology.\n"
        "  system_itp/               Martini force-field, protein, Go-model, solvent "
        "and free-molecule ITP files.\n\n"
        "1_mdp/\n"
        "  minimization.mdp          Build validation/minimisation input.\n"
        "  nvt.mdp                   Short NVT template.\n"
        "  npt.mdp                   Short NPT template.\n"
        "  production.mdp            Short production template; extend nsteps for a "
        "scientific production run.\n\n"
        "2_system/\n"
        "  system.gro                Final solvated/packed starting coordinates.\n"
        "  protein_cg.pdb            Coarse-grained protein, when available.\n"
        "  protein_clean.pdb         Cleaned atomistic input, when available.\n\n"
        + short_md_text +
        "metadata/\n"
        "  manifest.json             Exact MartiniSolv configuration and composition.\n"
        "  build.log                 Build commands and tool output.\n\n"
        "Example workflow from inside Simulation_Files:\n"
        "  gmx grompp -f 1_mdp/minimization.mdp -c 2_system/system.gro "
        "-r 2_system/system.gro -p 0_topology/system.top -o em.tpr\n"
        "  gmx mdrun -deffnm em\n"
        "  gmx grompp -f 1_mdp/nvt.mdp -c em.gro -r em.gro "
        "-p 0_topology/system.top -o nvt.tpr\n"
        "  gmx mdrun -deffnm nvt\n\n"
        "The supplied NVT/NPT/Production MDP files are short validation templates. "
        "Review timestep, duration, coupling groups and scientific protocol before "
        "using them for production research.\n"
    )


def archive(system: Path) -> bytes:
    """Create a MartiniSurf-style, simulation-ready Simulation_Files package."""
    system = Path(system).resolve()
    required = ("system.gro", "system.top", "minimization.mdp")
    missing = [name for name in required if not (system / name).is_file()]
    if missing:
        raise FileNotFoundError("Cannot package system; missing: " + ", ".join(missing))

    root = Path("Simulation_Files")
    data = io.BytesIO()
    short_md = system / "short_md"
    has_short_md = short_md.is_dir()

    with zipfile.ZipFile(data, "w", zipfile.ZIP_DEFLATED) as zipped:
        # 0_topology: main TOP plus every local include required by that topology.
        zipped.writestr(str(root / "0_topology" / "system.top"), _packaged_topology_text(system))
        for path in sorted(system.glob("*.itp")):
            zipped.write(path, str(root / "0_topology" / "system_itp" / path.name))

        # 1_mdp: always provide a complete minimal continuation protocol.
        zipped.write(
            system / "minimization.mdp",
            str(root / "1_mdp" / "minimization.mdp"),
        )
        for stage_name in STAGE_ORDER:
            defaults = DEFAULT_STAGE_SETTINGS[stage_name]
            stage = StageSettings(
                stage_name,
                True,
                float(defaults["dt_ps"]),
                float(defaults["time_ns"]),
            )
            zipped.writestr(
                str(root / "1_mdp" / f"{stage_name}.mdp"),
                _stage_mdp(stage, DEFAULT_XTC_WRITE_EVERY_PS),
            )

        # 2_system: canonical coordinates plus useful structure references.
        zipped.write(system / "system.gro", str(root / "2_system" / "system.gro"))
        for name in ("protein_cg.pdb", "protein_clean.pdb", "input.pdb"):
            path = system / name
            if path.is_file():
                target = "input_atomistic.pdb" if name == "input.pdb" else name
                zipped.write(path, str(root / "2_system" / target))

        # Optional Short MD results remain separate from reusable setup files.
        if has_short_md:
            for path in sorted(short_md.rglob("*")):
                if path.is_file():
                    zipped.write(
                        path,
                        str(root / "3_short_md" / path.relative_to(short_md)),
                    )

        # Reproducibility metadata and logs.
        for name in ("manifest.json", "build.log"):
            path = system / name
            if path.is_file():
                zipped.write(path, str(root / "metadata" / name))

        zipped.writestr(str(root / "README.txt"), _package_readme(has_short_md))

    return data.getvalue()
