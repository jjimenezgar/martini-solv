# MartiniSolv

An independent, small protein system builder inspired by [MartiniSurf](https://github.com/jjimenezgar/MartiniSurf). Its Streamlit interface keeps MartiniSurf's dark cards and clear step navigation, with a teal accent. The scientific core also runs without a browser.

**PDB → martinize2 (Martini 3) → cubic box → water or reline → optional free molecules → GROMACS files.** The output includes `system.gro`, `system.top`, force-field and molecular `.itp` files, `minimization.mdp`, a provenance manifest and a command log. A successful build runs `gmx grompp` with zero tolerated warnings; it does **not** claim the structure is equilibrated.

## Install and run

Linux/Conda, GROMACS, martinize2, DSSP and INSANE are required. The official Martinize2 and INSANE projects install their CLI via the `vermouth` and `insane` Python packages, respectively; `environment.yml` uses pip for those two tools inside Conda. This is an initial environment specification rather than a tested cross-platform lockfile:

```bash
mamba env create -f environment.yml
conda activate martini-solv
streamlit run app.py
```

Download [1UBQ from RCSB](https://files.rcsb.org/download/1UBQ.pdb), save it as `1ubq.pdb`, then use the upload control or the CLI:

```bash
python -m martini_solv.cli --pdb 1ubq.pdb --out builds/ubiquitin --solvent water --salt 0.15 --box-distance 1.2
```

For reline (ChCl:urea 1:2) at requested water mole fraction 0.10:

```bash
python -m martini_solv.cli --pdb 1ubq.pdb --out builds/ubiquitin_reline --solvent reline --salt 0 --water-fraction 0.10
```

For optional freely dissolved small molecules: add rows in Streamlit or use `--solute NAME:SMILES:COUNT` multiple times. This path requires a separate `martini_mapper` executable and its compatible dependencies; [MartiniSurf's integration](https://github.com/jjimenezgar/MartiniSurf/blob/master/streamlit_app/linker_generator.py) documents the same CLI invocation. Only neutral additional molecules are currently accepted. A generated topology must be scientifically reviewed for chemical accuracy and applicability to the chosen solvent. The core verifies its `.gro`/`.itp` outputs and fails explicitly if unavailable.

## Solvent models and limitations

The builder retrieves Martini 3 force-field files from the [official force-fields repository](https://github.com/marrink-lab/martini-forcefields) and the published [Vainikka et al. DES models](https://github.com/vainikanpete/martini3-DES-models), each at an immutable Git commit recorded in `martini_solv/builder.py` and every build's `manifest.json`. Model files retain their authorship and original licenses. The ChCl and urea models are associated with [Vainikka et al., ACS Sustainable Chemistry & Engineering (2021)](https://doi.org/10.1021/acssuschemeng.1c06521). The published DES repository is Apache-2.0 licensed.

For reline, the program requests one `CHOL`, one `CL`, and two `UREA` molecules per formula unit. The requested `x_H2O` counts **real water molecules**; a Martini W bead represents four waters. The actual fraction after integer rounding is saved. `des_pairs_per_nm3` is currently an initial packing target, **not** a density validated for each temperature, hydration level or protein. If insertion fails, enlarge the box or adjust packing in a later release. Water fraction beyond published benchmarks is exploratory and must be validated by the user. Added NaCl is supported only for the water path. For the DES path, only an electrically neutral protein is currently supported until ion accounting is extended.

Input PDBs must already have a sensible backbone and chain assignment. The current build does not repair missing loops, guarantee chirality or conduct production MD. DSSP and the elastic network provide the first protein model, whose appropriateness should be checked for the target protein. Results are prepared for *subsequent minimization and equilibration*, not ready for immediate production simulation.

## Layout and reproducibility

- `app.py`: Streamlit presentation; no molecular logic.
- `martini_solv/models.py`: validated settings and solvent stoichiometry.
- `martini_solv/builder.py`: headless tool orchestration, exact counts, topologies, GROMACS check.
- `martini_solv/cli.py`: same builder on a workstation or HPC login node.
- `tests/`: quick validation of composition and parameters (`python -m unittest discover -s tests -v`).

Each build writes its input PDB, actual settings, pinned upstream URLs and `build.log`. A failed build remains in a `.incomplete` folder for diagnosis. The Streamlit process currently runs the build synchronously; deploy it on a machine with the scientific tools and enough CPU/memory. A persistent job queue and optional short MD are follow-up work, not part of this first usable prototype.

The project deliberately excludes surfaces, immobilization, linkers, deposition and production simulations. This project extracts a small, self-contained protein-preparation workflow inspired by MartiniSurf and reimplements the orchestration layer without its surface workflow.

## Test status

The pure-Python tests are runnable without GROMACS. An end-to-end scientific smoke run requires the Conda environment and upstream model downloads. Consult the build log and check `gmx grompp` before using any generated files.
