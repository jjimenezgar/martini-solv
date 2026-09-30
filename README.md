<div align="center">

<img src="assets/martinisolv_logo.png" alt="MartiniSolv logo" width="220">

# MartiniSolv

### Martini 3 protein-in-solution system builder

Prepare, solvate, inspect and validate coarse-grained protein systems through a clean Streamlit workflow or from the command line.

[![Open MartiniSolv](https://img.shields.io/badge/Launch-MartiniSolv-42C7D5?style=for-the-badge&logo=streamlit&logoColor=white)](https://martinisolv.streamlit.app/)
[![Python](https://img.shields.io/badge/Python-3.11%2B-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![Martini 3](https://img.shields.io/badge/Martini-3-8FEAF2?style=flat-square)](https://cgmartini.nl/)
[![License: MIT](https://img.shields.io/badge/License-MIT-F4F5F7?style=flat-square)](LICENSE)

**Developed by [Juan Carlos Jiménez-García](https://github.com/jjimenezgar)**

</div>

---

## Overview

**MartiniSolv** is a compact, reproducible workflow for preparing **Martini 3 protein systems in solution**. It is inspired by [MartiniSurf](https://github.com/jjimenezgar/MartiniSurf), but deliberately removes surfaces, linkers and immobilization logic.

The application guides the user through:

```text
Home → Structure → Model → Environment → Review & Build → Short MD
```

The scientific workflow is:

```text
PDB
 ↓
martinize2 / Martini 3
 ↓
simulation box
 ↓
water or DES environment
 ↓
optional free molecules
 ↓
GROMACS-ready topology + coordinates
 ↓
optional short MD validation
```

A successful build is checked with **GROMACS `grompp`** before it is presented as ready for download.

---

## Launch the app

<div align="center">

### [▶ Open MartiniSolv on Streamlit](https://martinisolv.streamlit.app/)

No local installation is required for the hosted interface.

</div>

---

## Main features

### Protein preparation

- Upload a PDB file or enter an RCSB PDB ID.
- Automatic protein-only cleanup.
- Martini 3 coarse-graining with `martinize2`.
- DSSP secondary-structure assignment.
- Optional **GōMartini** model.
- Configurable Gō interaction strength.
- Optional elastic network.
- Position restraints.
- Chain merging controls.
- Interactive coarse-grained structure visualization.

### Solvent environments

MartiniSolv currently supports:

| Environment | Composition | Water |
|---|---|---|
| Martini water | Water + optional NaCl | native |
| Reline | ChCl : urea = **1 : 2** | dry or user-defined wet fraction |
| ChCl:sorbitol | ChCl : sorbitol = **1 : 1** | dry or user-defined wet fraction |

For the DES systems, the builder preserves the intrinsic ChCl stoichiometry and adds only the extra counterions required to neutralize the complete system.

### Free molecules from SMILES

Additional dissolved molecules can be generated directly in the **Environment** step.

The workflow:

1. Enter a neutral SMILES.
2. Assign a molecule name and number of copies.
3. Generate the coarse-grained model with **Martini Mapper**.
4. Inspect the generated beads and Martini types.
5. Optionally assign a molecular charge.
6. Select which mapped bead carries that charge.
7. Insert the requested number of copies into the final system.

Supported manual net charges are currently:

```text
-2  -1   0   +1   +2
```

The default is **0**.

Charged free molecules are automatically counterbalanced during the build:

- positive free-molecule charge → additional **Cl⁻**
- negative free-molecule charge → additional **Na⁺**

The total protein + free-molecule charge is considered when neutralizing the system.

> Manual charge assignment changes the numerical charge in the generated ITP. MartiniSolv does not automatically change the Martini bead type. The selected bead and charge should therefore be consistent with a scientifically validated CG model.

---

## Review & Build

Before running MD, MartiniSolv provides an interactive inspection of the generated system.

The viewer can independently show or hide:

- Protein
- Free molecules
- Solvent
- Ions

Free-molecule copy counts are explicitly verified against the final topology. The interface reports both:

```text
Copies requested
Copies included
```

A mismatch causes the build to fail instead of silently generating an inconsistent system.

---

## Short MD validation

MartiniSolv includes an optional short validation workflow:

```text
Energy minimization
        ↓
       NVT
        ↓
       NPT
        ↓
   Production
```

The NVT, NPT and short Production stages expose:

- timestep
- simulation length
- XTC output frequency
- GROMACS warning tolerance
- CPU-thread count

Completed trajectories can be animated directly in the app.

### Analysis

For the selected MD stage, MartiniSolv can calculate:

- **RMSD** of protein backbone beads
- **RMSF** of protein backbone beads
- **system density**

The analyses are stage-specific: selecting Production analyzes the full Production stage rather than mixing data from NVT or NPT.

The Short MD workflow is intended as a **technical and structural validation step**, not as evidence of an equilibrated or converged production simulation.

---

## Simulation-ready download

The generated download is organized as:

```text
Simulation_Files/
├── 0_topology/
│   ├── system.top
│   └── system_itp/
│       ├── Protein.itp
│       ├── Martini force-field files
│       ├── Gō files when enabled
│       ├── solvent ITPs
│       └── free-molecule ITPs
│
├── 1_mdp/
│   ├── minimization.mdp
│   ├── nvt.mdp
│   ├── npt.mdp
│   └── production.mdp
│
├── 2_system/
│   ├── system.gro
│   ├── protein_cg.pdb
│   └── cleaned/input structures
│
├── 3_short_md/          # present after Short MD
│   ├── TPR / GRO / XTC
│   ├── EDR / CPT
│   ├── logs
│   └── analyses
│
├── metadata/
│   ├── manifest.json
│   └── build.log
│
└── README.txt
```

The packaged `system.top` is rewritten so local includes point correctly to `0_topology/system_itp/`, making the extracted package self-contained for continued GROMACS work.

---

## Local installation

The hosted Streamlit app is the easiest way to use MartiniSolv, but the project can also run locally or on an HPC system.

### Conda / Mamba

```bash
mamba env create -f conda/environment.yml
conda activate martini-solv
streamlit run app.py
```

The environment includes the main external tools used by the workflow:

- GROMACS
- martinize2 / vermouth
- INSANE
- Packmol
- MDTraj
- RDKit
- Martini Mapper

---

## Command-line use

MartiniSolv also exposes the headless builder through the CLI.

Example with water:

```bash
python -m martini_solv.cli \
  --pdb 1ubq.pdb \
  --out builds/ubiquitin \
  --solvent water \
  --salt 0.15 \
  --box-distance 1.0
```

Example with wet Reline:

```bash
python -m martini_solv.cli \
  --pdb 1ubq.pdb \
  --out builds/ubiquitin_reline \
  --solvent reline \
  --salt 0 \
  --water-fraction 0.10
```

---

## Reproducibility

Each build records:

- exact user configuration
- solvent composition
- protein net charge
- free-molecule charge metadata
- requested free-molecule copy counts
- generated composition
- pinned upstream model sources
- `martinize2` command
- build log
- GROMACS validation status

The main implementation is separated into reusable modules:

```text
app.py
└── Streamlit interface

martini_solv/
├── models.py
│   └── validated configuration and solvent composition
├── builder.py
│   └── system generation and GROMACS orchestration
├── molecular_viewer.py
│   └── structure and trajectory visualization
├── short_md.py
│   └── validation MD and analysis
└── cli.py
    └── command-line interface

tests/
└── regression and scientific consistency checks
```

---

## Scientific models

MartiniSolv retrieves Martini 3 force-field files from the official [Martini force-field repository](https://github.com/marrink-lab/martini-forcefields).

Reline components use the published Martini 3 DES models from [Vainikka et al.](https://github.com/vainikanpete/martini3-DES-models):

> Vainikka et al., *ACS Sustainable Chemistry & Engineering* (2021)  
> [https://doi.org/10.1021/acssuschemeng.1c06521](https://doi.org/10.1021/acssuschemeng.1c06521)

The ChCl:sorbitol implementation uses the ChCl Martini model together with the three-bead Martini 3 sorbitol model included in this project.

Upstream files are pinned to immutable revisions in the builder for reproducibility.

---

## Scope

MartiniSolv is intended to automate **system preparation and technical validation**.

It does not claim that generated systems are automatically scientifically validated for every molecule, solvent composition, temperature or research question.

Users should review:

- molecular mappings
- bead types
- manually assigned charges
- solvent composition
- equilibration protocol
- production MD settings

before using generated systems for scientific conclusions.

---

## Relationship to MartiniSurf

[MartiniSurf](https://github.com/jjimenezgar/MartiniSurf) focuses on biomolecular systems interacting with surfaces and immobilization environments.

**MartiniSolv** intentionally removes that complexity and focuses on proteins and dissolved molecules in bulk solution.

```text
MartiniSurf  → surfaces · linkers · immobilization
MartiniSolv  → proteins · solvents · free molecules
```

---

## Development

Run the lightweight regression suite with:

```bash
python -m unittest discover -s tests -v
```

GitHub Actions also checks Python compilation and the repository's automated tests.

---

<div align="center">

### MartiniSolv

**Protein systems in solution · Martini 3**

[Launch app](https://martinisolv.streamlit.app/) · [GitHub profile](https://github.com/jjimenezgar) · [MartiniSurf](https://github.com/jjimenezgar/MartiniSurf)

Developed by **Juan Carlos Jiménez-García**

</div>
