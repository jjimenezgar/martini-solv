from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import re


@dataclass(frozen=True)
class Solute:
    name: str
    smiles: str
    count: int
    net_charge: int = 0
    charged_bead: int | None = None

    def validate(self) -> None:
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]{0,15}", self.name):
            raise ValueError(f"Invalid molecule name: {self.name!r}")
        if not self.smiles.strip() or self.count < 1:
            raise ValueError("Each molecule needs SMILES and a positive count")
        try:
            from rdkit import Chem
        except ImportError as exc:
            raise RuntimeError("RDKit is required for SMILES validation") from exc
        if Chem.MolFromSmiles(self.smiles) is None:
            raise ValueError(f"Invalid SMILES for {self.name}")
        if not -2 <= int(self.net_charge) <= 2:
            raise ValueError("Free-molecule net charge must be between -2 and +2")
        if int(self.net_charge) != 0 and (self.charged_bead is None or int(self.charged_bead) < 1):
            raise ValueError("Choose which Martini bead carries the requested free-molecule charge")
        if int(self.net_charge) == 0 and self.charged_bead is not None:
            raise ValueError("Neutral free molecules must not define a charged bead")


@dataclass(frozen=True)
class BuildConfig:
    solvent: str = "water"
    water_fraction: float = 0.0
    salt_m: float = 0.15
    box_distance_nm: float = 1.0
    reline_density_g_cm3: float = 1.20
    chcl_sorbitol_density_g_cm3: float = 1.20
    seed: int = 2026

    # Protein model controls mirror the useful MartiniSurf protein controls.
    molecule_name: str = "Protein"
    merge_chains: str = "A"
    dssp: bool = True
    go: bool = True
    go_eps: float = 9.414
    elastic: bool = False
    elastic_force: int = 700
    position_restraints: str = "backbone"
    position_restraint_force: float = 1000.0
    maxwarn: int = 1
    martinize_extra_args: list[str] = field(default_factory=list)

    solutes: list[Solute] = field(default_factory=list)

    def validate(self) -> None:
        if self.solvent not in {"water", "reline", "chcl_sorbitol"}:
            raise ValueError("Solvent must be water, reline or chcl_sorbitol")
        if not 0 <= self.water_fraction < 0.5:
            raise ValueError("Water mole fraction must be between 0 and 0.5")
        if self.solvent == "water" and self.water_fraction:
            raise ValueError("Water fraction applies only to DES solvent modes")
        if self.solvent in {"reline", "chcl_sorbitol"} and self.salt_m:
            raise ValueError("Added NaCl is currently supported only for water")
        if self.salt_m < 0 or not 0.5 <= self.box_distance_nm <= 5:
            raise ValueError("Check salt concentration and box distance")
        if not 0.5 <= self.reline_density_g_cm3 <= 2.0:
            raise ValueError("Reline density must be between 0.5 and 2.0 g/cm³")
        if not 0.5 <= self.chcl_sorbitol_density_g_cm3 <= 2.0:
            raise ValueError("ChCl:sorbitol density must be between 0.5 and 2.0 g/cm³")
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]{0,15}", self.molecule_name):
            raise ValueError("Molecule name must start with a letter and contain only letters, numbers or underscores")
        if self.position_restraints not in {"backbone", "all", "none"}:
            raise ValueError("Position restraints must be backbone, all or none")
        if not 0 <= self.maxwarn <= 20:
            raise ValueError("martinize2 max warnings must be between 0 and 20")
        if self.go_eps < 0:
            raise ValueError("Go epsilon must be non-negative")
        if not 100 <= self.elastic_force <= 1500:
            raise ValueError("Elastic network force constant must be 100–1500 kJ/mol/nm²")
        if self.position_restraint_force <= 0:
            raise ValueError("Position-restraint force constant must be positive")
        if len({s.name.upper() for s in self.solutes}) != len(self.solutes):
            raise ValueError("Molecule names must be unique")
        reserved = {"W", "NA", "CL", "CHOL", "UREA", "SOR", "MOLECULE_0", self.molecule_name.upper()}
        if any(s.name.upper() in reserved for s in self.solutes):
            raise ValueError("Additional molecule name conflicts with a solvent, ion or protein type")
        for s in self.solutes:
            s.validate()

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True)


RELINE_DENSITY_G_CM3 = 1.20
CHCL_MOLAR_MASS_G_MOL = 139.62
UREA_MOLAR_MASS_G_MOL = 60.06
AVOGADRO = 6.02214076e23


def reline_counts(
    box_nm: float,
    water_fraction: float,
    density_g_cm3: float = RELINE_DENSITY_G_CM3,
) -> dict[str, float | int]:
    """Estimate a ChCl:urea 1:2 composition from the experimental bulk density.

    The dry DES formula unit is 1 choline chloride + 2 urea molecules.
    1 nm³ = 1e-21 cm³. Water fraction counts molecular waters; one Martini W
    bead represents four water molecules.
    """
    formula_mass = CHCL_MOLAR_MASS_G_MOL + 2 * UREA_MOLAR_MASS_G_MOL
    pairs_per_nm3 = density_g_cm3 * 1e-21 * AVOGADRO / formula_mass
    pairs = max(1, round(box_nm**3 * pairs_per_nm3))
    urea = 2 * pairs
    water_beads = round(water_fraction * (pairs + urea) / (4 * (1 - water_fraction)))
    actual = 4 * water_beads / (3 * pairs + 4 * water_beads) if (pairs or water_beads) else 0.0
    return {
        "CHOL": pairs,
        "CL": pairs,
        "UREA": urea,
        "W": water_beads,
        "x_water_actual": actual,
        "target_density_g_cm3": density_g_cm3,
        "dry_formula_units_per_nm3": pairs_per_nm3,
    }


SORBITOL_MOLAR_MASS_G_MOL = 182.17


def chcl_sorbitol_counts(
    box_nm: float,
    water_fraction: float,
    density_g_cm3: float,
) -> dict[str, float | int]:
    """Estimate ChCl:sorbitol 1:1 counts from a user-supplied dry-mixture density."""
    formula_mass = CHCL_MOLAR_MASS_G_MOL + SORBITOL_MOLAR_MASS_G_MOL
    units_per_nm3 = density_g_cm3 * 1e-21 * AVOGADRO / formula_mass
    pairs = max(1, round(box_nm**3 * units_per_nm3))
    # Water fraction counts molecular water; one Martini W bead represents four waters.
    water_beads = round(water_fraction * (2 * pairs) / (4 * (1 - water_fraction)))
    actual = 4 * water_beads / (2 * pairs + 4 * water_beads) if (pairs or water_beads) else 0.0
    return {
        "CHOL": pairs,
        "CL": pairs,
        "SOR": pairs,
        "W": water_beads,
        "x_water_actual": actual,
        "target_density_g_cm3": density_g_cm3,
        "dry_formula_units_per_nm3": units_per_nm3,
    }
