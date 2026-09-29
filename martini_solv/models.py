from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import re


@dataclass(frozen=True)
class Solute:
    name: str
    smiles: str
    count: int

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


@dataclass(frozen=True)
class BuildConfig:
    solvent: str = "water"
    water_fraction: float = 0.0
    salt_m: float = 0.15
    box_distance_nm: float = 1.0
    des_pairs_per_nm3: float = 3.2
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
        if self.solvent not in {"water", "reline"}:
            raise ValueError("Solvent must be water or reline")
        if not 0 <= self.water_fraction < 0.5:
            raise ValueError("Water mole fraction must be between 0 and 0.5")
        if self.solvent == "water" and self.water_fraction:
            raise ValueError("Water fraction applies only to reline")
        if self.solvent == "reline" and self.salt_m:
            raise ValueError("Added NaCl is currently supported only for water")
        if self.salt_m < 0 or not 0.5 <= self.box_distance_nm <= 5:
            raise ValueError("Check salt concentration and box distance")
        if not 0 < self.des_pairs_per_nm3 <= 4:
            raise ValueError("DES pair density is outside the supported initial packing range")
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
        reserved = {"W", "NA", "CL", "CHOL", "UREA", "MOLECULE_0", self.molecule_name.upper()}
        if any(s.name.upper() in reserved for s in self.solutes):
            raise ValueError("Additional molecule name conflicts with a solvent, ion or protein type")
        for s in self.solutes:
            s.validate()

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True)


def reline_counts(box_nm: float, water_fraction: float, pairs_per_nm3: float) -> dict[str, float | int]:
    """x_H2O counts molecular waters; one Martini W bead represents four waters."""
    pairs = max(1, round(box_nm**3 * pairs_per_nm3))
    urea = 2 * pairs
    water_beads = round(water_fraction * (pairs + urea) / (4 * (1 - water_fraction)))
    actual = 4 * water_beads / (3 * pairs + 4 * water_beads)
    return {"CHOL": pairs, "CL": pairs, "UREA": urea, "W": water_beads, "x_water_actual": actual}
