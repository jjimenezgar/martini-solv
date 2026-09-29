from __future__ import annotations

import argparse
import json
from pathlib import Path

from .builder import build
from .models import BuildConfig, Solute


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a Martini 3 protein solvent box")
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--solvent", choices=["water", "reline"], default="water")
    parser.add_argument("--water-fraction", type=float, default=0)
    parser.add_argument("--salt", type=float, default=0.15)
    parser.add_argument("--box-distance", type=float, default=1.2)
    parser.add_argument("--solute", action="append", default=[], metavar="NAME:SMILES:COUNT")
    args = parser.parse_args()
    solutes = []
    for entry in args.solute:
        name, separator, rest = entry.partition(":")
        if not separator or ":" not in rest:
            parser.error("--solute requires NAME:SMILES:COUNT")
        smiles, _, count = rest.rpartition(":")
        solutes.append(Solute(name, smiles, int(count)))
    config = BuildConfig(args.solvent, args.water_fraction, args.salt,
                         args.box_distance, solutes=solutes)
    result = build(args.pdb, args.out, config)
    print(json.dumps({"output": str(result), "manifest": str(result / "manifest.json")}))


if __name__ == "__main__":
    main()
