"""Keep GōMartini virtual CA sites out of user-facing trajectory renderings."""
from __future__ import annotations

from . import molecular_viewer as mv


_ORIGINAL_COMPONENT_MAP = mv._gif_component_map


def _component_map_without_go_virtuals(gro_path, reline_chloride_count: int = 0):
    components = list(_ORIGINAL_COMPONENT_MAP(gro_path, reline_chloride_count))
    atoms = mv._parse_gro_atoms(gro_path)
    for index, (atom, component) in enumerate(zip(atoms, components)):
        # GōMartini adds one virtual site per residue, written as a CA bead after
        # the physical BB/SC protein beads.  It is part of the force-field
        # machinery, not a bead that should be shown to the user.
        if component == "protein" and str(atom["name"]).strip().upper() == "CA":
            components[index] = "hidden"
    return components


def apply() -> None:
    mv._gif_component_map = _component_map_without_go_virtuals
