from __future__ import annotations

import unittest

from martini_solv.gif_backbone_only import _backbone_only_visible_indices


class GifBackboneOnlyTests(unittest.TestCase):
    def test_hides_protein_sidechains_but_keeps_other_components(self) -> None:
        atoms = [
            {"name": "BB"},
            {"name": "SC1"},
            {"name": "BB2"},
            {"name": "W"},
            {"name": "NA"},
        ]
        components = ["protein", "protein", "protein", "solvent", "ions"]
        self.assertEqual(
            _backbone_only_visible_indices(atoms, components, [0, 1, 2, 3, 4]),
            [0, 2, 3, 4],
        )

    def test_nonprotein_selection_is_unchanged(self) -> None:
        atoms = [{"name": "W"}, {"name": "CL"}]
        components = ["solvent", "ions"]
        self.assertEqual(_backbone_only_visible_indices(atoms, components, [0, 1]), [0, 1])


if __name__ == "__main__":
    unittest.main()
