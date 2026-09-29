import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from martini_solv.builder import _clean_protein_pdb, _normalize_insane_ions
from martini_solv.models import BuildConfig, reline_counts


class TestModels(unittest.TestCase):
    def test_reline_stoichiometry_and_actual_water_fraction(self):
        counts = reline_counts(5.0, 0.10, 3.2)
        self.assertEqual(counts["UREA"], 2 * counts["CHOL"])
        self.assertEqual(counts["CL"], counts["CHOL"])
        self.assertLess(abs(counts["x_water_actual"] - 0.10), 0.002)


    def test_reline_disallows_added_salt(self):
        with self.assertRaisesRegex(ValueError, "only for water"):
            BuildConfig(solvent="reline", salt_m=0.15).validate()


    def test_dry_reline_has_no_water(self):
        self.assertEqual(reline_counts(5.0, 0, 3.2)["W"], 0)

    def test_protein_pdb_filters_water_and_second_model(self):
        with TemporaryDirectory() as directory:
            source, target = Path(directory) / "input.pdb", Path(directory) / "protein.pdb"
            atom = "ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00 34.00           N"
            water = "HETATM    2  O   HOH A 100      27.340  24.430   2.614  1.00 34.00           O"
            source.write_text("MODEL        1\n" + atom + "\n" + water + "\nENDMDL\nMODEL        2\n" + atom + "\nENDMDL\n")
            _clean_protein_pdb(source, target)
            self.assertEqual(sum(line.startswith("ATOM  ") for line in target.read_text().splitlines()), 1)
            self.assertNotIn("HOH", target.read_text())

    def test_insane_ion_labels_match_martini3(self):
        with TemporaryDirectory() as directory:
            gro = Path(directory) / "system.gro"
            gro.write_text("ions\n2\n    1NA+   NA+    1   0.100   0.200   0.300\n"
                           "    2CL-   CL-    2   0.400   0.500   0.600\n"
                           "   1.00000   1.00000   1.00000\n")
            _normalize_insane_ions(gro)
            lines = gro.read_text().splitlines()
            self.assertEqual([(row[5:10].strip(), row[10:15].strip()) for row in lines[2:4]],
                             [("NA", "NA"), ("CL", "CL")])
