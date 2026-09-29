import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from martini_solv.builder import _clean_protein_pdb, _normalize_insane_ions, _name_molecule_type, _topology
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

    def test_mapper_molecule_type_is_renamed(self):
        with TemporaryDirectory() as directory:
            source, target = Path(directory) / "original.itp", Path(directory) / "BENZ.itp"
            source.write_text("[ moleculetype ]\n; name nrexcl\nres 1\n[ atoms ]\n1 C1 1 res C1 1 0\n")
            _name_molecule_type(source, target, "BENZ")
            self.assertIn("BENZ 1", target.read_text())
            self.assertIn("1 C1 1 res C1 1 0", target.read_text())

    def test_go_topology_preserves_martinize_directives_before_protein(self):
        with TemporaryDirectory() as directory:
            work = Path(directory)
            (work / "martini_v3.0.0.itp").write_text("[ defaults ]\\n1 1 yes 1.0 1.0\\n")
            (work / "martini_v3.0.0_ions_v1.itp").write_text("; ions\\n")
            (work / "martini_v3.0.0_solvents_v1.itp").write_text("; solvents\\n")
            (work / "go_atomtypes.itp").write_text("[ atomtypes ]\\nProtein_1 72 0 A 0.0 0.0\\n")
            (work / "Protein.itp").write_text("[ moleculetype ]\\nProtein 1\\n[ atoms ]\\n1 Protein_1 1 MET BB 1 0\\n")
            (work / "protein.top").write_text(
                "#define GO_VIRT\\n"
                "#include \\"martini_v3.0.0.itp\\"\\n"
                "#include \\"go_atomtypes.itp\\"\\n"
                "#include \\"Protein.itp\\"\\n\\n"
                "[ system ]\\nProtein\\n\\n"
                "[ molecules ]\\nProtein 1\\n"
            )
            _topology(work, [("Protein", 1)], [("W", 20)], go_enabled=True)
            text = (work / "system.top").read_text()
            self.assertLess(text.index("#define GO_VIRT"), text.index('#include "Protein.itp"'))
            self.assertLess(text.index('#include "go_atomtypes.itp"'), text.index('#include "Protein.itp"'))
            self.assertEqual(text.count('#include "Protein.itp"'), 1)
            self.assertIn("Protein          1", text)
            self.assertIn("W                20", text)
