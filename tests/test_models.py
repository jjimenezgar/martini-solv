import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from martini_solv.builder import _assign_itp_net_charge, _clean_protein_pdb, _gro_from_itp, _itp_net_charge, _normalize_insane_ions, _name_molecule_type, _prepare_uploaded_solute, _topology, _verify_existing_coordinates_preserved, _verify_free_molecule_topology_counts, _verify_reline_composition
from martini_solv.models import BuildConfig, Solute, chcl_sorbitol_counts, reline_counts
from martini_solv.molecular_viewer import _gif_component_map, _gif_frame_indices


class TestModels(unittest.TestCase):
    def test_reline_stoichiometry_and_actual_water_fraction(self):
        counts = reline_counts(5.0, 0.10)
        self.assertEqual(counts["UREA"], 2 * counts["CHOL"])
        self.assertEqual(counts["CL"], counts["CHOL"])
        self.assertLess(abs(counts["x_water_actual"] - 0.10), 0.002)


    def test_reline_disallows_added_salt(self):
        with self.assertRaisesRegex(ValueError, "only for water"):
            BuildConfig(solvent="reline", salt_m=0.15).validate()


    def test_dry_reline_has_no_water(self):
        self.assertEqual(reline_counts(5.0, 0)["W"], 0)

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
            (work / "martini_v3.0.0.itp").write_text("[ defaults ]\n1 1 yes 1.0 1.0\n")
            (work / "martini_v3.0.0_ions_v1.itp").write_text("; ions\n")
            (work / "martini_v3.0.0_solvents_v1.itp").write_text("; solvents\n")
            (work / "go_atomtypes.itp").write_text("[ atomtypes ]\nProtein_1 72 0 A 0.0 0.0\n")
            (work / "go_nbparams.itp").write_text("[ nonbond_params ]\nProtein_1 Protein_1 1 0.47 5.0\n")
            (work / "Protein.itp").write_text("[ moleculetype ]\nProtein 1\n[ atoms ]\n1 Protein_1 1 MET BB 1 0\n")
            (work / "protein.top").write_text(
                '#define GO_VIRT\n'
                '#include "martini_v3.0.0.itp"\n'
                '#include "go_atomtypes.itp"\n'
                '#include "Protein.itp"\n\n'
                '[ system ]\nProtein\n\n'
                '[ molecules ]\nProtein 1\n'
            )
            _topology(work, [("Protein", 1)], [("W", 20)], go_enabled=True)
            text = (work / "system.top").read_text()
            self.assertLess(text.index("#define GO_VIRT"), text.index('#include "Protein.itp"'))
            self.assertLess(text.index('#include "go_atomtypes.itp"'), text.index('#include "go_nbparams.itp"'))
            self.assertLess(text.index('#include "go_nbparams.itp"'), text.index('#include "martini_v3.0.0_ions_v1.itp"'))
            self.assertLess(text.index('#include "go_nbparams.itp"'), text.index('#include "martini_v3.0.0_solvents_v1.itp"'))
            self.assertLess(text.index('#include "go_nbparams.itp"'), text.index('#include "Protein.itp"'))
            self.assertEqual(text.count('#include "go_atomtypes.itp"'), 1)
            self.assertEqual(text.count('#include "go_nbparams.itp"'), 1)
            self.assertEqual(text.count('#include "Protein.itp"'), 1)
            self.assertIn("Protein          1", text)
            self.assertIn("W                20", text)

    def test_generic_martini_include_is_normalized(self):
        with TemporaryDirectory() as directory:
            work = Path(directory)
            (work / "martini_v3.0.0.itp").write_text("[ defaults ]\n1 1 yes 1.0 1.0\n")
            (work / "martini_v3.0.0_ions_v1.itp").write_text("; ions\n")
            (work / "martini_v3.0.0_solvents_v1.itp").write_text("; solvents\n")
            (work / "Protein.itp").write_text("[ moleculetype ]\nProtein 1\n[ atoms ]\n1 P5 1 MET BB 1 0\n")
            (work / "protein.top").write_text(
                '#include "martini.itp"\n'
                '#include "Protein.itp"\n\n'
                '[ system ]\nProtein\n\n'
                '[ molecules ]\nProtein 1\n'
            )
            _topology(work, [("Protein", 1)], [("W", 20)], go_enabled=False)
            text = (work / "system.top").read_text()
            self.assertNotIn('#include "martini.itp"', text)
            self.assertIn('#include "martini_v3.0.0.itp"', text)
            self.assertEqual(text.count('#include "martini_v3.0.0.itp"'), 1)


    def test_solvation_preserves_inserted_coordinate_block(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            before = root / "before.gro"
            after = root / "after.gro"
            before.write_text(
                "before\n2\n"
                "    1PROT   BB    1   0.100   0.100   0.100\n"
                "    2LIG    C1    2   0.200   0.200   0.200\n"
                "   1.00000   1.00000   1.00000\n"
            )
            after.write_text(
                "after\n3\n"
                "    1PROT   BB    1   0.100   0.100   0.100\n"
                "    2LIG    C1    2   0.200   0.200   0.200\n"
                "    3W      W    3   0.300   0.300   0.300\n"
                "   1.00000   1.00000   1.00000\n"
            )
            _verify_existing_coordinates_preserved(before, after)

    def test_solvation_detects_missing_inserted_molecule(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            before = root / "before.gro"
            after = root / "after.gro"
            before.write_text(
                "before\n2\n"
                "    1PROT   BB    1   0.100   0.100   0.100\n"
                "    2LIG    C1    2   0.200   0.200   0.200\n"
                "   1.00000   1.00000   1.00000\n"
            )
            after.write_text(
                "after\n2\n"
                "    1PROT   BB    1   0.100   0.100   0.100\n"
                "    2W      W    2   0.300   0.300   0.300\n"
                "   1.00000   1.00000   1.00000\n"
            )
            with self.assertRaisesRegex(RuntimeError, "did not preserve"):
                _verify_existing_coordinates_preserved(before, after)


    def test_reline_coordinate_composition_accepts_packmol_truncated_names(self):
        with TemporaryDirectory() as directory:
            gro = Path(directory) / "system.gro"
            gro.write_text(
                "reline\n8\n"
                "    1CHO    N1    1   0.100   0.100   0.100\n"
                "    1CHO    OH    2   0.200   0.100   0.100\n"
                "    2URE    N1    3   0.300   0.100   0.100\n"
                "    2URE    UP    4   0.400   0.100   0.100\n"
                "    2URE    UN    5   0.500   0.100   0.100\n"
                "    3URE    N1    6   0.600   0.100   0.100\n"
                "    4CL     CL    7   0.700   0.100   0.100\n"
                "    5W       W    8   0.800   0.100   0.100\n"
                "   2.00000   2.00000   2.00000\n"
            )
            expected = {"CHOL": 1, "UREA": 2, "SOR": 0, "CL": 1, "NA": 0, "W": 1}
            actual = _verify_reline_composition(gro, expected)
            self.assertEqual(actual, expected)

    def test_reline_uses_experimental_density(self):
        counts = reline_counts(1.0, 0.0)
        self.assertAlmostEqual(float(counts["target_density_g_cm3"]), 1.20, places=6)
        self.assertAlmostEqual(float(counts["dry_formula_units_per_nm3"]), 2.782, places=3)



    def test_reline_composition_accepts_extra_counterions(self):
        with TemporaryDirectory() as directory:
            gro = Path(directory) / "system.gro"
            gro.write_text(
                "charged reline\n10\n"
                "    1CHO    N1    1   0.100   0.100   0.100\n"
                "    1CHO    OH    2   0.200   0.100   0.100\n"
                "    2URE    N1    3   0.300   0.100   0.100\n"
                "    2URE    UP    4   0.400   0.100   0.100\n"
                "    2URE    UN    5   0.500   0.100   0.100\n"
                "    3URE    N1    6   0.600   0.100   0.100\n"
                "    4CL     CL    7   0.700   0.100   0.100\n"
                "    5CL     CL    8   0.800   0.100   0.100\n"
                "    6CL     CL    9   0.900   0.100   0.100\n"
                "    7W       W   10   1.000   0.100   0.100\n"
                "   2.00000   2.00000   2.00000\n"
            )
            expected = {"CHOL": 1, "UREA": 2, "SOR": 0, "CL": 3, "NA": 0, "W": 1}
            actual = _verify_reline_composition(gro, expected)
            self.assertEqual(actual, expected)


    def test_chcl_sorbitol_stoichiometry_and_wet_fraction(self):
        counts = chcl_sorbitol_counts(5.0, 0.10, 1.20)
        self.assertEqual(counts["SOR"], counts["CHOL"])
        self.assertEqual(counts["CL"], counts["CHOL"])
        self.assertLess(abs(counts["x_water_actual"] - 0.10), 0.003)

    def test_dry_chcl_sorbitol_has_no_water(self):
        self.assertEqual(chcl_sorbitol_counts(5.0, 0.0, 1.20)["W"], 0)

    def test_chcl_sorbitol_config_accepts_des_mode(self):
        BuildConfig(
            solvent="chcl_sorbitol",
            salt_m=0.0,
            water_fraction=0.10,
            chcl_sorbitol_density_g_cm3=1.20,
        ).validate()


    def test_manual_free_molecule_charge_assignment(self):
        with TemporaryDirectory() as directory:
            itp = Path(directory) / "AMP.itp"
            itp.write_text(
                "[ moleculetype ]\nAMP 1\n\n"
                "[ atoms ]\n"
                "1 P4 1 AMP B1 1 0.0 72.0\n"
                "2 P4 1 AMP B2 2 0.0 72.0\n"
                "3 P4 1 AMP B3 3 0.0 72.0\n"
            )
            _assign_itp_net_charge(itp, 1, 3)
            self.assertAlmostEqual(_itp_net_charge(itp), 1.0)
            atom_lines = [
                line.split()
                for line in itp.read_text().splitlines()
                if line.strip() and line.strip()[0].isdigit()
            ]
            self.assertEqual(float(atom_lines[2][6]), 1.0)

    def test_negative_free_molecule_charge_assignment(self):
        with TemporaryDirectory() as directory:
            itp = Path(directory) / "LIG.itp"
            itp.write_text(
                "[ moleculetype ]\nLIG 1\n\n"
                "[ atoms ]\n"
                "1 P4 1 LIG B1 1 0.0 72.0\n"
                "2 P4 1 LIG B2 2 0.0 72.0\n"
            )
            _assign_itp_net_charge(itp, -1, 1)
            self.assertAlmostEqual(_itp_net_charge(itp), -1.0)

    def test_solute_charge_metadata_defaults(self):
        solute = Solute("AMP", "CCO", 1)
        self.assertEqual(solute.net_charge, 0)
        self.assertIsNone(solute.charged_bead)

    def test_solute_charge_metadata_records_assignment(self):
        solute = Solute("AMP", "CCO", 3, 1, 2)
        self.assertEqual(solute.net_charge, 1)
        self.assertEqual(solute.charged_bead, 2)


    def test_free_molecule_copy_count_is_verified_in_system_top(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            gro = root / "AMP.gro"
            gro.write_text(
                "AMP\n2\n"
                "    1AMP    B1    1   0.100   0.100   0.100\n"
                "    1AMP    B2    2   0.200   0.100   0.100\n"
                "   1.00000   1.00000   1.00000\n"
            )
            itp = root / "AMP.itp"
            itp.write_text(
                "[ moleculetype ]\nAMP 1\n\n"
                "[ atoms ]\n1 P4 1 AMP B1 1 0.0 72.0\n2 P4 1 AMP B2 2 0.0 72.0\n"
            )
            top = root / "system.top"
            top.write_text("[ molecules ]\nProtein 1\nAMP 5\n")
            spec = Solute("AMP", "CCO", 5)
            rows = _verify_free_molecule_topology_counts(top, [(spec, gro, itp)])
            self.assertEqual(rows[0]["requested"], 5)
            self.assertEqual(rows[0]["included"], 5)
            self.assertEqual(rows[0]["beads_per_molecule"], 2)

    def test_free_molecule_copy_count_mismatch_fails(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            gro = root / "AMP.gro"
            gro.write_text(
                "AMP\n1\n"
                "    1AMP    B1    1   0.100   0.100   0.100\n"
                "   1.00000   1.00000   1.00000\n"
            )
            itp = root / "AMP.itp"
            itp.write_text(
                "[ moleculetype ]\nAMP 1\n\n"
                "[ atoms ]\n1 P4 1 AMP B1 1 0.0 72.0\n"
            )
            top = root / "system.top"
            top.write_text("[ molecules ]\nProtein 1\nAMP 1\n")
            spec = Solute("AMP", "CCO", 5)
            with self.assertRaisesRegex(RuntimeError, "requested 5 copies"):
                _verify_free_molecule_topology_counts(top, [(spec, gro, itp)])


    def test_uploaded_itp_generates_gro_from_connectivity(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.itp"
            source.write_text(
                "[ moleculetype ]\nOLD 1\n\n"
                "[ atoms ]\n"
                "1 P4 1 LIG B1 1 0.0 72.0\n"
                "2 P4 1 LIG B2 2 0.0 72.0\n"
                "3 P4 1 LIG B3 3 0.0 72.0\n\n"
                "[ constraints ]\n"
                "1 2 1 0.33\n"
                "2 3 1 0.34\n"
            )
            work = root / "work"
            work.mkdir()
            gro, itp = _prepare_uploaded_solute(work, "LIG", source)
            self.assertTrue(gro.is_file())
            self.assertTrue(itp.is_file())
            self.assertEqual(int(gro.read_text().splitlines()[1]), 3)
            self.assertIn("LIG 1", itp.read_text())

    def test_uploaded_itp_preserves_existing_charge(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "charged.itp"
            source.write_text(
                "[ moleculetype ]\nION 1\n\n"
                "[ atoms ]\n"
                "1 Q1 1 ION Q1 1 1.0 72.0\n"
            )
            work = root / "work"
            work.mkdir()
            _gro, itp = _prepare_uploaded_solute(work, "ION", source)
            self.assertAlmostEqual(_itp_net_charge(itp), 1.0)

    def test_itp_only_geometry_rejects_disconnected_multibead_model(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            itp = root / "disconnected.itp"
            itp.write_text(
                "[ moleculetype ]\nLIG 1\n\n"
                "[ atoms ]\n"
                "1 P4 1 LIG B1 1 0.0 72.0\n"
                "2 P4 1 LIG B2 2 0.0 72.0\n"
            )
            with self.assertRaisesRegex(ValueError, "no \[ bonds \] or \[ constraints \]"):
                _gro_from_itp(itp, root / "LIG.gro", "LIG")

    def test_uploaded_solute_requires_only_itp(self):
        with TemporaryDirectory() as directory:
            itp = Path(directory) / "LIG.itp"
            itp.write_text(
                "[ moleculetype ]\nLIG 1\n"
                "[ atoms ]\n1 P4 1 LIG B1 1 0.0 72.0\n"
            )
            Solute(
                name="LIG",
                smiles="",
                count=5,
                source="upload",
                template_itp=str(itp),
            ).validate()


    def test_gif_frame_indices_are_evenly_sampled(self):
        indices = _gif_frame_indices(101, max_frames=5)
        self.assertEqual(indices, [0, 25, 50, 75, 100])

    def test_gif_component_map_matches_viewer_categories(self):
        with TemporaryDirectory() as directory:
            gro = Path(directory) / "system.gro"
            gro.write_text(
                "gif categories\n8\n"
                "    1ALA    BB    1   0.100   0.100   0.100\n"
                "    1ALA   SC1    2   0.200   0.100   0.100\n"
                "    2LIG    C1    3   0.300   0.100   0.100\n"
                "    3W       W    4   0.400   0.100   0.100\n"
                "    4CHO    N1    5   0.500   0.100   0.100\n"
                "    5CL     CL    6   0.600   0.100   0.100\n"
                "    6CL     CL    7   0.700   0.100   0.100\n"
                "    7NA     NA    8   0.800   0.100   0.100\n"
                "   2.00000   2.00000   2.00000\n"
            )
            groups = _gif_component_map(gro, reline_chloride_count=1)
            self.assertEqual(
                groups,
                ["protein", "protein", "solute", "solvent", "solvent", "solvent", "ions", "ions"],
            )
