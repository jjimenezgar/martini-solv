import unittest

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
