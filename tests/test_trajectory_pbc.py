import unittest

try:
    import numpy as np
except ImportError:  # Lightweight CI intentionally does not install runtime scientific deps.
    np = None

from martini_solv.trajectory_pbc import _make_protein_whole, _minimum_image, _rigid_align


def _atoms():
    return [
        {"resid": 1, "resn": "ALA", "name": "BB"},
        {"resid": 1, "resn": "ALA", "name": "SC1"},
        {"resid": 2, "resn": "GLY", "name": "BB"},
        {"resid": 3, "resn": "LEU", "name": "BB"},
        {"resid": 3, "resn": "LEU", "name": "SC1"},
        {"resid": 4, "resn": "W", "name": "W"},
    ]


@unittest.skipUnless(np is not None, "NumPy is a runtime dependency, not installed in lightweight CI")
class TestTrajectoryPBC(unittest.TestCase):
    def test_minimum_image_wraps_across_orthorhombic_box(self):
        box = np.diag([2.0, 2.0, 2.0])
        delta = np.array([-1.7, 0.0, 0.0])
        wrapped = _minimum_image(delta, box)
        self.assertTrue(np.allclose(wrapped, [0.3, 0.0, 0.0]))

    def test_make_protein_whole_repairs_backbone_and_sidechain_pbc_split(self):
        atoms = _atoms()
        components = ["protein", "protein", "protein", "protein", "protein", "solvent"]
        xyz = np.array(
            [[
                [1.80, 1.00, 1.00],
                [1.95, 1.00, 1.00],
                [0.10, 1.00, 1.00],
                [0.45, 1.00, 1.00],
                [0.55, 1.10, 1.00],
                [0.05, 0.20, 0.20],
            ]],
            dtype=float,
        )
        boxes = np.array([np.diag([2.0, 2.0, 2.0])])

        whole = _make_protein_whole(xyz, boxes, atoms, components)

        self.assertTrue(np.isclose(np.linalg.norm(whole[0, 2] - whole[0, 0]), 0.30))
        self.assertTrue(np.isclose(np.linalg.norm(whole[0, 3] - whole[0, 2]), 0.35))
        self.assertTrue(np.isclose(np.linalg.norm(whole[0, 1] - whole[0, 0]), 0.15))
        self.assertTrue(np.allclose(whole[0, 5], xyz[0, 5]))

    def test_rigid_align_removes_translation_without_changing_internal_distances(self):
        frame0 = np.array([[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.4, 0.3, 0.0]])
        frame1 = frame0 + np.array([4.0, -2.0, 1.0])
        xyz = np.stack([frame0, frame1])

        aligned = _rigid_align(xyz, [0, 1, 2])

        self.assertTrue(np.allclose(aligned[1], frame0, atol=1e-7))
        self.assertTrue(
            np.isclose(
                np.linalg.norm(aligned[1, 2] - aligned[1, 0]),
                np.linalg.norm(frame1[2] - frame1[0]),
            )
        )


if __name__ == "__main__":
    unittest.main()
