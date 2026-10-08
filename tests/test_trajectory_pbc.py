import numpy as np

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


def test_minimum_image_wraps_across_orthorhombic_box():
    box = np.diag([2.0, 2.0, 2.0])
    delta = np.array([-1.7, 0.0, 0.0])
    wrapped = _minimum_image(delta, box)
    assert np.allclose(wrapped, [0.3, 0.0, 0.0])


def test_make_protein_whole_repairs_backbone_and_sidechain_pbc_split():
    atoms = _atoms()
    components = ["protein", "protein", "protein", "protein", "protein", "solvent"]
    xyz = np.array(
        [[
            [1.80, 1.00, 1.00],  # residue 1 BB
            [1.95, 1.00, 1.00],  # residue 1 SC
            [0.10, 1.00, 1.00],  # residue 2 BB across PBC
            [0.45, 1.00, 1.00],  # residue 3 BB
            [0.55, 1.10, 1.00],  # residue 3 SC
            [0.05, 0.20, 0.20],  # water: should remain untouched
        ]],
        dtype=float,
    )
    boxes = np.array([np.diag([2.0, 2.0, 2.0])])

    whole = _make_protein_whole(xyz, boxes, atoms, components)

    assert np.isclose(np.linalg.norm(whole[0, 2] - whole[0, 0]), 0.30)
    assert np.isclose(np.linalg.norm(whole[0, 3] - whole[0, 2]), 0.35)
    assert np.isclose(np.linalg.norm(whole[0, 1] - whole[0, 0]), 0.15)
    assert np.allclose(whole[0, 5], xyz[0, 5])


def test_rigid_align_removes_translation_without_changing_internal_distances():
    frame0 = np.array([[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.4, 0.3, 0.0]])
    frame1 = frame0 + np.array([4.0, -2.0, 1.0])
    xyz = np.stack([frame0, frame1])

    aligned = _rigid_align(xyz, [0, 1, 2])

    assert np.allclose(aligned[1], frame0, atol=1e-7)
    assert np.isclose(
        np.linalg.norm(aligned[1, 2] - aligned[1, 0]),
        np.linalg.norm(frame1[2] - frame1[0]),
    )
