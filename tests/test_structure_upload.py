from __future__ import annotations

import io

from martini_solv.structure_upload import (
    SUPPORTED_STRUCTURE_EXTENSIONS,
    _structure_types,
    normalize_structure_bytes,
)


PDB_TEXT = """\
ATOM      1  N   ALA A   1      -0.525   1.362   0.000  1.00 20.00           N\n
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n
ATOM      3  C   ALA A   1       1.525   0.000   0.000  1.00 20.00           C\n
ATOM      4  O   ALA A   1       2.050  -1.120   0.000  1.00 20.00           O\n
TER\nEND\n
"""


def test_supported_structure_extensions_include_pdb_and_cif() -> None:
    assert SUPPORTED_STRUCTURE_EXTENSIONS == ("pdb", "cif", "mmcif")
    assert _structure_types(["pdb"]) == ["pdb", "cif", "mmcif"]


def test_pdb_upload_is_preserved_verbatim() -> None:
    data = PDB_TEXT.encode("utf-8")
    assert normalize_structure_bytes("protein.pdb", data) == data


def test_mmcif_upload_is_converted_to_pdb() -> None:
    from openmm.app import PDBFile, PDBxFile

    pdb = PDBFile(io.StringIO(PDB_TEXT))
    cif = io.StringIO()
    PDBxFile.writeFile(pdb.topology, pdb.positions, cif, keepIds=True)

    converted = normalize_structure_bytes("alphafold_model.cif", cif.getvalue().encode("utf-8"))
    text = converted.decode("utf-8")

    assert "ATOM" in text
    assert " ALA " in text
    assert " A   1" in text
    assert text.rstrip().endswith("END")


def test_empty_cif_is_rejected() -> None:
    try:
        normalize_structure_bytes("model.cif", b"")
    except ValueError as exc:
        assert "empty" in str(exc).lower()
    else:
        raise AssertionError("Expected an empty CIF to be rejected")
