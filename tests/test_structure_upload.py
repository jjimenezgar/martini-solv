from __future__ import annotations

import importlib.util
import io
import unittest

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


class StructureUploadTests(unittest.TestCase):
    def test_supported_structure_extensions_include_pdb_and_cif(self) -> None:
        self.assertEqual(SUPPORTED_STRUCTURE_EXTENSIONS, ("pdb", "cif", "mmcif"))
        self.assertEqual(_structure_types(["pdb"]), ["pdb", "cif", "mmcif"])

    def test_pdb_upload_is_preserved_verbatim(self) -> None:
        data = PDB_TEXT.encode("utf-8")
        self.assertEqual(normalize_structure_bytes("protein.pdb", data), data)

    @unittest.skipUnless(importlib.util.find_spec("openmm"), "OpenMM not installed in lightweight CI")
    def test_mmcif_upload_is_converted_to_pdb(self) -> None:
        from openmm.app import PDBFile, PDBxFile

        pdb = PDBFile(io.StringIO(PDB_TEXT))
        cif = io.StringIO()
        PDBxFile.writeFile(pdb.topology, pdb.positions, cif, keepIds=True)

        converted = normalize_structure_bytes("alphafold_model.cif", cif.getvalue().encode("utf-8"))
        text = converted.decode("utf-8")

        self.assertIn("ATOM", text)
        self.assertIn(" ALA ", text)
        self.assertIn(" A   1", text)
        self.assertTrue(text.rstrip().endswith("END"))

    def test_empty_cif_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "empty"):
            normalize_structure_bytes("model.cif", b"")


if __name__ == "__main__":
    unittest.main()
