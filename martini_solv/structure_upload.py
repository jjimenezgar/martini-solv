"""Broaden MartiniSolv protein structure uploads to PDB and mmCIF/CIF.

The Streamlit app historically accepted only ``.pdb`` files and downstream
code expects PDB-formatted bytes.  This module keeps that downstream contract
stable: uploaded mmCIF/CIF structures are converted to PDB in memory before
the existing app sees them.
"""
from __future__ import annotations

from dataclasses import dataclass
import io
from pathlib import Path
import tempfile
from typing import Any


SUPPORTED_STRUCTURE_EXTENSIONS = ("pdb", "cif", "mmcif")
_CIF_SUFFIXES = {".cif", ".mmcif"}


def normalize_structure_bytes(filename: str, data: bytes) -> bytes:
    """Return PDB-formatted bytes for a supported uploaded structure.

    PDB input is passed through unchanged.  mmCIF/CIF input (including
    AlphaFold downloads) is parsed with OpenMM's PDBx reader and serialized as
    PDB so the existing MartiniSolv preparation/build pipeline does not need a
    second structure code path.
    """
    suffix = Path(filename or "").suffix.lower()
    if suffix == ".pdb":
        if not data.strip():
            raise ValueError("The uploaded PDB file is empty.")
        return data
    if suffix not in _CIF_SUFFIXES:
        raise ValueError("Protein structure must be a .pdb, .cif, or .mmcif file.")
    if not data.strip():
        raise ValueError("The uploaded mmCIF/CIF file is empty.")

    try:
        from openmm.app import PDBFile, PDBxFile
    except ImportError as exc:  # pragma: no cover - pdbfixer installs OpenMM in production
        raise RuntimeError("mmCIF/CIF upload support requires OpenMM.") from exc

    # PDBxFile is most robust when given a real path. The temporary input is
    # removed immediately after parsing; no user structure is retained.
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=True) as handle:
        handle.write(data)
        handle.flush()
        try:
            structure = PDBxFile(handle.name)
        except Exception as exc:
            raise ValueError(f"Could not parse uploaded mmCIF/CIF structure: {exc}") from exc

    output = io.StringIO()
    try:
        PDBFile.writeFile(structure.topology, structure.positions, output, keepIds=True)
    except Exception:
        # PDB has stricter identifier limits than mmCIF. For unusual multi-
        # character chain/residue identifiers, let OpenMM generate safe PDB
        # identifiers rather than rejecting an otherwise valid structure.
        output = io.StringIO()
        try:
            PDBFile.writeFile(structure.topology, structure.positions, output, keepIds=False)
        except Exception as exc:
            raise ValueError(f"Could not convert uploaded mmCIF/CIF structure to PDB: {exc}") from exc

    pdb = output.getvalue().encode("utf-8")
    if not any(line.startswith((b"ATOM  ", b"HETATM")) for line in pdb.splitlines()):
        raise ValueError("The uploaded mmCIF/CIF file contains no atom records.")
    return pdb


def _structure_types(existing: Any) -> list[str]:
    """Preserve existing allowed types while adding protein structure formats."""
    if existing is None:
        values: list[str] = []
    elif isinstance(existing, str):
        values = [existing]
    else:
        try:
            values = [str(value) for value in existing]
        except TypeError:
            values = []
    seen = {value.lower().lstrip(".") for value in values}
    for extension in SUPPORTED_STRUCTURE_EXTENSIONS:
        if extension not in seen:
            values.append(extension)
            seen.add(extension)
    return values


@dataclass
class _NormalizedUpload:
    """Small proxy that keeps Streamlit's UploadedFile interface intact."""

    _uploaded: Any
    _data: bytes

    @property
    def name(self) -> str:
        return str(self._uploaded.name)

    @property
    def type(self) -> str:
        return str(getattr(self._uploaded, "type", ""))

    @property
    def size(self) -> int:
        return len(self._data)

    def getvalue(self) -> bytes:
        return self._data

    def read(self, size: int = -1) -> bytes:
        return self._data if size is None or size < 0 else self._data[:size]

    def __getattr__(self, name: str) -> Any:
        return getattr(self._uploaded, name)


def _normalize_uploaded(uploaded: Any) -> Any:
    if uploaded is None:
        return None
    if isinstance(uploaded, (list, tuple)):
        normalized = [_normalize_uploaded(item) for item in uploaded]
        return type(uploaded)(normalized) if isinstance(uploaded, tuple) else normalized
    name = str(getattr(uploaded, "name", ""))
    if Path(name).suffix.lower() not in _CIF_SUFFIXES:
        return uploaded
    return _NormalizedUpload(uploaded, normalize_structure_bytes(name, uploaded.getvalue()))


def install_streamlit_structure_upload_support() -> None:
    """Extend only the app's protein-structure uploader, once per process."""
    import streamlit as st

    current = st.file_uploader
    if getattr(current, "_martinisolv_structure_formats", False):
        return

    def file_uploader(label: str, *args: Any, **kwargs: Any) -> Any:
        if str(label).strip().lower() == "protein structure":
            kwargs["type"] = _structure_types(kwargs.get("type"))
            uploaded = current(label, *args, **kwargs)
            try:
                return _normalize_uploaded(uploaded)
            except (ValueError, RuntimeError) as exc:
                st.error(str(exc))
                return None
        return current(label, *args, **kwargs)

    file_uploader._martinisolv_structure_formats = True  # type: ignore[attr-defined]
    st.file_uploader = file_uploader
