"""Resolve Diamond databases: protein FASTA or ``.dmnd`` files (not mixed)."""

from __future__ import annotations

import logging
from pathlib import Path

from biocol.blast.databases import infer_database_label, list_blast_databases
from biocol.exceptions import DatabaseError, DiamondError, MixedDatabaseTypeError
from biocol.sequence.validator import FASTA_EXTENSIONS

logger = logging.getLogger(__name__)

DMND_EXTENSION = ".dmnd"

DiamondEntry = tuple[Path, str, str]  # path, kind (fasta|dmnd), label


def _is_dmnd(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() == DMND_EXTENSION


def _is_fasta(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS


def _iter_in_directory(directory: Path) -> tuple[list[Path], list[Path]]:
    fasta_files = sorted(path for path in directory.rglob("*") if _is_fasta(path))
    dmnd_files = sorted(path for path in directory.rglob("*") if _is_dmnd(path))
    return fasta_files, dmnd_files


def _fasta_entries(paths: list[Path]) -> list[DiamondEntry]:
    entries: list[DiamondEntry] = []
    types: set[str] = set()
    for path in paths:
        listed = list_blast_databases(path)
        fasta_path, db_type = listed[0]
        types.add(db_type)
        if db_type != "protein":
            raise DiamondError(
                "diamond requires a protein query and a protein database"
            )
        entries.append((fasta_path, "fasta", infer_database_label(fasta_path)))
    if len(types) > 1:
        raise MixedDatabaseTypeError(
            "Input mixes nucleotide and protein databases"
        )
    return entries


def list_diamond_databases(source: str | Path) -> list[DiamondEntry]:
    """Return ``(path, kind, label)`` for a FASTA, a ``.dmnd``, or a folder.

    A folder must contain only protein FASTA files or only ``.dmnd`` files.
    """
    raw = Path(source)
    if raw.is_dir():
        logger.info("Scanning Diamond databases in folder: %s", raw)
        fasta_files, dmnd_files = _iter_in_directory(raw)
        if fasta_files and dmnd_files:
            raise DatabaseError(
                "Input mixes FASTA and Diamond (.dmnd) databases"
            )
        if dmnd_files:
            logger.info("Found %s Diamond database file(s)", len(dmnd_files))
            return [(path, "dmnd", path.stem) for path in dmnd_files]
        if fasta_files:
            return _fasta_entries(fasta_files)
        raise DatabaseError(
            f"Folder contains no FASTA or .dmnd files: {raw}"
        )

    if _is_dmnd(raw):
        logger.info("Using Diamond database: %s", raw)
        return [(raw, "dmnd", raw.stem)]

    if _is_fasta(raw):
        return _fasta_entries([raw])

    if raw.is_file():
        raise DatabaseError(
            f"Diamond database must be FASTA or .dmnd, not {raw.suffix}: {raw}"
        )
    raise DatabaseError(
        f"Unrecognized Diamond database: {source}. "
        "Pass a protein FASTA, a .dmnd file, or a folder of one type."
    )
