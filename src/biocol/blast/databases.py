from __future__ import annotations

import logging
import re
from collections import Counter
from pathlib import Path

from Bio import SeqIO

from biocol.exceptions import DatabaseError, MixedDatabaseTypeError
from biocol.sequence.classifier import detect_query_type
from biocol.sequence.reader import read_fasta
from biocol.sequence.validator import FASTA_EXTENSIONS

logger = logging.getLogger(__name__)

SequenceType = str
_ORGANISM_BRACKET = re.compile(r"\[([^\]]+)\]\s*$")


def infer_database_label(fasta_path: Path, sample: int = 40) -> str:
    """Species/database name: organism from NCBI headers ``[Genus species]``.

    If there is no majority organism, use the FASTA file stem.
    """
    organisms: list[str] = []
    with fasta_path.open(encoding="utf-8", errors="replace") as handle:
        for index, record in enumerate(SeqIO.parse(handle, "fasta")):
            if index >= sample:
                break
            match = _ORGANISM_BRACKET.search(record.description)
            if match:
                organisms.append(match.group(1).strip())
    if organisms:
        name, count = Counter(organisms).most_common(1)[0]
        if count >= max(1, (len(organisms) + 1) // 2):
            logger.debug("infer_database_label: %s → %s", fasta_path.name, name)
            return name
    return fasta_path.stem


def _is_fasta_file(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS


def _type_from_fasta(path: Path) -> SequenceType:
    records = read_fasta(path)
    return detect_query_type(records)


def _iter_fasta_in_directory(directory: Path) -> list[Path]:
    files = [
        path
        for path in directory.rglob("*")
        if _is_fasta_file(path)
    ]
    return sorted(files)


def list_blast_databases(source: str | Path) -> list[tuple[Path, SequenceType]]:
    """Resolve a FASTA file or a folder (including subfolders).

    Returns ``(path, type)`` pairs.
    """
    raw = Path(source)

    if raw.is_dir():
        logger.info("Scanning FASTA databases in folder: %s", raw)
        fasta_files = _iter_fasta_in_directory(raw)
        if not fasta_files:
            raise DatabaseError(
                f"Folder contains no FASTA files ({', '.join(sorted(FASTA_EXTENSIONS))}): {raw}"
            )
        logger.info("Found %s FASTA database file(s)", len(fasta_files))
        return [(path, _type_from_fasta(path)) for path in fasta_files]

    if _is_fasta_file(raw):
        logger.info("Using FASTA database: %s", raw)
        return [(raw, _type_from_fasta(raw))]

    if raw.is_file():
        raise DatabaseError(
            f"Database must be FASTA ({', '.join(sorted(FASTA_EXTENSIONS))}), not {raw.suffix}: {raw}"
        )

    raise DatabaseError(
        f"Unrecognized database: {source}. "
        "Pass a FASTA file or a folder of FASTA files."
    )


def detect_database_type(source: str | Path) -> SequenceType:
    """Type of a FASTA database: ``nucleotide`` or ``protein``.

    Accepts a FASTA file or a folder (and subfolders) of FASTA files
    of the same type.
    """
    entries = list_blast_databases(source)
    types = [db_type for _, db_type in entries]
    unique = set(types)
    if len(unique) > 1:
        logger.debug(
            "detect_database_type: mixed databases %s",
            [(str(path), db_type) for path, db_type in entries],
        )
        raise MixedDatabaseTypeError(
            "Input mixes nucleotide and protein databases"
        )
    database_type = types[0]
    logger.info(
        "Database type: %s (%s FASTA file(s))",
        database_type,
        len(entries),
    )
    return database_type
