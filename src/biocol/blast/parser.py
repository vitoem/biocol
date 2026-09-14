"""Columns and parsing of BLAST tabular ``outfmt 6``."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

OUTFMT6_STANDARD_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]
OUTFMT6_COLUMNS = [*OUTFMT6_STANDARD_COLUMNS, "nident", "qseq", "sseq"]
BLAST_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen "
    "qstart qend sstart send evalue bitscore nident qseq sseq"
)

logger = logging.getLogger(__name__)


def parse_blast_results(path: str | Path) -> pd.DataFrame:
    """Read a BLAST tabular file (outfmt 6) into a DataFrame.

    An empty file (no hits) returns a DataFrame with no rows.
    Accepts classic 12-column outfmt 6 or 15 columns including nident/qseq/sseq.
    """
    blast_path = Path(path)
    if not blast_path.exists():
        raise FileNotFoundError(f"BLAST file not found: {blast_path}")
    logger.info("Reading BLAST results: %s", blast_path.name)
    if blast_path.stat().st_size == 0:
        logger.debug("BLAST tabular is empty: %s", blast_path)
        return pd.DataFrame(columns=OUTFMT6_COLUMNS)

    frame = pd.read_csv(
        blast_path,
        sep="\t",
        header=None,
        comment="#",
    )
    width = frame.shape[1]
    if width >= len(OUTFMT6_COLUMNS):
        frame = frame.iloc[:, : len(OUTFMT6_COLUMNS)]
        frame.columns = OUTFMT6_COLUMNS
    elif width == len(OUTFMT6_STANDARD_COLUMNS):
        frame.columns = OUTFMT6_STANDARD_COLUMNS
        frame["nident"] = pd.NA
        frame["qseq"] = pd.NA
        frame["sseq"] = pd.NA
    else:
        raise ValueError(
            f"BLAST tabular must have 12 or 15 columns, found {width}: {blast_path}"
        )
    logger.debug("Parsed BLAST tabular %s (%s hit rows)", blast_path, len(frame))
    return frame


def fill_missing_hits(
    results: pd.DataFrame,
    query_ids: list[str],
    database_name: str,
) -> pd.DataFrame:
    """Ensure one row per query; if there was no hit, BLAST fields are empty."""
    frame = results.copy()
    if "database" not in frame.columns:
        frame["database"] = database_name
    present = set(frame["qseqid"].astype(str)) if not frame.empty else set()
    missing = [query_id for query_id in query_ids if query_id not in present]
    if not missing:
        return frame
    empty = pd.DataFrame(
        {
            "qseqid": missing,
            "database": database_name,
        }
    )
    for column in OUTFMT6_COLUMNS:
        if column != "qseqid" and column not in empty.columns:
            empty[column] = pd.NA
    return pd.concat([frame, empty], ignore_index=True)
