"""Parse hmmscan ``--tblout`` (space-delimited, ``#`` comments)."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

HMMSCAN_COLUMNS = [
    "target_name",
    "target_accession",
    "query_name",
    "query_accession",
    "evalue",
    "score",
    "bias",
    "dom_evalue",
    "dom_score",
    "dom_bias",
    "exp",
    "reg",
    "clu",
    "ov",
    "env",
    "dom",
    "rep",
    "inc",
    "description",
]


def parse_hmmscan_tblout(path: str | Path) -> pd.DataFrame:
    """Read an hmmscan tblout file into a DataFrame."""
    tbl_path = Path(path)
    if not tbl_path.exists():
        raise FileNotFoundError(f"hmmscan tblout not found: {tbl_path}")
    logger.info("Reading hmmscan tblout: %s", tbl_path.name)
    rows: list[list[str]] = []
    with tbl_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split(None, 18)
            if len(parts) < 18:
                continue
            if len(parts) < 19:
                parts.append("")
            rows.append(parts[:19])
    if not rows:
        return pd.DataFrame(columns=HMMSCAN_COLUMNS)
    frame = pd.DataFrame(rows, columns=HMMSCAN_COLUMNS)
    frame["evalue"] = pd.to_numeric(frame["evalue"], errors="coerce")
    frame["score"] = pd.to_numeric(frame["score"], errors="coerce")
    return frame


def best_hmmscan_hits(hits: pd.DataFrame) -> pd.DataFrame:
    """One row per query: lowest full-sequence e-value, then highest score."""
    if hits.empty:
        return hits.copy()
    ordered = hits.sort_values(
        ["query_name", "evalue", "score"],
        ascending=[True, True, False],
        na_position="last",
    )
    return ordered.drop_duplicates(subset=["query_name"], keep="first").reset_index(
        drop=True
    )
