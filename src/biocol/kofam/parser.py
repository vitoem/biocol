"""Parse KofamScan ``detail-tsv`` output."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

KOFAM_DETAIL_COLUMNS = [
    "mark",
    "gene_name",
    "ko",
    "threshold",
    "score",
    "evalue",
    "definition",
]


def parse_kofam_detail_tsv(path: str | Path) -> pd.DataFrame:
    """Read ``detail-tsv`` and keep hits that passed the KO threshold (``*``)."""
    tsv_path = Path(path)
    if not tsv_path.exists():
        raise FileNotFoundError(f"KofamScan TSV not found: {tsv_path}")
    logger.info("Reading KofamScan detail-tsv: %s", tsv_path.name)
    rows: list[list[str]] = []
    with tsv_path.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for parts in reader:
            if not parts or (parts[0].startswith("#") and parts[0] != "*"):
                continue
            if parts[0].strip() != "*":
                continue
            padded = parts + [""] * (len(KOFAM_DETAIL_COLUMNS) - len(parts))
            rows.append(padded[: len(KOFAM_DETAIL_COLUMNS)])
    if not rows:
        return pd.DataFrame(columns=KOFAM_DETAIL_COLUMNS)
    frame = pd.DataFrame(rows, columns=KOFAM_DETAIL_COLUMNS)
    frame["definition"] = frame["definition"].astype(str).str.strip().str.strip('"')
    frame["score"] = pd.to_numeric(frame["score"], errors="coerce")
    frame["evalue"] = pd.to_numeric(frame["evalue"], errors="coerce")
    logger.info("KofamScan hits above threshold: %s row(s)", len(frame))
    return frame
