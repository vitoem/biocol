from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

from biocol.exceptions import MetadataError

logger = logging.getLogger(__name__)

_DB_PREFIXES = {
    "ref",
    "gb",
    "emb",
    "dbj",
    "sp",
    "pdb",
    "tpg",
    "tpe",
    "tpd",
    "lcl",
    "gi",
}


def normalize_accession(value: object) -> str:
    """Strip prefixes such as ``ref|XP_123.1|`` and keep the accession."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    if "|" not in text:
        return text.split()[0]
    tokens = [token for token in text.split("|") if token]
    for token in reversed(tokens):
        lowered = token.lower()
        if lowered in _DB_PREFIXES:
            continue
        if re.fullmatch(r"\d+", token):
            continue
        return token.split()[0]
    return tokens[-1]


def load_accessions(path: str | Path) -> pd.DataFrame:
    """Read ``accession<TAB>descriptor`` (no header)."""
    acc_path = Path(path)
    if not acc_path.exists():
        raise FileNotFoundError(f"Accessions file not found: {acc_path}")
    try:
        frame = pd.read_csv(
            acc_path,
            sep="\t",
            header=None,
            names=["accession", "description"],
            dtype=str,
            comment="#",
        )
    except Exception as exc:
        raise MetadataError(f"Could not read accessions: {acc_path}") from exc

    frame = frame.dropna(how="all")
    if frame.empty:
        raise MetadataError(f"Accessions file has no rows: {acc_path}")
    if frame["accession"].isna().any() or (frame["accession"].str.strip() == "").any():
        raise MetadataError("Accessions file has rows without an accession")

    frame["accession"] = frame["accession"].str.strip()
    frame["description"] = frame["description"].fillna("").str.strip()
    frame["accession_norm"] = frame["accession"].map(normalize_accession)
    logger.info("Loaded %s accession descriptor(s) from %s", len(frame), acc_path)
    return frame
