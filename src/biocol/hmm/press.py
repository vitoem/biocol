"""Press an HMMER profile database with hmmpress if indexes are missing."""

from __future__ import annotations

import logging
from pathlib import Path

from biocol.hmm.execute import run_hmmer

logger = logging.getLogger(__name__)

PRESSED_SUFFIXES = (".h3m", ".h3i", ".h3f", ".h3p")


def pressed_index_paths(hmm_path: Path) -> list[Path]:
    return [Path(str(hmm_path) + suffix) for suffix in PRESSED_SUFFIXES]


def hmm_database_is_pressed(hmm_path: Path) -> bool:
    return all(path.exists() for path in pressed_index_paths(hmm_path))


def hmm_has_gathering_threshold(hmm_path: Path) -> bool:
    """True if the HMM file defines Pfam-style ``GA`` cutoffs."""
    with hmm_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("GA "):
                return True
    return False


def ensure_pressed_hmm(hmm_path: str | Path) -> Path:
    """Run ``hmmpress`` when binary indexes are absent. Returns the HMM path."""
    path = Path(hmm_path)
    if not path.is_file():
        raise FileNotFoundError(f"HMM database not found: {path}")
    if hmm_database_is_pressed(path):
        logger.info("HMM database already pressed: %s", path.name)
        return path
    logger.info("Running hmmpress on %s", path)
    run_hmmer(["hmmpress", str(path)])
    return path
