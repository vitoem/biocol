"""Run hmmscan of a protein FASTA against a pressed HMM database."""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path

import pandas as pd

from biocol.exceptions import HmmError
from biocol.hmm.execute import run_hmmer
from biocol.hmm.parser import parse_hmmscan_tblout
from biocol.hmm.press import ensure_pressed_hmm, hmm_has_gathering_threshold
from biocol.sequence.classifier import detect_query_type
from biocol.sequence.reader import read_fasta

logger = logging.getLogger(__name__)

DEFAULT_HMM_DIR = "hmm"
DEFAULT_NUM_THREADS = 1
DEFAULT_EVALUE = 10


def build_hmmscan_command(
    query: Path,
    hmm_db: Path,
    tblout: Path,
    *,
    use_cut_ga: bool,
    evalue: float = DEFAULT_EVALUE,
    num_threads: int = DEFAULT_NUM_THREADS,
) -> list[str]:
    command = [
        "hmmscan",
        "--tblout",
        str(tblout),
        "--noali",
        "--cpu",
        str(num_threads),
    ]
    if use_cut_ga:
        command.append("--cut_ga")
    else:
        command.extend(["-E", str(evalue)])
    command.extend([str(hmm_db), str(query)])
    return command


def run_hmmscan(
    query: str | Path,
    hmm_db: str | Path,
    *,
    hmm_dir: str | Path | None = None,
    num_threads: int = DEFAULT_NUM_THREADS,
    evalue: float = DEFAULT_EVALUE,
) -> pd.DataFrame:
    """Search protein queries against an HMM database (Pfam-style).

    Requires a protein FASTA. Presses the HMM with hmmpress if needed.
    Uses ``--cut_ga`` when the HMM defines gathering thresholds, otherwise
    ``-E`` (default 10). Returns every tblout hit that passed that cutoff
    (several domains per protein are kept).
    """
    query_path = Path(query)
    hmm_path = Path(hmm_db)
    records = read_fasta(query_path)
    if detect_query_type(records) != "protein":
        raise HmmError("hmmscan requires a protein sequence")
    logger.info("Query has %s protein sequence(s) for hmmscan", len(records))
    ensure_pressed_hmm(hmm_path)
    use_cut_ga = hmm_has_gathering_threshold(hmm_path)
    if use_cut_ga:
        logger.info("Using hmmscan --cut_ga")
    else:
        logger.info("HMM has no GA line; using hmmscan -E %s", evalue)

    def _scan(tblout: Path) -> pd.DataFrame:
        run_hmmer(
            build_hmmscan_command(
                query_path,
                hmm_path,
                tblout,
                use_cut_ga=use_cut_ga,
                evalue=evalue,
                num_threads=num_threads,
            )
        )
        parsed = parse_hmmscan_tblout(tblout)
        logger.info("hmmscan hits passing cutoff: %s row(s)", len(parsed))
        return parsed

    if hmm_dir is not None:
        dest = Path(hmm_dir)
        dest.mkdir(parents=True, exist_ok=True)
        tblout = dest / "hmmscan.tbl"
        logger.info("Saving hmmscan tblout in %s", tblout)
        return _scan(tblout)

    with tempfile.TemporaryDirectory(prefix="biocol_hmm_") as tmp:
        return _scan(Path(tmp) / "hmmscan.tbl")
