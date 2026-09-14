"""Run KofamScan ``exec_annotation`` on a protein FASTA."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from biocol.exceptions import KofamError
from biocol.kofam.execute import run_exec_annotation
from biocol.kofam.parser import parse_kofam_detail_tsv
from biocol.sequence.classifier import detect_query_type
from biocol.sequence.reader import read_fasta

logger = logging.getLogger(__name__)

DEFAULT_KOFAM_DIR = "kofam"
DEFAULT_NUM_THREADS = 1
OUTPUT_NAME = "kofam.tsv"
TMP_NAME = "tmp"
_PROFILE_FILE_SUFFIXES = {".hal", ".hmm"}


def build_exec_annotation_command(
    query: Path,
    profile: Path,
    ko_list: Path,
    out_file: Path,
    tmp_dir: Path,
    *,
    num_threads: int = DEFAULT_NUM_THREADS,
) -> list[str]:
    return [
        "exec_annotation",
        "-o",
        str(out_file),
        "-f",
        "detail-tsv",
        "--profile",
        str(profile),
        "--ko-list",
        str(ko_list),
        "--cpu",
        str(num_threads),
        "--tmp-dir",
        str(tmp_dir),
        str(query),
    ]


def _resolve_profile(profile: Path) -> Path:
    if not profile.exists():
        raise KofamError(f"KOfam profile not found: {profile}")
    if profile.is_dir():
        return profile
    if profile.is_file() and profile.suffix.lower() in _PROFILE_FILE_SUFFIXES:
        return profile
    raise KofamError(
        "KOfam profile must be a directory of .hmm files, a .hmm file, or a .hal file"
    )


def _resolve_ko_list(ko_list: Path) -> Path:
    if not ko_list.is_file():
        raise KofamError(f"KO list file not found: {ko_list}")
    return ko_list


def run_kofamscan(
    query: str | Path,
    profile: str | Path,
    ko_list: str | Path,
    *,
    kofam_dir: str | Path | None = None,
    num_threads: int = DEFAULT_NUM_THREADS,
) -> pd.DataFrame:
    """Annotate protein queries with KofamScan (KOfam + ko_list).

    Requires a protein FASTA. ``profile`` is a directory of HMM profiles, a
    ``.hmm`` file, or a ``.hal`` list. Output is ``detail-tsv``; only hits
    above the KO threshold (``*``) are returned. The annotation TSV is saved
    in ``kofam_dir`` (default ``kofam/``) and hmmsearch intermediates in
    ``kofam_dir/tmp`` (kept).
    """
    query_path = Path(query)
    records = read_fasta(query_path)
    if detect_query_type(records) != "protein":
        raise KofamError("kofamscan requires a protein sequence")
    logger.info("Query has %s protein sequence(s) for KofamScan", len(records))

    profile_path = _resolve_profile(Path(profile))
    ko_list_path = _resolve_ko_list(Path(ko_list))
    dest = Path(kofam_dir) if kofam_dir is not None else Path(DEFAULT_KOFAM_DIR)
    dest.mkdir(parents=True, exist_ok=True)
    tmp_dir = dest / TMP_NAME
    tmp_dir.mkdir(parents=True, exist_ok=True)
    out_file = dest / OUTPUT_NAME
    logger.info("Saving KofamScan detail-tsv in %s", out_file)
    logger.info("KofamScan tmp-dir (kept): %s", tmp_dir)

    run_exec_annotation(
        build_exec_annotation_command(
            query_path,
            profile_path,
            ko_list_path,
            out_file,
            tmp_dir,
            num_threads=num_threads,
        )
    )
    return parse_kofam_detail_tsv(out_file)
