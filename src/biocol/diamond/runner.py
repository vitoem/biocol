"""Diamond blastp (protein vs protein) as an optional BLASTP substitute."""

from __future__ import annotations

import logging
import re
import tempfile
from pathlib import Path

import pandas as pd

from biocol.blast.parser import BLAST_OUTFMT, fill_missing_hits, parse_blast_results
from biocol.blast.reciprocal import keep_reciprocal_best_hits
from biocol.diamond.databases import list_diamond_databases
from biocol.diamond.execute import run_diamond_command
from biocol.exceptions import DiamondError
from biocol.processing.hsp_filter import filter_hits_by_pident
from biocol.sequence.classifier import detect_query_type
from biocol.sequence.reader import read_fasta

logger = logging.getLogger(__name__)

DEFAULT_DIAMOND_DIR = "diamond"
_OUTFMT_FIELDS = BLAST_OUTFMT.split()


def _tabular_filename(stem: str, used: set[str]) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("._") or "db"
    name = f"{cleaned}.txt"
    index = 2
    while name in used:
        name = f"{cleaned}_{index}.txt"
        index += 1
    used.add(name)
    return name


def build_makedb_command(fasta: Path, db_prefix: Path) -> list[str]:
    return [
        "diamond",
        "makedb",
        "--in",
        str(fasta),
        "--db",
        str(db_prefix),
    ]


def build_blastp_command(
    query: Path,
    db_prefix: Path,
    out_file: Path,
    *,
    num_threads: int,
) -> list[str]:
    return [
        "diamond",
        "blastp",
        "--query",
        str(query),
        "--db",
        str(db_prefix),
        "--out",
        str(out_file),
        "--outfmt",
        *_OUTFMT_FIELDS,
        "--threads",
        str(num_threads),
    ]


def _out_path(
    stem: str,
    hits_dir: Path | None,
    tmp_path: Path,
    used_names: set[str],
) -> Path:
    if hits_dir is not None:
        return hits_dir / _tabular_filename(stem, used_names)
    return tmp_path / f"{stem}.txt"


def run_diamond(
    query: str | Path,
    database: str | Path,
    *,
    num_threads: int = 1,
    diamond_dir: str | Path | None = None,
    min_identity: float | None = None,
    reciprocal: bool = False,
) -> pd.DataFrame:
    """Run Diamond blastp (one run per protein FASTA or ``.dmnd``).

    ``evalue`` and ``max-target-seqs`` use Diamond defaults. Threads default
    to 1. Indexes built from FASTA are temporary. If ``diamond_dir`` is set,
    tabular files are kept there. Reciprocal search requires FASTA databases
    (not ``.dmnd``).
    """
    query_path = Path(query)
    logger.info("Reading query FASTA: %s", query_path)
    records = read_fasta(query_path)
    query_ids = [record.id for record in records]
    logger.info("Query has %s sequence(s)", len(query_ids))
    query_type = detect_query_type(records)
    if query_type != "protein":
        raise DiamondError(
            "diamond requires a protein query and a protein database"
        )

    db_entries = list_diamond_databases(database)
    kinds = {kind for _, kind, _ in db_entries}
    if "fasta" in kinds and "dmnd" in kinds:
        raise DiamondError("Input mixes FASTA and Diamond (.dmnd) databases")
    if reciprocal and "dmnd" in kinds:
        raise DiamondError(
            "Reciprocal Diamond search requires protein FASTA databases, not .dmnd"
        )

    hits_dir: Path | None = None
    used_names: set[str] = set()
    if diamond_dir is not None:
        hits_dir = Path(diamond_dir)
        hits_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Saving Diamond tabular files in %s", hits_dir)

    logger.info(
        "Selected diamond blastp (query=protein, database=protein, threads=%s)",
        num_threads,
    )
    forward = _run_diamond_blastp(
        query_path,
        query_ids,
        db_entries,
        num_threads=num_threads,
        hits_dir=hits_dir,
        used_names=used_names,
        min_identity=min_identity,
    )
    if not reciprocal:
        logger.info("Diamond finished (%s row(s) including empty hits)", len(forward))
        return forward

    reverse = _run_reverse_diamond(
        query_path,
        query_ids,
        db_entries,
        num_threads=num_threads,
        hits_dir=hits_dir,
        used_names=used_names,
        min_identity=min_identity,
    )
    combined = keep_reciprocal_best_hits(forward, reverse)
    logger.info("Diamond finished (%s row(s) including empty hits)", len(combined))
    return combined


def _run_diamond_blastp(
    query_path: Path,
    query_ids: list[str],
    db_entries: list[tuple[Path, str, str]],
    *,
    num_threads: int,
    hits_dir: Path | None,
    used_names: set[str],
    min_identity: float | None,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    with tempfile.TemporaryDirectory(prefix="biocol_diamond_") as tmp:
        tmp_path = Path(tmp)
        total = len(db_entries)
        for index, (db_path, kind, db_label) in enumerate(db_entries, start=1):
            file_stem = db_path.stem
            out_file = _out_path(file_stem, hits_dir, tmp_path, used_names)
            if kind == "fasta":
                prefix = tmp_path / f"db_{index}"
                logger.info(
                    "[%s/%s] Building Diamond database for %s (%s)",
                    index,
                    total,
                    db_label,
                    db_path.name,
                )
                run_diamond_command(build_makedb_command(db_path, prefix))
                db_arg = prefix
            else:
                db_arg = db_path
                logger.info("[%s/%s] Using existing Diamond database %s", index, total, db_path.name)
            logger.info("[%s/%s] Running diamond blastp against %s", index, total, db_label)
            run_diamond_command(
                build_blastp_command(
                    query_path,
                    db_arg,
                    out_file,
                    num_threads=num_threads,
                )
            )
            parsed = parse_blast_results(out_file)
            parsed["database"] = db_label
            filled = fill_missing_hits(parsed, query_ids, db_label)
            frames.append(filled)
            logger.info(
                "[%s/%s] %s: %s Diamond hit(s)%s",
                index,
                total,
                db_label,
                0 if parsed.empty else len(parsed),
                f" → {out_file}" if hits_dir is not None else "",
            )
    combined = pd.concat(frames, ignore_index=True)
    return filter_hits_by_pident(combined, min_identity, query_ids=query_ids)


def _run_reverse_diamond(
    query_path: Path,
    query_ids: list[str],
    db_entries: list[tuple[Path, str, str]],
    *,
    num_threads: int,
    hits_dir: Path | None,
    used_names: set[str],
    min_identity: float | None,
) -> pd.DataFrame:
    logger.info("Reciprocal Diamond: blastp with database FASTA as query")
    frames: list[pd.DataFrame] = []
    with tempfile.TemporaryDirectory(prefix="biocol_diamond_rev_") as tmp:
        tmp_path = Path(tmp)
        query_prefix = tmp_path / "query_db"
        logger.info("Building reverse Diamond database from query FASTA")
        run_diamond_command(build_makedb_command(query_path, query_prefix))
        total = len(db_entries)
        for index, (fasta_path, kind, db_label) in enumerate(db_entries, start=1):
            if kind != "fasta":
                raise DiamondError(
                    "Reciprocal Diamond search requires protein FASTA databases, not .dmnd"
                )
            species_ids = [record.id for record in read_fasta(fasta_path)]
            out_file = _out_path(
                f"reverse_{fasta_path.stem}",
                hits_dir,
                tmp_path,
                used_names,
            )
            logger.info(
                "[%s/%s] Reverse diamond blastp: %s as query against original FASTA",
                index,
                total,
                db_label,
            )
            run_diamond_command(
                build_blastp_command(
                    fasta_path,
                    query_prefix,
                    out_file,
                    num_threads=num_threads,
                )
            )
            parsed = parse_blast_results(out_file)
            parsed["database"] = db_label
            filled = fill_missing_hits(parsed, species_ids, db_label)
            filled = filter_hits_by_pident(filled, min_identity, query_ids=species_ids)
            frames.append(filled)
            logger.info(
                "[%s/%s] reverse %s: %s Diamond hit(s)%s",
                index,
                total,
                db_label,
                0 if parsed.empty else len(parsed),
                f" → {out_file}" if hits_dir is not None else "",
            )
    return pd.concat(frames, ignore_index=True)
