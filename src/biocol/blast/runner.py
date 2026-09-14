from __future__ import annotations

import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from biocol.blast.databases import infer_database_label, list_blast_databases
from biocol.blast.parser import BLAST_OUTFMT, fill_missing_hits, parse_blast_results
from biocol.blast.reciprocal import keep_reciprocal_best_hits
from biocol.blast.selection import select_blast_program
from biocol.exceptions import BlastExecutionError, MixedDatabaseTypeError
from biocol.processing.hsp_filter import filter_hits_by_pident
from biocol.sequence.classifier import detect_query_type
from biocol.sequence.reader import read_fasta

logger = logging.getLogger(__name__)

_DBTYPE = {"nucleotide": "nucl", "protein": "prot"}
DEFAULT_MAX_TARGET_SEQS = 3
DEFAULT_NUM_THREADS = 1
DEFAULT_BLAST_DIR = "blast"


def _tabular_filename(stem: str, used: set[str]) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("._") or "db"
    name = f"{cleaned}.txt"
    index = 2
    while name in used:
        name = f"{cleaned}_{index}.txt"
        index += 1
    used.add(name)
    return name


def build_makeblastdb_command(
    fasta: Path,
    out_prefix: Path,
    dbtype: str,
) -> list[str]:
    return [
        "makeblastdb",
        "-in",
        str(fasta),
        "-dbtype",
        dbtype,
        "-out",
        str(out_prefix),
        "-parse_seqids",
    ]


def build_blast_command(
    program: str,
    query: Path,
    db_prefix: Path,
    out_file: Path,
    *,
    evalue: float = 10,
    max_target_seqs: int = DEFAULT_MAX_TARGET_SEQS,
    num_threads: int = DEFAULT_NUM_THREADS,
) -> list[str]:
    return [
        program,
        "-query",
        str(query),
        "-db",
        str(db_prefix),
        "-out",
        str(out_file),
        "-outfmt",
        BLAST_OUTFMT,
        "-evalue",
        str(evalue),
        "-max_target_seqs",
        str(max_target_seqs),
        "-num_threads",
        str(num_threads),
    ]


def _run_command(command: list[str]) -> None:
    executable = command[0]
    if shutil.which(executable) is None:
        raise BlastExecutionError(
            f"'{executable}' was not found in PATH. Activate the conda environment inecol."
        )
    logger.debug("BLAST+ command: %s", " ".join(command))
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise BlastExecutionError(
            f"{executable} failed (exit code {completed.returncode}): {completed.stderr.strip()}"
        )


def _run_forward_blast(
    query_path: Path,
    query_ids: list[str],
    program: str,
    db_entries: list[tuple[Path, str]],
    dbtype: str,
    *,
    evalue: float,
    max_target_seqs: int,
    num_threads: int,
    hits_dir: Path | None,
    used_names: set[str],
    min_identity: float | None,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    with tempfile.TemporaryDirectory(prefix="biocol_blast_") as tmp:
        tmp_path = Path(tmp)
        total = len(db_entries)
        for index, (fasta_path, _db_type) in enumerate(db_entries, start=1):
            file_stem = fasta_path.stem
            db_label = infer_database_label(fasta_path)
            prefix = tmp_path / file_stem
            if hits_dir is not None:
                out_file = hits_dir / _tabular_filename(file_stem, used_names)
            else:
                out_file = tmp_path / f"{file_stem}.txt"
            logger.info(
                "[%s/%s] Building BLAST database for %s (%s)",
                index,
                total,
                db_label,
                fasta_path.name,
            )
            _run_command(build_makeblastdb_command(fasta_path, prefix, dbtype))
            logger.info("[%s/%s] Running %s against %s", index, total, program, db_label)
            _run_command(
                build_blast_command(
                    program,
                    query_path,
                    prefix,
                    out_file,
                    evalue=evalue,
                    max_target_seqs=max_target_seqs,
                    num_threads=num_threads,
                )
            )
            parsed = parse_blast_results(out_file)
            parsed["database"] = db_label
            filled = fill_missing_hits(parsed, query_ids, db_label)
            frames.append(filled)
            logger.info(
                "[%s/%s] %s: %s BLAST hit(s)%s",
                index,
                total,
                db_label,
                0 if parsed.empty else len(parsed),
                f" → {out_file}" if hits_dir is not None else "",
            )

    combined = pd.concat(frames, ignore_index=True)
    return filter_hits_by_pident(combined, min_identity, query_ids=query_ids)


def _run_reverse_blast(
    query_path: Path,
    query_ids: list[str],
    query_type: str,
    db_entries: list[tuple[Path, str]],
    database_type: str,
    *,
    translated: bool,
    evalue: float,
    num_threads: int,
    hits_dir: Path | None,
    used_names: set[str],
    min_identity: float | None,
) -> pd.DataFrame:
    program = select_blast_program(
        database_type, query_type, translated=translated
    )
    max_target_seqs = max(1, len(query_ids))
    query_dbtype = _DBTYPE[query_type]
    logger.info(
        "Reciprocal BLAST: selected %s (query=%s, database=%s, max_target_seqs=%s)",
        program,
        database_type,
        query_type,
        max_target_seqs,
    )
    frames: list[pd.DataFrame] = []
    with tempfile.TemporaryDirectory(prefix="biocol_blast_rev_") as tmp:
        tmp_path = Path(tmp)
        query_prefix = tmp_path / "query_db"
        logger.info("Building reverse BLAST database from query FASTA")
        _run_command(build_makeblastdb_command(query_path, query_prefix, query_dbtype))
        total = len(db_entries)
        for index, (fasta_path, _db_type) in enumerate(db_entries, start=1):
            file_stem = fasta_path.stem
            db_label = infer_database_label(fasta_path)
            species_ids = [record.id for record in read_fasta(fasta_path)]
            tabular_stem = f"reverse_{file_stem}"
            if hits_dir is not None:
                out_file = hits_dir / _tabular_filename(tabular_stem, used_names)
            else:
                out_file = tmp_path / f"{tabular_stem}.txt"
            logger.info(
                "[%s/%s] Reverse %s: %s as query against original FASTA",
                index,
                total,
                program,
                db_label,
            )
            _run_command(
                build_blast_command(
                    program,
                    fasta_path,
                    query_prefix,
                    out_file,
                    evalue=evalue,
                    max_target_seqs=max_target_seqs,
                    num_threads=num_threads,
                )
            )
            parsed = parse_blast_results(out_file)
            parsed["database"] = db_label
            filled = fill_missing_hits(parsed, species_ids, db_label)
            filled = filter_hits_by_pident(filled, min_identity, query_ids=species_ids)
            frames.append(filled)
            logger.info(
                "[%s/%s] reverse %s: %s BLAST hit(s)%s",
                index,
                total,
                db_label,
                0 if parsed.empty else len(parsed),
                f" → {out_file}" if hits_dir is not None else "",
            )
    return pd.concat(frames, ignore_index=True)


def run_blast(
    query: str | Path,
    database: str | Path,
    *,
    translated: bool = False,
    evalue: float = 10,
    max_target_seqs: int = DEFAULT_MAX_TARGET_SEQS,
    num_threads: int = DEFAULT_NUM_THREADS,
    blast_dir: str | Path | None = None,
    min_identity: float | None = None,
    reciprocal: bool = False,
    diamond: bool = False,
    diamond_dir: str | Path | None = None,
) -> pd.DataFrame:
    """Run BLAST+ (one run per database FASTA) and parse outfmt 6.

    BLAST indexes are built in a temporary directory and removed afterwards.
    If ``blast_dir`` is set, each BLAST tabular file is kept there (one
    ``.txt`` per database FASTA, unfiltered). Otherwise the tabular files
    are discarded. ``min_identity`` drops HSPs whose ``pident`` is below
    that cutoff after parsing. If a query has no hit in a database, an
    empty row is included.

    If ``reciprocal`` is True, each database FASTA is also searched against
    the original query (program types swapped). Only pairs that are the
    top hit in both directions are kept.

    If ``diamond`` is True, run Diamond blastp instead of BLAST+ (protein
    query and protein FASTA or ``.dmnd`` only). Tabular files go to
    ``diamond_dir``. ``evalue`` / ``max_target_seqs`` / ``translated`` /
    ``blast_dir`` are ignored; Diamond defaults apply except threads.
    """
    if diamond:
        from biocol.diamond.runner import run_diamond

        return run_diamond(
            query,
            database,
            num_threads=num_threads,
            diamond_dir=diamond_dir,
            min_identity=min_identity,
            reciprocal=reciprocal,
        )
    query_path = Path(query)
    logger.info("Reading query FASTA: %s", query_path)
    records = read_fasta(query_path)
    query_ids = [record.id for record in records]
    logger.info("Query has %s sequence(s)", len(query_ids))
    query_type = detect_query_type(records)
    db_entries = list_blast_databases(database)
    db_types = {db_type for _, db_type in db_entries}
    if len(db_types) > 1:
        raise MixedDatabaseTypeError(
            "Input mixes nucleotide and protein databases"
        )
    database_type = next(iter(db_types))
    logger.info(
        "Database type: %s (%s FASTA file(s))",
        database_type,
        len(db_entries),
    )
    program = select_blast_program(
        query_type, database_type, translated=translated
    )
    logger.info(
        "Selected %s (query=%s, database=%s, evalue=%s, max_target_seqs=%s, threads=%s)",
        program,
        query_type,
        database_type,
        evalue,
        max_target_seqs,
        num_threads,
    )

    hits_dir: Path | None = None
    used_names: set[str] = set()
    if blast_dir is not None:
        hits_dir = Path(blast_dir)
        hits_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Saving BLAST tabular files in %s", hits_dir)

    forward = _run_forward_blast(
        query_path,
        query_ids,
        program,
        db_entries,
        _DBTYPE[database_type],
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        num_threads=num_threads,
        hits_dir=hits_dir,
        used_names=used_names,
        min_identity=min_identity,
    )
    if not reciprocal:
        logger.info("BLAST finished (%s row(s) including empty hits)", len(forward))
        return forward

    reverse = _run_reverse_blast(
        query_path,
        query_ids,
        query_type,
        db_entries,
        database_type,
        translated=translated,
        evalue=evalue,
        num_threads=num_threads,
        hits_dir=hits_dir,
        used_names=used_names,
        min_identity=min_identity,
    )
    combined = keep_reciprocal_best_hits(forward, reverse)
    logger.info("BLAST finished (%s row(s) including empty hits)", len(combined))
    return combined
