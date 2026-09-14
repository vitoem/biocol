"""CLI commands. Calls only the public biocol API."""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path

from biocol import (
    DEFAULT_BLAST_DIR,
    DEFAULT_DIAMOND_DIR,
    DEFAULT_HMM_DIR,
    DEFAULT_KOFAM_DIR,
    DEFAULT_OUTPUT,
    HmmError,
    KofamError,
    build_result_table,
    filter_hits_by_pident,
    parse_blast_results,
    run_blast,
    run_hmmscan,
    run_kofamscan,
    write_results_csv,
)


def _dir_next_to_output(output: str | None, explicit: str | None, default_name: str) -> Path:
    if explicit:
        return Path(explicit)
    parent = Path(output).parent if output else Path(DEFAULT_OUTPUT).parent
    return parent / default_name


def _maybe_hmmscan(args: Namespace, protein_fasta: str | None):
    if not args.hmm_db:
        return None
    if not protein_fasta:
        raise HmmError("hmmscan requires a protein sequence")
    return run_hmmscan(
        protein_fasta,
        args.hmm_db,
        hmm_dir=_dir_next_to_output(args.output, args.hmm_dir, DEFAULT_HMM_DIR),
        num_threads=getattr(args, "threads", 1),
    )


def _maybe_kofamscan(args: Namespace, protein_fasta: str | None):
    profile = args.kofam_profile
    ko_list = args.ko_list
    if not profile and not ko_list:
        return None
    if not profile or not ko_list:
        raise KofamError("kofamscan requires --kofam-profile and --ko-list")
    if not protein_fasta:
        raise KofamError("kofamscan requires a protein sequence")
    return run_kofamscan(
        protein_fasta,
        profile,
        ko_list,
        kofam_dir=_dir_next_to_output(args.output, args.kofam_dir, DEFAULT_KOFAM_DIR),
        num_threads=getattr(args, "threads", 1),
    )


def run_from_fasta(args: Namespace) -> Path:
    use_diamond = args.diamond
    hits = run_blast(
        args.query,
        args.db,
        translated=args.tblastx,
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        num_threads=args.threads,
        blast_dir=(
            None
            if use_diamond
            else _dir_next_to_output(args.output, args.blast_dir, DEFAULT_BLAST_DIR)
        ),
        min_identity=args.min_identity,
        reciprocal=args.reciprocal,
        diamond=use_diamond,
        diamond_dir=(
            _dir_next_to_output(args.output, args.diamond_dir, DEFAULT_DIAMOND_DIR)
            if use_diamond
            else None
        ),
    )
    hmm_hits = _maybe_hmmscan(args, args.query)
    kofam_hits = _maybe_kofamscan(args, args.query)
    table = build_result_table(
        hits,
        args.accessions,
        query_fasta=args.query,
        cdna_fasta=args.cdna,
        protein_fasta=args.protein_fasta,
        hmm_hits=hmm_hits,
        kofam_hits=kofam_hits,
    )
    return write_results_csv(table, args.output)


def run_from_blast(args: Namespace) -> Path:
    hits = parse_blast_results(args.blast)
    if "database" not in hits.columns:
        hits["database"] = Path(args.accessions).stem
    hits = filter_hits_by_pident(hits, args.min_identity)
    hmm_hits = _maybe_hmmscan(args, args.protein_fasta)
    kofam_hits = _maybe_kofamscan(args, args.protein_fasta)
    table = build_result_table(
        hits,
        args.accessions,
        query_fasta=None,
        protein_fasta=args.protein_fasta,
        hmm_hits=hmm_hits,
        kofam_hits=kofam_hits,
    )
    return write_results_csv(table, args.output)
