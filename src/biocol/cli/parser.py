"""Argument parser for the biocol CLI. No backend logic here."""

from __future__ import annotations

import argparse
import sys

from biocol import (
    DEFAULT_BLAST_DIR,
    DEFAULT_DIAMOND_DIR,
    DEFAULT_HMM_DIR,
    DEFAULT_KOFAM_DIR,
    DEFAULT_MAX_TARGET_SEQS,
    DEFAULT_NUM_THREADS,
    DEFAULT_OUTPUT,
)

from biocol.cli.helptext import render_from_blast_help, render_run_help, render_top_help
from biocol.cli.style import BOLD, RED, paint


def _min_identity(value: str) -> float:
    try:
        percent = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("min-identity must be a number") from exc
    if percent < 0 or percent > 100:
        raise argparse.ArgumentTypeError("min-identity must be between 0 and 100")
    return percent


class _HelpParser(argparse.ArgumentParser):
    """argparse parser that prints the custom BIOCOL help screens."""

    def __init__(self, *args, help_renderer=None, **kwargs) -> None:
        self._help_renderer = help_renderer
        super().__init__(*args, **kwargs)

    def format_help(self) -> str:
        if self._help_renderer is not None:
            return self._help_renderer()
        return super().format_help()

    def error(self, message: str) -> None:
        sys.stderr.write(self.format_help())
        sys.stderr.write(paint(f"\nerror: {message}\n", BOLD, RED, stream=sys.stderr))
        self.exit(2)


def build_parser() -> argparse.ArgumentParser:
    parser = _HelpParser(
        prog="biocol",
        add_help=True,
        help_renderer=render_top_help,
    )
    parser.add_argument(
        "--no-color",
        action="store_true",
        help="Disable ANSI colors in help, errors, and status lines",
    )
    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
        metavar="COMMAND",
        parser_class=_HelpParser,
    )

    run = subparsers.add_parser(
        "run",
        help="FASTA + databases → BLAST+ → TSV",
        help_renderer=render_run_help,
    )
    run.add_argument(
        "--query",
        required=True,
        metavar="FASTA",
        help="Query FASTA / multifasta used for BLAST (all records same type)",
    )
    run.add_argument(
        "--cdna",
        default=None,
        metavar="FASTA",
        help="Optional query CDS FASTA: Length (nt) and cDNA columns",
    )
    run.add_argument(
        "--protein",
        default=None,
        dest="protein_fasta",
        metavar="FASTA",
        help="Optional protein FASTA (gene models): fills Length(aa) and protein columns",
    )
    run.add_argument(
        "--db",
        required=True,
        metavar="PATH",
        help="Subject FASTA file, folder of FASTA files, or with --diamond a .dmnd file/folder",
    )
    run.add_argument(
        "--accessions",
        required=True,
        metavar="FILE",
        help="TSV of SUBJECT hits: accession<TAB>descriptor (no header). Missing ids → ---",
    )
    run.add_argument(
        "--output",
        default=None,
        metavar="TSV",
        help=f"Output TSV path (default: {DEFAULT_OUTPUT})",
    )
    run.add_argument(
        "--blast-dir",
        default=None,
        dest="blast_dir",
        metavar="DIR",
        help=(
            "Directory for BLAST tabular files (one .txt per database FASTA). "
            f"Default: '{DEFAULT_BLAST_DIR}' next to the TSV"
        ),
    )
    run.add_argument(
        "--tblastx",
        action="store_true",
        help="If query and database are both nucleotide, use tblastx instead of blastn",
    )
    run.add_argument(
        "--evalue",
        type=float,
        default=10,
        metavar="N",
        help="BLAST E-value threshold (default: 10)",
    )
    run.add_argument(
        "--max-target-seqs",
        type=int,
        default=DEFAULT_MAX_TARGET_SEQS,
        dest="max_target_seqs",
        metavar="N",
        help=f"Maximum aligned sequences to keep per query (default: {DEFAULT_MAX_TARGET_SEQS})",
    )
    run.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_NUM_THREADS,
        metavar="N",
        help=f"BLAST+, Diamond, hmmscan, or KofamScan CPU threads (default: {DEFAULT_NUM_THREADS})",
    )
    run.add_argument(
        "--min-identity",
        type=_min_identity,
        default=None,
        dest="min_identity",
        metavar="N",
        help="Keep HSPs with pident >= N (0-100). Omit for no identity cutoff",
    )
    run.add_argument(
        "--reciprocal",
        action="store_true",
        dest="reciprocal",
        help=(
            "Also BLAST each database FASTA against the query; keep only pairs "
            "that are the top hit in both directions"
        ),
    )
    run.add_argument(
        "--diamond",
        action="store_true",
        dest="diamond",
        help=(
            "Use Diamond blastp instead of BLASTP (protein query and protein "
            "FASTA or .dmnd). evalue and max-target-seqs use Diamond defaults"
        ),
    )
    run.add_argument(
        "--diamond-dir",
        default=None,
        dest="diamond_dir",
        metavar="DIR",
        help=(
            "Directory for Diamond tabular files (only with --diamond). "
            f"Default: '{DEFAULT_DIAMOND_DIR}' next to the TSV"
        ),
    )
    run.add_argument(
        "--hmm-db",
        default=None,
        dest="hmm_db",
        metavar="HMM",
        help=(
            "Optional HMMER HMM database (e.g. Pfam-A.hmm). Runs hmmscan on a "
            "protein query; hmmpress if indexes are missing"
        ),
    )
    run.add_argument(
        "--hmm-dir",
        default=None,
        dest="hmm_dir",
        metavar="DIR",
        help=(
            "Directory for hmmscan tblout (hmmscan.tbl). "
            f"Default: '{DEFAULT_HMM_DIR}' next to the TSV (only if --hmm-db is set)"
        ),
    )
    run.add_argument(
        "--kofam-profile",
        default=None,
        dest="kofam_profile",
        metavar="PATH",
        help=(
            "Optional KOfam profiles: directory of .hmm files, a .hmm file, or a "
            ".hal list (e.g. eukaryote.hal). Requires --ko-list; protein query"
        ),
    )
    run.add_argument(
        "--ko-list",
        default=None,
        dest="ko_list",
        metavar="FILE",
        help="KOfam ko_list file (required with --kofam-profile)",
    )
    run.add_argument(
        "--kofam-dir",
        default=None,
        dest="kofam_dir",
        metavar="DIR",
        help=(
            "Directory for exec_annotation detail-tsv (kofam.tsv) and tmp/. "
            f"Default: '{DEFAULT_KOFAM_DIR}' next to the TSV "
            "(only if --kofam-profile and --ko-list are set)"
        ),
    )

    from_blast = subparsers.add_parser(
        "from-blast",
        help="Existing BLAST tabular + accessions → TSV (no BLAST+)",
        help_renderer=render_from_blast_help,
    )
    from_blast.add_argument(
        "--blast",
        required=True,
        metavar="FILE",
        help="BLAST tabular file (outfmt 6, typically .txt)",
    )
    from_blast.add_argument(
        "--accessions",
        required=True,
        metavar="FILE",
        help="TSV of SUBJECT hits: accession<TAB>descriptor (no header)",
    )
    from_blast.add_argument(
        "--output",
        default=None,
        metavar="TSV",
        help=f"Output TSV path (default: {DEFAULT_OUTPUT})",
    )
    from_blast.add_argument(
        "--min-identity",
        type=_min_identity,
        default=None,
        dest="min_identity",
        metavar="N",
        help="Keep HSPs with pident >= N (0-100). Omit for no identity cutoff",
    )
    from_blast.add_argument(
        "--hmm-db",
        default=None,
        dest="hmm_db",
        metavar="HMM",
        help="Optional HMMER HMM database; requires --protein (amino-acid FASTA)",
    )
    from_blast.add_argument(
        "--protein",
        default=None,
        dest="protein_fasta",
        metavar="FASTA",
        help="Protein FASTA for hmmscan and/or KofamScan",
    )
    from_blast.add_argument(
        "--hmm-dir",
        default=None,
        dest="hmm_dir",
        metavar="DIR",
        help=(
            f"Directory for hmmscan tblout. Default: '{DEFAULT_HMM_DIR}' next to the TSV"
        ),
    )
    from_blast.add_argument(
        "--kofam-profile",
        default=None,
        dest="kofam_profile",
        metavar="PATH",
        help=(
            "Optional KOfam profiles (directory, .hmm, or .hal). "
            "Requires --ko-list and --protein"
        ),
    )
    from_blast.add_argument(
        "--ko-list",
        default=None,
        dest="ko_list",
        metavar="FILE",
        help="KOfam ko_list file (required with --kofam-profile)",
    )
    from_blast.add_argument(
        "--kofam-dir",
        default=None,
        dest="kofam_dir",
        metavar="DIR",
        help=(
            "Directory for exec_annotation detail-tsv. "
            f"Default: '{DEFAULT_KOFAM_DIR}' next to the TSV"
        ),
    )
    return parser
