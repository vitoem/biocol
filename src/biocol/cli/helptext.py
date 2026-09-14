"""Colored help screens for the biocol CLI."""

from __future__ import annotations

import os
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

from biocol.cli.style import BOLD, CYAN, DIM, GREEN, MAGENTA, YELLOW, paint

WIDTH = 76


def _box_ok() -> bool:
    if os.environ.get("BIOCOL_ASCII"):
        return False
    encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
    try:
        "─│╭╮╰╯".encode(encoding)
    except Exception:
        return False
    return True


def _center(text: str, width: int) -> str:
    pad = max(width - len(text), 0)
    left = pad // 2
    return (" " * left) + text + (" " * (pad - left))


def banner(subtitle: str, stream=None) -> str:
    inner = WIDTH - 2
    if _box_ok():
        top, bot, side, bar = "╭", "╰", "│", "─"
        top_r, bot_r = "╮", "╯"
    else:
        top = bot = side = "+"
        bar = "-"
        top_r = bot_r = "+"
    line1 = _center("BIOCOL", inner)
    line2 = _center(subtitle, inner)
    return "\n".join(
        [
            paint(f"{top}{bar * inner}{top_r}", CYAN, BOLD, stream=stream),
            paint(side, CYAN, BOLD, stream=stream)
            + paint(line1, BOLD, stream=stream)
            + paint(side, CYAN, BOLD, stream=stream),
            paint(side, CYAN, BOLD, stream=stream)
            + paint(line2, DIM, stream=stream)
            + paint(side, CYAN, BOLD, stream=stream),
            paint(f"{bot}{bar * inner}{bot_r}", CYAN, BOLD, stream=stream),
        ]
    )


def heading(title: str, stream=None) -> str:
    return paint(title, BOLD, CYAN, stream=stream)


def cmd(name: str, stream=None) -> str:
    return paint(name, BOLD, MAGENTA, stream=stream)


def opt(name: str, stream=None) -> str:
    return paint(name, YELLOW, stream=stream)


def example(block: str, stream=None) -> str:
    return paint(block, GREEN, stream=stream)


def dim(text: str, stream=None) -> str:
    return paint(text, DIM, stream=stream)


def render_top_help() -> str:
    s = sys.stdout
    parts = [
        banner("BLAST Annotation & Integration Tool", stream=s),
        "",
        "BIOCOL processes biological sequence data and generates integrated",
        "BLAST annotation tables in Dataset S2-style TSV format.",
        "",
        "It supports two workflows:",
        "",
        "  1. Run BLAST+ locally from FASTA sequences.",
        "  2. Build the result table from an existing BLAST outfmt 6 file.",
        "",
        f"{heading('USAGE', stream=s)}",
        "",
        f"  {cmd('biocol', stream=s)} [{opt('--no-color', stream=s)}] {cmd('COMMAND', stream=s)} ...",
        "",
        f"{heading('WORKFLOWS', stream=s)}",
        "",
        f"  {cmd('run', stream=s)}",
        "    FASTA query → BLAST+ → annotations → TSV",
        "",
        "    Runs BLAST+ locally against one FASTA file or a directory of FASTA",
        "    databases. Sequence types are detected automatically and the matching",
        "    BLAST program is selected (blastp, blastn, blastx, tblastn, tblastx).",
        "",
        "    Optional --cdna / --protein fill gene-model columns.",
        "",
        f"  {cmd('from-blast', stream=s)}",
        "    BLAST outfmt 6 → annotations → TSV",
        "",
        "    Parses an existing BLAST tabular file and joins accession descriptors",
        "    without running BLAST+ again. Query sequence columns stay empty.",
        "",
        f"{heading('QUICK START', stream=s)}",
        "",
        "  Complete workflow:",
        "",
        example(
            "    biocol run \\\n"
            "      --query query.faa \\\n"
            "      --db species.faa \\\n"
            "      --accessions species.txt \\\n"
            "      --output results.tsv",
            stream=s,
        ),
        "",
        "  Include query cDNA columns:",
        "",
        example(
            "    biocol run \\\n"
            "      --query query.faa \\\n"
            "      --cdna query_cds.fna \\\n"
            "      --db species.faa \\\n"
            "      --accessions species.txt \\\n"
            "      --output results.tsv",
            stream=s,
        ),
        "",
        "  Use existing BLAST results:",
        "",
        example(
            "    biocol from-blast \\\n"
            "      --blast hits.txt \\\n"
            "      --accessions species.txt \\\n"
            "      --output results.tsv",
            stream=s,
        ),
        "",
        f"{heading('LEARN MORE', stream=s)}",
        "",
        f"  {cmd('biocol run --help', stream=s)}",
        "      Inputs, BLAST options, program selection and examples.",
        "",
        f"  {cmd('biocol from-blast --help', stream=s)}",
        "      Options for processing an existing BLAST tabular file.",
        "",
        f"{heading('GLOBAL OPTIONS', stream=s)}",
        "",
        f"  {opt('-h, --help', stream=s)}",
        "      Show this help message and exit.",
        "",
        f"  {opt('--no-color', stream=s)}",
        "      Disable ANSI colors. Also honored via NO_COLOR=1.",
        "",
        f"{heading('NOTES', stream=s)}",
        "",
        "  --accessions is a .txt file: accession<TAB>descriptor.",
        "  It must list SUBJECT ids from --db. Unmatched hits → ---.",
        "  --cdna must be CDS/transcripts, not a whole-genome FASTA.",
        "  The TSV shows only the best hit per query and species.",
        "  Optional --reciprocal keeps that hit only if it is also the top",
        "  hit in the reverse BLAST (not available in from-blast).",
        f"  Defaults: e-value 10, max-target-seqs {DEFAULT_MAX_TARGET_SEQS},",
        f"  threads {DEFAULT_NUM_THREADS}, output {DEFAULT_OUTPUT}.",
        "",
        dim("BIOCOL — bioinformatics sequence annotation from the command line.", stream=s),
        "",
    ]
    return "\n".join(parts)


def render_run_help() -> str:
    s = sys.stdout
    parts = [
        banner("run — FASTA → BLAST+ → TSV", stream=s),
        "",
        "Run BLAST+ on a FASTA query against one FASTA file or a folder of FASTA",
        "files (subfolders included). Join descriptors and write a TSV file.",
        "",
        f"{heading('USAGE', stream=s)}",
        "",
        f"  {cmd('biocol run', stream=s)} {opt('--query', stream=s)} FASTA {opt('--db', stream=s)} PATH "
        f"{opt('--accessions', stream=s)} FILE [options]",
        "",
        f"{heading('REQUIRED', stream=s)}",
        "",
        f"  {opt('--query FASTA', stream=s)}",
        "      Query FASTA / multifasta. All records must be the same type.",
        "      A protein query also fills Length(aa) and protein sequence columns.",
        "",
        f"  {opt('--db PATH', stream=s)}",
        "      Subject FASTA, or a directory of FASTA files (same molecule type).",
        "      With --diamond: protein FASTA/.dmnd, or a folder of one type only.",
        "",
        f"  {opt('--accessions FILE', stream=s)}",
        "      TSV of SUBJECT hits: accession<TAB>descriptor (no header).",
        "",
        f"{heading('OPTIONAL GENE MODELS', stream=s)}  {dim('(query organism)', stream=s)}",
        "",
        f"  {opt('--cdna FASTA', stream=s)}",
        "      Fill Length (nt) and cDNA columns. Must be CDS, not genomic DNA.",
        "",
        f"  {opt('--protein FASTA', stream=s)}",
        "      Fill Length(aa) and protein columns (gene models).",
        "",
        f"{heading('BLAST OPTIONS', stream=s)}",
        "",
        f"  {opt('--tblastx', stream=s)}",
        "      If query and database are both nucleotide transcripts.",
        "",
        f"  {opt('--evalue N', stream=s)}",
        "      BLAST E-value threshold (default: 10). Ignored with --diamond.",
        "",
        f"  {opt('--max-target-seqs N', stream=s)}",
        f"      Hits kept by BLAST per query (default: {DEFAULT_MAX_TARGET_SEQS}).",
        "      The TSV still shows only the best hit per species. Ignored with --diamond.",
        "",
        f"  {opt('--threads N', stream=s)}",
        f"      BLAST+, Diamond, hmmscan, or KofamScan CPU threads "
        f"(default: {DEFAULT_NUM_THREADS}).",
        "",
        f"  {opt('--min-identity N', stream=s)}",
        "      Keep HSPs with BLAST pident >= N (0-100, decimals allowed).",
        "      Omit for no identity cutoff. Same for blastn/blastp/blastx/tblastn/tblastx.",
        "",
        f"  {opt('--reciprocal', stream=s)}",
        "      After the query vs database search, BLAST each database FASTA",
        "      against the original query (program types swapped). Keep only",
        "      pairs that are the top hit both ways (lowest e-value, then",
        "      highest bitscore). Rank 2 is not tried. Reverse max-target-seqs",
        "      is the number of query sequences. Same --min-identity on both",
        "      directions. Saves reverse_*.txt under --blast-dir.",
        "",
        f"  {opt('--output TSV', stream=s)}",
        f"      Output path (default: {DEFAULT_OUTPUT}).",
        "",
        f"  {opt('--blast-dir DIR', stream=s)}",
        "      Keep BLAST tabular files (one .txt per database FASTA).",
        f"      Default: '{DEFAULT_BLAST_DIR}/' next to the TSV. Unused with --diamond.",
        "",
        f"{heading('DIAMOND', stream=s)}  {dim('(optional; blastp substitute)', stream=s)}",
        "",
        f"  {opt('--diamond', stream=s)}",
        "      Run Diamond blastp instead of BLASTP. Protein query and protein",
        "      FASTA or .dmnd required (error otherwise). One database type per",
        "      run (all FASTA or all .dmnd). evalue and max-target-seqs use",
        "      Diamond defaults. FASTA indexes (makedb) are temporary.",
        "      With --reciprocal, reverse is also diamond blastp (FASTA dbs only;",
        "      .dmnd cannot be used as the reverse query).",
        "",
        f"  {opt('--diamond-dir DIR', stream=s)}",
        "      Keep Diamond tabular files (one .txt per database; reverse_*.txt",
        "      if --reciprocal).",
        f"      Default: '{DEFAULT_DIAMOND_DIR}/' next to the TSV (only with --diamond).",
        "",
        f"{heading('HMMER / PFAM', stream=s)}  {dim('(optional)', stream=s)}",
        "",
        f"  {opt('--hmm-db HMM', stream=s)}",
        "      HMMER profile database (e.g. Pfam-A.hmm). Runs hmmscan on the",
        "      protein --query. Presses with hmmpress if .h3m/.h3i/.h3f/.h3p",
        "      are missing. Uses --cut_ga when the HMM has GA lines, else -E 10.",
        "      Query must be protein (error otherwise).",
        "",
        f"  {opt('--hmm-dir DIR', stream=s)}",
        "      Save hmmscan tblout as hmmscan.tbl.",
        f"      Default: '{DEFAULT_HMM_DIR}/' next to the TSV (only with --hmm-db).",
        "",
        f"{heading('KOFAMSCAN', stream=s)}  {dim('(optional; metabolic KO assignment)', stream=s)}",
        "",
        f"  {opt('--kofam-profile PATH', stream=s)}",
        "      KOfam HMM profiles: a directory of .hmm files, one .hmm file, or",
        "      a .hal list (e.g. profiles/eukaryote.hal). Requires --ko-list.",
        "      Runs exec_annotation (detail-tsv) on a protein --query.",
        "      Only hits above the adaptive KO threshold (*) are kept.",
        "",
        f"  {opt('--ko-list FILE', stream=s)}",
        "      KOfam ko_list (KEGG ortholog thresholds). Required with",
        "      --kofam-profile.",
        "",
        f"  {opt('--kofam-dir DIR', stream=s)}",
        "      Save kofam.tsv and keep hmmsearch intermediates in tmp/.",
        f"      Default: '{DEFAULT_KOFAM_DIR}/' next to the TSV",
        "      (only with --kofam-profile and --ko-list).",
        "",
        f"{heading('PROGRAM SELECTION', stream=s)}",
        "",
        "  Query vs database decides the BLAST program:",
        "",
        "    protein    × protein       blastp",
        "    nucleotide × nucleotide    blastn   (or tblastx with --tblastx)",
        "    nucleotide × protein       blastx",
        "    protein    × nucleotide    tblastn",
        "",
        f"{heading('EXAMPLES', stream=s)}",
        "",
        "  Protein query vs one proteome:",
        "",
        example(
            "    biocol run --query query.faa --db species.faa \\\n"
            "      --accessions species.txt --output results.tsv",
            stream=s,
        ),
        "",
        "  Several species (folder of FASTA files, same type):",
        "",
        example(
            "    biocol run --query query.faa --db databases/ \\\n"
            "      --accessions all_species.txt --output results.tsv",
            stream=s,
        ),
        "",
        "  Nucleotide vs nucleotide:",
        "",
        example(
            "    biocol run --query query.fna --db other.fna --accessions other.txt\n"
            "    biocol run --query query.fna --db other.fna --accessions other.txt --tblastx",
            stream=s,
        ),
        "",
        "  Protein query plus Pfam (hmmscan):",
        "",
        example(
            "    biocol run --query query.faa --db species.faa \\\n"
            "      --accessions species.txt --hmm-db Pfam-A.hmm",
            stream=s,
        ),
        "",
        "  Reciprocal BLAST (RBH; top hit in both directions):",
        "",
        example(
            "    biocol run --query query.faa --db species.faa \\\n"
            "      --accessions species.txt --reciprocal",
            stream=s,
        ),
        "",
        "  Diamond blastp (protein vs protein FASTA or .dmnd):",
        "",
        example(
            "    biocol run --query query.faa --db species.faa \\\n"
            "      --accessions species.txt --diamond\n"
            "    biocol run --query query.faa --db species.dmnd \\\n"
            "      --accessions species.txt --diamond --diamond-dir diamond",
            stream=s,
        ),
        "",
        "  Protein query plus KofamScan (KO assignment):",
        "",
        example(
            "    biocol run --query query.faa --db species.faa \\\n"
            "      --accessions species.txt \\\n"
            "      --kofam-profile profiles/eukaryote.hal --ko-list ko_list",
            stream=s,
        ),
        "",
        f"{heading('NOTES', stream=s)}",
        "",
        "  --accessions must match SUBJECT ids in --db, not the query.",
        "  Missing descriptors are written as ---.",
        "  Empty query columns (e.g. cDNA on a protein-only run) are omitted.",
        "  --hmm-db is optional; omit it to skip hmmscan.",
        "  --reciprocal is optional; omit it for a one-way BLAST.",
        "  --diamond is optional; omit it to use BLAST+.",
        "  --kofam-profile and --ko-list are optional together; omit both to",
        "  skip KofamScan. Passing only one of them is an error.",
        "",
        f"{heading('SEE ALSO', stream=s)}",
        "",
        f"  {cmd('biocol --help', stream=s)}           Overview and quick start.",
        f"  {cmd('biocol from-blast --help', stream=s)}  TSV from an existing BLAST file.",
        "",
    ]
    return "\n".join(parts)


def render_from_blast_help() -> str:
    s = sys.stdout
    parts = [
        banner("from-blast — BLAST tabular → TSV", stream=s),
        "",
        "Parse a BLAST outfmt 6 file, join accession descriptors, and write the",
        "same TSV as path 1. Does not run BLAST+. Sequence columns stay empty.",
        "",
        f"{heading('USAGE', stream=s)}",
        "",
        f"  {cmd('biocol from-blast', stream=s)} {opt('--blast', stream=s)} FILE "
        f"{opt('--accessions', stream=s)} FILE [{opt('--output', stream=s)} TSV] "
        f"[{opt('--min-identity', stream=s)} N] [{opt('--hmm-db', stream=s)} HMM] "
        f"[{opt('--kofam-profile', stream=s)} PATH]",
        "",
        f"{heading('REQUIRED', stream=s)}",
        "",
        f"  {opt('--blast FILE', stream=s)}",
        "      BLAST tabular file (outfmt 6, typically .txt).",
        "",
        f"  {opt('--accessions FILE', stream=s)}",
        "      .txt file: accession<TAB>descriptor (no header).",
        "",
        f"{heading('OPTIONAL', stream=s)}",
        "",
        f"  {opt('--output TSV', stream=s)}",
        f"      Output path (default: {DEFAULT_OUTPUT}).",
        "",
        f"  {opt('--min-identity N', stream=s)}",
        "      Keep HSPs with BLAST pident >= N (0-100, decimals allowed).",
        "      Omit for no identity cutoff.",
        "",
        f"  {opt('--hmm-db HMM', stream=s)}",
        "      Run hmmscan; requires --protein (amino-acid FASTA).",
        "",
        f"  {opt('--protein FASTA', stream=s)}",
        "      Protein query for hmmscan and/or KofamScan.",
        "",
        f"  {opt('--hmm-dir DIR', stream=s)}",
        f"      hmmscan tblout directory (default: {DEFAULT_HMM_DIR}/ next to the TSV).",
        "",
        f"  {opt('--kofam-profile PATH', stream=s)}",
        "      Run KofamScan; requires --ko-list and --protein.",
        "",
        f"  {opt('--ko-list FILE', stream=s)}",
        "      KOfam ko_list (required with --kofam-profile).",
        "",
        f"  {opt('--kofam-dir DIR', stream=s)}",
        f"      exec_annotation output directory (default: {DEFAULT_KOFAM_DIR}/ next to the TSV).",
        "",
        f"{heading('EXAMPLES', stream=s)}",
        "",
        example(
            "    biocol from-blast --blast hits.txt \\\n"
            "      --accessions species.txt --output results.tsv",
            stream=s,
        ),
        "",
        "  With KofamScan:",
        "",
        example(
            "    biocol from-blast --blast hits.txt --accessions species.txt \\\n"
            "      --protein query.faa --kofam-profile profiles/eukaryote.hal \\\n"
            "      --ko-list ko_list",
            stream=s,
        ),
        "",
        f"{heading('NOTES', stream=s)}",
        "",
        "  No --query / --cdna: sequence columns stay empty.",
        "  Reciprocal BLAST is only available in 'biocol run', not from-blast.",
        "  Diamond search is only available in 'biocol run'; a Diamond outfmt 6",
        "  file can still be passed as --blast.",
        "  KofamScan needs --kofam-profile, --ko-list, and --protein together.",
        "  The species block name is the accessions file stem",
        "  (Benincasa_hispida_gd.txt → Benincasa_hispida_gd).",
        "",
        f"{heading('SEE ALSO', stream=s)}",
        "",
        f"  {cmd('biocol --help', stream=s)}     Overview and quick start.",
        f"  {cmd('biocol run --help', stream=s)}  Full FASTA → BLAST+ workflow.",
        "",
    ]
    return "\n".join(parts)
