# biocol

Backend for INECOL's BLAST tool. It takes a FASTA query (or an existing tabular BLAST result), compares it against local FASTA databases, and writes a **plain TSV** annotation file.


## Requirements

### Conda environment

The tests were performed using the Conda environment provided in `environment.yml`.

To create the environment, run:

```bash
conda env create -f environment.yml
```

Activate the environment:
```bash
conda activate inecol
```

- Python 3.10+
- BLAST+ in PATH (`conda activate inecol` on Ubuntu/WSL)
- Diamond in PATH only if you use `diamond=True`
- KofamScan (`exec_annotation`), HMMER, and GNU Parallel in PATH if you use `run_kofamscan`
- Dependencies: Biopython and pandas

```bash
cd biocol
pip install -e ".[dev]"
```

`-e` installs the package in editable mode: backend changes are visible without reinstalling.

## Inputs

| Input | Format |
|---------|---------|
| Query | FASTA / multifasta: `.fa`, `.fasta`, `.fna`, `.faa`, `.fas`. The entire file must be the same type (DNA, RNA, or protein). `U` counts as a nucleotide. |
| Databases | A FASTA file or a **folder** (including subfolders). Same extensions. All databases must be the same type. One BLAST per FASTA file. |
| Accessions | Text `accession<TAB>descriptor`, without a header. |
| Tabular BLAST (path 2) | `outfmt 6`: 12 NCBI columns, or 15 if it includes `nident`, `qseq`, and `sseq` (as written by `run_blast`). |

BLAST parameters (modifiable): `evalue` default **10**, `max_target_seqs` default **3**, `threads` default **1**. `min_identity` is optional: if it is not provided, there is **no cutoff** on identity percentage (`pident`). The TSV shows only the **best hit** per query and species.

## How the BLAST program is selected

`run_blast` detects the type of the query and databases and calls `select_blast_program`. If **both** are nucleotide, the default is **blastn**. `translated=True` (CLI flag `--tblastx`) selects **tblastx**. In all other combinations, `translated` is ignored.

| Query | Database | `translated` | Program | What it compares |
|-------|---------------|--------------|----------|-------------|
| Nucleotide | Nucleotide | not provided / `False` (default) | `blastn` | Nucleotide against nucleotide |
| Nucleotide | Nucleotide | `True` (explicit) | `tblastx` | Query and database translated in six frames (protein) |
| Nucleotide | Protein | ignored | `blastx` | Translated nucleotide query against proteins |
| Protein | Protein | ignored | `blastp` | Protein against protein |
| Protein | Nucleotide | ignored | `tblastn` | Protein against translations of the nucleotide database |

## Usage

There are two paths to the same TSV.

### Path 1 — FASTA + databases

```python
from biocol import run_blast, build_result_table, write_results_csv

hits = run_blast(
    "query.fa",
    "bases/",                 # un FASTA o carpeta
    translated=False,         # True → tblastx si query y base son nucleótido
    evalue=10,
    max_target_seqs=3,
    blast_dir="blast",        # tabular outfmt 6; si se omite, no se guarda
    min_identity=80,          # opcional; si se omite, no se filtra por pident
    reciprocal=False,         # True → BLAST de vuelta y solo pares top en ambos sentidos
    diamond=False,            # True → diamond blastp (query y base proteína; FASTA o .dmnd)
    diamond_dir="diamond",    # tabular Diamond; si se omite, no se guarda
)
table = build_result_table(hits, "accessions.txt", query_fasta="query.fa")
# optional: hmmscan (protein query + Pfam-A.hmm or another HMMER3 database)
# hmm_hits = run_hmmscan("query.faa", "Pfam-A.hmm", hmm_dir="hmm")
# table = build_result_table(hits, "accessions.txt", query_fasta="query.faa", hmm_hits=hmm_hits)
# kofam_hits = run_kofamscan("query.faa", "profiles/eukaryote.hal", "ko_list", kofam_dir="kofam")
# table = build_result_table(hits, "accessions.txt", query_fasta="query.faa", kofam_hits=kofam_hits)
write_results_csv(table, "results.tsv")  # si se omite, usa results.tsv
```

`run_blast` creates temporary databases with `makeblastdb` (they are deleted when finished), runs one BLAST per FASTA file, and parses `outfmt 6`. If you pass `blast_dir`, it stores one `.txt` file per database there **without filtering**. `min_identity` discards HSPs with a `pident` below the threshold (the same applies to blastn and blastp). If a query has no hit in a database, an empty row is kept.

`reciprocal=True` (default `False`) also performs the reverse comparison: each FASTA from `--db` becomes the query, and the original FASTA becomes the database. The program is selected using the swapped types (`blastx` ↔ `tblastn`; `tblastx` ↔ `tblastx`). In the reverse run, `max_target_seqs` is the number of sequences in the original query. The same `min_identity` is applied in both directions. Only the HSPs from pairs that are the **top hit in both directions** are retained (lowest e-value, then highest bitscore; hit #2 is not tested). If there is no reciprocal hit, the cell is set to `---`. With `blast_dir`, `reverse_<stem>.txt` is also saved. `from-blast` does not support this mode.

`diamond=True` (default `False`, backend only) replaces **blastp**: the query and database must both be protein. Two database input types are supported, with only one type per run: a **`.dmnd`** file/folder, or a protein FASTA (in which case `diamond makedb` is run on a temporary database that is deleted afterward). `evalue` and `max-target-seqs` use Diamond's defaults; `threads` defaults to **1**. The tabular output is saved in `diamond_dir` (not in `blast/`). If `reciprocal=True` and the databases are protein FASTA files, the reverse run also uses `diamond blastp`. A `.dmnd` file cannot be used for the reverse run (there is no FASTA query for the species) → error. Nucleotide query or database → `diamond requires a protein query and a protein database`. A Diamond tabular file using the same `outfmt 6` format can be read with `from-blast` / `parse_blast_results`.

`run_kofamscan` calls `exec_annotation` on a **protein** FASTA. `profile` can be a folder containing `.hmm` files, a `.hmm` file, or a `.hal` file (e.g. `profiles/eukaryote.hal`); `ko_list` is the KOfam `ko_list` file. The format is `detail-tsv`; only rows containing `*` are retained. The KofamScan TSV is stored in `kofam_dir` (`kofam/kofam.tsv`), and the temporary hmmsearch files are stored in `kofam_dir/tmp` (they are not deleted). `num_threads` defaults to **1**. Nucleotide query → `kofamscan requires a protein sequence`. Pass it to `build_result_table(..., kofam_hits=...)`. In the CLI: `--kofam-profile` and `--ko-list` must be provided together (providing only one → error).

### Path 2 — Existing tabular BLAST

Without a query FASTA: `qseqid` + hits + descriptors are populated; sequence columns remain empty.

```python
from biocol import parse_blast_results, filter_hits_by_pident, build_result_table, write_results_csv

hits = parse_blast_results("blast_outfmt6.txt")
hits = filter_hits_by_pident(hits, 80)  # opcional
table = build_result_table(hits, "accessions.txt", query_fasta=None)
write_results_csv(table)
```

If the tabular result does not contain a `database` column, the `hit` name is used.

### CLI

After `pip install -e ".[dev]`:

```bash
biocol run --query query.fa --db bases/ --accessions accessions.txt
biocol run --query query.fa --db bases/ --accessions accessions.txt --min-identity 80
biocol run --query query.fa --db bases/ --accessions accessions.txt --tblastx --evalue 1e-5 --max-target-seqs 50 --threads 4 --output my_results.tsv --blast-dir blast_hits
biocol run --query genes.faa --db bases/ --accessions accessions.txt --cdna genes.fna --protein genes.faa
biocol run --query genes.faa --db bases/ --accessions accessions.txt --hmm-db Pfam-A.hmm --hmm-dir hmm
biocol run --query query.fa --db bases/ --accessions accessions.txt --reciprocal
biocol run --query query.faa --db species.faa --accessions accessions.txt --diamond
biocol run --query query.faa --db species.dmnd --accessions accessions.txt --diamond --diamond-dir diamond
biocol run --query query.faa --db species.faa --accessions accessions.txt --kofam-profile profiles/eukaryote.hal --ko-list ko_list

biocol from-blast --blast hits.txt --accessions Benincasa_hispida_gd.txt
biocol from-blast --blast hits.txt --accessions Benincasa_hispida_gd.txt --min-identity 80.5 --output my_results.tsv
biocol from-blast --blast hits.txt --accessions Benincasa_hispida_gd.txt --protein query.faa --kofam-profile profiles/eukaryote.hal --ko-list ko_list
```

`--output` is optional (default: `results.tsv`). `biocol run` also writes BLAST tabular files to `--blast-dir` (default: `blast/` next to the TSV). Optional `--hmm-db Pfam-A.hmm` runs hmmscan (protein query; `hmmpress` if needed) and writes `hmm/hmmscan.tbl`. Optional `--reciprocal` runs reverse BLAST and keeps only rank-1 pairs in both directions. Optional `--diamond` runs Diamond blastp (protein vs protein FASTA or `.dmnd`; tabular output in `--diamond-dir`, default `diamond/` next to the TSV). Optional `--kofam-profile` + `--ko-list` run KofamScan `exec_annotation` (protein query; detail-tsv in `--kofam-dir`, default `kofam/` next to the TSV, tmp retained). Help text and errors are in English.

`from-blast` does not take a query FASTA and does not support `--reciprocal` or `--diamond` (a Diamond `outfmt 6` file can still be passed as `--blast`). KofamScan in `from-blast` requires `--protein`. The species name in the TSV header is the accessions file stem (`Benincasa_hispida_gd.txt` → `Benincasa hispida gd`).

## Output TSV

TSV with **three header rows**, like Dataset S2 (Pfam and KOfam optional; no GO). There are no merged cells: the section name and species name are repeated or placed in the first column of each block.

1. Section: empty for query; `Annotation based on top-BLAST-hit method` in each BLAST block; `Pfam domains` if hmmscan was run; `KOfamScan` if KofamScan was run.
2. Species: `stem` of the database FASTA (e.g. `protein`, `amborella`). Empty for Pfam and KOfam.
3. Column names: `Gene ID`, `Length (nt)`, `cDNA Sequences (nt)`, `Length(aa)`, `Protein Sequences (aa)`, and per species `Accesion No.`, `Description`, `Identity %`, `Identity % (full query)`, `Alignment length`, `e-value`, `Score`. If hmmscan is present: `# of Pfam domain identified`, `e-value`, `score`, `Accesion`, `Name`, `Description of target` (all domains that passed the cutoff, separated by `; `; no domain → `---`). If KofamScan is present: `KO`, `score`, `e-value`, `KO definition` (only assignments above the KO threshold, marked with `*`; multiple KOs from the same protein joined with `; `; no KO → `---`).

`Identity %` is the BLAST `pident` (with respect to the alignment). `Identity % (full query)` uses the aligned `qseq`/`sseq` sequences: **unique** identical query positions ÷ full query length × 100 (aa if the query is protein, nt if it is nucleotide). Multiple HSPs from the same subject are merged without counting an overlap twice. In blastx/tblastx, each identical amino acid covers 3 nt of the query. `--cdna` is not included in the denominator. No hit, no `qseq`/`sseq` (12-column tabular output), or length 0: `---`.

`Length (aa)` and hit sequence are not included. No hit or no descriptor: `---`.

Empty query columns are **not written** (protein query → no cDNA; nucleotide query → no protein). If gene models are available, both FASTA files (`--cdna` and `--protein`) can be provided and the first block is complete.

One row per query (only the best hit per species). In Excel, import the TSV and optionally merge cells in the first two rows.

## Public API

```python
from biocol import (
    # FASTA
    validate_fasta_file,
    read_fasta,
    detect_sequence_type,
    detect_query_type,
    # BLAST
    detect_database_type,
    list_blast_databases,
    select_blast_program,
    run_blast,
    parse_blast_results,
    run_hmmscan,
    parse_hmmscan_tblout,
    run_kofamscan,
    parse_kofam_detail_tsv,
    # Resultados
    load_accessions,
    filter_hits_by_pident,
    build_result_table,
    write_results_csv,
    QUERY_COLUMNS,
    DEFAULT_OUTPUT,
)
```

`detect_sequence_type()` accepts `str`, `Bio.Seq.Seq`, or `SeqRecord`. `detect_query_type()` classifies a complete FASTA by the majority type (protein in case of a tie).

Errores: `FastaError`, `EmptyFastaError`, `InvalidFastaError`, `MixedSequenceTypeError`, `BlastError`, `DatabaseError`, `MixedDatabaseTypeError`, `BlastExecutionError`, `DiamondError`, `DiamondExecutionError`, `HmmError`, `HmmExecutionError`, `KofamError`, `KofamExecutionError`, `MetadataError`.

## Tests (Dana)

```bash
pytest -q
pytest tests/test_detect_sequence_type.py -q
```

Running Pytest displays INFO logs. Fixtures are in `tests/fixtures/` (FASTA, BLAST `outfmt 6`, accessions). Import only from `biocol`.

`run_blast` against real BLAST+ only in the `conda` `inecol` environment.

| Input | Expected result |
|---------|--------------------|
| DNA FASTA | `nucleotide` |
| RNA FASTA | `nucleotide` |
| Protein FASTA | `protein` |
| Multifasta of the same type | single type |
| Mixed DNA + protein multifasta | majority type (tie → protein) |
| Nonexistent path | `FileNotFoundError` |
| Non-FASTA extension (e.g. `.txt`) | `InvalidFastaError` |
| Empty file / no sequences | `EmptyFastaError` |
| Content that is not FASTA | `InvalidFastaError` |
| Empty sequence or invalid characters | `InvalidFastaError` |

## Team and structure

| Person  | Role    | Works on |
|----------|--------|----------------|
| Alondra  | logic | `src/biocol/`  |
| Emiliano | CLI    | `src/biocol/cli/` (only the public `biocol` API) |
| Dana     | pruebas | `tests/` |

```
src/biocol/           backend
  sequence/           reading, validation, and query type
  blast/              database type, selection, execution, tabular parsing, and reciprocal BLAST
  diamond/            optional makedb / blastp (protein vs protein)
  hmm/                hmmpress, hmmscan, and tblout parsing
  kofam/              exec_annotation (KofamScan) and detail-tsv parsing
  metadata/           accessions and descriptors
  processing/         wide table, identity (full query), and HSP filtering
  output/             TSV writing
  cli/                CLI (argparse + run / from-blast commands)
tests/                tests (pytest)
tests/fixtures/       FASTA, BLAST outfmt 6, and example accessions
```

FASTA reading uses `Bio.SeqIO`. The nucleotide/protein alphabets come from `Bio.Data.IUPACData`.

Developer reference: [`docs/backend.md`](docs/backend.md) (backend), [`docs/cli.md`](docs/cli.md) (CLI).
