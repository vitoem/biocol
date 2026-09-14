# biocol backend — module and function reference

This document describes how the backend is organized and what each module and function does. Use it when changing or extending the package. 

---

## Pipeline (what the backend does)

1. Read and type a query FASTA (`sequence/`).
2. Search subjects:
  - BLAST+ (`blast/`), or
  - Diamond blastp (`diamond/`, via `run_blast(..., diamond=True)`).
3. Optionally filter HSPs by `pident` (`processing/hsp_filter.py`).
4. Optionally keep reciprocal best hits (`blast/reciprocal.py`).
5. Optionally run hmmscan (`hmm/`) and/or KofamScan (`kofam/`).
6. Join accession descriptors and build one wide row per query (`processing/table.py`).
7. Write a Dataset S2-style TSV with three header rows (`output/`).

One TSV row per query. BLAST/Diamond show **one best hit per species** (lowest e-value, then highest bitscore; ties keep the first row). Pfam and KOfam concatenate all passing hits with `"; "`. Missing values are `---`.

---



## Package layout

```
src/biocol/
  __init__.py          Public re-exports
  exceptions.py        Error types
  sequence/            FASTA I/O, alphabets, nucleotide vs protein
  blast/               BLAST+ databases, program choice, run, parse, RBH
  diamond/             Diamond makedb / blastp
  hmm/                 hmmpress + hmmscan
  kofam/               exec_annotation (KofamScan)
  metadata/            accession → descriptor
  processing/          identity, HSP filter, wide table
  output/              TSV writer
  cli/                 argparse + run / from-blast
```

Private helpers are prefixed with `_`. Callers outside the same submodule should not depend on them.

---



## `biocol/__init__.py`

Re-exports the public API. When you add a user-facing function or exception, export it here **and** add it to `__all__`.

---



## `exceptions.py`

All domain errors inherit from `ValueError` so the CLI can catch a single family.


| Class                    | When it is raised                                                      |
| ------------------------ | ---------------------------------------------------------------------- |
| `FastaError`             | Base for FASTA problems                                                |
| `EmptyFastaError`        | File exists but has no sequences                                       |
| `InvalidFastaError`      | Bad extension, unparseable FASTA, empty seq, illegal residues          |
| `MixedSequenceTypeError` | Defined; mixed FASTA currently does **not** abort (majority type wins) |
| `BlastError`             | BLAST selection / execution base                                       |
| `DatabaseError`          | Missing or unrecognized FASTA / `.dmnd` source                         |
| `MixedDatabaseTypeError` | Folder mixes nucleotide and protein (or FASTA + `.dmnd` for Diamond)   |
| `BlastExecutionError`    | `makeblastdb` or BLAST+ missing or non-zero exit                       |
| `MetadataError`          | Accessions file unreadable or empty / missing ids                      |
| `HmmError`               | HMMER / hmmscan (e.g. nucleotide query)                                |
| `HmmExecutionError`      | `hmmpress` / `hmmscan` missing or failed                               |
| `DiamondError`           | Diamond usage (non-protein, `.dmnd` + reciprocal, etc.)                |
| `DiamondExecutionError`  | `diamond` missing or failed                                            |
| `KofamError`             | KofamScan usage (non-protein, bad profile, one of two CLI flags)       |
| `KofamExecutionError`    | `exec_annotation` missing or failed                                    |


Messages are English strings. Keep them stable if tests match them.

---



## `sequence/`



### `alphabets.py`

IUPAC sets from `Bio.Data.IUPACData`.


| Name                             | Role                                                                          |
| -------------------------------- | ----------------------------------------------------------------------------- |
| `NUCLEOTIDE_LETTERS`             | Ambiguous DNA + RNA (includes K, R, Y, … which also appear in proteins)       |
| `UNAMBIGUOUS_NUCLEOTIDE_LETTERS` | A, C, G, T, U — used so a peptide of only A/C/G is not called DNA             |
| `PROTEIN_LETTERS`                | Extended protein alphabet                                                     |
| `GAP_LETTERS`                    | `-` and `.` (ignored for typing)                                              |
| `VALID_RESIDUES`                 | Letters allowed in a FASTA sequence                                           |
| `PROTEIN_ONLY_LETTERS`           | Letters that force protein classification (E, F, I, L, P, Q, X, O, Z, J, `*`) |


If you change alphabets, re-check `detect_sequence_type` and `validate_seq_record`.

### `validator.py`


| Function              | Role                                                       |
| --------------------- | ---------------------------------------------------------- |
| `check_fasta_path`    | Exists, is a file, suffix in `{.fa,.fasta,.fna,.faa,.fas}` |
| `validate_seq_record` | Non-empty sequence; every character in `VALID_RESIDUES`    |
| `validate_fasta_file` | Path check + full parse via `read_fasta`                   |




### `reader.py`


| Function     | Role                                                         |
| ------------ | ------------------------------------------------------------ |
| `read_fasta` | `Bio.SeqIO.parse`, validate each record, uppercase sequences |




### `classifier.py`


| Function               | Role                                                                                                                                                                          |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `_residues`            | Uppercase letters, gaps stripped                                                                                                                                              |
| `_sequence_id`         | `SeqRecord.id` or `"<no id>"` for logs                                                                                                                                        |
| `detect_sequence_type` | One sequence → `"nucleotide"` or `"protein"`: any protein-only letter → protein; else all nucleotide letters **and** at least one unambiguous base → nucleotide; else protein |
| `detect_query_type`    | Whole FASTA: majority type; **ties → protein**; mixed files do not raise                                                                                                      |


---



## `blast/`



### `databases.py`


| Function                   | Role                                                                                                          |
| -------------------------- | ------------------------------------------------------------------------------------------------------------- |
| `infer_database_label`     | Majority `[Organism]` in NCBI headers (sample 40 records); else file stem. Becomes the TSV species block name |
| `_is_fasta_file`           | File with a FASTA suffix                                                                                      |
| `_type_from_fasta`         | `detect_query_type` on that file                                                                              |
| `_iter_fasta_in_directory` | Recursive FASTA list, sorted                                                                                  |
| `list_blast_databases`     | File or folder → `list[(Path, type)]`                                                                         |
| `detect_database_type`     | Single type for the whole source; mixed → `MixedDatabaseTypeError`                                            |




### `selection.py`


| Function               | Role                                                                                                                                                       |
| ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `select_blast_program` | `(query_type, database_type, translated)` → `blastn` / `tblastx` / `blastx` / `blastp` / `tblastn`. `translated` only applies when **both** are nucleotide |




### `parser.py`


| Name                       | Role                                                                            |
| -------------------------- | ------------------------------------------------------------------------------- |
| `OUTFMT6_STANDARD_COLUMNS` | 12 NCBI columns                                                                 |
| `OUTFMT6_COLUMNS`          | 15 columns including `nident`, `qseq`, `sseq`                                   |
| `BLAST_OUTFMT`             | BLAST `-outfmt` string used by `run_blast` / Diamond                            |
| `parse_blast_results`      | Read tabular; 12 or 15 columns; empty file → empty frame with 15 columns        |
| `fill_missing_hits`        | Ensure every `query_id` has a row for `database_name`; no hit → NA BLAST fields |




### `ranking.py`


| Function          | Role                                                                                                                                            |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `assign_hit_rank` | Per `(qseqid, database)`: sort e-value ascending, bitscore descending; `hit_rank` 1, 2, … Empty `sseqid` stays rank 1 so the TSV can show `---` |


This ranking **must** stay the same in the result table and in reciprocal BLAST.

### `reciprocal.py`

RBH: keep the forward pair only if it is rank-1 **both ways**. Rank 2 is never tried. All HSPs of an accepted pair are kept.


| Function                         | Role                                                        |
| -------------------------------- | ----------------------------------------------------------- |
| `_strip_version`                 | Drop `.digits` suffix                                       |
| `sequence_ids_match`             | Same sequence after `normalize_accession` and version strip |
| `_has_subject`                   | Non-empty `sseqid`                                          |
| `_rank1_rows` / `_top_subject`   | Rank-1 subject for a query in one database                  |
| `_reverse_query_ids_for_subject` | Reverse-search query ids that match a forward subject       |
| `_drop_rank_helpers`             | Drop `_evalue_sort`, `_score_sort`, `hit_rank`              |
| `keep_reciprocal_best_hits`      | Filter forward HSPs using reverse rank-1                    |




### `runner.py`


| Name                        | Role                                                                                                                                                                              |
| --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `DEFAULT_MAX_TARGET_SEQS`   | `3`                                                                                                                                                                               |
| `DEFAULT_NUM_THREADS`       | `1`                                                                                                                                                                               |
| `DEFAULT_BLAST_DIR`         | `"blast"`                                                                                                                                                                         |
| `_tabular_filename`         | Unique `stem.txt` names                                                                                                                                                           |
| `build_makeblastdb_command` | `makeblastdb -parse_seqids`                                                                                                                                                       |
| `build_blast_command`       | BLAST+ with `BLAST_OUTFMT`, evalue, max_target_seqs, threads                                                                                                                      |
| `_run_command`              | PATH check + subprocess; failure → `BlastExecutionError`                                                                                                                          |
| `_run_forward_blast`        | Temp `makeblastdb`, one BLAST per FASTA, parse, optional save to `blast_dir`, `min_identity`, `fill_missing_hits`                                                                 |
| `_run_reverse_blast`        | Each DB FASTA as query vs original query; reverse program via swapped types; `max_target_seqs` = number of original query sequences; writes `reverse_*.txt` if `blast_dir` is set |
| `run_blast`                 | Public BLAST+ driver. `diamond=True` delegates to `run_diamond` and ignores evalue / max_target_seqs / translated / blast_dir                                                     |


Indexes live in a temp directory and are deleted. Tabular in `blast_dir` is **unfiltered** (filter is applied after parse).

---



## `diamond/`

Protein vs protein only. Not mixed FASTA and `.dmnd` in one run.

### `databases.py`


| Function                 | Role                                                        |
| ------------------------ | ----------------------------------------------------------- |
| `_is_dmnd` / `_is_fasta` | Suffix checks                                               |
| `_iter_in_directory`     | FASTA vs `.dmnd` lists                                      |
| `_fasta_entries`         | Protein FASTA only; label via `infer_database_label`        |
| `list_diamond_databases` | `(path, kind, label)` where `kind` is `"fasta"` or `"dmnd"` |




### `execute.py`


| Function              | Role                                                                     |
| --------------------- | ------------------------------------------------------------------------ |
| `run_diamond_command` | Run `diamond`; missing binary or non-zero exit → `DiamondExecutionError` |




### `runner.py`


| Function               | Role                                                                                                               |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `build_makedb_command` | `diamond makedb` (temp DB, deleted)                                                                                |
| `build_blastp_command` | `diamond blastp` with same 15-column outfmt; **no** evalue / max-target-seqs (Diamond defaults); threads default 1 |
| `_out_path`            | Keep tabular under `diamond_dir` or a temp file                                                                    |
| `run_diamond`          | Public Diamond driver (also used from `run_blast`)                                                                 |
| `_run_diamond_blastp`  | One blastp per database                                                                                            |
| `_run_reverse_diamond` | Reciprocal Diamond; **FASTA databases only** (`.dmnd` cannot be the reverse query)                                 |


---



## `hmm/`



### `execute.py`


| Function    | Role                       |
| ----------- | -------------------------- |
| `run_hmmer` | Run `hmmpress` / `hmmscan` |




### `press.py`


| Function                      | Role                                        |
| ----------------------------- | ------------------------------------------- |
| `pressed_index_paths`         | `.h3m` `.h3i` `.h3f` `.h3p` next to the HMM |
| `hmm_database_is_pressed`     | All four indexes exist                      |
| `hmm_has_gathering_threshold` | File has a `GA` line                        |
| `ensure_pressed_hmm`          | `hmmpress` if indexes missing               |




### `parser.py`


| Function               | Role                                                                                                          |
| ---------------------- | ------------------------------------------------------------------------------------------------------------- |
| `parse_hmmscan_tblout` | Space-delimited `--tblout`; skip `#`; 19 fields                                                               |
| `best_hmmscan_hits`    | One row per query (lowest e-value). The **TSV does not use this**; the table concatenates all passing domains |




### `runner.py`


| Function                | Role                                                              |
| ----------------------- | ----------------------------------------------------------------- |
| `build_hmmscan_command` | `--tblout --noali --cpu`; `--cut_ga` if GA present else `-E 10`   |
| `run_hmmscan`           | Protein query only; writes `hmm_dir/hmmscan.tbl` (default `hmm/`) |


---



## `kofam/`

Wraps KofamScan `exec_annotation`. Does **not** call hmmsearch/hmmscan itself.

### `execute.py`


| Function              | Role                                    |
| --------------------- | --------------------------------------- |
| `run_exec_annotation` | PATH + subprocess for `exec_annotation` |




### `parser.py`


| Function                 | Role                                                                                                        |
| ------------------------ | ----------------------------------------------------------------------------------------------------------- |
| `parse_kofam_detail_tsv` | Keep only rows whose first field is `*`; columns: mark, gene_name, ko, threshold, score, evalue, definition |




### `runner.py`


| Function                        | Role                                                                                            |
| ------------------------------- | ----------------------------------------------------------------------------------------------- |
| `build_exec_annotation_command` | `-f detail-tsv --profile --ko-list --cpu --tmp-dir`; no `-E` / `-T`                             |
| `_resolve_profile`              | Directory of `.hmm`, or a `.hmm` / `.hal` file                                                  |
| `_resolve_ko_list`              | Must be a file                                                                                  |
| `run_kofamscan`                 | Protein FASTA; output `kofam_dir/kofam.tsv`; intermediates in `kofam_dir/tmp` (**not deleted**) |


---



## `metadata/`



### `accessions.py`


| Function              | Role                                                         |
| --------------------- | ------------------------------------------------------------ |
| `normalize_accession` | Strip BLAST-style prefixes (`ref|XP_1.1|` → `XP_1.1`)        |
| `load_accessions`     | `accession<TAB>descriptor`, no header; adds `accession_norm` |


---



## `processing/`



### `hsp_filter.py`


| Function                | Role                                                                                              |
| ----------------------- | ------------------------------------------------------------------------------------------------- |
| `filter_hits_by_pident` | Drop HSPs with `pident < min_identity`; `None` = no cutoff; restore empty rows per query/database |




### `identity.py`

Full-query identity (not BLAST `pident`).


| Function                      | Role                                                                                                                    |
| ----------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `_query_step`                 | 3 nt per aligned residue for blastx/tblastx when `qend-qstart+1 == 3 * residues`; else 1                                |
| `identical_query_positions`   | 1-based query coords identical in one HSP (`qseq`/`sseq`)                                                               |
| `full_query_identity_percent` | Union of those coords across HSPs of **one** query–subject pair, divided by full query length × 100. `None` → TSV `---` |


Overlapping HSPs must not double-count positions.

### `table.py`

Wide Dataset S2 table.


| Name                             | Role                                                           |
| -------------------------------- | -------------------------------------------------------------- |
| `QUERY_COLUMNS`                  | gene_id, length_nt, cdna_sequence, length_aa, protein_sequence |
| `HIT_FIELDS`                     | Maps TSV fields to BLAST columns + full-query identity         |
| `PFAM_COLUMNS` / `KOFAM_COLUMNS` | Optional trailing blocks                                       |
| `_safe_name`                     | Species prefix from database label                             |
| `_empty_query_row`               | NA query fields                                                |
| `_strip_version`                 | Same version strip as reciprocal                               |
| `_record_lookup_keys`            | Join NCBI CDS FASTA to protein ids (`[protein_id=…]`, `_cds_`) |
| `_apply_fasta_to_rows`           | Fill nt or aa columns from a FASTA                             |
| `_query_metadata`                | Combine `--query` / `--cdna` / `--protein`                     |
| `_query_length_info`             | Lengths for identity (gaps/dots ignored)                       |
| `_lookup_description`            | Match `sseqid` to accessions                                   |
| `_annotation_lookup`             | Group Pfam/KOfam rows by query, sort by e-value                |
| `_pfam_*` / `_kofam_*`           | Cells joined with `"; "`; empty → `---`                        |
| `build_result_table`             | Public join: BLAST + accessions + optional hmm/kofam           |


`build_result_table(..., hmm_hits=, kofam_hits=)` only adds those blocks if the DataFrame is passed (not `None`).

---



## `output/`



### `csv_writer.py`


| Function                   | Role                                                                                                                                                                                                                      |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `drop_empty_query_columns` | Drop unused query columns (e.g. cDNA on a protein-only run); keep `gene_id`                                                                                                                                               |
| `_hit_prefix_and_field`    | Split `species_accession` into prefix + field                                                                                                                                                                             |
| `format_s2_csv`            | Three header rows: section, species, column labels. BLAST section = `Annotation based on top-BLAST-hit method`; Pfam = `Pfam domains`; KOfam = `KOfamScan`. Species name only on the accession column of each BLAST block |
| `write_results_csv`        | Default `results.tsv`; forces `.tsv` suffix                                                                                                                                                                               |


BLAST hit prefixes must exclude `QUERY_COLUMNS`, `PFAM_COLUMNS`, and `KOFAM_COLUMNS` so `kofam_score` is not treated as a BLAST field.

---



## `cli/`

Uses **only** `from biocol import ...`.


| Module        | Role                                                               |
| ------------- | ------------------------------------------------------------------ |
| `parser.py`   | argparse; custom help via `helptext.py`                            |
| `helptext.py` | Colored English help for top / run / from-blast                    |
| `commands.py` | `run_from_fasta`, `run_from_blast`; optional hmmscan and KofamScan |
| `main.py`     | Logging, dispatch, catch domain errors → exit 1                    |
| `style.py`    | ANSI colors; `--no-color` / `NO_COLOR`                             |
| `__main__.py` | `python -m biocol.cli`                                             |




### `commands.py` helpers


| Function              | Role                                                                     |
| --------------------- | ------------------------------------------------------------------------ |
| `_dir_next_to_output` | Default `blast/`, `hmm/`, `kofam/`, `diamond/` next to the TSV           |
| `_maybe_hmmscan`      | If `--hmm-db`; protein FASTA required                                    |
| `_maybe_kofamscan`    | Needs **both** `--kofam-profile` and `--ko-list`; protein FASTA required |
| `run_from_fasta`      | BLAST/Diamond then optional Pfam/KOfam then TSV                          |
| `run_from_blast`      | Parse existing tabular; hmm/kofam need `--protein`                       |


`from-blast` has no `--reciprocal` or `--diamond` (a Diamond outfmt 6 file can still be `--blast`).

---



# Traducción al español

# Backend de biocol — referencia de módulos y funciones

Este documento describe cómo está organizado el backend y qué hace cada módulo y cada función. Úsalo al cambiar o ampliar el paquete.

## Flujo (qué hace el backend)

1. Lee y clasifica un FASTA de query (`sequence/`).
2. Busca sujetos:
  - BLAST+ (`blast/`), o
  - Diamond blastp (`diamond/`, vía `run_blast(..., diamond=True)`).
3. Opcionalmente filtra HSP por `pident` (`processing/hsp_filter.py`).
4. Opcionalmente conserva reciprocal best hits (`blast/reciprocal.py`).
5. Opcionalmente corre hmmscan (`hmm/`) y/o KofamScan (`kofam/`).
6. Une descriptores de accesiones y arma una fila ancha por query (`processing/table.py`).
7. Escribe un TSV estilo Dataset S2 con tres filas de cabecera (`output/`).

Una fila de TSV por query. BLAST/Diamond muestran **un mejor hit por especie** (menor e-value, luego mayor bitscore; empates conservan la primera fila). Pfam y KOfam concatenan todos los hits que pasaron el corte con `"; "`. Los valores faltantes son `---`.

---



## Estructura del paquete

```
src/biocol/
  __init__.py          Reexportación pública
  exceptions.py        Tipos de error
  sequence/            E/S FASTA, alfabetos, nucleótido vs proteína
  blast/               Bases BLAST+, elección de programa, ejecución, parseo, RBH
  diamond/             Diamond makedb / blastp
  hmm/                 hmmpress + hmmscan
  kofam/               exec_annotation (KofamScan)
  metadata/            accesión → descriptor
  processing/          identidad, filtro de HSP, tabla ancha
  output/              escritura del TSV
  cli/                 argparse + run / from-blast
```

Las funciones auxiliares privadas llevan prefijo `_`. Quienes estén fuera del mismo submódulo no deben depender de ellas.

---



## `biocol/__init__.py`

Reexporta la API pública. Si agregas una función o excepción visible para el usuario, expórtala aquí **y** añádela a `__all__`.

---



## `exceptions.py`

Todos los errores de dominio heredan de `ValueError` para que la CLI pueda atrapar una sola familia.


| Clase                    | Cuándo se lanza                                                                |
| ------------------------ | ------------------------------------------------------------------------------ |
| `FastaError`             | Base de problemas de FASTA                                                     |
| `EmptyFastaError`        | El archivo existe pero no tiene secuencias                                     |
| `InvalidFastaError`      | Extensión inválida, FASTA no parseable, secuencia vacía, residuos ilegales     |
| `MixedSequenceTypeError` | Definida; un FASTA mixto **no** aborta (gana el tipo mayoritario)              |
| `BlastError`             | Base de selección / ejecución BLAST                                            |
| `DatabaseError`          | Fuente FASTA / `.dmnd` ausente o no reconocida                                 |
| `MixedDatabaseTypeError` | La carpeta mezcla nucleótido y proteína (o FASTA + `.dmnd` en Diamond)         |
| `BlastExecutionError`    | Falta `makeblastdb` o BLAST+, o el proceso sale distinto de cero               |
| `MetadataError`          | Archivo de accesiones ilegible, vacío o sin ids                                |
| `HmmError`               | HMMER / hmmscan (p. ej. query nucleótido)                                      |
| `HmmExecutionError`      | Falta o falla `hmmpress` / `hmmscan`                                           |
| `DiamondError`           | Uso de Diamond (no proteína, `.dmnd` + recíproco, etc.)                        |
| `DiamondExecutionError`  | Falta o falla `diamond`                                                        |
| `KofamError`             | Uso de KofamScan (no proteína, perfil inválido, solo una de las dos flags CLI) |
| `KofamExecutionError`    | Falta o falla `exec_annotation`                                                |


Los mensajes son cadenas en inglés. Manténlos estables si hay tests que los comparan.

---



## `sequence/`



### `alphabets.py`

Conjuntos IUPAC de `Bio.Data.IUPACData`.


| Nombre                           | Rol                                                                                |
| -------------------------------- | ---------------------------------------------------------------------------------- |
| `NUCLEOTIDE_LETTERS`             | ADN + ARN ambiguos (incluye K, R, Y, … que también aparecen en proteínas)          |
| `UNAMBIGUOUS_NUCLEOTIDE_LETTERS` | A, C, G, T, U — para que un péptido solo con A/C/G no se clasifique como ADN       |
| `PROTEIN_LETTERS`                | Alfabeto proteico extendido                                                        |
| `GAP_LETTERS`                    | `-` y `.` (se ignoran al clasificar)                                               |
| `VALID_RESIDUES`                 | Letras permitidas en un FASTA                                                      |
| `PROTEIN_ONLY_LETTERS`           | Letras que fuerzan clasificación como proteína (E, F, I, L, P, Q, X, O, Z, J, `*`) |


Si cambias los alfabetos, vuelve a revisar `detect_sequence_type` y `validate_seq_record`.

### `validator.py`


| Función               | Rol                                                         |
| --------------------- | ----------------------------------------------------------- |
| `check_fasta_path`    | Existe, es archivo, sufijo en `{.fa,.fasta,.fna,.faa,.fas}` |
| `validate_seq_record` | Secuencia no vacía; cada carácter en `VALID_RESIDUES`       |
| `validate_fasta_file` | Chequeo de ruta + parseo completo vía `read_fasta`          |




### `reader.py`


| Función      | Rol                                                                   |
| ------------ | --------------------------------------------------------------------- |
| `read_fasta` | `Bio.SeqIO.parse`, valida cada registro, pasa secuencias a mayúsculas |




### `classifier.py`


| Función                | Rol                                                                                                                                                                              |
| ---------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `_residues`            | Letras en mayúsculas, gaps quitados                                                                                                                                              |
| `_sequence_id`         | `SeqRecord.id` o `"<no id>"` para logs                                                                                                                                           |
| `detect_sequence_type` | Una secuencia → `"nucleotide"` o `"protein"`: cualquier letra solo-proteína → proteína; si no, todas nucleótido **y** al menos una base no ambigua → nucleótido; si no, proteína |
| `detect_query_type`    | FASTA completo: tipo mayoritario; **empate → proteína**; archivos mixtos no lanzan error                                                                                         |


---



## `blast/`



### `databases.py`


| Función                    | Rol                                                                                                                                          |
| -------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `infer_database_label`     | Mayoría `[Organismo]` en cabeceras NCBI (muestra 40 registros); si no, stem del archivo. Pasa a ser el nombre de bloque de especie en el TSV |
| `_is_fasta_file`           | Archivo con sufijo FASTA                                                                                                                     |
| `_type_from_fasta`         | `detect_query_type` sobre ese archivo                                                                                                        |
| `_iter_fasta_in_directory` | Lista recursiva de FASTA, ordenada                                                                                                           |
| `list_blast_databases`     | Archivo o carpeta → `list[(Path, tipo)]`                                                                                                     |
| `detect_database_type`     | Un solo tipo para toda la fuente; mixto → `MixedDatabaseTypeError`                                                                           |




### `selection.py`


| Función                | Rol                                                                                                                                                          |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `select_blast_program` | `(query_type, database_type, translated)` → `blastn` / `tblastx` / `blastx` / `blastp` / `tblastn`. `translated` solo aplica cuando **ambos** son nucleótido |




### `parser.py`


| Nombre                     | Rol                                                                                   |
| -------------------------- | ------------------------------------------------------------------------------------- |
| `OUTFMT6_STANDARD_COLUMNS` | 12 columnas NCBI                                                                      |
| `OUTFMT6_COLUMNS`          | 15 columnas incluyendo `nident`, `qseq`, `sseq`                                       |
| `BLAST_OUTFMT`             | Cadena `-outfmt` de BLAST usada por `run_blast` / Diamond                             |
| `parse_blast_results`      | Lee tabular; 12 o 15 columnas; archivo vacío → frame vacío con 15 columnas            |
| `fill_missing_hits`        | Asegura una fila por cada `query_id` en `database_name`; sin hit → campos BLAST en NA |




### `ranking.py`


| Función           | Rol                                                                                                                                                                     |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `assign_hit_rank` | Por `(qseqid, database)`: ordena e-value ascendente, bitscore descendente; `hit_rank` 1, 2, … Un `sseqid` vacío se queda en rango 1 para que el TSV pueda mostrar `---` |


Este ranking **debe** ser el mismo en la tabla de resultados y en el BLAST recíproco.

### `reciprocal.py`

RBH: se conserva el par de ida solo si es rango 1 **en ambos sentidos**. Nunca se prueba el rango 2. Se conservan todos los HSP del par aceptado.


| Función                          | Rol                                                                    |
| -------------------------------- | ---------------------------------------------------------------------- |
| `_strip_version`                 | Quita el sufijo `.dígitos`                                             |
| `sequence_ids_match`             | Misma secuencia tras `normalize_accession` y quitar versión            |
| `_has_subject`                   | `sseqid` no vacío                                                      |
| `_rank1_rows` / `_top_subject`   | Sujeto de rango 1 de una query en una base                             |
| `_reverse_query_ids_for_subject` | Ids de query de la búsqueda inversa que coinciden con un sujeto de ida |
| `_drop_rank_helpers`             | Quita `_evalue_sort`, `_score_sort`, `hit_rank`                        |
| `keep_reciprocal_best_hits`      | Filtra HSP de ida usando el rango 1 de la vuelta                       |




### `runner.py`


| Nombre                      | Rol                                                                                                                                                                                                          |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `DEFAULT_MAX_TARGET_SEQS`   | `3`                                                                                                                                                                                                          |
| `DEFAULT_NUM_THREADS`       | `1`                                                                                                                                                                                                          |
| `DEFAULT_BLAST_DIR`         | `"blast"`                                                                                                                                                                                                    |
| `_tabular_filename`         | Nombres únicos `stem.txt`                                                                                                                                                                                    |
| `build_makeblastdb_command` | `makeblastdb -parse_seqids`                                                                                                                                                                                  |
| `build_blast_command`       | BLAST+ con `BLAST_OUTFMT`, evalue, max_target_seqs, threads                                                                                                                                                  |
| `_run_command`              | Comprueba PATH + subprocess; fallo → `BlastExecutionError`                                                                                                                                                   |
| `_run_forward_blast`        | `makeblastdb` temporal, un BLAST por FASTA, parseo, guardado opcional en `blast_dir`, `min_identity`, `fill_missing_hits`                                                                                    |
| `_run_reverse_blast`        | Cada FASTA de base como query contra la query original; programa inverso con tipos intercambiados; `max_target_seqs` = número de secuencias de la query original; escribe `reverse_*.txt` si hay `blast_dir` |
| `run_blast`                 | Driver público de BLAST+. `diamond=True` delega en `run_diamond` e ignora evalue / max_target_seqs / translated / blast_dir                                                                                  |


Los índices viven en un directorio temporal y se borran. El tabular en `blast_dir` va **sin filtrar** (el filtro se aplica después del parseo).

---



## `diamond/`

Solo proteína contra proteína. No se mezclan FASTA y `.dmnd` en una misma corrida.

### `databases.py`


| Función                  | Rol                                                        |
| ------------------------ | ---------------------------------------------------------- |
| `_is_dmnd` / `_is_fasta` | Comprueban sufijo                                          |
| `_iter_in_directory`     | Listas de FASTA vs `.dmnd`                                 |
| `_fasta_entries`         | Solo FASTA proteico; etiqueta vía `infer_database_label`   |
| `list_diamond_databases` | `(path, kind, label)` donde `kind` es `"fasta"` o `"dmnd"` |




### `execute.py`


| Función               | Rol                                                                                    |
| --------------------- | -------------------------------------------------------------------------------------- |
| `run_diamond_command` | Ejecuta `diamond`; binario ausente o salida distinta de cero → `DiamondExecutionError` |




### `runner.py`


| Función                | Rol                                                                                                                            |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `build_makedb_command` | `diamond makedb` (base temporal, se borra)                                                                                     |
| `build_blastp_command` | `diamond blastp` con el mismo outfmt de 15 columnas; **sin** evalue / max-target-seqs (defaults de Diamond); threads default 1 |
| `_out_path`            | Guarda el tabular en `diamond_dir` o en un archivo temporal                                                                    |
| `run_diamond`          | Driver público de Diamond (también lo usa `run_blast`)                                                                         |
| `_run_diamond_blastp`  | Un blastp por base                                                                                                             |
| `_run_reverse_diamond` | Diamond recíproco; **solo bases FASTA** (un `.dmnd` no puede ser la query de vuelta)                                           |


---



## `hmm/`



### `execute.py`


| Función     | Rol                            |
| ----------- | ------------------------------ |
| `run_hmmer` | Ejecuta `hmmpress` / `hmmscan` |




### `press.py`


| Función                       | Rol                                      |
| ----------------------------- | ---------------------------------------- |
| `pressed_index_paths`         | `.h3m` `.h3i` `.h3f` `.h3p` junto al HMM |
| `hmm_database_is_pressed`     | Existen los cuatro índices               |
| `hmm_has_gathering_threshold` | El archivo tiene una línea `GA`          |
| `ensure_pressed_hmm`          | `hmmpress` si faltan índices             |




### `parser.py`


| Función                | Rol                                                                                                           |
| ---------------------- | ------------------------------------------------------------------------------------------------------------- |
| `parse_hmmscan_tblout` | `--tblout` delimitado por espacios; ignora `#`; 19 campos                                                     |
| `best_hmmscan_hits`    | Una fila por query (menor e-value). El **TSV no usa esto**; la tabla concatena todos los dominios que pasaron |




### `runner.py`


| Función                 | Rol                                                                 |
| ----------------------- | ------------------------------------------------------------------- |
| `build_hmmscan_command` | `--tblout --noali --cpu`; `--cut_ga` si hay GA, si no `-E 10`       |
| `run_hmmscan`           | Solo query proteína; escribe `hmm_dir/hmmscan.tbl` (default `hmm/`) |


---



## `kofam/`

Envuelve `exec_annotation` de KofamScan. **No** llama hmmsearch/hmmscan por su cuenta.

### `execute.py`


| Función               | Rol                                    |
| --------------------- | -------------------------------------- |
| `run_exec_annotation` | PATH + subprocess de `exec_annotation` |




### `parser.py`


| Función                  | Rol                                                                                                               |
| ------------------------ | ----------------------------------------------------------------------------------------------------------------- |
| `parse_kofam_detail_tsv` | Conserva solo filas cuyo primer campo es `*`; columnas: mark, gene_name, ko, threshold, score, evalue, definition |




### `runner.py`


| Función                         | Rol                                                                                             |
| ------------------------------- | ----------------------------------------------------------------------------------------------- |
| `build_exec_annotation_command` | `-f detail-tsv --profile --ko-list --cpu --tmp-dir`; sin `-E` / `-T`                            |
| `_resolve_profile`              | Directorio de `.hmm`, o un archivo `.hmm` / `.hal`                                              |
| `_resolve_ko_list`              | Debe ser un archivo                                                                             |
| `run_kofamscan`                 | FASTA proteico; salida `kofam_dir/kofam.tsv`; intermedios en `kofam_dir/tmp` (**no se borran**) |


---



## `metadata/`



### `accessions.py`


| Función               | Rol                                                             |
| --------------------- | --------------------------------------------------------------- |
| `normalize_accession` | Quita prefijos estilo BLAST (`ref|XP_1.1|` → `XP_1.1`)          |
| `load_accessions`     | `accesión<TAB>descriptor`, sin cabecera; añade `accession_norm` |


---



## `processing/`



### `hsp_filter.py`


| Función                 | Rol                                                                                                |
| ----------------------- | -------------------------------------------------------------------------------------------------- |
| `filter_hits_by_pident` | Descarta HSP con `pident < min_identity`; `None` = sin corte; restaura filas vacías por query/base |




### `identity.py`

Identidad de query completa (no el `pident` de BLAST).


| Función                       | Rol                                                                                                                                            |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `_query_step`                 | 3 nt por residuo alineado en blastx/tblastx cuando `qend-qstart+1 == 3 * residuos`; si no, 1                                                   |
| `identical_query_positions`   | Coordenadas 1-based de la query idénticas en un HSP (`qseq`/`sseq`)                                                                            |
| `full_query_identity_percent` | Unión de esas coordenadas en los HSP de **un** par query–sujeto, dividida por la longitud completa de la query × 100. `None` → `---` en el TSV |


Los HSP que se solapan no deben contar dos veces la misma posición.

### `table.py`

Tabla ancha estilo Dataset S2.


| Nombre                           | Rol                                                                   |
| -------------------------------- | --------------------------------------------------------------------- |
| `QUERY_COLUMNS`                  | gene_id, length_nt, cdna_sequence, length_aa, protein_sequence        |
| `HIT_FIELDS`                     | Mapea campos del TSV a columnas BLAST + identidad de query completa   |
| `PFAM_COLUMNS` / `KOFAM_COLUMNS` | Bloques finales opcionales                                            |
| `_safe_name`                     | Prefijo de especie a partir de la etiqueta de la base                 |
| `_empty_query_row`               | Campos de query en NA                                                 |
| `_strip_version`                 | El mismo recorte de versión que en recíproco                          |
| `_record_lookup_keys`            | Une FASTA CDS de NCBI con ids de proteína (`[protein_id=…]`, `_cds_`) |
| `_apply_fasta_to_rows`           | Rellena columnas nt o aa desde un FASTA                               |
| `_query_metadata`                | Combina `--query` / `--cdna` / `--protein`                            |
| `_query_length_info`             | Longitudes para identidad (se ignoran gaps/puntos)                    |
| `_lookup_description`            | Empata `sseqid` con accesiones                                        |
| `_annotation_lookup`             | Agrupa filas Pfam/KOfam por query, ordena por e-value                 |
| `_pfam_*` / `_kofam_*`           | Celdas unidas con `"; "`; vacío → `---`                               |
| `build_result_table`             | Unión pública: BLAST + accesiones + hmm/kofam opcionales              |


`build_result_table(..., hmm_hits=, kofam_hits=)` solo agrega esos bloques si se pasa el DataFrame (no `None`).

---



## `output/`



### `csv_writer.py`


| Función                    | Rol                                                                                                                                                                                                                                            |
| -------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `drop_empty_query_columns` | Quita columnas de query sin datos (p. ej. cDNA en una corrida solo proteína); conserva `gene_id`                                                                                                                                               |
| `_hit_prefix_and_field`    | Parte `species_accession` en prefijo + campo                                                                                                                                                                                                   |
| `format_s2_csv`            | Tres filas de cabecera: sección, especie, nombres de columna. Sección BLAST = `Annotation based on top-BLAST-hit method`; Pfam = `Pfam domains`; KOfam = `KOfamScan`. El nombre de especie solo en la columna de accesión de cada bloque BLAST |
| `write_results_csv`        | Default `results.tsv`; fuerza sufijo `.tsv`                                                                                                                                                                                                    |


Los prefijos de hit BLAST deben excluir `QUERY_COLUMNS`, `PFAM_COLUMNS` y `KOFAM_COLUMNS` para que `kofam_score` no se trate como campo BLAST.

---



## `cli/`

Usa **solo** `from biocol import ...`.


| Módulo        | Rol                                                                |
| ------------- | ------------------------------------------------------------------ |
| `parser.py`   | argparse; ayuda personalizada vía `helptext.py`                    |
| `helptext.py` | Ayuda en inglés con color para top / run / from-blast              |
| `commands.py` | `run_from_fasta`, `run_from_blast`; hmmscan y KofamScan opcionales |
| `main.py`     | Logging, despacho, captura errores de dominio → exit 1             |
| `style.py`    | Colores ANSI; `--no-color` / `NO_COLOR`                            |
| `__main__.py` | `python -m biocol.cli`                                             |




### Auxiliares de `commands.py`


| Función               | Rol                                                                           |
| --------------------- | ----------------------------------------------------------------------------- |
| `_dir_next_to_output` | Default `blast/`, `hmm/`, `kofam/`, `diamond/` junto al TSV                   |
| `_maybe_hmmscan`      | Si hay `--hmm-db`; se exige FASTA proteico                                    |
| `_maybe_kofamscan`    | Hace falta **ambos** `--kofam-profile` y `--ko-list`; se exige FASTA proteico |
| `run_from_fasta`      | BLAST/Diamond, luego Pfam/KOfam opcionales, luego TSV                         |
| `run_from_blast`      | Parsea tabular existente; hmm/kofam necesitan `--protein`                     |


`from-blast` no tiene `--reciprocal` ni `--diamond` (un outfmt 6 de Diamond sí puede pasarse como `--blast`).

---



## Invariantes (no romperlos sin actualizar tests y este archivo)

1. Importaciones públicas: `from biocol import ...` solo para CLI y notebooks.
2. Ranking BLAST: menor e-value, luego mayor bitscore; la misma función para TSV y RBH.
3. Recíproco: rango 1 en ambos sentidos; no se prueba el hit #2.
4. Identity % = `pident` de BLAST. Identity % (full query) = posiciones idénticas únicas de la query / longitud de la query.
5. KofamScan: solo proteína; `detail-tsv`; conservar filas `*`; no borrar tmp.
6. Diamond: proteína vs proteína; no mezclar FASTA y `.dmnd`; el recíproco necesita bases FASTA.
7. hmmscan: solo proteína; en el TSV van todos los dominios que pasaron, no solo el mejor.
8. Celdas de anotación vacías: `---`.
9. Texto de excepciones y CLI: inglés.

---





