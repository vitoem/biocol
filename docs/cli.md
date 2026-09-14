# biocol CLI — module and function reference

This document describes how the command-line interface is organized and what each module and function does. Use it when changing or extending the CLI.  Backend internals are in `[docs/backend.md](backend.md)`.

---

## What the CLI does

Two commands, both write the same TSV:

1. `biocol run` — FASTA query + FASTA (or Diamond) databases → search → optional Pfam/KOfam → TSV.
2. `biocol from-blast` — existing BLAST/Diamond `outfmt 6` + accessions → optional Pfam/KOfam → TSV (no BLAST+ / Diamond search).

Entry point (after `pip install -e .`): `biocol` → `biocol.cli.main:console`. Also `python -m biocol.cli`.

---



## Package layout

```
src/biocol/cli/
  __init__.py     Exports main, console
  __main__.py     python -m biocol.cli
  parser.py       argparse (no backend logic)
  helptext.py     Colored help screens
  commands.py     Dispatch to public API
  main.py         Logging, errors, exit codes
  style.py        ANSI colors
```

Private helpers are prefixed with `_`.

---



## Commands and flags

Global (before the subcommand): `--no-color`, `-h` / `--help`.

`--no-color` is stripped before `parse_args` so argparse does not treat it as unknown on subcommands.

### `biocol run`

Required: `--query`, `--db`, `--accessions`.


| Flag                   | Dest              | Default                    | Role                                                                                      |
| ---------------------- | ----------------- | -------------------------- | ----------------------------------------------------------------------------------------- |
| `--query FASTA`        | `query`           | required                   | Query FASTA (same molecule type throughout). Also the protein input for hmmscan/KofamScan |
| `--db PATH`            | `db`              | required                   | Subject FASTA, folder of FASTA, or with `--diamond` protein FASTA/`.dmnd`                 |
| `--accessions FILE`    | `accessions`      | required                   | `accession<TAB>descriptor` for **subject** ids                                            |
| `--cdna FASTA`         | `cdna`            | `None`                     | CDS FASTA → Length (nt) and cDNA columns                                                  |
| `--protein FASTA`      | `protein_fasta`   | `None`                     | Gene-model protein FASTA → Length(aa) and protein columns                                 |
| `--output TSV`         | `output`          | `results.tsv`              | Result TSV                                                                                |
| `--blast-dir DIR`      | `blast_dir`       | `blast/` next to the TSV   | Unfiltered BLAST tabular (unused with `--diamond`)                                        |
| `--tblastx`            | `tblastx`         | `False`                    | Both nucleotide → `tblastx` instead of `blastn`                                           |
| `--evalue N`           | `evalue`          | `10`                       | BLAST E-value (ignored with `--diamond`)                                                  |
| `--max-target-seqs N`  | `max_target_seqs` | `3`                        | BLAST hits kept per query (ignored with `--diamond`)                                      |
| `--threads N`          | `threads`         | `1`                        | BLAST+, Diamond, hmmscan, and KofamScan CPUs                                              |
| `--min-identity N`     | `min_identity`    | `None`                     | Keep HSP with `pident >= N` (0–100)                                                       |
| `--reciprocal`         | `reciprocal`      | `False`                    | Reverse search; keep rank-1 pairs both ways                                               |
| `--diamond`            | `diamond`         | `False`                    | Diamond blastp instead of BLASTP                                                          |
| `--diamond-dir DIR`    | `diamond_dir`     | `diamond/` next to the TSV | Diamond tabular (only with `--diamond`)                                                   |
| `--hmm-db HMM`         | `hmm_db`          | `None`                     | hmmscan on protein `--query`                                                              |
| `--hmm-dir DIR`        | `hmm_dir`         | `hmm/` next to the TSV     | `hmmscan.tbl`                                                                             |
| `--kofam-profile PATH` | `kofam_profile`   | `None`                     | KOfam dir / `.hmm` / `.hal`; requires `--ko-list`                                         |
| `--ko-list FILE`       | `ko_list`         | `None`                     | KOfam `ko_list`; requires `--kofam-profile`                                               |
| `--kofam-dir DIR`      | `kofam_dir`       | `kofam/` next to the TSV   | `kofam.tsv` + `tmp/`                                                                      |


Passing only one of `--kofam-profile` / `--ko-list` is an error.

### `biocol from-blast`

Required: `--blast`, `--accessions`. No `--query`, `--cdna`, `--reciprocal`, `--diamond`, `--threads`.


| Flag                                            | Dest            | Default       | Role                                                     |
| ----------------------------------------------- | --------------- | ------------- | -------------------------------------------------------- |
| `--blast FILE`                                  | `blast`         | required      | `outfmt 6` (12 or 15 columns; Diamond tabular OK)        |
| `--accessions FILE`                             | `accessions`    | required      | Subject descriptors; stem becomes the species block name |
| `--output TSV`                                  | `output`        | `results.tsv` | Result TSV                                               |
| `--min-identity N`                              | `min_identity`  | `None`        | Same `pident` cutoff as `run`                            |
| `--protein FASTA`                               | `protein_fasta` | `None`        | Protein FASTA for hmmscan and/or KofamScan               |
| `--hmm-db` / `--hmm-dir`                        |                 |               | Same as `run`; need `--protein`                          |
| `--kofam-profile` / `--ko-list` / `--kofam-dir` |                 |               | Same as `run`; need `--protein`                          |


If the tabular has no `database` column, the accessions file stem is used.

---



## Exit codes


| Code | Meaning                                                                                                   |
| ---- | --------------------------------------------------------------------------------------------------------- |
| `0`  | Success; TSV path printed on stdout; `Done.` on stderr                                                    |
| `1`  | Domain / I/O error (`FileNotFoundError`, FASTA/BLAST/Diamond/Hmm/Kofam/Metadata, `OSError`, `ValueError`) |
| `2`  | argparse error or no command (`run` / `from-blast` missing)                                               |


---



## `parser.py`

No search or table logic — flags only.


| Name                      | Role                                                                                          |
| ------------------------- | --------------------------------------------------------------------------------------------- |
| `_min_identity`           | Parse float 0–100; else `ArgumentTypeError`                                                   |
| `_HelpParser`             | argparse subclass; `--help` uses `help_renderer`; errors print help + `error: ...` and exit 2 |
| `_HelpParser.format_help` | Custom screen if `help_renderer` is set                                                       |
| `_HelpParser.error`       | Help + red error on stderr                                                                    |
| `build_parser`            | Root parser + `run` + `from-blast`                                                            |


---



## `helptext.py`

Colored English screens. Box drawing falls back to ASCII if `BIOCOL_ASCII` is set or the encoding cannot print box characters.


| Name                                          | Role                                   |
| --------------------------------------------- | -------------------------------------- |
| `WIDTH`                                       | Inner help width (76)                  |
| `_box_ok`                                     | Whether Unicode box drawing is allowed |
| `_center`                                     | Center a string in a width             |
| `banner`                                      | `BIOCOL` + subtitle box                |
| `heading` / `cmd` / `opt` / `example` / `dim` | Paint helpers                          |
| `render_top_help`                             | `biocol --help`                        |
| `render_run_help`                             | `biocol run --help`                    |
| `render_from_blast_help`                      | `biocol from-blast --help`             |


When you add a flag, add it here as well as in `parser.py`.

---



## `commands.py`

Calls only `run_blast`, `parse_blast_results`, `filter_hits_by_pident`, `run_hmmscan`, `run_kofamscan`, `build_result_table`, `write_results_csv`.


| Function              | Role                                                                                                            |
| --------------------- | --------------------------------------------------------------------------------------------------------------- |
| `_dir_next_to_output` | Explicit dir, or `<parent of output>/<default_name>` (`results.tsv` → current directory)                        |
| `_maybe_hmmscan`      | Skip if no `--hmm-db`; else protein FASTA required; `run_hmmscan`                                               |
| `_maybe_kofamscan`    | Skip if both profile and ko-list absent; error if only one; protein FASTA required; `run_kofamscan`             |
| `run_from_fasta`      | `run_blast` (Diamond path sets `blast_dir=None` and fills `diamond_dir`); hmmscan/Kofam on `--query`; write TSV |
| `run_from_blast`      | Parse tabular, optional `min_identity`, hmmscan/Kofam on `--protein`, write TSV                                 |


`--threads` is not defined on `from-blast`; hmmscan/Kofam then use `getattr(args, "threads", 1)` → **1**.

---



## `main.py`


| Name                     | Role                                                          |
| ------------------------ | ------------------------------------------------------------- |
| `_COMMANDS`              | `"run"` → `run_from_fasta`, `"from-blast"` → `run_from_blast` |
| `_CliFormatter`          | `biocol <message>` on stderr; errors red, warnings yellow     |
| `_configure_cli_logging` | `logging.getLogger("biocol")` at INFO, no propagate           |
| `main`                   | Color setup, parse, dispatch, catch errors, print TSV path    |
| `console`                | `raise SystemExit(main())` for the `biocol` script            |


`main(argv)` is what tests call.

---



## `style.py`

No extra dependencies. Honors `--no-color`, `NO_COLOR`, `TERM=dumb`, and non-TTY streams.


| Name                                                                | Role                                     |
| ------------------------------------------------------------------- | ---------------------------------------- |
| `RESET`, `BOLD`, `DIM`, `RED`, `GREEN`, `YELLOW`, `CYAN`, `MAGENTA` | ANSI codes                               |
| `reset_color` / `disable_color`                                     | Session color on/off                     |
| `use_color`                                                         | Whether to emit ANSI                     |
| `paint`                                                             | Wrap text in codes if color is on        |
| `enable_windows_vt`                                                 | Enable VT processing on Windows consoles |


---



## `__init__.py` / `__main__.py`


| Name          | Role                                   |
| ------------- | -------------------------------------- |
| `__init__.py` | `from biocol.cli import main, console` |
| `__main__.py` | `console()` when run as a module       |


---

# Traducción al español



# CLI de biocol — referencia de módulos y funciones

Este documento describe cómo está organizada la interfaz de línea de comandos y qué hace cada módulo y cada función. Úsalo al cambiar o ampliar la CLI.

---

## Qué hace la CLI

Dos comandos; ambos escriben el mismo TSV estilo Dataset S2:

1. `biocol run` — FASTA de query + bases FASTA (o Diamond) → búsqueda → Pfam/KOfam opcionales → TSV.
2. `biocol from-blast` — `outfmt 6` BLAST/Diamond ya existente + accesiones → Pfam/KOfam opcionales → TSV (sin BLAST+ / Diamond).

Punto de entrada (tras `pip install -e .`): `biocol` → `biocol.cli.main:console`. También `python -m biocol.cli`.

---



## Estructura del paquete

```
src/biocol/cli/
  __init__.py     Exporta main, console
  __main__.py     python -m biocol.cli
  parser.py       argparse (sin lógica de backend)
  helptext.py     Pantallas de ayuda con color
  commands.py     Despacho a la API pública
  main.py         Logging, errores, códigos de salida
  style.py        Colores ANSI
```

Las funciones auxiliares privadas llevan prefijo `_`.

---



## Comandos y flags

Globales (antes del subcomando): `--no-color`, `-h` / `--help`.

`--no-color` se quita antes de `parse_args` para que argparse no lo trate como desconocido en los subcomandos.

### `biocol run`

Obligatorios: `--query`, `--db`, `--accessions`.


| Flag                   | Dest              | Default                 | Rol                                                                                           |
| ---------------------- | ----------------- | ----------------------- | --------------------------------------------------------------------------------------------- |
| `--query FASTA`        | `query`           | obligatorio             | FASTA de query (mismo tipo en todo el archivo). También es la proteína para hmmscan/KofamScan |
| `--db PATH`            | `db`              | obligatorio             | FASTA sujeto, carpeta de FASTA, o con `--diamond` FASTA/`.dmnd` proteico                      |
| `--accessions FILE`    | `accessions`      | obligatorio             | `accesión<TAB>descriptor` de los ids del **sujeto**                                           |
| `--cdna FASTA`         | `cdna`            | `None`                  | FASTA CDS → columnas Length (nt) y cDNA                                                       |
| `--protein FASTA`      | `protein_fasta`   | `None`                  | FASTA de proteína (modelos de gen) → Length(aa) y proteína                                    |
| `--output TSV`         | `output`          | `results.tsv`           | TSV de resultado                                                                              |
| `--blast-dir DIR`      | `blast_dir`       | `blast/` junto al TSV   | Tabular BLAST sin filtrar (no se usa con `--diamond`)                                         |
| `--tblastx`            | `tblastx`         | `False`                 | Ambos nucleótido → `tblastx` en lugar de `blastn`                                             |
| `--evalue N`           | `evalue`          | `10`                    | E-value de BLAST (se ignora con `--diamond`)                                                  |
| `--max-target-seqs N`  | `max_target_seqs` | `3`                     | Hits de BLAST por query (se ignora con `--diamond`)                                           |
| `--threads N`          | `threads`         | `1`                     | CPUs de BLAST+, Diamond, hmmscan y KofamScan                                                  |
| `--min-identity N`     | `min_identity`    | `None`                  | Conservar HSP con `pident >= N` (0–100)                                                       |
| `--reciprocal`         | `reciprocal`      | `False`                 | Búsqueda inversa; solo pares rango 1 en ambos sentidos                                        |
| `--diamond`            | `diamond`         | `False`                 | Diamond blastp en lugar de BLASTP                                                             |
| `--diamond-dir DIR`    | `diamond_dir`     | `diamond/` junto al TSV | Tabular Diamond (solo con `--diamond`)                                                        |
| `--hmm-db HMM`         | `hmm_db`          | `None`                  | hmmscan sobre `--query` proteico                                                              |
| `--hmm-dir DIR`        | `hmm_dir`         | `hmm/` junto al TSV     | `hmmscan.tbl`                                                                                 |
| `--kofam-profile PATH` | `kofam_profile`   | `None`                  | Perfiles KOfam (dir / `.hmm` / `.hal`); exige `--ko-list`                                     |
| `--ko-list FILE`       | `ko_list`         | `None`                  | `ko_list` de KOfam; exige `--kofam-profile`                                                   |
| `--kofam-dir DIR`      | `kofam_dir`       | `kofam/` junto al TSV   | `kofam.tsv` + `tmp/`                                                                          |


Pasar solo uno de `--kofam-profile` / `--ko-list` es un error.

### `biocol from-blast`

Obligatorios: `--blast`, `--accessions`. No hay `--query`, `--cdna`, `--reciprocal`, `--diamond`, `--threads`.


| Flag                                            | Dest            | Default       | Rol                                                                |
| ----------------------------------------------- | --------------- | ------------- | ------------------------------------------------------------------ |
| `--blast FILE`                                  | `blast`         | obligatorio   | `outfmt 6` (12 o 15 columnas; vale tabular Diamond)                |
| `--accessions FILE`                             | `accessions`    | obligatorio   | Descriptores del sujeto; el stem es el nombre de bloque de especie |
| `--output TSV`                                  | `output`        | `results.tsv` | TSV de resultado                                                   |
| `--min-identity N`                              | `min_identity`  | `None`        | El mismo corte de `pident` que en `run`                            |
| `--protein FASTA`                               | `protein_fasta` | `None`        | FASTA proteico para hmmscan y/o KofamScan                          |
| `--hmm-db` / `--hmm-dir`                        |                 |               | Igual que en `run`; hace falta `--protein`                         |
| `--kofam-profile` / `--ko-list` / `--kofam-dir` |                 |               | Igual que en `run`; hace falta `--protein`                         |


Si el tabular no trae columna `database`, se usa el stem del archivo de accesiones.

---



## Códigos de salida


| Código | Significado                                                                                                   |
| ------ | ------------------------------------------------------------------------------------------------------------- |
| `0`    | Éxito; ruta del TSV en stdout; `Done.` en stderr                                                              |
| `1`    | Error de dominio / E/S (`FileNotFoundError`, FASTA/BLAST/Diamond/Hmm/Kofam/Metadata, `OSError`, `ValueError`) |
| `2`    | Error de argparse o falta el comando (`run` / `from-blast`)                                                   |


---



## `parser.py`

Sin lógica de búsqueda ni de tabla: solo flags.


| Nombre                    | Rol                                                                                                            |
| ------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `_min_identity`           | Parsea float 0–100; si no, `ArgumentTypeError`                                                                 |
| `_HelpParser`             | Subclase de argparse; `--help` usa `help_renderer`; los errores imprimen la ayuda + `error: ...` y salen con 2 |
| `_HelpParser.format_help` | Pantalla propia si hay `help_renderer`                                                                         |
| `_HelpParser.error`       | Ayuda + error en rojo en stderr                                                                                |
| `build_parser`            | Parser raíz + `run` + `from-blast`                                                                             |


---



## `helptext.py`

Pantallas en inglés con color. El recuadro pasa a ASCII si está `BIOCOL_ASCII` o si la codificación no puede pintar caracteres de caja.


| Nombre                                        | Rol                                           |
| --------------------------------------------- | --------------------------------------------- |
| `WIDTH`                                       | Ancho interior de la ayuda (76)               |
| `_box_ok`                                     | Si se permiten caracteres Unicode de recuadro |
| `_center`                                     | Centra un texto en un ancho                   |
| `banner`                                      | Caja `BIOCOL` + subtítulo                     |
| `heading` / `cmd` / `opt` / `example` / `dim` | Ayudas de color                               |
| `render_top_help`                             | `biocol --help`                               |
| `render_run_help`                             | `biocol run --help`                           |
| `render_from_blast_help`                      | `biocol from-blast --help`                    |


Si agregas una flag, agrégala aquí y en `parser.py`.

---



## `commands.py`

Solo llama a `run_blast`, `parse_blast_results`, `filter_hits_by_pident`, `run_hmmscan`, `run_kofamscan`, `build_result_table`, `write_results_csv`.


| Función               | Rol                                                                                                                 |
| --------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `_dir_next_to_output` | Directorio explícito, o `<padre del output>/<nombre_default>` (`results.tsv` → directorio actual)                   |
| `_maybe_hmmscan`      | No hace nada si no hay `--hmm-db`; si no, exige FASTA proteico; `run_hmmscan`                                       |
| `_maybe_kofamscan`    | No hace nada si faltan perfil y ko-list; error si solo hay uno; exige FASTA proteico; `run_kofamscan`               |
| `run_from_fasta`      | `run_blast` (con Diamond: `blast_dir=None` y se llena `diamond_dir`); hmmscan/Kofam sobre `--query`; escribe el TSV |
| `run_from_blast`      | Parsea tabular, `min_identity` opcional, hmmscan/Kofam sobre `--protein`, escribe el TSV                            |


`--threads` no existe en `from-blast`; hmmscan/Kofam usan `getattr(args, "threads", 1)` → **1**.

---



## `main.py`


| Nombre                   | Rol                                                                  |
| ------------------------ | -------------------------------------------------------------------- |
| `_COMMANDS`              | `"run"` → `run_from_fasta`, `"from-blast"` → `run_from_blast`        |
| `_CliFormatter`          | `biocol <mensaje>` en stderr; errores en rojo, avisos en amarillo    |
| `_configure_cli_logging` | `logging.getLogger("biocol")` en INFO, sin propagate                 |
| `main`                   | Color, parseo, despacho, captura de errores, imprime la ruta del TSV |
| `console`                | `raise SystemExit(main())` para el script `biocol`                   |


Los tests llaman `main(argv)`.

---



## `style.py`

Sin dependencias extra. Respeta `--no-color`, `NO_COLOR`, `TERM=dumb` y flujos que no son TTY.


| Nombre                                                              | Rol                                                  |
| ------------------------------------------------------------------- | ---------------------------------------------------- |
| `RESET`, `BOLD`, `DIM`, `RED`, `GREEN`, `YELLOW`, `CYAN`, `MAGENTA` | Códigos ANSI                                         |
| `reset_color` / `disable_color`                                     | Color de la sesión encendido/apagado                 |
| `use_color`                                                         | Si se emiten ANSI                                    |
| `paint`                                                             | Envuelve el texto en códigos si el color está activo |
| `enable_windows_vt`                                                 | Activa procesamiento VT en consolas de Windows       |


---



## `__init__.py` / `__main__.py`


| Nombre        | Rol                                    |
| ------------- | -------------------------------------- |
| `__init__.py` | `from biocol.cli import main, console` |
| `__main__.py` | `console()` al ejecutar el módulo      |


---



