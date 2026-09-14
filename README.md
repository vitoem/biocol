# biocol

Backend de la herramienta BLAST de INECOL. Toma una query FASTA (o un BLAST tabular ya existente), la compara contra bases FASTA locales y escribe un **TSV plano** de anotación.


## Requisitos

- Python 3.10+
- BLAST+ en el PATH (`conda activate inecol` en Ubuntu/WSL)
- Diamond en el PATH solo si usas `diamond=True`
- KofamScan (`exec_annotation`), HMMER y GNU Parallel en el PATH si usas `run_kofamscan`
- Dependencias: Biopython y pandas

```bash
cd biocol
pip install -e ".[dev]"
```

`-e` instala el paquete en modo editable: los cambios del backend se ven sin reinstalar.

## Entradas

| Entrada | Formato |
|---------|---------|
| Query | FASTA / multifasta: `.fa`, `.fasta`, `.fna`, `.faa`, `.fas`. Todo el archivo debe ser del mismo tipo (ADN, ARN o proteína). `U` cuenta como nucleótido. |
| Bases | Un FASTA o una **carpeta** (incluye subcarpetas). Mismas extensiones. Todas las bases del mismo tipo. Un BLAST por archivo FASTA. |
| Accesiones | Texto `accession<TAB>descriptor`, sin encabezado. |
| BLAST tabular (camino 2) | `outfmt 6`: 12 columnas NCBI, o 15 si incluye `nident`, `qseq` y `sseq` (lo que escribe `run_blast`). |

Parámetros de BLAST (modificables): `evalue` default **10**, `max_target_seqs` default **3**, `threads` default **1**. `min_identity` es opcional: si no se pasa, **no hay corte** por % de identidad (`pident`). El TSV muestra solo el **mejor hit** por query y especie.

## Cómo se elige el programa BLAST

`run_blast` detecta el tipo de la query y de las bases y llama a `select_blast_program`. Si **ambos** son nucleótido, el default es **blastn**. `translated=True` (flag CLI `--tblastx`) elige **tblastx**. En las demás combinaciones `translated` se ignora.

| Query | Base de datos | `translated` | Programa | Qué compara |
|-------|---------------|--------------|----------|-------------|
| Nucleótido | Nucleótido | no se pasa / `False` (default) | `blastn` | Nucleótido contra nucleótido |
| Nucleótido | Nucleótido | `True` (explícito) | `tblastx` | Query y base traducidas en seis marcos (proteína) |
| Nucleótido | Proteína | se ignora | `blastx` | Query nucleotídica traducida contra proteínas |
| Proteína | Proteína | se ignora | `blastp` | Proteína contra proteína |
| Proteína | Nucleótido | se ignora | `tblastn` | Proteína contra traducciones de la base nucleotídica |

## Uso

Hay dos caminos al mismo TSV.

### Camino 1 — FASTA + bases

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
# opcional: hmmscan (query proteína + Pfam-A.hmm u otra base HMMER3)
# hmm_hits = run_hmmscan("query.faa", "Pfam-A.hmm", hmm_dir="hmm")
# table = build_result_table(hits, "accessions.txt", query_fasta="query.faa", hmm_hits=hmm_hits)
# kofam_hits = run_kofamscan("query.faa", "profiles/eukaryote.hal", "ko_list", kofam_dir="kofam")
# table = build_result_table(hits, "accessions.txt", query_fasta="query.faa", kofam_hits=kofam_hits)
write_results_csv(table, "results.tsv")  # si se omite, usa results.tsv
```

`run_blast` crea bases temporales con `makeblastdb` (se borran al terminar), lanza un BLAST por FASTA y parsea `outfmt 6`. Si pasas `blast_dir`, deja ahí un `.txt` por base **sin filtrar**. `min_identity` descarta HSP con `pident` menor al umbral (igual en blastn y blastp). Si una query no tiene hit en una base, queda una fila vacía.

`reciprocal=True` (default `False`) hace además el contraste inverso: cada FASTA de `--db` pasa a ser query y el FASTA original es la base. El programa se elige con los tipos intercambiados (`blastx` ↔ `tblastn`; `tblastx` ↔ `tblastx`). En la vuelta, `max_target_seqs` es el número de secuencias del query original. El mismo `min_identity` se aplica en ida y vuelta. Solo se conservan los HSP del par que es **top hit en ambos sentidos** (menor e-value, luego mayor bitscore; no se prueba el hit #2). Si no hay recíproco, la celda queda `---`. Con `blast_dir` también se guarda `reverse_<stem>.txt`. `from-blast` no aplica este modo.

`diamond=True` (default `False`, solo backend) sustituye **blastp**: query y base tienen que ser proteína. Dos entradas de base, un solo tipo por corrida: archivo/carpeta **`.dmnd`**, o FASTA proteico (entonces `diamond makedb` en un temporal que se borra). `evalue` y `max-target-seqs` son los default de Diamond; `threads` default **1**. El tabular se guarda en `diamond_dir` (no en `blast/`). Si `reciprocal=True` y las bases son FASTA proteína, la vuelta también es `diamond blastp`. Un `.dmnd` no sirve para la vuelta (no hay FASTA query de la especie) → error. Query o base nucleótido → `diamond requires a protein query and a protein database`. Un tabular Diamond con el mismo `outfmt 6` se puede leer con `from-blast` / `parse_blast_results`.

`run_kofamscan` llama a `exec_annotation` sobre un FASTA **proteico**. `profile` es una carpeta de `.hmm`, un archivo `.hmm`, o un `.hal` (p. ej. `profiles/eukaryote.hal`); `ko_list` es el archivo `ko_list` de KOfam. Formato `detail-tsv`; se conservan solo las filas con `*`. El TSV de KofamScan queda en `kofam_dir` (`kofam/kofam.tsv`) y los temporales de hmmsearch en `kofam_dir/tmp` (no se borran). `num_threads` default **1**. Query nucleótido → `kofamscan requires a protein sequence`. Pásalo a `build_result_table(..., kofam_hits=...)`. En CLI: `--kofam-profile` y `--ko-list` juntos (uno solo → error).

### Camino 2 — BLAST tabular ya existente

Sin FASTA de query: se rellenan `qseqid` + hits + descriptores; las columnas de secuencia quedan vacías.

```python
from biocol import parse_blast_results, filter_hits_by_pident, build_result_table, write_results_csv

hits = parse_blast_results("blast_outfmt6.txt")
hits = filter_hits_by_pident(hits, 80)  # opcional
table = build_result_table(hits, "accessions.txt", query_fasta=None)
write_results_csv(table)
```

Si el tabular no trae columna `database`, se usa el nombre `hit`.

### CLI

After `pip install -e ".[dev]"`:

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

`--output` is optional (default: `results.tsv`). `biocol run` also writes BLAST tabular files to `--blast-dir` (default: `blast/` next to the TSV). Optional `--hmm-db Pfam-A.hmm` runs hmmscan (protein query; `hmmpress` if needed) and writes `hmm/hmmscan.tbl`. Optional `--reciprocal` runs reverse BLAST and keeps only rank-1 pairs in both directions. Optional `--diamond` runs Diamond blastp (protein vs protein FASTA or `.dmnd`; tabular in `--diamond-dir`, default `diamond/` next to the TSV). Optional `--kofam-profile` + `--ko-list` run KofamScan `exec_annotation` (protein query; detail-tsv in `--kofam-dir`, default `kofam/` next to the TSV, tmp kept). Help text and errors are in English.

`from-blast` does not take a query FASTA and does not support `--reciprocal` or `--diamond` (a Diamond `outfmt 6` file can still be passed as `--blast`). KofamScan in `from-blast` needs `--protein`. The species name in the TSV header is the accessions file stem (`Benincasa_hispida_gd.txt` → `Benincasa hispida gd`).

## TSV de salida

TSV con **tres filas de cabecera**, como Dataset S2 (Pfam y KOfam opcionales; sin GO). No hay celdas combinadas: el nombre de sección y el de la especie se repiten o quedan en la primera columna de cada bloque.

1. Sección: vacío en query; `Annotation based on top-BLAST-hit method` en cada bloque BLAST; `Pfam domains` si hubo hmmscan; `KOfamScan` si hubo KofamScan.
2. Especie: `stem` del FASTA de base (p. ej. `protein`, `amborella`). Vacío en Pfam y KOfam.
3. Nombres de columna: `Gene ID`, `Length (nt)`, `cDNA Sequences (nt)`, `Length(aa)`, `Protein Sequences (aa)`, y por especie `Accesion No.`, `Description`, `Identity %`, `Identity % (full query)`, `Alignment length`, `e-value`, `Score`. Si hay hmmscan: `# of Pfam domain identified`, `e-value`, `score`, `Accesion`, `Name`, `Description of target` (todos los dominios que pasaron el corte, separados por `; `; sin dominio → `---`). Si hay KofamScan: `KO`, `score`, `e-value`, `KO definition` (solo asignaciones sobre el umbral del KO, con `*`; varios KO de la misma proteína unidos con `; `; sin KO → `---`).

`Identity %` es el `pident` de BLAST (respecto al alineamiento). `Identity % (full query)` usa las secuencias alineadas `qseq`/`sseq`: posiciones idénticas **únicas** de la query ÷ longitud completa de la query × 100 (aa si la query es proteína, nt si es nucleótido). Varios HSP del mismo sujeto se unen sin contar dos veces un solape. En blastx/tblastx cada aminoácido idéntico cubre 3 nt de la query. `--cdna` no entra en el denominador. Sin hit, sin `qseq`/`sseq` (tabular de 12 columnas) o longitud 0: `---`.

No se incluyen Length (aa) ni secuencia del hit. Sin hit o sin descriptor: `---`.

Las columnas de query que vayan vacías **no se escriben** (query proteína → sin cDNA; query nucleótido → sin proteína). Si hay modelos de gen, se pueden pasar ambos FASTA (`--cdna` y `--protein`) y el primer bloque queda completo.

Una fila por query (solo el mejor hit por especie). En Excel, importar el TSV y opcionalmente combinar celdas de las dos primeras filas.

## API pública

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

`detect_sequence_type()` acepta `str`, `Bio.Seq.Seq` o `SeqRecord`. `detect_query_type()` clasifica un FASTA completo por el tipo mayoritario (en empate, proteína).

Errores: `FastaError`, `EmptyFastaError`, `InvalidFastaError`, `MixedSequenceTypeError`, `BlastError`, `DatabaseError`, `MixedDatabaseTypeError`, `BlastExecutionError`, `DiamondError`, `DiamondExecutionError`, `HmmError`, `HmmExecutionError`, `KofamError`, `KofamExecutionError`, `MetadataError`.

## Pruebas (Dana)

```bash
pytest -q
pytest tests/test_detect_sequence_type.py -q
```

Al correr Pytest se muestran logs INFO. Fixtures en `tests/fixtures/` (FASTA, BLAST `outfmt 6`, accesiones). Importar solo desde `biocol`.

`run_blast` contra BLAST+ real solo en el entorno `conda` `inecol`.

| Entrada | Resultado esperado |
|---------|--------------------|
| FASTA de ADN | `nucleotide` |
| FASTA de ARN | `nucleotide` |
| FASTA de proteína | `protein` |
| Multifasta del mismo tipo | tipo único |
| Multifasta mixto ADN + proteína | tipo mayoritario (empate → proteína) |
| Ruta inexistente | `FileNotFoundError` |
| Extensión no FASTA (p. ej. `.txt`) | `InvalidFastaError` |
| Archivo vacío / sin secuencias | `EmptyFastaError` |
| Contenido que no es FASTA | `InvalidFastaError` |
| Secuencia vacía o caracteres inválidos | `InvalidFastaError` |

## Equipo y estructura

| Persona  | Rol    | Trabaja sobre |
|----------|--------|----------------|
| Alondra  | lógica | `src/biocol/`  |
| Emiliano | CLI    | `src/biocol/cli/` (solo API pública de `biocol`) |
| Dana     | pruebas | `tests/` |

```
src/biocol/           backend
  sequence/           lectura, validación y tipo de query
  blast/              tipo de base, selección, ejecución, parseo tabular y BLAST recíproco
  diamond/            makedb / blastp opcional (proteína vs proteína)
  hmm/                hmmpress, hmmscan y parseo tblout
  kofam/              exec_annotation (KofamScan) y parseo detail-tsv
  metadata/           accesiones y descriptores
  processing/         tabla ancha, identidad (full query) y filtro de HSP
  output/             escritura del TSV
  cli/                CLI (argparse + comandos run / from-blast)
tests/                pruebas (pytest)
tests/fixtures/       FASTA, BLAST outfmt 6 y accesiones de ejemplo
```

La lectura de FASTA usa `Bio.SeqIO`. Los alfabetos nucleótido/proteína salen de `Bio.Data.IUPACData`.

Developer reference: [`docs/backend.md`](docs/backend.md) (backend), [`docs/cli.md`](docs/cli.md) (CLI).
