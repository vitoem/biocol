from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from biocol.processing.table import KOFAM_COLUMNS, PFAM_COLUMNS, QUERY_COLUMNS

logger = logging.getLogger(__name__)

DEFAULT_OUTPUT = "results.tsv"

QUERY_LABELS = {
    "gene_id": "Gene ID",
    "length_nt": "Length (nt)",
    "cdna_sequence": "cDNA Sequences (nt)",
    "length_aa": "Length(aa)",
    "protein_sequence": "Protein Sequences (aa)",
}

HIT_LABELS = {
    "accession": "Accesion No.",
    "description": "Description",
    "identity_pct": "Identity %",
    "identity_full_query": "Identity % (full query)",
    "alignment_length": "Alignment length",
    "evalue": "e-value",
    "score": "Score",
}

BLAST_SECTION = "Annotation based on top-BLAST-hit method"
PFAM_SECTION = "Pfam domains"
KOFAM_SECTION = "KOfamScan"

PFAM_LABELS = {
    "pfam_n_domains": "# of Pfam domain identified",
    "pfam_evalue": "e-value",
    "pfam_score": "score",
    "pfam_accession": "Accesion",
    "pfam_name": "Name",
    "pfam_description": "Description of target",
}

KOFAM_LABELS = {
    "kofam_ko": "KO",
    "kofam_score": "score",
    "kofam_evalue": "e-value",
    "kofam_definition": "KO definition",
}


def _column_has_values(series: pd.Series) -> bool:
    if series.isna().all():
        return False
    text = series.astype(str).str.strip()
    return text.replace({"nan": "", "<NA>": "", "None": ""}).ne("").any()


def drop_empty_query_columns(table: pd.DataFrame) -> pd.DataFrame:
    """Drop query columns that have no data (e.g. cDNA when the query is protein)."""
    keep = []
    for column in table.columns:
        if column in QUERY_COLUMNS and column != "gene_id" and not _column_has_values(table[column]):
            continue
        keep.append(column)
    return table.loc[:, keep]


def _hit_prefix_and_field(column: str) -> tuple[str, str] | None:
    for field in sorted(HIT_LABELS, key=len, reverse=True):
        suffix = f"_{field}"
        if column.endswith(suffix):
            return column[: -len(suffix)], field
    return None


def format_s2_csv(table: pd.DataFrame) -> pd.DataFrame:
    """Table with a 3-row header and per-species blocks, Dataset S2 style."""
    frame = drop_empty_query_columns(table)
    query_cols = [column for column in QUERY_COLUMNS if column in frame.columns]
    pfam_cols = [column for column in PFAM_COLUMNS if column in frame.columns]
    kofam_cols = [column for column in KOFAM_COLUMNS if column in frame.columns]
    hit_cols = [
        column
        for column in frame.columns
        if column not in QUERY_COLUMNS
        and column not in PFAM_COLUMNS
        and column not in KOFAM_COLUMNS
    ]

    prefixes: list[str] = []
    for column in hit_cols:
        parsed = _hit_prefix_and_field(column)
        if parsed and parsed[0] not in prefixes:
            prefixes.append(parsed[0])

    ordered = list(query_cols)
    for prefix in prefixes:
        for field in HIT_LABELS:
            name = f"{prefix}_{field}"
            if name in frame.columns:
                ordered.append(name)
    ordered.extend(pfam_cols)
    ordered.extend(kofam_cols)
    frame = frame.loc[:, ordered]

    section_row: list[str] = []
    species_row: list[str] = []
    label_row: list[str] = []
    for column in ordered:
        if column in QUERY_LABELS:
            section_row.append("")
            species_row.append("")
            label_row.append(QUERY_LABELS[column])
            continue
        if column in PFAM_LABELS:
            is_block_start = column == pfam_cols[0] if pfam_cols else False
            section_row.append(PFAM_SECTION if is_block_start else "")
            species_row.append("")
            label_row.append(PFAM_LABELS[column])
            continue
        if column in KOFAM_LABELS:
            is_block_start = column == kofam_cols[0] if kofam_cols else False
            section_row.append(KOFAM_SECTION if is_block_start else "")
            species_row.append("")
            label_row.append(KOFAM_LABELS[column])
            continue
        parsed = _hit_prefix_and_field(column)
        prefix, field = parsed if parsed else (column, "")
        is_block_start = field == "accession"
        section_row.append(BLAST_SECTION if is_block_start else "")
        species_row.append(prefix.replace("_", " ") if is_block_start else "")
        label_row.append(HIT_LABELS.get(field, column))

    header = pd.DataFrame([section_row, species_row, label_row], columns=ordered)
    body = frame.copy()
    body.columns = ordered
    return pd.concat([header, body], ignore_index=True)


def write_results_csv(
    table: pd.DataFrame,
    output: str | Path | None = None,
) -> Path:
    """Write the Dataset S2-style TSV. Defaults to ``results.tsv``."""
    path = Path(output) if output else Path(DEFAULT_OUTPUT)
    if path.suffix.lower() != ".tsv":
        path = path.with_suffix(".tsv")
    format_s2_csv(table).to_csv(path, index=False, header=False, sep="\t")
    logger.info("Wrote TSV: %s", path)
    return path
