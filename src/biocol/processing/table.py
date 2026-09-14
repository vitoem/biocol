from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

from biocol.blast.parser import OUTFMT6_COLUMNS
from biocol.blast.ranking import assign_hit_rank
from biocol.metadata.accessions import load_accessions, normalize_accession
from biocol.processing.identity import full_query_identity_percent
from biocol.sequence.classifier import detect_query_type
from biocol.sequence.reader import read_fasta

logger = logging.getLogger(__name__)

QUERY_COLUMNS = [
    "gene_id",
    "length_nt",
    "cdna_sequence",
    "length_aa",
    "protein_sequence",
]

HIT_FIELDS = [
    ("accession", "sseqid"),
    ("description", "description"),
    ("identity_pct", "pident"),
    ("identity_full_query", None),
    ("alignment_length", "length"),
    ("evalue", "evalue"),
    ("score", "bitscore"),
]

PFAM_COLUMNS = [
    "pfam_n_domains",
    "pfam_evalue",
    "pfam_score",
    "pfam_accession",
    "pfam_name",
    "pfam_description",
]

KOFAM_COLUMNS = [
    "kofam_ko",
    "kofam_score",
    "kofam_evalue",
    "kofam_definition",
]


def _safe_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
    return cleaned or "db"


def _empty_query_row(query_id: str) -> dict:
    return {
        "gene_id": query_id,
        "length_nt": pd.NA,
        "cdna_sequence": pd.NA,
        "length_aa": pd.NA,
        "protein_sequence": pd.NA,
    }


def _strip_version(accession: str) -> str:
    if "." in accession:
        head, tail = accession.rsplit(".", 1)
        if tail.isdigit():
            return head
    return accession


_PROTEIN_ID_TAG = re.compile(r"\[protein_id=([^\]]+)\]", re.IGNORECASE)
_CDS_ACCESSION = re.compile(r"_cds_([A-Z]{1,3}_?\d+(?:\.\d+)?)", re.IGNORECASE)


def _record_lookup_keys(record) -> list[str]:
    """IDs that can join an NCBI CDS FASTA with a .faa file."""
    keys = [record.id, normalize_accession(record.id)]
    description = record.description or ""
    tagged = _PROTEIN_ID_TAG.search(description)
    if tagged:
        keys.append(tagged.group(1).strip())
    cds = _CDS_ACCESSION.search(record.id)
    if cds:
        keys.append(cds.group(1))
    expanded: list[str] = []
    for key in keys:
        if not key:
            continue
        expanded.append(key)
        stripped = _strip_version(key)
        if stripped != key:
            expanded.append(stripped)
    return list(dict.fromkeys(expanded))


def _apply_fasta_to_rows(
    rows: dict[str, dict],
    fasta_path: str | Path,
    *,
    force_type: str | None = None,
) -> None:
    records = read_fasta(fasta_path)
    seq_type = force_type or detect_query_type(records)
    by_key: dict[str, str] = {}
    for record in records:
        sequence = str(record.seq)
        for key in _record_lookup_keys(record):
            by_key.setdefault(key, sequence)
    for query_id, row in rows.items():
        sequence = by_key.get(query_id) or by_key.get(_strip_version(query_id), "")
        if not sequence:
            continue
        if seq_type == "nucleotide":
            row["length_nt"] = len(sequence)
            row["cdna_sequence"] = sequence
        else:
            row["length_aa"] = len(sequence)
            row["protein_sequence"] = sequence


def _query_metadata(
    query_ids: list[str],
    query_fasta: str | Path | None = None,
    cdna_fasta: str | Path | None = None,
    protein_fasta: str | Path | None = None,
) -> pd.DataFrame:
    if not query_ids:
        return pd.DataFrame(columns=QUERY_COLUMNS)
    rows = {query_id: _empty_query_row(query_id) for query_id in query_ids}
    if query_fasta is not None:
        _apply_fasta_to_rows(rows, query_fasta)
    if cdna_fasta is not None:
        _apply_fasta_to_rows(rows, cdna_fasta, force_type="nucleotide")
    if protein_fasta is not None:
        _apply_fasta_to_rows(rows, protein_fasta, force_type="protein")
    return pd.DataFrame(list(rows.values()))


def _query_length_info(
    query_fasta: str | Path | None,
) -> tuple[dict[str, int], bool]:
    if query_fasta is None:
        return {}, False
    records = read_fasta(query_fasta)
    query_is_nucleotide = detect_query_type(records) == "nucleotide"
    lengths: dict[str, int] = {}
    for record in records:
        length = len("".join(char for char in str(record.seq) if char not in "-."))
        for key in _record_lookup_keys(record):
            lengths.setdefault(key, length)
    return lengths, query_is_nucleotide


def _lookup_description(accession: object, accessions: pd.DataFrame) -> str:
    if accession is None or (isinstance(accession, float) and pd.isna(accession)):
        return "---"
    raw = str(accession).strip()
    candidates = [raw, normalize_accession(raw)]
    candidates.append(_strip_version(candidates[-1]))
    for key in dict.fromkeys(candidates):
        exact = accessions.loc[accessions["accession"] == key, "description"]
        if not exact.empty:
            return exact.iloc[0] or "---"
        by_norm = accessions.loc[accessions["accession_norm"] == key, "description"]
        if not by_norm.empty:
            return by_norm.iloc[0] or "---"
    return "---"


def _annotation_lookup(
    hits: pd.DataFrame,
    query_column: str,
) -> dict[str, list[dict]]:
    by_query: dict[str, list[dict]] = {}
    if hits.empty:
        return by_query
    ordered = hits.copy()
    if "evalue" in ordered.columns:
        ordered["_e"] = pd.to_numeric(ordered["evalue"], errors="coerce")
        ordered = ordered.sort_values([query_column, "_e"], na_position="last")
    for record in ordered.to_dict("records"):
        query_name = str(record.get(query_column, "")).strip()
        if not query_name:
            continue
        by_query.setdefault(query_name, []).append(record)
        stripped = _strip_version(query_name)
        if stripped != query_name:
            by_query.setdefault(stripped, []).append(record)
    return by_query


def _empty_pfam_row() -> dict:
    return {column: "---" for column in PFAM_COLUMNS}


def _pfam_lookup(hmm_hits: pd.DataFrame) -> dict[str, list[dict]]:
    return _annotation_lookup(hmm_hits, "query_name")


def _pfam_cell(records: list[dict], key: str) -> str:
    values: list[str] = []
    for record in records:
        raw = record.get(key)
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        text = str(raw).strip()
        if text in {"", "-", "nan", "<NA>"}:
            continue
        values.append(text)
    return "; ".join(values) if values else "---"


def _pfam_fields_for_query(query_id: str, by_query: dict[str, list[dict]]) -> dict:
    records = by_query.get(query_id) or by_query.get(_strip_version(query_id)) or []
    if not records:
        return _empty_pfam_row()
    return {
        "pfam_n_domains": len(records),
        "pfam_evalue": _pfam_cell(records, "evalue"),
        "pfam_score": _pfam_cell(records, "score"),
        "pfam_accession": _pfam_cell(records, "target_accession"),
        "pfam_name": _pfam_cell(records, "target_name"),
        "pfam_description": _pfam_cell(records, "description"),
    }


def _empty_kofam_row() -> dict:
    return {column: "---" for column in KOFAM_COLUMNS}


def _kofam_lookup(kofam_hits: pd.DataFrame) -> dict[str, list[dict]]:
    return _annotation_lookup(kofam_hits, "gene_name")


def _kofam_fields_for_query(query_id: str, by_query: dict[str, list[dict]]) -> dict:
    records = by_query.get(query_id) or by_query.get(_strip_version(query_id)) or []
    if not records:
        return _empty_kofam_row()
    return {
        "kofam_ko": _pfam_cell(records, "ko"),
        "kofam_score": _pfam_cell(records, "score"),
        "kofam_evalue": _pfam_cell(records, "evalue"),
        "kofam_definition": _pfam_cell(records, "definition"),
    }


def build_result_table(
    blast_hits: pd.DataFrame,
    accessions_path: str | Path,
    query_fasta: str | Path | None = None,
    cdna_fasta: str | Path | None = None,
    protein_fasta: str | Path | None = None,
    hmm_hits: pd.DataFrame | None = None,
    kofam_hits: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Wide Dataset S2-style table (Pfam and KOfam optional).

    One row per query with the best hit per database (lowest e-value,
    highest score). BLAST may return more subjects; only the top hit is shown.
    If ``hmm_hits`` is provided, a Pfam block is appended with every
    hmmscan hit that passed the cutoff (several domains per protein).
    If ``kofam_hits`` is provided, a KOfamScan block is appended with every
    KO assignment that passed the adaptive threshold.
    """
    logger.info("Joining BLAST hits with accession descriptors")
    accessions = load_accessions(accessions_path)
    hits = blast_hits.copy()
    if hits.empty:
        hits = pd.DataFrame(columns=[*OUTFMT6_COLUMNS, "database"])

    if "database" not in hits.columns:
        hits["database"] = "hit"

    hits["qseqid"] = hits["qseqid"].astype(str)
    query_ids = list(dict.fromkeys(hits["qseqid"].tolist()))
    if query_fasta is not None:
        fasta_ids = [record.id for record in read_fasta(query_fasta)]
        query_ids = list(dict.fromkeys([*fasta_ids, *query_ids]))
    if protein_fasta is not None:
        fasta_ids = [record.id for record in read_fasta(protein_fasta)]
        query_ids = list(dict.fromkeys([*fasta_ids, *query_ids]))

    hits["description"] = hits["sseqid"].map(
        lambda value: _lookup_description(value, accessions)
    )
    hits = assign_hit_rank(hits)

    databases = list(dict.fromkeys(hits["database"].astype(str).tolist()))
    query_meta = _query_metadata(
        query_ids,
        query_fasta=query_fasta,
        cdna_fasta=cdna_fasta,
        protein_fasta=protein_fasta,
    )
    query_lengths, query_is_nucleotide = _query_length_info(query_fasta)
    pfam_by_query = _pfam_lookup(hmm_hits) if hmm_hits is not None else None
    kofam_by_query = _kofam_lookup(kofam_hits) if kofam_hits is not None else None

    table_rows: list[dict] = []
    dash_fields = {"accession", "description", "identity_full_query"}
    for query_id in query_ids:
        base = query_meta.loc[query_meta["gene_id"] == query_id]
        query_row = base.iloc[0].to_dict() if not base.empty else {
            "gene_id": query_id,
            "length_nt": pd.NA,
            "cdna_sequence": pd.NA,
            "length_aa": pd.NA,
            "protein_sequence": pd.NA,
        }
        row = {column: query_row[column] for column in QUERY_COLUMNS}
        query_length = query_lengths.get(query_id) or query_lengths.get(
            _strip_version(query_id), 0
        )
        for database in databases:
            prefix = _safe_name(str(database))
            match = hits[
                (hits["qseqid"] == query_id)
                & (hits["database"].astype(str) == str(database))
                & (hits["hit_rank"] == 1)
            ]
            for field, source in HIT_FIELDS:
                column = f"{prefix}_{field}"
                if match.empty:
                    row[column] = "---" if field in dash_fields else pd.NA
                    continue
                if field == "identity_full_query":
                    subject = match.iloc[0]["sseqid"]
                    if pd.isna(subject) or str(subject) == "":
                        row[column] = "---"
                        continue
                    hsps = hits[
                        (hits["qseqid"] == query_id)
                        & (hits["database"].astype(str) == str(database))
                        & (hits["sseqid"] == subject)
                    ]
                    percent = full_query_identity_percent(
                        hsps.to_dict("records"),
                        query_length,
                        query_is_nucleotide=query_is_nucleotide,
                    )
                    row[column] = "---" if percent is None else percent
                    continue
                value = match.iloc[0][source]
                if field == "accession":
                    if pd.isna(value) or str(value) == "":
                        row[column] = "---"
                    else:
                        row[column] = normalize_accession(value) or str(value)
                elif field == "description":
                    row[column] = value if value else "---"
                else:
                    row[column] = value
        if pfam_by_query is not None:
            row.update(_pfam_fields_for_query(query_id, pfam_by_query))
        if kofam_by_query is not None:
            row.update(_kofam_fields_for_query(query_id, kofam_by_query))
        table_rows.append(row)

    table = pd.DataFrame(table_rows)
    logger.info(
        "Built result table: %s gene(s), %s column(s) (best hit per species)",
        len(table),
        table.shape[1],
    )
    return table
