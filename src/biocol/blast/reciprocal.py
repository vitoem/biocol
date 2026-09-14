"""Keep forward HSPs whose top subject is also the reverse top hit."""

from __future__ import annotations

import logging

import pandas as pd

from biocol.blast.parser import fill_missing_hits
from biocol.blast.ranking import assign_hit_rank
from biocol.metadata.accessions import normalize_accession

logger = logging.getLogger(__name__)

_RANK_HELPERS = ("_evalue_sort", "_score_sort", "hit_rank")


def _strip_version(accession: str) -> str:
    if "." in accession:
        head, tail = accession.rsplit(".", 1)
        if tail.isdigit():
            return head
    return accession


def sequence_ids_match(left: object, right: object) -> bool:
    """True if BLAST/FASTA identifiers refer to the same sequence."""
    first = normalize_accession(left)
    second = normalize_accession(right)
    if not first or not second:
        return False
    if first == second:
        return True
    return _strip_version(first) == _strip_version(second)


def _has_subject(value: object) -> bool:
    return value is not None and not (isinstance(value, float) and pd.isna(value)) and str(value) != ""


def _rank1_rows(ranked: pd.DataFrame, query_id: str, database: str) -> pd.DataFrame:
    return ranked[
        (ranked["qseqid"].astype(str) == str(query_id))
        & (ranked["database"].astype(str) == str(database))
        & (ranked["hit_rank"] == 1)
    ]


def _top_subject(ranked: pd.DataFrame, query_id: str, database: str) -> object | None:
    match = _rank1_rows(ranked, query_id, database)
    if match.empty:
        return None
    subject = match.iloc[0]["sseqid"]
    return subject if _has_subject(subject) else None


def _reverse_query_ids_for_subject(ranked: pd.DataFrame, subject: object, database: str) -> list[str]:
    if ranked.empty or "database" not in ranked.columns:
        return []
    db_rows = ranked[ranked["database"].astype(str) == str(database)]
    if db_rows.empty:
        return []
    seen: list[str] = []
    for query_id in db_rows["qseqid"].astype(str):
        if query_id in seen:
            continue
        if sequence_ids_match(query_id, subject):
            seen.append(query_id)
    return seen


def _drop_rank_helpers(frame: pd.DataFrame) -> pd.DataFrame:
    drop = [column for column in _RANK_HELPERS if column in frame.columns]
    return frame.drop(columns=drop) if drop else frame


def keep_reciprocal_best_hits(
    forward: pd.DataFrame,
    reverse: pd.DataFrame,
) -> pd.DataFrame:
    """Keep forward HSPs of pairs that are rank-1 in both directions.

    For each query and species, the forward top subject S is kept only if
    the reverse top subject of S is that query. Other forward subjects
    (rank 2+) are not tried. All HSPs of an accepted pair are retained.
    """
    if forward.empty:
        return forward.copy()

    fwd = assign_hit_rank(forward)
    rev = assign_hit_rank(reverse) if not reverse.empty else reverse
    databases = list(dict.fromkeys(fwd["database"].astype(str).tolist()))
    query_ids = list(dict.fromkeys(fwd["qseqid"].astype(str).tolist()))

    kept_parts: list[pd.DataFrame] = []
    accepted = 0
    for database in databases:
        fwd_db = fwd[fwd["database"].astype(str) == database]
        rows: list[pd.DataFrame] = []
        for query_id in query_ids:
            subject = _top_subject(fwd, query_id, database)
            if subject is None:
                continue
            reverse_queries = _reverse_query_ids_for_subject(rev, subject, database)
            reciprocal = False
            for reverse_query in reverse_queries:
                reverse_top = _top_subject(rev, reverse_query, database)
                if reverse_top is not None and sequence_ids_match(reverse_top, query_id):
                    reciprocal = True
                    break
            if not reciprocal:
                continue
            pair = fwd_db[
                (fwd_db["qseqid"].astype(str) == query_id)
                & fwd_db["sseqid"].map(lambda value, target=subject: sequence_ids_match(value, target))
            ]
            if pair.empty:
                continue
            rows.append(_drop_rank_helpers(pair))
            accepted += 1
        if rows:
            kept_parts.append(fill_missing_hits(pd.concat(rows, ignore_index=True), query_ids, database))
        else:
            empty = pd.DataFrame(columns=_drop_rank_helpers(fwd).columns)
            kept_parts.append(fill_missing_hits(empty, query_ids, database))

    combined = pd.concat(kept_parts, ignore_index=True)
    logger.info(
        "Reciprocal BLAST kept %s query-species pair(s)",
        accepted,
    )
    return combined
