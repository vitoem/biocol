"""Keep HSPs whose BLAST ``pident`` meets a minimum identity cutoff."""

from __future__ import annotations

import logging

import pandas as pd

from biocol.blast.parser import fill_missing_hits

logger = logging.getLogger(__name__)


def filter_hits_by_pident(
    hits: pd.DataFrame,
    min_identity: float | None,
    query_ids: list[str] | None = None,
) -> pd.DataFrame:
    """Drop HSP rows with ``pident`` below ``min_identity``.

    ``None`` leaves the table unchanged (no cutoff). Empty BLAST rows (no
    subject) are kept. Queries that lose every HSP in a database get an
    empty row again so the TSV still shows ``---``.
    """
    if min_identity is None:
        return hits
    if min_identity < 0 or min_identity > 100:
        raise ValueError("min_identity must be between 0 and 100")
    if hits.empty:
        return hits

    logger.info("Keeping HSPs with pident >= %s", min_identity)
    empty = hits["sseqid"].isna() | (hits["sseqid"].astype(str) == "")
    pident = pd.to_numeric(hits["pident"], errors="coerce")
    kept = hits.loc[empty | (pident >= min_identity)].copy()
    dropped = int((~empty & (pident < min_identity)).sum())
    logger.info("Dropped %s HSP(s) below identity cutoff", dropped)

    ids = query_ids or list(dict.fromkeys(hits["qseqid"].astype(str).tolist()))
    if "database" not in hits.columns:
        return fill_missing_hits(kept, ids, "hit")

    frames = [
        fill_missing_hits(
            kept.loc[kept["database"].astype(str) == str(database)],
            ids,
            str(database),
        )
        for database in dict.fromkeys(hits["database"].astype(str).tolist())
    ]
    return pd.concat(frames, ignore_index=True)
