"""Rank BLAST HSPs the same way the result table picks the top hit."""

from __future__ import annotations

import pandas as pd


def assign_hit_rank(hits: pd.DataFrame) -> pd.DataFrame:
    """Sort by lowest e-value, then highest bitscore, and number ranks.

    Rank 1 is the first remaining HSP per query and database after that
    sort. Empty subject rows are left at rank 1 so the TSV can show
    ``---``. Ties keep the first row after sorting.
    """
    ranked = hits.copy()
    if ranked.empty:
        ranked["hit_rank"] = pd.Series(dtype="int64")
        return ranked
    if "database" not in ranked.columns:
        ranked["database"] = "hit"
    ranked["_evalue_sort"] = pd.to_numeric(ranked["evalue"], errors="coerce")
    ranked["_score_sort"] = pd.to_numeric(ranked["bitscore"], errors="coerce")
    ranked = ranked.sort_values(
        ["qseqid", "database", "_evalue_sort", "_score_sort"],
        ascending=[True, True, True, False],
        na_position="last",
    )
    ranked["hit_rank"] = ranked.groupby(["qseqid", "database"], dropna=False).cumcount() + 1
    empty_mask = ranked["sseqid"].isna() | (ranked["sseqid"].astype(str) == "")
    ranked.loc[empty_mask, "hit_rank"] = 1
    return ranked
