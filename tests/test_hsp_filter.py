import pandas as pd
import pytest

from biocol.processing.hsp_filter import filter_hits_by_pident


def _hits() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "qseqid": ["Q1", "Q1", "Q1", "Q2", "Q2"],
            "sseqid": ["S1", "S1", "S2", "S3", "S3"],
            "pident": [95.0, 78.0, 85.0, 62.0, 91.0],
            "database": ["lyrata"] * 5,
        }
    )


def test_no_cutoff_returns_same_table() -> None:
    hits = _hits()
    filtered = filter_hits_by_pident(hits, None)
    assert filtered is hits


def test_keeps_hsps_at_or_above_cutoff() -> None:
    filtered = filter_hits_by_pident(_hits(), 80)
    pidents = list(filtered.loc[filtered["sseqid"].notna(), "pident"])
    assert pidents == [95.0, 85.0, 91.0]


def test_decimal_cutoff() -> None:
    filtered = filter_hits_by_pident(_hits(), 90.5)
    pidents = list(filtered.loc[filtered["sseqid"].notna(), "pident"])
    assert pidents == [95.0, 91.0]


def test_query_with_no_remaining_hsp_becomes_empty_row() -> None:
    hits = pd.DataFrame(
        {
            "qseqid": ["Q1", "Q2"],
            "sseqid": ["S1", "S2"],
            "pident": [95.0, 62.0],
            "database": ["lyrata", "lyrata"],
        }
    )
    filtered = filter_hits_by_pident(hits, 80)
    q2 = filtered[filtered["qseqid"] == "Q2"].iloc[0]
    assert pd.isna(q2["sseqid"])


def test_invalid_cutoff_raises() -> None:
    with pytest.raises(ValueError, match="between 0 and 100"):
        filter_hits_by_pident(_hits(), 120)
