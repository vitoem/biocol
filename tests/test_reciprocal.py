import pandas as pd

from biocol.blast.reciprocal import keep_reciprocal_best_hits, sequence_ids_match


def _hits(rows: list[tuple]) -> pd.DataFrame:
    qseqid, sseqid, evalue, bitscore, pident, database = zip(*rows)
    return pd.DataFrame(
        {
            "qseqid": list(qseqid),
            "sseqid": list(sseqid),
            "evalue": list(evalue),
            "bitscore": list(bitscore),
            "pident": list(pident),
            "database": list(database),
        }
    )


def test_sequence_ids_match_normalizes_prefixes() -> None:
    assert sequence_ids_match("ref|XP_001.1|", "XP_001.1")
    assert sequence_ids_match("XP_001.1", "XP_001")
    assert not sequence_ids_match("XP_001", "XP_002")


def test_keeps_pair_that_is_rank1_both_ways() -> None:
    forward = _hits(
        [
            ("Q1", "S1", 1e-50, 200.0, 99.0, "lyrata"),
            ("Q1", "S2", 1e-10, 150.0, 80.0, "lyrata"),
        ]
    )
    reverse = _hits(
        [
            ("S1", "Q1", 1e-48, 190.0, 98.0, "lyrata"),
            ("S1", "Q2", 1e-5, 40.0, 70.0, "lyrata"),
        ]
    )
    kept = keep_reciprocal_best_hits(forward, reverse)
    subjects = list(kept.loc[kept["sseqid"].notna(), "sseqid"])
    assert subjects == ["S1"]
    assert "Q1" in set(kept["qseqid"].astype(str))


def test_rejects_when_reverse_top_is_another_query() -> None:
    forward = _hits([("Q1", "S1", 1e-50, 200.0, 99.0, "lyrata")])
    reverse = _hits(
        [
            ("S1", "Q2", 1e-60, 210.0, 99.0, "lyrata"),
            ("S1", "Q1", 1e-10, 80.0, 80.0, "lyrata"),
        ]
    )
    kept = keep_reciprocal_best_hits(forward, reverse)
    q1 = kept[kept["qseqid"] == "Q1"].iloc[0]
    assert pd.isna(q1["sseqid"])


def test_does_not_fall_back_to_forward_rank2() -> None:
    forward = _hits(
        [
            ("Q1", "S1", 1e-50, 200.0, 99.0, "lyrata"),
            ("Q1", "S2", 1e-10, 150.0, 90.0, "lyrata"),
        ]
    )
    reverse = _hits(
        [
            ("S1", "Q2", 1e-40, 180.0, 95.0, "lyrata"),
            ("S2", "Q1", 1e-30, 160.0, 92.0, "lyrata"),
        ]
    )
    kept = keep_reciprocal_best_hits(forward, reverse)
    q1 = kept[kept["qseqid"] == "Q1"].iloc[0]
    assert pd.isna(q1["sseqid"])


def test_keeps_all_hsps_of_accepted_pair() -> None:
    forward = _hits(
        [
            ("Q1", "S1", 1e-50, 200.0, 99.0, "lyrata"),
            ("Q1", "S1", 1e-20, 80.0, 70.0, "lyrata"),
            ("Q1", "S2", 1e-10, 60.0, 60.0, "lyrata"),
        ]
    )
    reverse = _hits([("S1", "Q1", 1e-48, 190.0, 98.0, "lyrata")])
    kept = keep_reciprocal_best_hits(forward, reverse)
    pair = kept[(kept["qseqid"] == "Q1") & (kept["sseqid"] == "S1")]
    assert len(pair) == 2
    assert "S2" not in set(kept["sseqid"].dropna().astype(str))


def test_empty_reverse_clears_hits() -> None:
    forward = _hits([("Q1", "S1", 1e-50, 200.0, 99.0, "lyrata")])
    kept = keep_reciprocal_best_hits(forward, pd.DataFrame())
    assert pd.isna(kept.iloc[0]["sseqid"])
