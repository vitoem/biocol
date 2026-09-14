from biocol.processing.identity import (
    full_query_identity_percent,
    identical_query_positions,
)


def test_identical_positions_plus_strand() -> None:
    matched = identical_query_positions(
        "ACGT",
        "ACGT",
        1,
        4,
        query_is_nucleotide=True,
    )
    assert matched == {1, 2, 3, 4}


def test_identical_positions_skips_gaps() -> None:
    matched = identical_query_positions(
        "A-CG",
        "ATCG",
        1,
        3,
        query_is_nucleotide=True,
    )
    assert matched == {1, 2, 3}


def test_identical_positions_skips_mismatches() -> None:
    matched = identical_query_positions(
        "ACGT",
        "ACCT",
        1,
        4,
        query_is_nucleotide=True,
    )
    assert matched == {1, 2, 4}


def test_protein_query_counts_one_position_per_residue() -> None:
    matched = identical_query_positions(
        "MVL",
        "MVL",
        1,
        3,
        query_is_nucleotide=False,
    )
    assert matched == {1, 2, 3}


def test_identical_positions_minus_strand() -> None:
    matched = identical_query_positions(
        "AT",
        "AT",
        10,
        9,
        query_is_nucleotide=True,
    )
    assert matched == {9, 10}


def test_blastx_marks_three_nucleotides_per_amino_acid() -> None:
    matched = identical_query_positions(
        "MVL",
        "MVL",
        1,
        9,
        query_is_nucleotide=True,
    )
    assert matched == {1, 2, 3, 4, 5, 6, 7, 8, 9}


def test_full_query_identity_unions_non_overlapping_hsps() -> None:
    hsps = [
        {"qseq": "AAA", "sseq": "AAA", "qstart": 1, "qend": 3},
        {"qseq": "TTT", "sseq": "TTT", "qstart": 4, "qend": 6},
    ]
    percent = full_query_identity_percent(hsps, 10, query_is_nucleotide=True)
    assert percent == 60.0


def test_full_query_identity_does_not_double_count_overlap() -> None:
    hsps = [
        {"qseq": "AAAA", "sseq": "AAAA", "qstart": 1, "qend": 4},
        {"qseq": "AAAA", "sseq": "AAAA", "qstart": 1, "qend": 4},
    ]
    percent = full_query_identity_percent(hsps, 10, query_is_nucleotide=True)
    assert percent == 40.0


def test_full_query_identity_none_without_sequences() -> None:
    hsps = [{"qseq": None, "sseq": None, "qstart": 1, "qend": 4}]
    assert full_query_identity_percent(hsps, 10, query_is_nucleotide=True) is None


def test_full_query_identity_none_when_length_zero() -> None:
    hsps = [{"qseq": "A", "sseq": "A", "qstart": 1, "qend": 1}]
    assert full_query_identity_percent(hsps, 0, query_is_nucleotide=True) is None
