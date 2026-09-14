from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biocol.sequence.classifier import detect_query_type, detect_sequence_type
from biocol.sequence.reader import read_fasta


def test_detect_dna() -> None:
    assert detect_sequence_type("ATGCGATCGTAGCTAG") == "nucleotide"


def test_detect_rna() -> None:
    assert detect_sequence_type("AUGCGAUCGAUCGAUC") == "nucleotide"


def test_detect_protein() -> None:
    assert detect_sequence_type("MVLSPADKTNVKAAWGKVGAHAGEYGAEALER") == "protein"


def test_detect_from_seq_and_seqrecord() -> None:
    assert detect_sequence_type(Seq("ATGCGATCGTAG")) == "nucleotide"
    record = SeqRecord(Seq("MVLSPADKTNVKAAWGKV"), id="p1")
    assert detect_sequence_type(record) == "protein"


def test_detect_from_dna_fasta(fixtures_dir: Path) -> None:
    records = read_fasta(fixtures_dir / "dna.fa")
    assert detect_query_type(records) == "nucleotide"


def test_detect_from_rna_fasta(fixtures_dir: Path) -> None:
    records = read_fasta(fixtures_dir / "rna.fa")
    assert detect_query_type(records) == "nucleotide"


def test_detect_from_protein_fasta(fixtures_dir: Path) -> None:
    records = read_fasta(fixtures_dir / "protein.fa")
    assert detect_query_type(records) == "protein"


def test_detect_multifasta_same_type(fixtures_dir: Path) -> None:
    records = read_fasta(fixtures_dir / "multi.fa")
    assert detect_query_type(records) == "nucleotide"


def test_mixed_multifasta_uses_majority(fixtures_dir: Path) -> None:
    records = read_fasta(fixtures_dir / "mixed.fa")
    assert detect_query_type(records) == "protein"


def test_mixed_fasta_nucleotide_majority() -> None:
    records = [
        SeqRecord(Seq("ATGCGATCGTAGCTAG"), id="n1"),
        SeqRecord(Seq("ATGCGATCGTAGCTAG"), id="n2"),
        SeqRecord(Seq("MVLSPADKTNVKAAWGKV"), id="p1"),
    ]
    assert detect_query_type(records) == "nucleotide"


def test_poly_lysine_is_protein_not_nucleotide() -> None:
    assert detect_sequence_type("K" * 169) == "protein"


def test_dna_with_iupac_ambiguity_is_nucleotide() -> None:
    assert detect_sequence_type("ATGCGATCGTAGCTAGNNNNKK") == "nucleotide"


def test_unclassifiable_sequence_raises() -> None:
    with pytest.raises(ValueError, match="classifiable"):
        detect_sequence_type("---")


def test_detect_sequence_type_logs_classification(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level("DEBUG", logger="biocol.sequence.classifier"):
        detect_sequence_type("ATGCGATCGTAGCTAG")
    assert "type=nucleotide" in caplog.text


def test_detect_query_type_logs_fasta_type(
    fixtures_dir: Path, caplog: pytest.LogCaptureFixture
) -> None:
    records = read_fasta(fixtures_dir / "protein.fa")
    with caplog.at_level("INFO", logger="biocol.sequence.classifier"):
        detect_query_type(records)
    assert "classified as protein" in caplog.text
