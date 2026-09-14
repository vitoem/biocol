from pathlib import Path

import pytest
from Bio.SeqRecord import SeqRecord

from biocol.sequence.reader import read_fasta


def test_read_single_dna_sequence(fixtures_dir: Path) -> None:
    records = read_fasta(fixtures_dir / "dna.fa")
    assert len(records) == 1
    assert isinstance(records[0], SeqRecord)
    assert records[0].id == "seq1"
    assert str(records[0].seq).startswith("ATGC")
    assert len(records[0].seq) == len(str(records[0].seq))


def test_read_multifasta(fixtures_dir: Path) -> None:
    records = read_fasta(fixtures_dir / "multi.fa")
    assert [record.id for record in records] == ["gene_a", "gene_b", "gene_c"]
    assert all(record.seq for record in records)


def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        read_fasta(tmp_path / "no_existe.fa")
