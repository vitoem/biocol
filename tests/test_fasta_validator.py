from pathlib import Path

import pytest

from biocol.exceptions import EmptyFastaError, InvalidFastaError
from biocol.sequence.validator import validate_fasta_file


def test_existing_valid_fasta(fixtures_dir: Path) -> None:
    path = validate_fasta_file(fixtures_dir / "dna.fa")
    assert path.exists()


def test_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        validate_fasta_file(tmp_path / "ausente.fa")


def test_wrong_extension(tmp_path: Path) -> None:
    bad = tmp_path / "secuencias.txt"
    bad.write_text(">s1\nATGC\n", encoding="utf-8")
    with pytest.raises(InvalidFastaError, match="Unsupported extension"):
        validate_fasta_file(bad)


def test_empty_fasta(tmp_path: Path) -> None:
    empty = tmp_path / "vacio.fa"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(EmptyFastaError):
        validate_fasta_file(empty)


def test_not_fasta_content(fixtures_dir: Path) -> None:
    with pytest.raises(InvalidFastaError, match="Could not parse"):
        validate_fasta_file(fixtures_dir / "not_fasta.fa")


def test_empty_sequence(fixtures_dir: Path) -> None:
    with pytest.raises(InvalidFastaError, match="empty"):
        validate_fasta_file(fixtures_dir / "empty_seq.fa")


def test_invalid_characters(fixtures_dir: Path) -> None:
    with pytest.raises(InvalidFastaError, match="invalid characters"):
        validate_fasta_file(fixtures_dir / "invalid_chars.fa")
