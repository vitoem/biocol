from pathlib import Path

import pytest

from biocol import (
    DatabaseError,
    MixedDatabaseTypeError,
    detect_database_type,
    list_blast_databases,
    select_blast_program,
)
from biocol.blast.databases import infer_database_label


def test_infer_database_label_from_ncbi_organism(tmp_path: Path) -> None:
    faa = tmp_path / "protein.faa"
    faa.write_text(
        ">NP_1 MLO-like protein 1 [Cucumis melo]\nMVLSPADKTNVKAAWGKV\n"
        ">NP_2 kinase [Cucumis melo]\nMVLSPADKTNVKAAWGKV\n",
        encoding="utf-8",
    )
    assert infer_database_label(faa) == "Cucumis melo"


def test_infer_database_label_falls_back_to_stem(tmp_path: Path) -> None:
    faa = tmp_path / "amborella.faa"
    faa.write_text(">a1\nMVLSPADKTNVKAAWGKV\n", encoding="utf-8")
    assert infer_database_label(faa) == "amborella"


def test_fasta_as_database(fixtures_dir: Path) -> None:
    assert detect_database_type(fixtures_dir / "dna.fa") == "nucleotide"
    assert detect_database_type(fixtures_dir / "protein.fa") == "protein"


def test_faa_extension(tmp_path: Path) -> None:
    faa = tmp_path / "hits.faa"
    faa.write_text(">p1\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALER\n", encoding="utf-8")
    assert detect_database_type(faa) == "protein"


def test_folder_of_same_type(tmp_path: Path) -> None:
    (tmp_path / "amborella.faa").write_text(
        ">a1\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALER\n", encoding="utf-8"
    )
    nested = tmp_path / "plantas"
    nested.mkdir()
    (nested / "vitis.fa").write_text(
        ">v1\nMVLSPADKTNVKAAWGKV\n", encoding="utf-8"
    )
    assert detect_database_type(tmp_path) == "protein"
    listed = list_blast_databases(tmp_path)
    assert len(listed) == 2
    assert {db_type for _, db_type in listed} == {"protein"}


def test_folder_mixed_types_raises(fixtures_dir: Path, tmp_path: Path) -> None:
    (tmp_path / "dna.fa").write_text(
        (fixtures_dir / "dna.fa").read_text(encoding="utf-8"), encoding="utf-8"
    )
    (tmp_path / "protein.fa").write_text(
        (fixtures_dir / "protein.fa").read_text(encoding="utf-8"), encoding="utf-8"
    )
    with pytest.raises(MixedDatabaseTypeError):
        detect_database_type(tmp_path)


def test_empty_folder_raises(tmp_path: Path) -> None:
    with pytest.raises(DatabaseError):
        detect_database_type(tmp_path)


def test_non_fasta_file_raises(tmp_path: Path) -> None:
    other = tmp_path / "db.txt"
    other.write_text("no fasta", encoding="utf-8")
    with pytest.raises(DatabaseError, match="must be FASTA"):
        detect_database_type(other)


def test_unknown_source_raises() -> None:
    with pytest.raises(DatabaseError):
        detect_database_type("no_existe_xyz")


@pytest.mark.parametrize(
    ("query", "database", "translated", "program"),
    [
        ("nucleotide", "protein", False, "blastx"),
        ("protein", "protein", False, "blastp"),
        ("protein", "nucleotide", False, "tblastn"),
        ("nucleotide", "nucleotide", False, "blastn"),
        ("nucleotide", "nucleotide", True, "tblastx"),
    ],
)
def test_select_blast_program(
    query: str, database: str, translated: bool, program: str
) -> None:
    assert select_blast_program(query, database, translated=translated) == program


def test_nuc_vs_nuc_defaults_to_blastn() -> None:
    assert select_blast_program("nucleotide", "nucleotide") == "blastn"


def test_translated_ignored_when_not_both_nucleotide() -> None:
    assert select_blast_program("protein", "protein", translated=True) == "blastp"
