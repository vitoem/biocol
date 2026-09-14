from pathlib import Path

import pandas as pd
import pytest

from biocol import (
    DEFAULT_OUTPUT,
    MetadataError,
    QUERY_COLUMNS,
    build_result_table,
    load_accessions,
    normalize_accession,
    write_results_csv,
)


def _hits() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "qseqid": ["seq1", "seq1", "seq1"],
            "sseqid": ["ref|XP_001.1|", "ref|XP_002.1|", "gb|NO_MATCH.1|"],
            "pident": [99.0, 80.0, 70.0],
            "length": [100, 90, 80],
            "mismatch": [0, 1, 2],
            "gapopen": [0, 0, 0],
            "qstart": [1, 1, 1],
            "qend": [100, 90, 80],
            "sstart": [1, 1, 1],
            "send": [100, 90, 80],
            "evalue": [1e-50, 1e-10, 1e-5],
            "bitscore": [200.0, 150.0, 80.0],
            "database": ["amborella", "vitis", "amborella"],
        }
    )


def test_normalize_accession_strips_ref_prefix() -> None:
    assert normalize_accession("ref|XP_001.1|") == "XP_001.1"
    assert normalize_accession("XP_001") == "XP_001"


def test_load_accessions(fixtures_dir: Path) -> None:
    frame = load_accessions(fixtures_dir / "accessions.txt")
    assert list(frame["accession"]) == ["XP_001", "XP_002", "XP_003"]
    assert "Amborella protein one" in list(frame["description"])


def test_load_accessions_empty_raises(tmp_path: Path) -> None:
    empty = tmp_path / "vacio.txt"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(MetadataError):
        load_accessions(empty)


def test_wide_table_best_hit_and_descriptions(fixtures_dir: Path) -> None:
    table = build_result_table(
        _hits(),
        fixtures_dir / "accessions.txt",
        query_fasta=fixtures_dir / "dna.fa",
    )
    assert list(table.columns[:5]) == QUERY_COLUMNS
    assert "amborella_accession" in table.columns
    assert "vitis_description" in table.columns
    assert len(table) == 1

    best = table.iloc[0]
    assert best["gene_id"] == "seq1"
    assert best["length_nt"] == 28
    assert isinstance(best["cdna_sequence"], str)
    assert pd.isna(best["length_aa"])
    assert pd.isna(best["protein_sequence"])
    assert best["amborella_accession"] == "XP_001.1"
    assert best["amborella_description"] == "Amborella protein one"
    assert best["vitis_accession"] == "XP_002.1"
    assert best["vitis_description"] == "Vitis protein two"
    assert best["amborella_identity_full_query"] == "---"


def test_protein_fasta_fills_aa_only(fixtures_dir: Path) -> None:
    hits = _hits()
    hits["qseqid"] = "prot1"
    table = build_result_table(
        hits,
        fixtures_dir / "accessions.txt",
        query_fasta=fixtures_dir / "protein.fa",
    )
    row = table.iloc[0]
    assert row["gene_id"] == "prot1"
    assert pd.isna(row["length_nt"])
    assert row["length_aa"] == 60
    assert isinstance(row["protein_sequence"], str)


def test_camino_2_without_fasta(fixtures_dir: Path) -> None:
    table = build_result_table(
        _hits(),
        fixtures_dir / "accessions.txt",
        query_fasta=None,
    )
    assert table.iloc[0]["gene_id"] == "seq1"
    assert pd.isna(table.iloc[0]["cdna_sequence"])
    assert pd.isna(table.iloc[0]["protein_sequence"])


def test_write_results_csv_default_name(tmp_path: Path, fixtures_dir: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    table = build_result_table(_hits(), fixtures_dir / "accessions.txt")
    path = write_results_csv(table)
    assert path.name == DEFAULT_OUTPUT
    assert path.exists()
    loaded = pd.read_csv(path, header=None, sep="\t")
    section = loaded.iloc[0].astype(str).tolist()
    assert section.count("Annotation based on top-BLAST-hit method") == 2
    assert "Accesion No." in set(loaded.iloc[2].astype(str))
    assert "e-value" in set(loaded.iloc[2].astype(str))
    assert "Score" in set(loaded.iloc[2].astype(str))
    assert "Identity % (full query)" in set(loaded.iloc[2].astype(str))
    assert "amborella" in set(loaded.iloc[1].astype(str))
    assert "vitis" in set(loaded.iloc[1].astype(str))


def test_write_results_csv_drops_empty_nucleotide_columns(tmp_path: Path, fixtures_dir: Path) -> None:
    hits = _hits()
    hits["qseqid"] = "prot1"
    table = build_result_table(
        hits,
        fixtures_dir / "accessions.txt",
        query_fasta=fixtures_dir / "protein.fa",
    )
    path = write_results_csv(table, tmp_path / "prot.tsv")
    labels = pd.read_csv(path, header=None, sep="\t").iloc[2].astype(str).tolist()
    assert "Gene ID" in labels
    assert "Length(aa)" in labels
    assert "Protein Sequences (aa)" in labels
    assert "Length (nt)" not in labels
    assert "cDNA Sequences (nt)" not in labels


def test_gene_model_fills_nt_and_aa(tmp_path: Path, fixtures_dir: Path) -> None:
    cdna = tmp_path / "gene.fna"
    pep = tmp_path / "gene.faa"
    cdna.write_text(">seq1\nATGCGATCGATCGATCGTAGCTAGCTAG\n", encoding="utf-8")
    pep.write_text(">seq1\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALER\n", encoding="utf-8")
    table = build_result_table(
        _hits(),
        fixtures_dir / "accessions.txt",
        query_fasta=pep,
        cdna_fasta=cdna,
        protein_fasta=pep,
    )
    row = table.iloc[0]
    assert row["length_nt"] == 28
    assert row["length_aa"] == 32
    assert isinstance(row["cdna_sequence"], str)
    assert isinstance(row["protein_sequence"], str)


def test_cds_matches_protein_id_in_header(tmp_path: Path, fixtures_dir: Path) -> None:
    pep = tmp_path / "prot.faa"
    cds = tmp_path / "cds.fna"
    pep.write_text(">NP_171609.1\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALER\n", encoding="utf-8")
    cds.write_text(
        ">lcl|NC_003070.9_cds_NP_171609.1_1 [protein_id=NP_171609.1] [gene=NAC001]\n"
        "ATGCGATCGATCGATCGTAGCTAGCTAG\n",
        encoding="utf-8",
    )
    hits = _hits()
    hits["qseqid"] = "NP_171609.1"
    table = build_result_table(
        hits,
        fixtures_dir / "accessions.txt",
        query_fasta=pep,
        cdna_fasta=cds,
        protein_fasta=pep,
    )
    row = table.iloc[0]
    assert row["gene_id"] == "NP_171609.1"
    assert row["length_nt"] == 28
    assert isinstance(row["cdna_sequence"], str)


def test_full_query_identity_column_from_aligned_sequences(fixtures_dir: Path) -> None:
    sequence = "ATGCGATCGATCGATCGTAGCTAGCTAG"
    hits = _hits()
    hits = hits[hits["sseqid"] == "ref|XP_001.1|"].copy()
    hits["qstart"] = 1
    hits["qend"] = 28
    hits["qseq"] = sequence
    hits["sseq"] = sequence
    table = build_result_table(
        hits,
        fixtures_dir / "accessions.txt",
        query_fasta=fixtures_dir / "dna.fa",
    )
    row = table.iloc[0]
    assert row["amborella_identity_full_query"] == 100.0


def test_pfam_block_after_blast(fixtures_dir: Path, tmp_path: Path) -> None:
    hmm_hits = pd.DataFrame(
        {
            "query_name": ["seq1", "seq1"],
            "target_name": ["Pkinase", "Pkinase_Tyr"],
            "target_accession": ["PF00069.25", "PF07714.17"],
            "evalue": [8.2e-26, 1.3e-14],
            "score": [90.9, 54.1],
            "description": [
                "Protein kinase domain",
                "Protein tyrosine kinase",
            ],
        }
    )
    table = build_result_table(
        _hits(),
        fixtures_dir / "accessions.txt",
        hmm_hits=hmm_hits,
    )
    row = table.iloc[0]
    assert row["pfam_n_domains"] == 2
    assert row["pfam_accession"] == "PF00069.25; PF07714.17"
    assert row["pfam_name"] == "Pkinase; Pkinase_Tyr"
    assert row["pfam_description"] == (
        "Protein kinase domain; Protein tyrosine kinase"
    )
    path = write_results_csv(table, tmp_path / "out.tsv")
    loaded = pd.read_csv(path, header=None, sep="\t")
    assert "Pfam domains" in set(loaded.iloc[0].astype(str))
    assert "# of Pfam domain identified" in set(loaded.iloc[2].astype(str))
    assert "Accesion" in set(loaded.iloc[2].astype(str))
    assert "Name" in set(loaded.iloc[2].astype(str))
    assert "Description of target" in set(loaded.iloc[2].astype(str))
