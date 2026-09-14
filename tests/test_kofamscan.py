from pathlib import Path

import pandas as pd
import pytest

from biocol import (
    KofamError,
    KofamExecutionError,
    build_result_table,
    parse_kofam_detail_tsv,
    run_kofamscan,
    write_results_csv,
)
from biocol.kofam.runner import build_exec_annotation_command


def test_parse_kofam_detail_tsv_keeps_star_hits_only(fixtures_dir: Path) -> None:
    hits = parse_kofam_detail_tsv(fixtures_dir / "kofam_detail.tsv")
    assert set(hits["gene_name"]) == {"gene2", "seq1"}
    assert "*" not in set(hits["ko"])
    gene2 = hits[hits["gene_name"] == "gene2"]
    assert list(gene2["ko"]) == ["K00002", "K00004"]


def test_exec_annotation_command(tmp_path: Path) -> None:
    command = build_exec_annotation_command(
        tmp_path / "q.faa",
        tmp_path / "eukaryote.hal",
        tmp_path / "ko_list",
        tmp_path / "kofam" / "kofam.tsv",
        tmp_path / "kofam" / "tmp",
        num_threads=40,
    )
    assert command[0] == "exec_annotation"
    assert command[command.index("-f") + 1] == "detail-tsv"
    assert command[command.index("--cpu") + 1] == "40"
    assert "--e-value" not in command
    assert "-E" not in command
    assert "-T" not in command


def test_run_kofamscan_rejects_nucleotide_query(
    tmp_path: Path, fixtures_dir: Path
) -> None:
    with pytest.raises(KofamError, match="protein sequence"):
        run_kofamscan(
            fixtures_dir / "dna.fa",
            tmp_path / "eukaryote.hal",
            tmp_path / "ko_list",
        )


def test_run_kofamscan_requires_executable(
    tmp_path: Path, fixtures_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    profile = tmp_path / "eukaryote.hal"
    profile.write_text("# list\n", encoding="utf-8")
    ko_list = tmp_path / "ko_list"
    ko_list.write_text("K00001\n", encoding="utf-8")
    monkeypatch.setattr("biocol.kofam.execute.shutil.which", lambda _name: None)
    with pytest.raises(KofamExecutionError, match="was not found"):
        run_kofamscan(
            fixtures_dir / "protein.fa",
            profile,
            ko_list,
            kofam_dir=tmp_path / "kofam",
        )


def test_run_kofamscan_writes_tsv_and_keeps_tmp(
    tmp_path: Path, fixtures_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    profile = tmp_path / "profiles"
    profile.mkdir()
    ko_list = tmp_path / "ko_list"
    ko_list.write_text("k\n", encoding="utf-8")
    dest = tmp_path / "kofam_out"

    def fake_run(command: list[str]) -> None:
        out = Path(command[command.index("-o") + 1])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            (fixtures_dir / "kofam_detail.tsv").read_text(encoding="utf-8"),
            encoding="utf-8",
        )

    monkeypatch.setattr("biocol.kofam.runner.run_exec_annotation", fake_run)
    hits = run_kofamscan(
        fixtures_dir / "protein.fa",
        profile,
        ko_list,
        kofam_dir=dest,
        num_threads=2,
    )
    assert (dest / "kofam.tsv").exists()
    assert (dest / "tmp").is_dir()
    assert not hits.empty
    assert set(hits["gene_name"]) <= {"gene2", "seq1"}


def test_run_kofamscan_rejects_non_profile_file(
    tmp_path: Path, fixtures_dir: Path
) -> None:
    bogus = tmp_path / "ko_list.txt"
    bogus.write_text("x\n", encoding="utf-8")
    ko_list = tmp_path / "ko_list"
    ko_list.write_text("k\n", encoding="utf-8")
    with pytest.raises(KofamError, match="directory of .hmm files"):
        run_kofamscan(fixtures_dir / "protein.fa", bogus, ko_list)


def test_kofam_block_in_result_table(fixtures_dir: Path, tmp_path: Path) -> None:
    kofam_hits = parse_kofam_detail_tsv(fixtures_dir / "kofam_detail.tsv")
    table = build_result_table(
        pd.DataFrame(
            {
                "qseqid": ["seq1"],
                "sseqid": ["XP_001"],
                "pident": [99.0],
                "length": [10],
                "mismatch": [0],
                "gapopen": [0],
                "qstart": [1],
                "qend": [10],
                "sstart": [1],
                "send": [10],
                "evalue": [1e-5],
                "bitscore": [50.0],
                "database": ["amborella"],
            }
        ),
        fixtures_dir / "accessions.txt",
        kofam_hits=kofam_hits,
    )
    row = table.iloc[0]
    assert row["kofam_ko"] == "K00003"
    assert row["kofam_definition"] != "---"
    path = write_results_csv(table, tmp_path / "out.tsv")
    loaded = pd.read_csv(path, header=None, sep="\t")
    assert "KOfamScan" in set(loaded.iloc[0].astype(str))
    assert "KO" in set(loaded.iloc[2].astype(str))
    assert "KO definition" in set(loaded.iloc[2].astype(str))
