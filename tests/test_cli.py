from pathlib import Path

import pandas as pd

from biocol.cli import main
from biocol.cli.parser import build_parser


def test_from_blast_writes_tsv(tmp_path: Path, fixtures_dir: Path) -> None:
    output = tmp_path / "out.tsv"
    code = main(
        [
            "from-blast",
            "--blast",
            str(fixtures_dir / "blast_outfmt6.txt"),
            "--accessions",
            str(fixtures_dir / "accessions.txt"),
            "--output",
            str(output),
        ]
    )
    assert code == 0
    assert output.exists()
    table = pd.read_csv(output, header=2, sep="\t")
    assert "Gene ID" in table.columns
    assert "Accesion No." in table.columns
    assert "q1" in set(table["Gene ID"].astype(str))


def test_from_blast_missing_file_returns_error(tmp_path: Path) -> None:
    code = main(
        [
            "from-blast",
            "--blast",
            str(tmp_path / "missing.txt"),
            "--accessions",
            str(tmp_path / "acc.txt"),
        ]
    )
    assert code == 1


def test_from_blast_hmm_db_without_protein_returns_error(
    tmp_path: Path, fixtures_dir: Path
) -> None:
    code = main(
        [
            "--no-color",
            "from-blast",
            "--blast",
            str(fixtures_dir / "blast_outfmt6.txt"),
            "--accessions",
            str(fixtures_dir / "accessions.txt"),
            "--hmm-db",
            str(tmp_path / "Pfam-A.hmm"),
            "--output",
            str(tmp_path / "out.tsv"),
        ]
    )
    assert code == 1


def test_run_help_is_english(capsys) -> None:
    try:
        main(["--no-color", "run", "--help"])
    except SystemExit as exc:
        assert exc.code == 0
    else:
        raise AssertionError("expected SystemExit from --help")
    text = capsys.readouterr().out
    assert "BIOCOL" in text
    assert "PROGRAM SELECTION" in text
    assert "--protein" in text
    assert "--blast-dir" in text
    assert "--min-identity" in text
    assert "--reciprocal" in text
    assert "--diamond" in text
    assert "--diamond-dir" in text
    assert "--hmm-db" in text
    assert "--hmm-dir" in text
    assert "--kofam-profile" in text
    assert "--ko-list" in text
    assert "--kofam-dir" in text
    assert "results.tsv" in text
    assert "\\033[" not in text


def test_cli_prints_stage_logs_on_stderr(
    tmp_path: Path, fixtures_dir: Path, capsys
) -> None:
    output = tmp_path / "out.tsv"
    code = main(
        [
            "--no-color",
            "from-blast",
            "--blast",
            str(fixtures_dir / "blast_outfmt6.txt"),
            "--accessions",
            str(fixtures_dir / "accessions.txt"),
            "--output",
            str(output),
        ]
    )
    captured = capsys.readouterr()
    assert code == 0
    assert str(output) in captured.out
    assert "biocol" in captured.err
    assert "Reading BLAST results" in captured.err
    assert "Joining BLAST hits" in captured.err
    assert "Wrote TSV" in captured.err
    assert "Done." in captured.err
    assert "\\033[" not in captured.err


def test_run_help_mentions_reciprocal_not_in_from_blast(capsys) -> None:
    try:
        main(["--no-color", "from-blast", "--help"])
    except SystemExit as exc:
        assert exc.code == 0
    else:
        raise AssertionError("expected SystemExit from --help")
    text = capsys.readouterr().out
    assert "--reciprocal" not in text
    assert "--diamond" not in text
    assert "--kofam-profile" in text
    assert "--ko-list" in text
    assert "only available in 'biocol run'" in text


def test_run_passes_reciprocal_to_run_blast(
    tmp_path: Path, fixtures_dir: Path, monkeypatch
) -> None:
    captured: dict = {}

    def fake_run_blast(*_args, **kwargs):
        captured.update(kwargs)
        return pd.DataFrame(
            {
                "qseqid": ["prot1"],
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
                "database": ["protein"],
            }
        )

    monkeypatch.setattr("biocol.cli.commands.run_blast", fake_run_blast)
    output = tmp_path / "out.tsv"
    code = main(
        [
            "--no-color",
            "run",
            "--query",
            str(fixtures_dir / "protein.fa"),
            "--db",
            str(fixtures_dir / "protein.fa"),
            "--accessions",
            str(fixtures_dir / "accessions.txt"),
            "--output",
            str(output),
            "--reciprocal",
        ]
    )
    assert code == 0
    assert captured.get("reciprocal") is True
    assert output.exists()


def test_run_reciprocal_defaults_to_false() -> None:
    args = build_parser().parse_args(
        [
            "run",
            "--query",
            "q.fa",
            "--db",
            "db.fa",
            "--accessions",
            "acc.txt",
        ]
    )
    assert args.reciprocal is False
    assert args.diamond is False
    assert args.kofam_profile is None
    assert args.ko_list is None


def test_run_passes_diamond_to_run_blast(
    tmp_path: Path, fixtures_dir: Path, monkeypatch
) -> None:
    captured: dict = {}

    def fake_run_blast(*_args, **kwargs):
        captured.update(kwargs)
        return pd.DataFrame(
            {
                "qseqid": ["prot1"],
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
                "database": ["protein"],
            }
        )

    monkeypatch.setattr("biocol.cli.commands.run_blast", fake_run_blast)
    output = tmp_path / "out.tsv"
    code = main(
        [
            "--no-color",
            "run",
            "--query",
            str(fixtures_dir / "protein.fa"),
            "--db",
            str(fixtures_dir / "protein.fa"),
            "--accessions",
            str(fixtures_dir / "accessions.txt"),
            "--output",
            str(output),
            "--diamond",
        ]
    )
    assert code == 0
    assert captured.get("diamond") is True
    assert captured.get("blast_dir") is None
    assert captured.get("diamond_dir") == tmp_path / "diamond"
    assert output.exists()


def _fake_blast_hits() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "qseqid": ["prot1"],
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
            "database": ["protein"],
        }
    )


def test_run_kofam_requires_profile_and_ko_list(
    tmp_path: Path, fixtures_dir: Path, monkeypatch
) -> None:
    monkeypatch.setattr("biocol.cli.commands.run_blast", lambda *_a, **_k: _fake_blast_hits())
    code = main(
        [
            "--no-color",
            "run",
            "--query",
            str(fixtures_dir / "protein.fa"),
            "--db",
            str(fixtures_dir / "protein.fa"),
            "--accessions",
            str(fixtures_dir / "accessions.txt"),
            "--output",
            str(tmp_path / "out.tsv"),
            "--kofam-profile",
            str(tmp_path / "eukaryote.hal"),
        ]
    )
    assert code == 1


def test_run_passes_kofam_to_run_kofamscan(
    tmp_path: Path, fixtures_dir: Path, monkeypatch
) -> None:
    captured: dict = {}
    monkeypatch.setattr("biocol.cli.commands.run_blast", lambda *_a, **_k: _fake_blast_hits())

    def fake_kofam(query, profile, ko_list, **kwargs):
        captured["query"] = query
        captured["profile"] = profile
        captured["ko_list"] = ko_list
        captured.update(kwargs)
        return pd.DataFrame(
            {
                "mark": ["*"],
                "gene_name": ["prot1"],
                "ko": ["K00003"],
                "threshold": [200.0],
                "score": [200.0],
                "evalue": [1.2e-10],
                "definition": ["homoserine dehydrogenase"],
            }
        )

    monkeypatch.setattr("biocol.cli.commands.run_kofamscan", fake_kofam)
    profile = tmp_path / "eukaryote.hal"
    profile.write_text("# list\n", encoding="utf-8")
    ko_list = tmp_path / "ko_list"
    ko_list.write_text("K00003\n", encoding="utf-8")
    output = tmp_path / "out.tsv"
    code = main(
        [
            "--no-color",
            "run",
            "--query",
            str(fixtures_dir / "protein.fa"),
            "--db",
            str(fixtures_dir / "protein.fa"),
            "--accessions",
            str(fixtures_dir / "accessions.txt"),
            "--output",
            str(output),
            "--kofam-profile",
            str(profile),
            "--ko-list",
            str(ko_list),
        ]
    )
    assert code == 0
    assert captured["profile"] == str(profile)
    assert captured["ko_list"] == str(ko_list)
    assert captured["kofam_dir"] == tmp_path / "kofam"
    table = pd.read_csv(output, header=None, sep="\t")
    assert "KOfamScan" in set(table.iloc[0].astype(str))
    assert "KO" in set(table.iloc[2].astype(str))


def test_from_blast_kofam_without_protein_returns_error(
    tmp_path: Path, fixtures_dir: Path
) -> None:
    code = main(
        [
            "--no-color",
            "from-blast",
            "--blast",
            str(fixtures_dir / "blast_outfmt6.txt"),
            "--accessions",
            str(fixtures_dir / "accessions.txt"),
            "--kofam-profile",
            str(tmp_path / "eukaryote.hal"),
            "--ko-list",
            str(tmp_path / "ko_list"),
            "--output",
            str(tmp_path / "out.tsv"),
        ]
    )
    assert code == 1
