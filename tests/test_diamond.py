from pathlib import Path

import pytest

from biocol import DatabaseError, DiamondError, DiamondExecutionError, run_blast
from biocol.diamond.runner import build_blastp_command, build_makedb_command


def test_blastp_command_uses_diamond_defaults_except_threads(tmp_path: Path) -> None:
    command = build_blastp_command(
        tmp_path / "q.faa",
        tmp_path / "db",
        tmp_path / "out.txt",
        num_threads=4,
    )
    assert command[:2] == ["diamond", "blastp"]
    assert "--evalue" not in command
    assert "--max-target-seqs" not in command
    assert "-k" not in command
    assert command[command.index("--threads") + 1] == "4"
    assert "nident" in command
    assert "qseq" in command


def test_makedb_command(tmp_path: Path) -> None:
    command = build_makedb_command(tmp_path / "db.faa", tmp_path / "db")
    assert command[:2] == ["diamond", "makedb"]
    assert command[command.index("--in") + 1].endswith("db.faa")


def test_diamond_requires_executable(
    fixtures_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("biocol.diamond.execute.shutil.which", lambda _name: None)
    with pytest.raises(DiamondExecutionError, match="was not found"):
        run_blast(
            fixtures_dir / "protein.fa",
            fixtures_dir / "protein.fa",
            diamond=True,
        )


def test_diamond_requires_protein_query(fixtures_dir: Path) -> None:
    with pytest.raises(DiamondError, match="protein query and a protein database"):
        run_blast(
            fixtures_dir / "dna.fa",
            fixtures_dir / "protein.fa",
            diamond=True,
        )


def test_diamond_requires_protein_database(fixtures_dir: Path) -> None:
    with pytest.raises(DiamondError, match="protein query and a protein database"):
        run_blast(
            fixtures_dir / "protein.fa",
            fixtures_dir / "dna.fa",
            diamond=True,
        )


def test_diamond_rejects_mixed_fasta_and_dmnd(tmp_path: Path, fixtures_dir: Path) -> None:
    folder = tmp_path / "mix"
    folder.mkdir()
    (folder / "a.faa").write_text(
        (fixtures_dir / "protein.fa").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    (folder / "b.dmnd").write_bytes(b"")
    with pytest.raises(DatabaseError, match="mixes FASTA and Diamond"):
        run_blast(fixtures_dir / "protein.fa", folder, diamond=True)


def _hsp_line(query: str = "prot1", subject: str = "XP_001") -> str:
    return (
        f"{query}\t{subject}\t99.0\t10\t0\t0\t1\t10\t1\t10\t1e-5\t50.0"
        "\t10\tAAAAAAAAAA\tAAAAAAAAAA\n"
    )


def test_diamond_makedb_then_blastp_from_fasta(
    tmp_path: Path, fixtures_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("biocol.diamond.execute.shutil.which", lambda name: f"/usr/bin/{name}")
    commands: list[list[str]] = []

    def fake_run(command: list[str]) -> None:
        commands.append(command)
        if command[1] == "makedb":
            return
        Path(command[command.index("--out") + 1]).write_text(_hsp_line(), encoding="utf-8")

    monkeypatch.setattr("biocol.diamond.runner.run_diamond_command", fake_run)
    dest = tmp_path / "diamond_hits"
    hits = run_blast(
        fixtures_dir / "protein.fa",
        fixtures_dir / "protein.fa",
        diamond=True,
        diamond_dir=dest,
        num_threads=2,
    )
    assert commands[0][:2] == ["diamond", "makedb"]
    assert commands[1][:2] == ["diamond", "blastp"]
    assert "--evalue" not in commands[1]
    assert commands[1][commands[1].index("--threads") + 1] == "2"
    assert list(dest.glob("*.txt"))
    assert not hits.empty


def test_diamond_uses_existing_dmnd(
    tmp_path: Path, fixtures_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("biocol.diamond.execute.shutil.which", lambda name: f"/usr/bin/{name}")
    dmnd = tmp_path / "lyrata.dmnd"
    dmnd.write_bytes(b"")
    commands: list[list[str]] = []

    def fake_run(command: list[str]) -> None:
        commands.append(command)
        Path(command[command.index("--out") + 1]).write_text(_hsp_line(), encoding="utf-8")

    monkeypatch.setattr("biocol.diamond.runner.run_diamond_command", fake_run)
    dest = tmp_path / "diamond_hits"
    run_blast(
        fixtures_dir / "protein.fa",
        dmnd,
        diamond=True,
        diamond_dir=dest,
    )
    assert all(cmd[1] != "makedb" for cmd in commands)
    assert commands[0][:2] == ["diamond", "blastp"]
    assert commands[0][commands[0].index("--db") + 1] == str(dmnd)


def test_diamond_reciprocal_requires_fasta_not_dmnd(
    tmp_path: Path, fixtures_dir: Path
) -> None:
    dmnd = tmp_path / "lyrata.dmnd"
    dmnd.write_bytes(b"")
    with pytest.raises(DiamondError, match="not .dmnd"):
        run_blast(
            fixtures_dir / "protein.fa",
            dmnd,
            diamond=True,
            reciprocal=True,
        )


def test_diamond_reciprocal_blastp_both_ways(
    tmp_path: Path, fixtures_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("biocol.diamond.execute.shutil.which", lambda name: f"/usr/bin/{name}")
    commands: list[list[str]] = []

    def fake_run(command: list[str]) -> None:
        commands.append(command)
        if command[1] == "makedb":
            return
        query = Path(command[command.index("--query") + 1])
        out = Path(command[command.index("--out") + 1])
        if query.name == "protein.fa":
            out.write_text(_hsp_line("prot1", "prot1"), encoding="utf-8")
        else:
            out.write_text(_hsp_line("prot1", "prot1"), encoding="utf-8")

    monkeypatch.setattr("biocol.diamond.runner.run_diamond_command", fake_run)
    dest = tmp_path / "diamond_hits"
    hits = run_blast(
        fixtures_dir / "protein.fa",
        fixtures_dir / "protein.fa",
        diamond=True,
        diamond_dir=dest,
        reciprocal=True,
    )
    blastp = [cmd for cmd in commands if cmd[1] == "blastp"]
    makedb = [cmd for cmd in commands if cmd[1] == "makedb"]
    assert len(makedb) == 2
    assert len(blastp) == 2
    names = sorted(path.name for path in dest.glob("*.txt"))
    assert any(name.startswith("reverse_") for name in names)
    assert not hits[hits["sseqid"].notna()].empty
