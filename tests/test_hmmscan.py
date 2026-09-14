from pathlib import Path

import pandas as pd
import pytest

from biocol.exceptions import HmmError
from biocol.hmm.parser import parse_hmmscan_tblout
from biocol.hmm.press import hmm_has_gathering_threshold
from biocol.hmm.runner import build_hmmscan_command, run_hmmscan


def test_parse_hmmscan_tblout_keeps_all_hits(fixtures_dir: Path) -> None:
    hits = parse_hmmscan_tblout(fixtures_dir / "hmmscan_tblout.txt")
    assert len(hits) == 3
    names = set(hits.loc[hits["query_name"] == "jg63.t1", "target_name"])
    assert names == {"Pkinase", "Pkinase_Tyr"}


def test_hmm_has_gathering_threshold(tmp_path: Path) -> None:
    hmm = tmp_path / "db.hmm"
    hmm.write_text("HMMER3/f\nNAME  x\nGA    21.10 21.10;\n", encoding="utf-8")
    assert hmm_has_gathering_threshold(hmm) is True
    hmm.write_text("HMMER3/f\nNAME  x\n", encoding="utf-8")
    assert hmm_has_gathering_threshold(hmm) is False


def test_hmmscan_command_cut_ga_or_evalue(tmp_path: Path) -> None:
    with_ga = build_hmmscan_command(
        tmp_path / "q.faa",
        tmp_path / "db.hmm",
        tmp_path / "out.tbl",
        use_cut_ga=True,
    )
    assert "--cut_ga" in with_ga
    without = build_hmmscan_command(
        tmp_path / "q.faa",
        tmp_path / "db.hmm",
        tmp_path / "out.tbl",
        use_cut_ga=False,
        evalue=10,
    )
    assert "-E" in without
    assert "10" in without


def test_run_hmmscan_rejects_nucleotide_query(fixtures_dir: Path) -> None:
    with pytest.raises(HmmError, match="protein sequence"):
        run_hmmscan(fixtures_dir / "dna.fa", fixtures_dir / "protein.fa")
