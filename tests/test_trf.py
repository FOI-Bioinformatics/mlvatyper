from pathlib import Path

import pytest

from mlvatyper.cli import main
from mlvatyper.models import Amplicon, Hit
from mlvatyper.trf import (
    TrfRepeat, discover_with_trf, parse_trf_ngs, select_repeat, trf_available, trf_command,
)
from tests.conftest import UNIT, make_amplicon, random_dna

DATA = Path(__file__).parent / "data"


def test_parse_trf_ngs_fixture():
    repeats = parse_trf_ngs((DATA / "trf_ngs.txt").read_text())
    assert len(repeats) == 6
    m1 = repeats["M1_9bp_131bp_4U"]
    assert [r.period for r in m1] == [9, 5, 4]
    assert m1[0] == TrfRepeat(46, 82, 9, 4.1, 9, 100, 0, 74, 1.98, "GATCGATTC", "GATCGATTCGATCGATTCGATCGATTCGATCGATTCG")
    assert repeats["M6_22bp_249bp_7U"][0].copies == 7.0


def test_parse_trf_ngs_errors():
    with pytest.raises(ValueError):
        parse_trf_ngs("1 2 3\n")
    with pytest.raises(ValueError):
        parse_trf_ngs("@x\n1 2 3\n")
    assert parse_trf_ngs("@x\n\n") == {"x": []}


def test_select_repeat_prefers_smallest_period_within_slack():
    repeats = parse_trf_ngs((DATA / "trf_ngs.txt").read_text())["M1_9bp_131bp_4U"]
    assert select_repeat(repeats).period == 9  # 4-mer and 5-mer score below 85 percent of 74
    assert select_repeat(repeats, period_slack=0.5).period == 4
    assert select_repeat([]) is None


def test_trf_command_pins_parameters(tmp_path):
    cmd = trf_command(tmp_path / "a.fasta", 30, 100)
    assert cmd[2:7] == ["2", "7", "7", "80", "10"] and cmd[7:9] == ["30", "100"] and cmd[-2:] == ["-h", "-ngs"]


@pytest.mark.skipif(not trf_available(), reason="trf not installed")
def test_discover_with_trf_on_synthetic():
    amps = [
        Amplicon(Hit("ref", "c", 1, 1, "+", "m1"), 1, make_amplicon(5)),
        Amplicon(Hit("ref", "c", 1, 1, "+", "m2"), 1, random_dna(200, 31)),
    ]
    found, missing, command = discover_with_trf(amps, "ref", min_score=20)
    assert list(found) == ["m1"] and missing == ["m2"]
    d = found["m1"][0]
    assert d.engine == "trf" and d.copies == 5 and d.unit_length == 6 and d.engine_copies >= 4.9
    assert d.primary and d.category == "A"
    assert set(d.unit) <= set(UNIT) and command[0] == "trf" and command[1] == "amplicons.fasta"
    assert command[7] == "20"  # TRF runs at the weak score floor


@pytest.mark.skipif(not trf_available(), reason="trf not installed")
def test_discover_cli_trf_engine(fixture_dir):
    out = fixture_dir / "out"
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "sampleA.fasta"),
                 "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"]) == 0
    assert main(["discover", "--outdir", str(out), "--engine", "trf", "--min-score", "20"]) == 0
    report = (out / "discover_report.tsv").read_text()
    assert "\ttrf\t" in report and "\tA\t1\t0\t" in report
    assert "VNTR1\t6\t70\t5\t\t" in (out / "marker_defs_discovered.tsv").read_text()
    assert "trf_command" in (out / "run_info.json").read_text()


def test_discover_cli_trf_missing(fixture_dir, monkeypatch):
    monkeypatch.setenv("PATH", "")
    out = fixture_dir / "out"
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "sampleA.fasta"),
                 "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"]) == 0
    with pytest.raises(FileNotFoundError):
        main(["discover", "--outdir", str(out), "--engine", "trf"])
    assert main(["discover", "--outdir", str(out)]) == 0  # auto falls back to builtin
    assert "\tbuiltin\t" in (out / "discover_report.tsv").read_text()
