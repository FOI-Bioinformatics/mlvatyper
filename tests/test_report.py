import json
import re

from mlvatyper.cli import main
from mlvatyper.report import STATUS_GLOSSARY, cell_colour, load_report_data, methods_text, render, status_counts
from tests.conftest import UNIT


def _run(fixture_dir, extra=()):
    out = fixture_dir / "out"
    (fixture_dir / "m.tsv").write_text(f"marker\tpatterns\tunit_length\nVNTR1\t{UNIT}\t6\n")
    assert main(["run", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
                 "--outdir", str(out), "--prefix", "t", "-N", "0", "-M", "50", "--skip-invalid", "--align", "none",
                 "--markers", str(fixture_dir / "m.tsv"), "--reference-sample", "sampleA", *extra]) == 0
    return out


def test_run_writes_report_with_all_sections(fixture_dir):
    out = _run(fixture_dir)
    html = (out / "t_report.html").read_text()
    for heading in ("Summary", "Allele profiles", "Figures", "Evidence per amplicon", "Markers without amplicon",
                    "Methods", "Appendix", "Input files"):
        assert heading in html
    for marker in ("VNTR1", "ABSENT"):
        assert marker in html
    assert "sampleA" in html and "sampleB" in html
    assert html.count("<svg") == 3
    # no external resources are loaded
    assert not re.search(r'(src|href)="https?://', html)
    assert "merPCR" in html and "mlvatyper" in html
    assert "NOT_DETECTED" in html and STATUS_GLOSSARY["NOT_DETECTED"] in html
    info = json.loads((out / "run_info.json").read_text())
    assert len(info["inputs"]) == 3 and all(len(i["md5"]) == 32 and i["size"] > 0 for i in info["inputs"])
    assert any(i["path"].endswith("primers.sts") for i in info["inputs"])


def test_report_subcommand_and_no_report(fixture_dir):
    out = _run(fixture_dir, ["--no-report"])
    assert not (out / "t_report.html").exists()
    assert main(["report", "--outdir", str(out), "--prefix", "t", "--title", "My run"]) == 0
    html = (out / "t_report.html").read_text()
    assert "<title>My run</title>" in html
    assert main(["report", "--outdir", str(out), "--prefix", "nothing"]) == 1


def test_report_helpers(fixture_dir):
    out = _run(fixture_dir)
    data = load_report_data(out, "t")
    counts = status_counts(data)
    assert counts.loc["sampleA", "OK"] == 1 and counts.loc["sampleA", "NOT_DETECTED"] == 1
    assert counts.loc["sampleB", "MULTIPLE"] == 1  # one marker with two amplicons on two contigs
    assert cell_colour("OK") == "#ffffff" and cell_colour("OK;MULTIPLE") != "#ffffff"
    text = methods_text(data.run_info)
    assert "0 mismatches per primer" in text and "sample sampleA" in text and "intermediate policy" in text
    assert "searched again with a hash word of 6 bases" in text
    assert render(data).startswith("<!DOCTYPE html>")
    assert methods_text({}).startswith("VNTR amplicons were located")
