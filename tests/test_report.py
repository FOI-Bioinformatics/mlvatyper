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
    assert 'class="legend"' in html and "Table 1." in html and "Figure 3." in html and "Software" in html
    assert "markers as rows" not in html  # two samples, two markers: no transposition
    info = json.loads((out / "run_info.json").read_text())
    assert len(info["inputs"]) == 4 and all(len(i["md5"]) == 32 and i["size"] > 0 for i in info["inputs"])
    assert any(i["path"].endswith("primers.sts") for i in info["inputs"])
    assert any(i["path"].endswith("m.tsv") for i in info["inputs"])  # the marker file is an input as well


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


def test_command_line_uses_the_negated_flags():
    from mlvatyper.report import _command_line
    on = _command_line({"search": {"parameters": {"dedup": True, "iupac": True, "fasta": ["a.fa"]}}})
    assert on == "mlvatyper run --fasta a.fa"
    off = _command_line({"search": {"parameters": {"dedup": False, "iupac": False, "report": False}}})
    assert off == "mlvatyper run --no-dedup --no-iupac --no-report"


def test_report_with_amplicons_of_undefined_markers(fixture_dir):
    """hits.tsv holds markers that the definitions used by summarize do not list."""
    out = _run(fixture_dir)
    (out / "primers_used.sts").unlink()
    (fixture_dir / "sub.tsv").write_text(f"marker\tpatterns\tunit_length\tref_size\tref_repeats\n"
                                         f"VNTR1\t{UNIT}\t6\t70\t5\n")
    hits = (out / "hits.tsv").read_text().splitlines()
    extra = hits[1].replace("\tVNTR1\t", "\tOTHER\t", 1)
    (out / "hits.tsv").write_text("\n".join([*hits, extra]) + "\n")
    assert main(["summarize", "--outdir", str(out), "--prefix", "u", "--markers", str(fixture_dir / "sub.tsv")]) == 0
    assert main(["report", "--outdir", str(out), "--prefix", "u"]) == 0
    assert (out / "u_report.html").exists()


def test_summary_counts_no_definition_as_found_and_methods_window(tmp_path):
    import pandas as pd

    from mlvatyper.report import ReportData, _summary_table
    alleles = pd.DataFrame({"m1": ["5"], "m2": ["NO_DEFINITION"], "m3": ["NOT_DETECTED"]}, index=["s"])
    status = pd.DataFrame({"m1": ["OK"], "m2": ["NO_DEFINITION"], "m3": ["NOT_DETECTED"]}, index=["s"])
    data = ReportData(tmp_path, "p", "t", {}, ["s"], ["m1", "m2", "m3"], alleles, status, alleles, pd.DataFrame(),
                      pd.DataFrame(), None, None)
    row = _summary_table(data).loc["s"]
    assert (row["markers typed"], row["missing"]) == (1, 1)
    text = methods_text({"search": {"parameters": {"markers": "m.tsv", "max_repeats": 40, "margin": 500}}})
    assert "0 to 40 repeats" in text
    assert "repeats instead" not in methods_text({"search": {"parameters": {"margin": 500}}})


def test_report_layout_conventions(fixture_dir):
    """Captions above tables and in sequence, short codes with a key, no favicon request, print styles."""
    out = _run(fixture_dir)
    html = (out / "t_report.html").read_text()
    assert '<link rel="icon" href="data:,">' in html and "@page" in html
    numbers = [int(n) for n in re.findall(r'<p class="tcap"><b>Table (\d+)\.', html)]
    assert numbers == list(range(1, len(numbers) + 1)) and len(numbers) >= 4
    # a table caption comes before its table, a figure caption after its figure
    assert html.index("Table 1.") < html.index("<table", html.index("Table 1."))
    assert html.index("<svg") < html.index("<figcaption>")
    assert 'scope="col"' in html and 'scope="row"' in html
    assert "In place of an allele:" in html and "ND</b> NOT_DETECTED" in html
    # the missing marker shows its short code and keeps the full one in the title
    assert '<td title="NOT_DETECTED"' in html or 'title="NOT_DETECTED">ND<' in html


def test_report_names_the_scheme_and_its_publication(tmp_path):
    from mlvatyper.report import methods_text, scheme_info
    run_info = {"search": {"parameters": {"scheme": "burkholderia_uren_2007", "mismatches": 2,
                                          "rescue_wordsize": 6, "rescue_indels": True}}}
    info = scheme_info(run_info)
    assert info["doi"] == "10.1186/1471-2180-7-23" and info["organism"] == "Burkholderia pseudomallei"
    assert len(info["interrupted"]) == 17
    text = methods_text(run_info)
    assert "burkholderia_uren_2007 scheme for Burkholderia pseudomallei" in text and "BMC Microbiol 2007;7:23" in text
    assert "up to 2 mismatches" in text and "17 markers whose arrays hold variant units" in text
    assert scheme_info({"search": {"parameters": {}}}) is None
    assert scheme_info({"search": {"parameters": {"scheme": "no-such-scheme"}}}) is None


def test_a_set_of_samples_gets_the_relationship_figures(tmp_path):
    import pandas as pd

    from mlvatyper.report import ReportData, figures
    markers = [f"m{i}" for i in range(10)]
    samples = [f"s{i}" for i in range(6)]
    alleles = pd.DataFrame([[str(3 + (i + j) % 4) for j in range(10)] for i in range(6)], index=samples,
                           columns=markers)
    status = pd.DataFrame("OK", index=samples, columns=markers)
    calls = pd.DataFrame(columns=["sample", "marker", "copy", "length", "expected_size", "residual_bp"])
    data = ReportData(tmp_path, "p", "t", {}, samples, markers, alleles, status, alleles, calls,
                      pd.DataFrame(columns=["sample", "marker", "reason"]), None, None)
    figs = figures(data)
    assert {"heatmap", "distance", "spectra", "callrate", "status", "mst", "mds"} <= set(figs)
    few = ReportData(tmp_path, "p", "t", {}, samples[:2], markers, alleles.iloc[:2], status.iloc[:2], alleles.iloc[:2],
                     calls, pd.DataFrame(columns=["sample", "marker", "reason"]), None, None)
    assert "heatmap" not in figures(few)


def test_a_figure_that_fails_does_not_stop_the_report(tmp_path, monkeypatch, caplog):
    import pandas as pd

    from mlvatyper import plots
    from mlvatyper.report import ReportData, figures
    markers = [f"m{i}" for i in range(6)]
    samples = [f"s{i}" for i in range(6)]
    alleles = pd.DataFrame([[str(3 + (i * j) % 4) for j in range(6)] for i in range(6)], index=samples,
                           columns=markers)
    status = pd.DataFrame("OK", index=samples, columns=markers)
    calls = pd.DataFrame(columns=["sample", "marker", "copy", "length", "expected_size", "residual_bp"])
    data = ReportData(tmp_path, "p", "t", {}, samples, markers, alleles, status, alleles, calls,
                      pd.DataFrame(columns=["sample", "marker", "reason"]), None, None)

    def broken(*args, **kwargs):
        raise RuntimeError("no convergence")

    monkeypatch.setattr(plots, "mds_plot", broken)
    figs = figures(data)
    assert figs["mds"] == "" and figs["heatmap"].startswith("<svg") and figs["status"].startswith("<svg")
    assert "The mds figure was left out of the report" in caplog.text


def test_command_line_reproduces_the_run_and_parses():
    """Only run options are written; internal values and global options are left out."""
    import shlex

    from mlvatyper.cli import build_parser
    from mlvatyper.report import _command_line
    params = {"scheme": "burkholderia_uren_2007", "fasta": ["a.fa"], "outdir": "o", "prefix": "p",
              "interrupted_arrays": ["1788k", "857k"], "log_json": "events.jsonl", "rescue_indels": False,
              "dedup": False, "jobs": 2}
    cmd = _command_line({"search": {"parameters": params}})
    assert "--interrupted-arrays" not in cmd and "--log-json" not in cmd
    assert "--no-rescue-indels" in cmd and "--no-dedup" in cmd and "--jobs 2" in cmd
    args = build_parser().parse_args(shlex.split(cmd)[1:])
    assert args.command == "run" and args.rescue_indels is False and args.dedup is False


def _report_data(tmp_path, alleles, status):
    import pandas as pd

    from mlvatyper.report import ReportData
    return ReportData(tmp_path, "p", "t", {}, list(alleles.index), list(alleles.columns), alleles, status,
                      alleles, pd.DataFrame(), pd.DataFrame(), None, None)


def test_summary_rows_stay_with_their_samples(tmp_path):
    """A sample without an OK call sorted last in the status counts, so typed/missing landed on another row."""
    import pandas as pd

    from mlvatyper.report import _summary_table
    alleles = pd.DataFrame({"m1": ["5", "NOT_DETECTED", "6"], "m2": ["4", "NOT_DETECTED", "4"]},
                           index=["a_mini", "b_empty", "c_mini"])
    status = alleles.where(alleles != "NOT_DETECTED", "NOT_DETECTED").mask(alleles != "NOT_DETECTED", "OK")
    table = _summary_table(_report_data(tmp_path, alleles, status))
    assert list(table.index) == ["a_mini", "b_empty", "c_mini"]
    assert table.loc["b_empty", ["markers typed", "missing"]].tolist() == [0, 2]
    assert table.loc["c_mini", ["markers typed", "missing"]].tolist() == [2, 0]


def test_summary_counts_an_exact_half_grid_allele_as_typed(tmp_path):
    import pandas as pd

    from mlvatyper.report import _summary_table
    alleles = pd.DataFrame({"m1": ["6.5"], "m2": ["3.1"], "m3": [">=50"], "m4": ["9/9"]}, index=["s"])
    status = pd.DataFrame({"m1": ["OK"], "m2": ["OK"], "m3": ["CENSORED"], "m4": ["MULTIPLE"]}, index=["s"])
    assert _summary_table(_report_data(tmp_path, alleles, status)).loc["s", "markers typed"] == 3


def test_one_locus_panel_still_gets_the_spectra_and_call_rates(tmp_path):
    import pandas as pd

    from mlvatyper.report import figures
    samples = [f"s{i}" for i in range(6)]
    alleles = pd.DataFrame({"m1": [str(2 + i % 3) for i in range(6)]}, index=samples)
    status = pd.DataFrame({"m1": ["OK"] * 6}, index=samples)
    data = _report_data(tmp_path, alleles, status)
    data.calls = pd.DataFrame({"sample": samples, "marker": ["m1"] * 6, "copy": ["1"] * 6, "length": ["100"] * 6,
                               "expected_size": ["100"] * 6, "residual_bp": ["0"] * 6, "allele": alleles["m1"],
                               "status": ["OK"] * 6})
    out = figures(data)
    assert out.get("spectra") and out.get("callrate"), "no clustering is needed for these two"
    assert "heatmap" not in out  # one locus has nothing to cluster


def test_summary_tables_of_an_earlier_larger_run_are_removed(fixture_dir):
    for i in range(5):
        (fixture_dir / f"extra{i}.fasta").write_text((fixture_dir / "sampleA.fasta").read_text())
    out = fixture_dir / "out"
    base = ["run", "--sts", str(fixture_dir / "primers.sts"), "--outdir", str(out), "--prefix", "p", "-N", "0",
            "-M", "50", "--skip-invalid", "--align", "none", "--markers", str(fixture_dir / "patterns.tsv"),
            "--reference-sample", "sampleA"]
    assert main([*base, "--fasta", str(fixture_dir / "*.fasta")]) == 0
    assert (out / "p_sample_summary.tsv").exists() and (out / "p_locus_summary.tsv").exists()
    assert main([*base, "--fasta", str(fixture_dir / "sampleA.fasta"), str(fixture_dir / "sampleB.fasta")]) == 0
    assert not (out / "p_sample_summary.tsv").exists() and not (out / "p_locus_summary.tsv").exists()
