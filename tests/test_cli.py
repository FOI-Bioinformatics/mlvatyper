import json

import pandas as pd

from mlvatyper.cli import main
from tests.conftest import FWD, REV


def read(path):
    return pd.read_csv(path, sep="\t", index_col=0, keep_default_na=False, dtype=str)


def test_run_end_to_end_with_reference_sample(fixture_dir):
    out = fixture_dir / "out"
    rc = main([
        "run", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
        "--outdir", str(out), "--prefix", "t", "-N", "0", "-M", "50", "--skip-invalid",
        "--markers", str(fixture_dir / "patterns.tsv"), "--reference-sample", "sampleA",
        "--clwclean", "--align", "none", "--legacy-tables",
    ])
    assert rc == 0
    assert (out / "invalid_primers.tsv").read_text().count("\n") == 2
    hits = pd.read_csv(out / "hits.tsv", sep="\t")
    assert set(hits["sample"]) == {"sampleA", "sampleB"}
    assert set(hits["fwd_mismatches"]) == {0} and set(hits["rev_mismatches"]) == {0}
    lengths = read(out / "t_lengths_summary.tsv")
    assert lengths.loc["sampleA", "VNTR1"] == "70"
    assert lengths.loc["sampleB", "VNTR1"] == "82/58"

    derived = read(out / "marker_defs_derived.tsv")
    assert derived.loc["VNTR1", "ref_size"] == "70" and derived.loc["VNTR1", "ref_repeats"] == "5"

    alleles = read(out / "t_alleles.tsv")
    assert alleles.loc["sampleA", "VNTR1"] == "5"
    assert alleles.loc["sampleB", "VNTR1"] == "7/3"
    assert alleles.loc["sampleA", "ABSENT"] == "NOT_DETECTED"
    status = read(out / "t_status.tsv")
    assert status.loc["sampleA", "VNTR1"] == "OK"
    assert "MULTIPLE" in status.loc["sampleB", "VNTR1"] and "HAS_N" in status.loc["sampleB", "VNTR1"]
    counts = read(out / "t_tandem_counts.tsv")
    assert counts.loc["sampleB", "VNTR1"] == "7;3"  # the unit with N is one mismatch from the pattern
    profile = read(out / "t_profile.tsv")
    assert profile.loc["sampleA", "profile"] == "5-NOT_DETECTED"
    long = pd.read_csv(out / "t_calls_long.tsv", sep="\t")
    assert len(long) == 3 + 2  # three amplicons plus ABSENT for two samples
    assert set(long.columns) >= {"raw_repeats", "allele", "residual_bp", "tandem_copies", "status"}

    # legacy tables still available on request
    legacy = read(out / "t_count_summary.tsv")
    assert legacy.loc["sampleA", "VNTR1"] == "5" and legacy.loc["sampleB", "VNTR1"] == "7;2"
    diff = read(out / "t_count_lengthdiff.tsv")
    assert diff.loc["sampleA", "VNTR1"] == "2.000"
    clean = read(out / "t_lengths_summary_clean.tsv")
    assert clean.loc["sampleB", "VNTR1"] == "82"
    assert (out / "sequences" / "VNTR1.fasta").exists()

    info = json.loads((out / "run_info.json").read_text())
    assert info["mlvatyper_version"] and "search" in info and "summarize" in info
    assert info["summarize"]["parameters"]["rounding"] == "intermediate"


def test_run_without_definitions_skips_calling(fixture_dir):
    out = fixture_dir / "out2"
    rc = main(["run", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "sampleA.fasta"),
               "--outdir", str(out), "--prefix", "t", "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"])
    assert rc == 0
    assert (out / "t_lengths_summary.tsv").exists()
    assert not (out / "t_alleles.tsv").exists()
    assert not (out / "t_count_summary.tsv").exists()


def test_summarize_reuses_derived_definitions(fixture_dir):
    out = fixture_dir / "out3"
    base = ["--outdir", str(out), "--prefix", "t"]
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
                 "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"]) == 0
    assert main(["summarize", *base, "--markers", str(fixture_dir / "patterns.tsv"),
                 "--reference-sample", "sampleA"]) == 0
    assert main(["summarize", *base, "--markers", str(out / "marker_defs_derived.tsv"), "--rounding", "up"]) == 0
    alleles = read(out / "t_alleles.tsv")
    assert alleles.loc["sampleB", "VNTR1"] == "7/3"


def test_summarize_unknown_reference_sample_fails(fixture_dir):
    out = fixture_dir / "out4"
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "sampleA.fasta"),
                 "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"]) == 0
    assert main(["summarize", "--outdir", str(out), "--prefix", "t", "--reference-sample", "nope"]) == 1


def test_invalid_sts_aborts_without_skip(fixture_dir):
    rc = main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "sampleA.fasta"),
               "--outdir", str(fixture_dir / "o2"), "--align", "none"])
    assert rc == 1


def test_summarize_requires_hits(tmp_path):
    assert main(["summarize", "--outdir", str(tmp_path), "--prefix", "x"]) == 1


def test_cli_passes_one_sided_and_wordsize(tmp_path):
    from tests.test_search import write_one_sided_fixture
    write_one_sided_fixture(tmp_path)
    base = ["search", "--sts", str(tmp_path / "p.sts"), "--fasta", str(tmp_path / "s.fasta"),
            "-N", "0", "-M", "50", "--align", "none"]
    assert main(base + ["--outdir", str(tmp_path / "a"), "--rescue-wordsize", "0"]) == 0
    assert pd.read_csv(tmp_path / "a" / "hits.tsv", sep="\t").empty
    assert (tmp_path / "a" / "missing_markers.tsv").read_text().splitlines()[1] == "s\tOneSided\tNOT_SEARCHED\t"
    assert main(base + ["--outdir", str(tmp_path / "b"), "--one-sided", "--rescue-wordsize", "0"]) == 0
    assert len(pd.read_csv(tmp_path / "b" / "hits.tsv", sep="\t")) == 1
    assert main(base + ["--outdir", str(tmp_path / "c"), "-W", "8", "--rescue-wordsize", "0"]) == 0
    assert len(pd.read_csv(tmp_path / "c" / "hits.tsv", sep="\t")) == 1
    assert main(base + ["--outdir", str(tmp_path / "d")]) == 0  # the default rescue pass recovers it
    hits = pd.read_csv(tmp_path / "d" / "hits.tsv", sep="\t")
    assert len(hits) == 1 and hits.search_pass[0] == "rescue" and hits.fwd_mismatches[0] == 0


def test_search_markers_writes_size_windows(fixture_dir):
    (fixture_dir / "m.tsv").write_text("marker\tpatterns\tunit_length\tref_size\tref_repeats\n"
                                       "VNTR1\tGATCGA\t6\t70\t5\n")
    out = fixture_dir / "out"
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
                 "--outdir", str(out), "--skip-invalid", "--align", "none", "-M", "5",
                 "--markers", str(fixture_dir / "m.tsv"), "--max-repeats", "20"]) == 0
    sts = (out / "primers_used.sts").read_text().splitlines()
    assert sts[0].split("\t")[3] == "40-160" and sts[1].split("\t")[3] == "100"  # ABSENT keeps its size
    hits = pd.read_csv(out / "hits.tsv", sep="\t")
    assert sorted(hits.length) == [58, 70, 82]  # -M 5 alone would have missed the 7 and 3 repeat amplicons
    assert set(hits.expected_size) == {70}  # the STS size is kept, not the midpoint of the window
    sts = (fixture_dir / "range.sts")
    sts.write_text(f"VNTR1\t{FWD}\t{REV}\t60-80\n")
    assert main(["search", "--sts", str(sts), "--fasta", str(fixture_dir / "sampleA.fasta"),
                 "--outdir", str(fixture_dir / "r"), "--align", "none", "--markers", str(fixture_dir / "m.tsv")]) == 0
    assert (fixture_dir / "r" / "primers_used.sts").read_text().split("\t")[3].strip() == "60-80"  # explicit range kept


def test_non_fasta_input_fails_only_that_sample(fixture_dir, caplog):
    (fixture_dir / "notes.fasta").write_text("<!doctype html>\n<html></html>\n")
    out = fixture_dir / "out"
    rc = main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
               "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"])
    assert rc == 2 and "1 samples failed" in caplog.text
    hits = pd.read_csv(out / "hits.tsv", sep="\t")
    assert set(hits["sample"]) == {"sampleA", "sampleB"}  # the valid samples were still processed and written


def test_run_types_the_good_samples_when_one_fails(fixture_dir, caplog):
    """A failed sample does not stop allele calling for the others, and is not reported as NOT_DETECTED."""
    (fixture_dir / "corrupt.fasta.gz").write_bytes(b"not gzip at all")
    out = fixture_dir / "out"
    rc = main(["run", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
               str(fixture_dir / "*.fasta.gz"), "--outdir", str(out), "--prefix", "t", "-N", "0", "-M", "50",
               "--skip-invalid", "--align", "none", "--markers", str(fixture_dir / "patterns.tsv"),
               "--reference-sample", "sampleA"])
    assert rc == 2 and "corrupt" in caplog.text
    assert (out / "samples.txt").read_text().split() == ["sampleA", "sampleB"]
    alleles = read(out / "t_alleles.tsv")
    assert list(alleles.index) == ["sampleA", "sampleB"] and alleles.loc["sampleA", "VNTR1"] == "5"
    assert (out / "t_report.html").exists()


def test_missing_fasta_path_stops_the_run(fixture_dir, caplog):
    rc = main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "sampleA.fasta"),
               str(fixture_dir / "sampleX.fasta"), "--outdir", str(fixture_dir / "o"), "--align", "none"])
    assert rc == 1 and "sampleX.fasta not found" in caplog.text
    rc = main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "sampleA.fasta"),
               str(fixture_dir / "nothing*.fasta"), "--outdir", str(fixture_dir / "o"), "--align", "none",
               "--skip-invalid", "-N", "0", "-M", "50"])
    assert rc == 0 and "matches no file" in caplog.text


def test_sample_name_with_spaces(fixture_dir, caplog):
    from Bio import SeqIO
    spaced = fixture_dir / "my sample A.fasta"
    spaced.write_text((fixture_dir / "sampleA.fasta").read_text())
    out = fixture_dir / "o"
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(spaced),
                 "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"]) == 0
    assert (out / "samples.txt").read_text().split() == ["my_sample_A"] and "replaced by _" in caplog.text
    assert [r.id for r in SeqIO.parse(out / "sequences" / "VNTR1.fasta", "fasta")] == ["my_sample_A_1"]


def test_masked_fasta_gives_the_same_amplicons(fixture_dir):
    text = (fixture_dir / "sampleA.fasta").read_text().splitlines()
    (fixture_dir / "sampleA.fasta").write_text("\n".join([text[0], "XXXXXXXXXX" + text[1].lower()]) + "\n")
    out = fixture_dir / "o"
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "sampleA.fasta"),
                 "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"]) == 0
    hits = pd.read_csv(out / "hits.tsv", sep="\t")
    assert hits.loc[0, "sequence"].startswith(FWD) and hits.loc[0, "length"] == 70


def test_summarize_errors_and_columns(fixture_dir, caplog):
    out = fixture_dir / "out"
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
                 "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"]) == 0
    # a marker without an amplicon in any sample still has a lengths column
    assert main(["summarize", "--outdir", str(out), "--prefix", "a",
                 "--markers", str(fixture_dir / "patterns.tsv")]) == 0
    assert "ABSENT" in read(out / "a_lengths_summary.tsv").columns
    # a duplicated marker is an error with a message, not two columns
    dup = fixture_dir / "dup.tsv"
    dup.write_text("marker\tunit_length\tref_size\tref_repeats\nVNTR1\t6\t70\t5\nVNTR1\t6\t76\t5\n")
    assert main(["summarize", "--outdir", str(out), "--prefix", "b", "--markers", str(dup)]) == 1
    assert "already defined on line 2" in caplog.text
    # a reference sample without amplicons is an error with a message, not a traceback
    (out / "samples.txt").write_text("sampleA\nsampleB\nempty\n")
    assert main(["summarize", "--outdir", str(out), "--prefix", "c", "--markers", str(fixture_dir / "patterns.tsv"),
                 "--reference-sample", "empty"]) == 1
    assert "has no amplicons" in caplog.text


def test_summarize_prefixes_do_not_mix(fixture_dir):
    from mlvatyper.report import load_report_data, methods_text
    out = fixture_dir / "out"
    assert main(["run", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
                 "--outdir", str(out), "--prefix", "p", "-N", "0", "-M", "50", "--skip-invalid", "--align", "none",
                 "--markers", str(fixture_dir / "patterns.tsv"), "--reference-sample", "sampleA"]) == 0
    assert (out / "p_definitions.tsv").exists()
    assert main(["summarize", "--outdir", str(out), "--prefix", "q", "--markers", str(fixture_dir / "patterns.tsv"),
                 "--reference-sample", "sampleA", "--rounding", "down", "--allele-tolerance", "0.1"]) == 0
    text = methods_text(json.loads((out / "run_info.json").read_text()), "p")
    assert "0.25 of an integer" in text and "intermediate policy" in text
    # a later summarize of p that cannot call removes p's tables instead of leaving them from the earlier run
    assert main(["summarize", "--outdir", str(out), "--prefix", "p"]) == 0
    assert not (out / "p_alleles.tsv").exists() and not (out / "p_report.html").exists()
    assert load_report_data(out, "q").definitions is not None


def test_summarize_without_mlvabank_removes_the_old_export(fixture_dir):
    out = fixture_dir / "out"
    base = ["--outdir", str(out), "--prefix", "t", "--markers", str(fixture_dir / "patterns.tsv"),
            "--reference-sample", "sampleA"]
    assert main(["run", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"), "-N", "0",
                 "-M", "50", "--skip-invalid", "--align", "none", "--no-report", "--mlvabank", *base]) == 0
    assert (out / "t_mlvabank.csv").exists()
    assert main(["summarize", *base]) == 0
    assert not (out / "t_mlvabank.csv").exists()


def test_jobs_parallel_matches_serial(fixture_dir):
    """--jobs N over several samples gives byte-identical hits, samples and missing-marker files."""
    for name in ("sampleC", "sampleD"):  # duplicate the two fixture samples so there are four
        (fixture_dir / f"{name}.fasta").write_text((fixture_dir / "sampleA.fasta").read_text())
    base = ["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
            "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"]
    assert main([*base, "--outdir", str(fixture_dir / "serial")]) == 0
    assert main([*base, "--outdir", str(fixture_dir / "par"), "--jobs", "4"]) == 0
    for name in ("hits.tsv", "samples.txt", "missing_markers.tsv", "near_contig_end.tsv"):
        assert (fixture_dir / "serial" / name).read_text() == (fixture_dir / "par" / name).read_text(), name


def test_jobs_parallel_isolates_a_failed_sample(fixture_dir, caplog):
    (fixture_dir / "broken.fasta").write_text("not a fasta file\n")
    rc = main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
               "--outdir", str(fixture_dir / "o"), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none",
               "--jobs", "3"])
    assert rc == 2 and "broken" in caplog.text
    assert "broken" not in (fixture_dir / "o" / "samples.txt").read_text()
    assert {"sampleA", "sampleB"} <= set((fixture_dir / "o" / "samples.txt").read_text().split())


def test_log_json_writes_json_lines(fixture_dir):
    log = fixture_dir / "events.jsonl"
    rc = main(["--log-json", str(log), "search", "--sts", str(fixture_dir / "primers.sts"),
               "--fasta", str(fixture_dir / "sampleA.fasta"), "--outdir", str(fixture_dir / "o"),
               "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"])
    assert rc == 0
    lines = [json.loads(line) for line in log.read_text().splitlines()]
    assert lines and all({"time", "level", "logger", "message"} <= set(e) for e in lines)
    assert any(e["logger"] == "mlvatyper" and "valid primer pairs" in e["message"] for e in lines)
    assert all(e["time"].endswith("+00:00") for e in lines)  # UTC


def test_log_json_records_a_failed_sample(fixture_dir):
    (fixture_dir / "broken.fasta").write_text("not a fasta\n")
    log = fixture_dir / "events.jsonl"
    rc = main(["--log-json", str(log), "search", "--sts", str(fixture_dir / "primers.sts"),
               "--fasta", str(fixture_dir / "*.fasta"), "--outdir", str(fixture_dir / "o"),
               "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"])
    assert rc == 2
    events = [json.loads(line) for line in log.read_text().splitlines()]
    assert any(e["level"] == "ERROR" and "broken" in e["message"] for e in events)


def _search_args(fixture_dir, out, *extra):
    return ["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "sample*.fasta"),
            "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none", *extra]


def test_search_writes_a_result_per_sample(fixture_dir):
    out = fixture_dir / "o"
    assert main(_search_args(fixture_dir, out)) == 0
    for sample in ("sampleA", "sampleB"):
        assert (out / "per_sample" / f"{sample}.hits.tsv").exists()
        assert (out / "per_sample" / f"{sample}.missing.tsv").exists()
        assert (out / "per_sample" / f"{sample}.done.json").exists()


def test_resume_reuses_finished_samples_and_gives_identical_files(fixture_dir, caplog):
    first, second = fixture_dir / "first", fixture_dir / "second"
    assert main(_search_args(fixture_dir, first)) == 0
    # a second run into a copy of the per-sample results takes them from there
    import shutil
    shutil.copytree(first / "per_sample", second / "per_sample")
    caplog.clear()
    assert main(_search_args(fixture_dir, second, "--resume")) == 0
    assert "Resuming: 2 of 2 samples" in caplog.text
    for name in ("hits.tsv", "missing_markers.tsv", "samples.txt", "near_contig_end.tsv"):
        assert (first / name).read_text() == (second / name).read_text(), name


def test_resume_searches_again_what_changed(fixture_dir, caplog):
    out = fixture_dir / "o"
    assert main(_search_args(fixture_dir, out)) == 0
    # a changed input file is searched again, the other sample is reused
    sample = fixture_dir / "sampleB.fasta"
    sample.write_text(sample.read_text() + ">extra\nACGTACGTAC\n")
    caplog.clear()
    assert main(_search_args(fixture_dir, out, "--resume")) == 0
    assert "Resuming: 1 of 2 samples" in caplog.text
    # a changed search setting invalidates every sample
    caplog.clear()
    assert main(_search_args(fixture_dir, out, "--resume", "--mismatches", "1")) == 0
    assert "Resuming" not in caplog.text


def test_resume_after_a_failed_sample_repeats_only_that_sample(fixture_dir, caplog):
    (fixture_dir / "sampleC.fasta").write_text("not a fasta file\n")
    out = fixture_dir / "o"
    assert main(_search_args(fixture_dir, out)) == 2  # sampleC fails, the others are kept
    (fixture_dir / "sampleC.fasta").write_text((fixture_dir / "sampleA.fasta").read_text())
    caplog.clear()
    assert main(_search_args(fixture_dir, out, "--resume")) == 0
    assert "Resuming: 2 of 3 samples" in caplog.text
    assert set((out / "samples.txt").read_text().split()) == {"sampleA", "sampleB", "sampleC"}


def test_resume_works_with_parallel_jobs(fixture_dir, caplog):
    out = fixture_dir / "o"
    assert main(_search_args(fixture_dir, out, "--jobs", "2")) == 0
    first = (out / "hits.tsv").read_text()
    caplog.clear()
    assert main(_search_args(fixture_dir, out, "--jobs", "2", "--resume")) == 0
    assert "Resuming: 2 of 2 samples" in caplog.text and (out / "hits.tsv").read_text() == first


def test_the_search_engine_is_importable_without_the_command_line():
    """runner holds the search and the per-sample cache, so a pipeline can call it without argparse."""
    from mlvatyper import runner
    assert {"search_sample", "search_samples", "save_sample", "load_sample", "search_key",
            "PER_SAMPLE_DIR"} <= set(dir(runner))
    assert "argparse" not in runner.search_samples.__doc__


def test_run_where_every_sample_fails_exits_2_without_a_traceback(fixture_dir, caplog):
    (fixture_dir / "bad.fasta").write_text("not a fasta\n")
    out = fixture_dir / "o"
    rc = main(["run", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "bad.fasta"),
               "--outdir", str(out), "--prefix", "p", "--skip-invalid", "--align", "none"])
    assert rc == 2 and "samples failed" in caplog.text
    # summarize on the empty sample list does not crash either
    assert main(["summarize", "--outdir", str(out), "--prefix", "p", "--markers",
                 str(fixture_dir / "patterns.tsv")]) == 0
