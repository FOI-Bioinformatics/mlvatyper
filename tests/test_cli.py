import json

import pandas as pd

from mlvatyper.cli import main


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
    assert (tmp_path / "a" / "missing_markers.tsv").read_text().splitlines()[1] == "s\tOneSided\tNOT_SEARCHED"
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
    assert set(hits.expected_size) == {100}
