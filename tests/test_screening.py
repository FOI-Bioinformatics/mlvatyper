"""The across-sample summary of a screen."""

import pandas as pd

from mlvatyper.report import screening_overview
from mlvatyper.screening import MIN_SAMPLES, locus_summary, overview, sample_summary


def _tables():
    samples = [f"s{i}" for i in range(6)]
    alleles = pd.DataFrame({
        "good": ["5", "5", "6", "6", "7", "7"],
        "flagged": ["4", "4", "4", "5", "4", "4"],
        "negative": ["-1", "2", "3", "3", "2", "2"],
        "gone": ["NOT_DETECTED", "NOT_DETECTED", "NEAR_CONTIG_END", "3", "ONE_PRIMER_ONLY", "NOT_DETECTED"],
    }, index=samples)
    status = pd.DataFrame({
        "good": ["OK"] * 6,
        "flagged": ["COUNT_DISAGREES"] * 5 + ["OK"],
        "negative": ["NEGATIVE", "OK", "OK", "OK", "OK", "OK"],
        "gone": ["NOT_DETECTED", "NOT_DETECTED", "NEAR_CONTIG_END", "OK", "ONE_PRIMER_ONLY", "NOT_DETECTED"],
    }, index=samples)
    return alleles, status


def test_sample_summary_counts():
    alleles, status = _tables()
    s = sample_summary(alleles, status)
    assert s.loc["s0", ["typed", "ok", "flagged", "missing"]].tolist() == [3, 1, 2, 1]
    assert s.loc["s3", ["typed", "ok", "flagged", "missing"]].tolist() == [4, 3, 1, 0]
    assert s.loc["s3", "ok_fraction"] == 0.75 and s.loc["s0", "ok_fraction"] == 0.25


def test_locus_summary_finds_the_loci_to_review():
    alleles, status = _tables()
    loci = locus_summary(alleles, status)
    assert loci.loc["good", "concern"] == "" and loci.loc["good", "distinct_alleles"] == 3
    assert loci.loc["flagged", "concern"] == "COUNT_DISAGREES in 83% of calls"
    assert loci.loc["negative", "concern"] == "allele -1" and loci.loc["negative", "lowest_allele"] == -1
    assert loci.loc["gone", "typed"] == 1 and loci.loc["gone", "missing"] == 5


def test_overview_lists_low_call_rate_samples_and_missing_reasons():
    alleles, status = _tables()
    missing = pd.DataFrame({"sample": ["s0", "s1", "s2"], "marker": ["gone"] * 3,
                            "reason": ["NOT_DETECTED", "NOT_DETECTED", "NEAR_CONTIG_END"]})
    ov = overview(alleles, status, missing)
    # a share of OK loci below 75%: s0, s1, s2 and s4 (s3 and s5 are at 75% exactly)
    assert set(ov["low_call_rate"]) == {"s0", "s1", "s2", "s4"}
    assert ov["missing_by_reason"].set_index("reason")["products"].to_dict() == {
        "NOT_DETECTED": 2, "NEAR_CONTIG_END": 1}
    assert list(ov["concerns"].index) == ["flagged", "negative"]


def test_small_sets_get_no_overview():
    class Data:
        samples = ["a", "b"]
    assert len(Data.samples) < MIN_SAMPLES and screening_overview(Data()) is None


def test_binned_counts_as_a_good_call():
    """A locus called from a bin table reports BINNED for every normal result."""
    samples = [f"s{i}" for i in range(5)]
    alleles = pd.DataFrame({"a": ["3", "4", "3", "5", "4"], "b": ["7", "7", "8", "8", "8"]}, index=samples)
    status = pd.DataFrame({"a": ["BINNED"] * 5, "b": ["BINNED", "BINNED", "BINNED", "OFF_TABLE", "BINNED;MULTIPLE"]},
                          index=samples)
    s = sample_summary(alleles, status)
    assert s["ok"].tolist() == [2, 2, 2, 1, 1]
    loci = locus_summary(alleles, status)
    assert loci.loc["a", "ok"] == 5 and loci.loc["a", "concern"] == ""
    assert loci.loc["b", "ok"] == 3


def test_the_summaries_and_the_figures_agree_on_what_is_typed():
    """One parser for allele text: a half or censored allele counts as typed in both."""
    from mlvatyper.plots import allele_matrix
    alleles = pd.DataFrame({"a": ["20.5"] * 5, "b": [">=50"] * 5, "c": ["9/9"] * 5,
                            "d": ["NOT_DETECTED"] * 5}, index=[f"s{i}" for i in range(5)])
    status = pd.DataFrame("OK", index=alleles.index, columns=alleles.columns)
    status["d"] = "NOT_DETECTED"
    typed_here = sample_summary(alleles, status)["typed"].tolist()
    import numpy as np
    typed_there = (~np.isnan(allele_matrix(alleles))).sum(axis=1).tolist()
    assert typed_here == typed_there == [2] * 5


def test_a_locus_the_scheme_explains_is_listed_apart_from_the_ones_to_review():
    samples = [f"s{i}" for i in range(6)]
    alleles = pd.DataFrame({"known": ["2.5"] * 6, "new": ["3.5"] * 6, "fine": [str(i) for i in range(6)]},
                           index=samples)
    status = pd.DataFrame({"known": ["INTERMEDIATE"] * 6, "new": ["INTERMEDIATE"] * 6, "fine": ["OK"] * 6},
                          index=samples)
    plain = overview(alleles, status, None)
    assert sorted(plain["concerns"].index) == ["known", "new"] and plain["explained"].empty
    marked = overview(alleles, status, None, ("known",))
    assert list(marked["concerns"].index) == ["new"], "only the undeclared locus needs review"
    assert list(marked["explained"].index) == ["known"], "the declared one is reported apart"
    # the finding itself is kept either way, so nothing is hidden
    assert marked["loci"].loc["known", "concern"] == plain["loci"].loc["known", "concern"] != ""
    assert marked["loci"].loc["known", "explained"] and not marked["loci"].loc["new", "explained"]
