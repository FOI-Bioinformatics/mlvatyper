from vntrfinder.models import Amplicon, Hit, PatternDef
from vntrfinder.summarize import (
    count_repeats, length_diff_from_shortest_table, lengths_table, pattern_counts_table,
)


def amp(sample, marker, copy, seq):
    return Amplicon(Hit(sample, "c", 1, len(seq), "+", marker), copy, seq)


AMPS = [
    amp("s1", "m1", 1, "GATCGA" * 5),
    amp("s2", "m1", 1, "GATCGA" * 7),
    amp("s2", "m1", 2, "GATCGA" * 2 + "GATNGA"),
    amp("s1", "m2", 1, "ACGT"),
]
SAMPLES = ["s1", "s2", "s3"]
DEFS = [PatternDef("m1", ("GATCGA",), 6), PatternDef("m2")]


def test_count_repeats_non_overlapping_case_insensitive():
    assert count_repeats("gatcgaGATCGAxx", ("GATCGA",)) == 2
    assert count_repeats("AAAA", ("AA",)) == 2
    assert count_repeats("ACGT", ("TTT", "ACG")) == 1


def test_lengths_table():
    df = lengths_table(AMPS, SAMPLES)
    assert list(df.index) == SAMPLES and list(df.columns) == ["m1", "m2"]
    assert df.loc["s1", "m1"] == "30"
    assert df.loc["s2", "m1"] == "42/18"
    assert df.loc["s1", "m2"] == "4"
    assert df.isna().loc["s3", "m1"] and df.isna().loc["s2", "m2"]


def test_lengths_table_clean_drops_N():
    df = lengths_table(AMPS, SAMPLES, clean=True)
    assert df.loc["s2", "m1"] == "42"


def test_pattern_counts_table():
    df = pattern_counts_table(AMPS, SAMPLES, DEFS)
    assert list(df.columns) == ["m1"]
    assert df.loc["s1", "m1"] == "5"
    assert df.loc["s2", "m1"] == "7;2"
    assert df.isna().loc["s3", "m1"]


def test_length_diff_from_shortest():
    df = length_diff_from_shortest_table(AMPS, SAMPLES, DEFS)
    assert df.loc["s1", "m1"] == "2.000"   # (30-18)/6
    assert df.loc["s2", "m1"] == "4.000;0.000"
    df_clean = length_diff_from_shortest_table(AMPS, SAMPLES, DEFS, clean=True)
    assert df_clean.loc["s1", "m1"] == "0.000" and df_clean.loc["s2", "m1"] == "2.000"
