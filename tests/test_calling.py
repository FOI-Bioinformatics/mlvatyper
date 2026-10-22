import pytest

from mlvatyper.calling import (
    allele_table, call_all, call_amplicon, derive_definitions, long_table, profile_table, read_bins,
    round_repeats, status_table, tandem_count_table,
)
from mlvatyper.models import Amplicon, Hit, MarkerDef

UNIT = "GATCGA"
FL = "ACGTACGGTACCAGTTGACA"
FR = "TGAAGGATGCTAGCCTGCAA"


def amp(sample, marker, copy, repeats, extra="", n_alt=0, has_n=False):
    seq = FL + UNIT * repeats + extra + FR
    if has_n:
        seq = seq.replace(UNIT, "GATNGA", 1)
    return Amplicon(Hit(sample, "c", 1, len(seq), "+", marker, n_alternatives=n_alt), copy, seq)


DEF = MarkerDef("m1", (UNIT,), 6, ref_size=len(FL + FR) + 5 * 6, ref_repeats=5)


@pytest.mark.parametrize("raw,policy,expected", [
    (5.0, "nearest", ("5", False)),
    (5.2, "nearest", ("5", False)),
    (4.8, "nearest", ("5", False)),
    (5.5, "nearest", ("6", True)),
    (5.25, "intermediate", ("5.5", True)),  # exactly at the tolerance is intermediate, as in MLVA_finder
    (4.75, "nearest", ("5", True)),
    (5.4, "up", ("6", True)),
    (5.6, "down", ("5", True)),
    (5.4, "intermediate", ("5.5", True)),
    (0.3, "intermediate", ("0.5", True)),
])
def test_round_repeats(raw, policy, expected):
    assert round_repeats(raw, 0.25, policy) == expected


def test_round_repeats_unknown_policy():
    with pytest.raises(ValueError):
        round_repeats(1.4, 0.25, "bogus")


def test_call_amplicon_ok_and_intermediate():
    c = call_amplicon(amp("s", "m1", 1, 7), DEF, 0.25, "intermediate")
    assert c.allele == "7" and c.raw == 7 and c.residual == 0 and c.status_text == "OK"
    assert c.array.copies == 7
    c = call_amplicon(amp("s", "m1", 1, 7, extra="ACG"), DEF, 0.25, "intermediate")
    assert c.allele == "7.5" and "INTERMEDIATE" in c.status and c.residual == 0.0
    c = call_amplicon(amp("s", "m1", 1, 7, extra="A"), DEF, 0.25, "nearest")
    assert c.allele == "7" and c.residual == 1 and c.status_text == "OK"


def test_call_amplicon_flags():
    c = call_amplicon(amp("s", "m1", 1, 7, n_alt=2, has_n=True), DEF, 0.25, "nearest")
    assert {"HAS_N", "AMBIGUOUS"} <= set(c.status) and c.allele == "7"
    # flanking deletion of one unit length: length says 6 but the array has 7 copies
    deleted = Amplicon(Hit("s", "c", 1, 1, "+", "m1"), 1, FL[6:] + UNIT * 7 + FR)
    c = call_amplicon(deleted, DEF, 0.25, "nearest")
    assert c.allele == "6" and "COUNT_DISAGREES" in c.status and c.array.copies == 7
    c = call_amplicon(amp("s", "m1", 1, 7), DEF, 0.25, "nearest", max_array_length=40)
    assert c.allele == ">=7" and "CENSORED" in c.status
    no_def = call_amplicon(amp("s", "m2", 1, 3), MarkerDef("m2", (UNIT,), 6), 0.25, "nearest")
    assert no_def.allele is None and no_def.status == ["NO_DEFINITION"] and no_def.array.copies == 3
    short = Amplicon(Hit("s", "c", 1, 1, "+", "m1"), 1, FL[:10])
    c = call_amplicon(short, DEF, 0.25, "down")
    assert "NEGATIVE" in c.status


def test_bins(tmp_path):
    p = tmp_path / "bins.tsv"
    p.write_text("marker\tsize\tallele\nm1\t82\t7\nm1\t60-71\t5\n")
    bins = read_bins(p)
    assert len(bins["m1"]) == 2
    c = call_amplicon(amp("s", "m1", 1, 7), DEF, 0.25, "nearest", bins=bins["m1"])
    assert c.allele == "7" and "BINNED" in c.status
    c = call_amplicon(amp("s", "m1", 1, 9), DEF, 0.25, "nearest", bins=bins["m1"])
    assert c.allele == "9" and "OFF_TABLE" in c.status


def test_call_all_tables():
    amps = [amp("s1", "m1", 1, 5), amp("s2", "m1", 1, 7), amp("s2", "m1", 2, 2)]
    defs = [DEF, MarkerDef("m2", (), 4, 50, 2)]
    calls = call_all(amps, ["s1", "s2", "s3"], defs, 0.25, "nearest", near_contig_end={("s3", "m1")})
    samples, markers = ["s1", "s2", "s3"], ["m1", "m2"]
    at = allele_table(calls, samples, markers)
    assert at.loc["s1", "m1"] == "5" and at.loc["s2", "m1"] == "7/2"
    assert at.loc["s3", "m1"] == "NEAR_CONTIG_END" and at.loc["s1", "m2"] == "NOT_DETECTED"
    st = status_table(calls, samples, markers)
    assert st.loc["s2", "m1"] == "MULTIPLE" and st.loc["s1", "m1"] == "OK"
    tc = tandem_count_table(calls, samples, markers)
    assert tc.loc["s2", "m1"] == "7;2" and tc.isna().loc["s3", "m1"]
    pt = profile_table(calls, samples, markers)
    assert pt.loc["s1", "profile"] == "5-NOT_DETECTED" and list(pt.columns) == ["profile", "m1", "m2"]
    lt = long_table(calls)
    assert len(lt) == 3 + 1 + 3  # three amplicons, s3/m1, and m2 for three samples
    row = lt[(lt["sample"] == "s1") & (lt["marker"] == "m1")].iloc[0]
    assert row["allele"] == "5" and row["tandem_copies"] == 5 and row["status"] == "OK"


def test_derive_definitions():
    amps = [amp("ref", "m1", 1, 5), amp("ref", "m1", 2, 6), amp("ref", "m3", 1, 4), amp("x", "m1", 1, 9)]
    defs = [MarkerDef("m1", (UNIT,), 6), MarkerDef("m2", (UNIT,), 6), MarkerDef("m3", (), 6),
            MarkerDef("m4", (UNIT,), 6, 10, 1)]
    derived = {d.primer_id: d for d in derive_definitions(amps, defs, "ref")}
    assert derived["m1"].ref_size == len(FL + FR) + 30 and derived["m1"].ref_repeats == 5
    assert not derived["m2"].can_call
    assert derived["m3"].ref_repeats == 0 and derived["m3"].ref_size == len(FL + FR) + 24
    assert derived["m4"] == defs[3]
    with pytest.raises(ValueError):
        derive_definitions(amps, defs, "absent")
