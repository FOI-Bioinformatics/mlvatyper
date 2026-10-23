from Bio.Seq import Seq

from mlvatyper.cli import main
from mlvatyper.discover import (
    Discovery, RawArray, Thresholds, array_score, builtin_arrays, canonical_rotation, canonical_unit, categorize,
    discover_units, discoveries_to_definitions, finish_candidates, is_mixed, make_discovery, merge_harmonics,
    rank_candidates,
    unit_entropy, write_candidates, write_report,
)
from mlvatyper.models import Amplicon, Hit
from mlvatyper.tandem import count_tandem_array, rotations
from tests.conftest import FWD, REV, UNIT, make_amplicon, random_dna

RC_REV = str(Seq(REV).reverse_complement())


def best_unit(seq, **kw):
    found, _ = discover_units([Amplicon(Hit("s", "c", 1, 1, "+", "m"), 1, seq)], "s", **kw)
    return found["m"][0] if "m" in found else None


def test_builtin_recovers_planted_unit_in_both_orientations():
    seq = make_amplicon(5)
    d = best_unit(seq)
    assert d.unit == canonical_rotation(UNIT) and d.copies == 5 and d.score == 2 * (30 - 6)
    rc = str(Seq(seq).reverse_complement())
    d_rc = best_unit(rc)
    assert d_rc.unit in rotations(str(Seq(UNIT).reverse_complement())) and d_rc.copies == 5
    assert d.unit_key == d_rc.unit_key == canonical_unit(UNIT)


def test_builtin_prefers_smallest_period_and_tolerates_variant():
    seq = random_dna(40, 21) + "GATCGATTC" * 3 + "GATCGTTTC" + "GATCGATTC" * 2 + random_dna(40, 22)
    d = best_unit(seq)
    assert d.unit_length == 9 and d.copies == 6 and d.array.variant == 1
    assert d.score == array_score(seq, d.array, 9)


def test_builtin_rejects_random_and_short_arrays():
    assert best_unit(random_dna(300, 23)) is None
    weak = random_dna(100, 24) + "ACGTAC" * 2 + random_dna(100, 25)  # score 12
    assert best_unit(weak) is None
    assert best_unit(weak, thresholds=Thresholds(weak_score=12)).category == "C"


def test_unit_entropy():
    assert unit_entropy("AAAA") == 0
    assert abs(unit_entropy("ACGT") - 2.0) < 1e-9
    assert 0.9 < unit_entropy("TAT") < 0.93


def test_merge_harmonics_keeps_smallest_period_within_slack():
    a6 = RawArray(20, 50, 6, "GATCGA", 48, 5.0)
    a12 = RawArray(20, 50, 12, "GATCGAGATCGA", 44, 2.5)
    other = RawArray(100, 130, 7, "ACGTTTG", 40, 4.3)
    merged = merge_harmonics([a12, a6, other])
    assert [m.period for m in merged] == [6, 7]
    assert merge_harmonics([a6, RawArray(20, 50, 3, "GAT", 10, 10.0)]) == [a6]  # 3-mer scores below slack


def test_two_arrays_in_one_amplicon_are_both_reported_and_marked_mixed():
    unit27 = "ACGGTTCAGATCCTGAAGTCCATGGTA"
    seq = FWD + random_dna(30, 41) + UNIT * 4 + random_dna(30, 42) + unit27 * 3 + random_dna(30, 43) + RC_REV
    amps = [Amplicon(Hit("s", "c", 1, 1, "+", "m"), 1, seq)]
    found, missing = discover_units(amps, "s", primer_lengths={"m": (len(FWD), len(REV))})
    cands = found["m"]
    assert sorted(c.unit_length for c in cands) == [6, 27] and missing == []
    assert cands[0].primary and cands[0].rank == 1 and not cands[1].primary and cands[1].rank == 2
    assert {c.category for c in cands} == {"A"} and is_mixed(cands)
    assert cands[0].unit_length == 6  # more integer copies wins within category A


def test_categories_at_the_boundaries():
    t = Thresholds()
    base = dict(marker="m", sample="s", amplicon_length=300, left_flank=50, right_flank=50)

    def disc(unit, copies, score, overlap=False):
        return Discovery(unit=unit, array=count_tandem_array(unit * copies, (unit,)), score=score,
                         primer_overlap=overlap, **base)

    assert categorize(disc("GATCGA", 3, 30), t) == "A"
    assert categorize(disc("GATCGA", 2, 30), t) == "B"
    assert categorize(disc("GATCGA", 3, 29), t) == "C"
    assert categorize(disc("GATC", 5, 40), t) == "C"  # short unit
    assert categorize(disc("AAAAAT", 4, 40), t) == "B"  # entropy 0.65 below 1.0
    assert categorize(disc("GATCGA", 4, 40, overlap=True), t) == "C"
    assert categorize(disc("GATCGA", 2, 30), Thresholds(strong_copies=2)) == "A"


def test_rank_candidates_orders_by_category_copies_score():
    base = dict(marker="m", sample="s", amplicon_length=300, left_flank=50, right_flank=50)
    a = Discovery(unit="GATCGA", array=count_tandem_array("GATCGA" * 3, ("GATCGA",)), score=30, category="B", **base)
    b = Discovery(unit="ACGTTTG", array=count_tandem_array("ACGTTTG" * 4, ("ACGTTTG",)), score=28, category="C",
                  **base)
    c = Discovery(unit="TTGACCA", array=count_tandem_array("TTGACCA" * 3, ("TTGACCA",)), score=40, category="B",
                  **base)
    ranked = rank_candidates([b, a, c])
    assert [d.unit for d in ranked] == ["TTGACCA", "GATCGA", "ACGTTTG"] and [d.rank for d in ranked] == [1, 2, 3]
    assert ranked[0].primary and not ranked[1].primary


def test_partial_second_copy_is_kept_as_category_b():
    unit = "ACGGTTCAGATC"
    seq = FWD + random_dna(40, 51) + unit * 2 + unit[:6] + "GGGGGG" + random_dna(40, 52) + RC_REV
    raw = RawArray(60, 60 + 30, 12, unit, 40, 2.5)
    d = make_discovery("m", "s", seq, raw, "trf")
    assert d.copies == 2 and categorize(d, Thresholds()) == "B"
    assert make_discovery("m", "s", seq, RawArray(0, 0, 12, "GGGGGGGGGGGG", 1, 1.0), "trf") is None


def test_primer_overlap_needs_a_whole_unit_inside_the_primer():
    seq = make_amplicon(5)
    raw = builtin_arrays(seq)[0]
    d = make_discovery("m", "s", seq, raw, "builtin", primer_lengths=(len(FWD), len(REV)))
    assert d.left_flank in (19, 20) and not d.primer_overlap  # phase may start one base early
    assert make_discovery("m", "s", seq, raw, "builtin", primer_lengths=(27, len(REV))).primer_overlap


def test_discover_units_and_definitions(tmp_path):
    amps = [
        Amplicon(Hit("ref", "c", 1, 1, "+", "m1"), 1, make_amplicon(5)),
        Amplicon(Hit("ref", "c", 1, 1, "+", "m1"), 2, make_amplicon(7)),
        Amplicon(Hit("ref", "c", 1, 1, "+", "m2"), 1, random_dna(200, 26)),
        Amplicon(Hit("other", "c", 1, 1, "+", "m3"), 1, make_amplicon(4)),
    ]
    found, missing = discover_units(amps, "ref")
    assert list(found) == ["m1"] and missing == ["m2"]
    d = found["m1"][0]
    assert d.copies == 5 and d.category == "A" and d.array.span == 30
    defs = {d.primer_id: d for d in discoveries_to_definitions(found, missing, amps, "ref")}
    assert defs["m1"].can_call and defs["m1"].ref_repeats == 5 and defs["m1"].ref_size == 70
    assert not defs["m2"].can_call and defs["m2"].ref_size == 200 and defs["m2"].unit_length == 0
    assert count_tandem_array(make_amplicon(7), defs["m1"].patterns).copies == 7
    # a category below --accept keeps only ref_size
    demoted = {m: [c.__class__(**{**c.__dict__, "category": "C"}) for c in cs] for m, cs in found.items()}
    assert all(not x.can_call for x in discoveries_to_definitions(demoted, missing, amps, "ref"))
    write_report(found, missing, tmp_path / "r.tsv")
    lines = (tmp_path / "r.tsv").read_text().splitlines()
    assert len(lines) == 3 and lines[1].startswith("m1\tref\t70\t6\t") and lines[2].startswith("m2\t")
    assert lines[0].split("\t")[-6:] == ["category", "n_candidates", "mixed", "n_samples", "copy_range", "explains"]
    write_candidates(found, tmp_path / "c.tsv")
    assert (tmp_path / "c.tsv").read_text().splitlines()[1].split("\t")[-7:-3] == ["0", "A", "1", "1"]


def test_discover_cli(fixture_dir):
    out = fixture_dir / "out"
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
                 "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"]) == 0
    assert main(["discover", "--outdir", str(out)]) == 1  # two samples, --sample required
    assert main(["discover", "--outdir", str(out), "--sample", "nobody"]) == 1
    assert main(["discover", "--outdir", str(out), "--sample", "sampleA", "--reference-sample", "sampleB"]) == 1
    assert main(["discover", "--outdir", str(out), "--sample", "sampleA"]) == 0
    text = (out / "marker_defs_discovered.tsv").read_text()
    assert "VNTR1\t6\t70\t5\t\t" in text
    assert (out / "discover_candidates.tsv").exists()
    assert main(["summarize", "--outdir", str(out), "--prefix", "t",
                 "--markers", str(out / "marker_defs_discovered.tsv")]) == 0
    assert "7/3" in (out / "t_alleles.tsv").read_text()
    assert main(["discover", "--outdir", str(fixture_dir / "empty")]) == 1


def test_mixed_marker_with_same_period_units_gets_both_patterns():
    unit_a, unit_b = "GATCGA", "ACGTTG"
    seq = FWD + random_dna(30, 44) + unit_a * 4 + random_dna(30, 45) + unit_b * 3 + random_dna(30, 46) + RC_REV
    amps = [Amplicon(Hit("s", "c", 1, 1, "+", "m"), 1, seq)]
    # the builtin engine keeps one array per period, so the two arrays are given as engine output here
    raw = [RawArray(50, 74, 6, unit_a, 36, 4.0), RawArray(104, 122, 6, unit_b, 24, 3.0)]
    cands = [make_discovery("m", "s", seq, r, "trf") for r in raw]
    found = {"m": finish_candidates(cands, Thresholds(min_score=20))}
    assert is_mixed(found["m"]) and [c.copies for c in found["m"]] == [4, 3]
    d = discoveries_to_definitions(found, [], amps, "s", Thresholds(min_score=20))[0]
    assert len(d.patterns) == 2 and set(d.patterns) == {canonical_rotation(unit_a), canonical_rotation(unit_b)}
    assert d.ref_repeats == 4  # the primary array alone anchors the reference count
    # a second unit below --accept, or of another period, is not combined
    strict = Thresholds(accept="A", min_score=20, strong_copies=4)  # the 3 copy array becomes B
    assert len(discoveries_to_definitions({"m": finish_candidates(cands, strict)}, [], amps, "s", strict)[0].patterns) == 1
    other = make_discovery("m", "s", seq + "ACGTTGA" * 3, RawArray(172, 193, 7, "ACGTTGA", 28, 3.0), "trf")
    mixed = {"m": finish_candidates([cands[0], other], Thresholds(min_score=20))}
    assert len(discoveries_to_definitions(mixed, [], amps, "s", Thresholds(min_score=20))[0].patterns) == 1
