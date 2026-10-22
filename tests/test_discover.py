from Bio.Seq import Seq

from mlvatyper.cli import main
from mlvatyper.discover import (
    array_score, discover_units, discoveries_to_definitions, find_repeat_unit, unit_entropy, write_report,
)
from mlvatyper.models import Amplicon, Hit
from mlvatyper.tandem import count_tandem_array, rotations
from tests.conftest import UNIT, make_amplicon, random_dna


def test_find_repeat_unit_recovers_planted_unit():
    seq = make_amplicon(5)
    unit, array, score = find_repeat_unit(seq)
    assert unit in rotations(UNIT) and array.copies == 5 and score == 2 * (30 - 6)
    rc = str(Seq(seq).reverse_complement())
    unit_rc, array_rc, _ = find_repeat_unit(rc)
    assert unit_rc in rotations(str(Seq(UNIT).reverse_complement())) and array_rc.copies == 5


def test_find_repeat_unit_prefers_smallest_period_and_tolerates_variant():
    seq = random_dna(40, 21) + "GATCGATTC" * 3 + "GATCGTTTC" + "GATCGATTC" * 2 + random_dna(40, 22)
    unit, array, score = find_repeat_unit(seq)
    assert len(unit) == 9 and array.copies == 6 and array.variant == 1
    assert score == array_score(seq, array, 9)


def test_find_repeat_unit_rejects_random_and_short_arrays():
    assert find_repeat_unit(random_dna(300, 23)) is None
    assert find_repeat_unit(random_dna(100, 24) + "ACGTAC" * 2 + random_dna(100, 25)) is None  # score 12
    assert find_repeat_unit(random_dna(100, 24) + "ACGTAC" * 2 + random_dna(100, 25), min_score=12) is not None


def test_unit_entropy():
    assert unit_entropy("AAAA") == 0
    assert abs(unit_entropy("ACGT") - 2.0) < 1e-9
    assert 0.9 < unit_entropy("TAT") < 0.93


def test_discover_units_and_definitions(tmp_path):
    amps = [
        Amplicon(Hit("ref", "c", 1, 1, "+", "m1"), 1, make_amplicon(5)),
        Amplicon(Hit("ref", "c", 1, 1, "+", "m1"), 2, make_amplicon(7)),
        Amplicon(Hit("ref", "c", 1, 1, "+", "m2"), 1, random_dna(200, 26)),
        Amplicon(Hit("other", "c", 1, 1, "+", "m3"), 1, make_amplicon(4)),
    ]
    found, missing = discover_units(amps, "ref")
    assert [d.marker for d in found] == ["m1"] and missing == ["m2"]
    # the primer ends in A, so the array may start one base early in a rotated phase
    assert found[0].copies == 5 and found[0].left_flank in (19, 20) and found[0].array.span == 30
    defs = {d.primer_id: d for d in discoveries_to_definitions(found, missing, amps, "ref")}
    assert defs["m1"].can_call and defs["m1"].ref_repeats == 5 and defs["m1"].ref_size == 70
    assert not defs["m2"].can_call and defs["m2"].ref_size == 200 and defs["m2"].unit_length == 0
    assert count_tandem_array(make_amplicon(7), defs["m1"].patterns).copies == 7
    write_report(found, missing, tmp_path / "r.tsv")
    lines = (tmp_path / "r.tsv").read_text().splitlines()
    assert len(lines) == 3 and lines[1].startswith("m1\tref\t70\t6\t") and lines[2].startswith("m2\t")


def test_discover_cli(fixture_dir):
    out = fixture_dir / "out"
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
                 "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"]) == 0
    assert main(["discover", "--outdir", str(out)]) == 1  # two samples, --sample required
    assert main(["discover", "--outdir", str(out), "--sample", "nobody"]) == 1
    assert main(["discover", "--outdir", str(out), "--sample", "sampleA"]) == 0
    text = (out / "marker_defs_discovered.tsv").read_text()
    assert "VNTR1\t6\t70\t5\t\t" in text
    assert main(["summarize", "--outdir", str(out), "--prefix", "t",
                 "--markers", str(out / "marker_defs_discovered.tsv")]) == 0
    assert "7/3" in (out / "t_alleles.tsv").read_text()
    assert main(["discover", "--outdir", str(fixture_dir / "empty")]) == 1
