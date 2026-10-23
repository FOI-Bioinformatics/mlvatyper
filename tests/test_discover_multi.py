from Bio.Seq import Seq

from mlvatyper.cli import main
from mlvatyper.discover import canonical_unit, discover_units
from mlvatyper.discover_multi import apply_support, explained_fraction, pool_units, reference_candidates
from mlvatyper.models import Amplicon, Hit
from tests.conftest import FWD, REV, UNIT, random_dna

RC_REV = str(Seq(REV).reverse_complement())
CHANCE = "ACGGTTCAGATCCTGAAGTCCATGGTA"  # 27 bp chance tandem present twice in every sample


def amplicon(sample, marker, seq, copy=1):
    return Amplicon(Hit(sample, "c", 1, 1, "+", marker), copy, seq)


def with_chance(repeats):
    return FWD + CHANCE * 2 + random_dna(20, 61) + UNIT * repeats + RC_REV


def test_explained_fraction():
    assert explained_fraction({"a": 100, "b": 106, "c": 118}, 6, 0.25) == 1.0
    assert explained_fraction({"a": 100, "b": 106, "c": 118}, 27, 0.25) == 0.0
    assert explained_fraction({"a": 100, "b": 100}, 6, 0.25) is None
    assert explained_fraction({"a": 100, "b": 103}, 6, 0.25) == 0.0  # 3 bp is half a unit


def test_cross_sample_support_promotes_the_varying_unit_and_demotes_the_constant_one():
    amps = [amplicon("ref", "m", with_chance(2)), amplicon("s2", "m", with_chance(5)),
            amplicon("s3", "m", with_chance(7))]
    by_sample = {s: discover_units(amps, s)[0] for s in ("ref", "s2", "s3")}
    assert [d.unit_length for d in by_sample["ref"]["m"]] == [27]  # two copies of a 6-mer score below the floor
    assert by_sample["s3"]["m"][0].unit_length == 6  # 7 copies of 6 bp are category A, the 2 x 27 bp tandem is B
    lengths = {"m": {a.hit.sample: len(a.sequence) for a in amps}}
    pooled = pool_units(by_sample, lengths)
    unit6 = pooled["m"][canonical_unit(UNIT)]
    assert unit6.samples == ("s2", "s3") and (unit6.copy_min, unit6.copy_max) == (5, 7) and unit6.explains == 1.0
    assert pooled["m"][canonical_unit(CHANCE)].explains < 0.5  # 30 bp of 18, 30 and 12 is near one 27 bp unit
    supported = apply_support(by_sample, pooled)
    assert supported["s3"]["m"][0].unit_length == 6 and supported["s3"]["m"][0].category == "A"
    chance = supported["ref"]["m"][0]
    assert chance.unit_length == 27 and chance.category == "C"  # demoted from B
    found, missing = reference_candidates(supported, "ref", amps)
    top = found["m"][0]
    assert top.unit_length == 6 and top.copies == 2 and top.engine == "pooled" and top.category == "A"
    assert (top.copy_min, top.copy_max, top.n_samples) == (5, 7, 2) and missing == []


def test_reference_without_the_array_gets_the_pooled_unit_counted_directly():
    amps = [amplicon("ref", "m", FWD + random_dna(30, 71) + UNIT + random_dna(30, 72) + RC_REV),
            amplicon("s2", "m", FWD + random_dna(30, 71) + UNIT * 4 + random_dna(30, 72) + RC_REV),
            amplicon("s3", "m", FWD + random_dna(30, 71) + UNIT * 6 + random_dna(30, 72) + RC_REV)]
    by_sample = {s: discover_units(amps, s)[0] for s in ("ref", "s2", "s3")}
    assert "m" not in by_sample["ref"]
    lengths = {"m": {a.hit.sample: len(a.sequence) for a in amps}}
    supported = apply_support(by_sample, pool_units(by_sample, lengths))
    found, missing = reference_candidates(supported, "ref", amps)
    d = found["m"][0]
    assert d.engine == "pooled" and d.copies == 1 and d.category == "A" and d.sample == "ref" and missing == []


def test_discover_cli_all_samples(fixture_dir):
    out = fixture_dir / "out"
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "*.fasta"),
                 "--outdir", str(out), "-N", "0", "-M", "50", "--skip-invalid", "--align", "none"]) == 0
    assert main(["discover", "--outdir", str(out), "--all-samples", "--engine", "builtin"]) == 0
    report = (out / "discover_report.tsv").read_text().splitlines()
    assert report[1].startswith("VNTR1\tsampleA\t70\t6\t") and report[1].rstrip("\n").endswith("\t2\t5-7\t1.00")
    cands = (out / "discover_candidates.tsv").read_text()
    assert "VNTR1\tsampleB\t" in cands and "VNTR1\tsampleA\t" in cands
    assert main(["discover", "--outdir", str(out), "--sample", "sampleB", "--sample", "sampleA",
                 "--reference-sample", "sampleA", "--engine", "builtin"]) == 0
    assert "VNTR1\t6\t70\t5\t\t" in (out / "marker_defs_discovered.tsv").read_text()
