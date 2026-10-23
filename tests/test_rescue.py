import pandas as pd
from Bio.Seq import Seq

from mlvatyper.cli import main
from mlvatyper.models import PrimerPair
from mlvatyper.rescue import align_primer, find_with_indels, missing_reasons, read_missing, write_missing
from mlvatyper.sts import has_hash_word
from tests.conftest import FWD, REV, UNIT, make_amplicon, random_dna


def test_has_hash_word():
    assert has_hash_word("ACGTACGTACGTA", 11) and not has_hash_word("ACGTACGTNACGTAC", 11)
    assert has_hash_word("ACGTACGTNACGTAC", 5) and not has_hash_word("", 1)


def test_align_primer_with_and_without_indel():
    primer = "ACGTTGCAGGATCCA"
    window = random_dna(30, 81) + primer + random_dna(30, 82)
    assert align_primer(primer, window) == (0, 0, 30, 45)
    deleted = random_dna(30, 81) + primer[:7] + primer[8:] + random_dna(30, 82)
    assert align_primer(primer, deleted) == (0, 1, 30, 44)
    inserted = random_dna(30, 81) + primer[:7] + "T" + primer[7:] + random_dna(30, 82)
    assert align_primer(primer, inserted)[:2] == (0, 1)
    assert align_primer("ACGTTGCAGGRTCCA", window)[:2] == (0, 0)  # IUPAC code matches
    substituted = random_dna(30, 81) + primer[:7] + "G" + primer[8:] + random_dna(30, 82)
    assert align_primer(primer, substituted) == (1, 0, 30, 45)


def _genome_with_indel():
    """VNTR1 amplicon whose reverse primer site lacks one base; the STS primer therefore has an insertion."""
    amp = make_amplicon(5)
    rev_site = amp[-len(REV):]
    mutated = amp[:-len(REV)] + rev_site[:10] + rev_site[11:]
    return random_dna(300, 91) + mutated + random_dna(300, 92), len(mutated)


def test_find_with_indels_on_both_strands():
    genome, size = _genome_with_indel()
    pair = PrimerPair("VNTR1", FWD, REV, len(make_amplicon(5)))
    hits = find_with_indels({"c": genome}, pair, "s", margin=50)
    assert len(hits) == 1
    h = hits[0]
    assert (h.strand, h.length, h.fwd_indels, h.rev_indels, h.search_pass) == ("+", size, 0, 1, "aligned")
    assert h.start == 301 and h.fwd_mismatches == 0 and h.rev_mismatches == 0
    rc = str(Seq(genome).reverse_complement())
    hits_rc = find_with_indels({"c": rc}, pair, "s", margin=50)
    assert len(hits_rc) == 1 and hits_rc[0].strand == "-" and hits_rc[0].length == size
    assert hits_rc[0].end == len(genome) - 300 and hits_rc[0].rev_indels == 1
    assert find_with_indels({"c": genome}, pair, "s", margin=50, max_mismatch=0) == hits


def test_missing_reasons(tmp_path):
    genome = random_dna(200, 93) + FWD + random_dna(100, 94)
    pairs = [
        PrimerPair("OnePrimer", FWD, REV, 70),
        PrimerPair("Absent", "GGGGGGGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCCCC", 100),
        PrimerPair("Unhashable", "ACGTNACGTNACGTNACGTN", "CCCCCCCCCCCCCCCCCCCC", 100),
        PrimerPair("NearEnd", "TTTTTTTTTTTTTTTTTTTT", "CCCCCCCCCCCCCCCCCCCC", 100),
    ]
    reasons = missing_reasons({"c": genome}, pairs, {"NearEnd"}, wordsize=8, one_sided=True)
    assert reasons == {"OnePrimer": "ONE_PRIMER_ONLY", "Absent": "NOT_FOUND", "NearEnd": "NEAR_CONTIG_END",
                       "Unhashable": "NOT_FOUND"}  # one sided: the hashable reverse primer is enough
    assert missing_reasons({"c": genome}, pairs[2:3], set(), 8, one_sided=False) == {"Unhashable": "NOT_SEARCHED"}
    write_missing({("s", "Absent"): "NOT_FOUND"}, tmp_path / "m.tsv")
    assert read_missing(tmp_path / "m.tsv") == {("s", "Absent"): "NOT_FOUND"}
    assert read_missing(tmp_path / "none.tsv") == {}


def test_run_with_indel_rescue_reports_primer_indel(tmp_path):
    genome, size = _genome_with_indel()
    (tmp_path / "s.fasta").write_text(">c\n" + genome + "\n")
    expected = len(make_amplicon(5))
    (tmp_path / "p.sts").write_text(f"VNTR1\t{FWD}\t{REV}\t{expected}\nABSENT\t{'G' * 20}\t{'C' * 20}\t100\n")
    (tmp_path / "m.tsv").write_text("marker\tpatterns\tunit_length\tref_size\tref_repeats\n"
                                    f"VNTR1\t{UNIT}\t6\t{expected}\t5\n")
    base = ["run", "--sts", str(tmp_path / "p.sts"), "--fasta", str(tmp_path / "s.fasta"), "--prefix", "t",
            "--markers", str(tmp_path / "m.tsv"), "-M", "50", "--align", "none"]
    out = tmp_path / "without"
    assert main(base + ["--outdir", str(out)]) == 0
    assert pd.read_csv(out / "hits.tsv", sep="\t").empty
    missing = pd.read_csv(out / "missing_markers.tsv", sep="\t")
    assert dict(zip(missing.marker, missing.reason, strict=True)) == {"VNTR1": "ONE_PRIMER_ONLY",
                                                                      "ABSENT": "NOT_FOUND"}
    status = pd.read_csv(out / "t_status.tsv", sep="\t", index_col=0, dtype=str)
    assert status.loc["s", "VNTR1"] == "ONE_PRIMER_ONLY" and status.loc["s", "ABSENT"] == "NOT_DETECTED"

    out = tmp_path / "with"
    assert main(base + ["--outdir", str(out), "--rescue-indels"]) == 0
    hits = pd.read_csv(out / "hits.tsv", sep="\t")
    assert list(hits.search_pass) == ["aligned"] and list(hits.rev_indels) == [1] and hits.length[0] == size
    status = pd.read_csv(out / "t_status.tsv", sep="\t", index_col=0, dtype=str)
    alleles = pd.read_csv(out / "t_alleles.tsv", sep="\t", index_col=0, dtype=str)
    assert "PRIMER_INDEL" in status.loc["s", "VNTR1"] and alleles.loc["s", "VNTR1"] == "5"  # 1 bp short of 5 units
    long = pd.read_csv(out / "t_calls_long.tsv", sep="\t")
    row = long[long.marker == "VNTR1"].iloc[0]
    assert row["search_pass"] == "aligned" and row["residual_bp"] == -1 and row["rev_indels"] == 1


def test_fix_sts_strips_stray_characters(fixture_dir):
    out = fixture_dir / "out"
    assert main(["search", "--sts", str(fixture_dir / "primers.sts"), "--fasta", str(fixture_dir / "sampleA.fasta"),
                 "--outdir", str(out), "--fix-sts", "--align", "none", "-M", "50"]) == 0
    assert "BAD1\t" + FWD in (out / "primers_used.sts").read_text()
    assert not (out / "invalid_primers.tsv").exists()
