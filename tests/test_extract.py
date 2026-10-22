from Bio.Seq import Seq

from mlvatyper.extract import (
    deduplicate_hits, extract_amplicons, load_fasta, primer_mismatches, primers_near_contig_ends, read_hits_table,
    write_hits_table, write_marker_fastas,
)
from mlvatyper.models import Hit, PrimerPair
from tests.conftest import FWD, REV


def test_extract_forward_and_reverse(fixture_dir):
    contigs = load_fasta(fixture_dir / "sampleA.fasta")
    assert list(contigs) == ["contig1"]
    seq = contigs["contig1"]
    hits = [
        Hit("sampleA", "contig1", 301, 370, "+", "VNTR1"),
        Hit("sampleA", "contig1", 301, 370, "-", "VNTR1"),
        Hit("sampleA", "missing", 1, 10, "+", "VNTR1"),
    ]
    amps = extract_amplicons(hits, contigs)
    assert len(amps) == 2
    assert amps[0].sequence == seq[300:370]
    assert amps[1].sequence == str(Seq(seq[300:370]).reverse_complement())
    assert [a.copy for a in amps] == [1, 2]
    assert amps[0].record_id == "sampleA_1"
    assert amps[0].hit.dist_to_contig_end == 300


def test_primer_mismatches_iupac():
    assert primer_mismatches("ACGT", "ACGT") == 0
    assert primer_mismatches("ACGT", "ACGA") == 1
    assert primer_mismatches("ACGT", "RCGW") == 0
    assert primer_mismatches("NCGT", "ACGT") == 1
    assert primer_mismatches("ACG", "ACGT") == 4


def test_extract_annotates_mismatches(fixture_dir):
    contigs = load_fasta(fixture_dir / "sampleA.fasta")
    pair = PrimerPair("VNTR1", FWD, REV, 70)
    amps = extract_amplicons([Hit("sampleA", "contig1", 301, 370, "+", "VNTR1")], contigs, {"VNTR1": pair})
    h = amps[0].hit
    assert (h.fwd_mismatches, h.rev_mismatches, h.expected_size) == (0, 0, 70)
    assert (h.total_mismatches, h.size_deviation) == (0, 0)
    bad = PrimerPair("VNTR1", "T" + FWD[1:], REV[:-1] + "G", 60)
    amps = extract_amplicons([Hit("sampleA", "contig1", 301, 370, "+", "VNTR1")], contigs, {"VNTR1": bad})
    assert (amps[0].hit.fwd_mismatches, amps[0].hit.rev_mismatches, amps[0].hit.size_deviation) == (1, 1, 10)


def test_deduplicate_overlapping_hits():
    base = dict(sample="s", contig="c", strand="+", primer_id="m")
    a = Hit(start=100, end=170, fwd_mismatches=0, rev_mismatches=1, expected_size=70, **base)
    b = Hit(start=100, end=190, fwd_mismatches=0, rev_mismatches=0, expected_size=70, **base)  # same forward position
    c = Hit(start=500, end=570, fwd_mismatches=0, rev_mismatches=0, expected_size=70, **base)  # separate locus
    d = Hit(start=175, end=245, fwd_mismatches=2, rev_mismatches=0, expected_size=70, **base)  # overlaps b by 16 bp
    kept = deduplicate_hits([a, b, c, d], primer_length=20)
    assert [(h.start, h.end, h.n_alternatives) for h in kept] == [(100, 190, 1), (175, 245, 0), (500, 570, 0)]
    kept = deduplicate_hits([a, b, c, d], primer_length=10)
    assert [(h.start, h.end, h.n_alternatives) for h in kept] == [(100, 190, 2), (500, 570, 0)]
    minus = [Hit(start=100, end=170, strand="-", sample="s", contig="c", primer_id="m"),
             Hit(start=120, end=170, strand="-", sample="s", contig="c", primer_id="m")]
    assert len(deduplicate_hits(minus)) == 1  # forward primer sits at the end on the minus strand


def test_extract_dedup_in_pipeline(fixture_dir):
    contigs = load_fasta(fixture_dir / "sampleA.fasta")
    hits = [Hit("sampleA", "contig1", 301, 370, "+", "VNTR1"), Hit("sampleA", "contig1", 301, 364, "+", "VNTR1")]
    assert len(extract_amplicons(hits, contigs)) == 1
    assert len(extract_amplicons(hits, contigs, deduplicate=False)) == 2


def test_primers_near_contig_ends():
    rc_rev = str(Seq(REV).reverse_complement())
    contigs = {"c1": "A" * 50 + FWD + "C" * 500, "c2": "G" * 500 + rc_rev + "T" * 5, "c3": "A" * 400}
    pairs = [PrimerPair("m", FWD, REV, 70), PrimerPair("x", "GGGGGGGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCCCC", 50)]
    assert primers_near_contig_ends(contigs, pairs, 100) == {"m", "x"}
    assert primers_near_contig_ends({"c3": contigs["c3"]}, pairs, 100) == set()
    assert primers_near_contig_ends({"c1": contigs["c1"]}, pairs[:1], 30) == set()


def test_hits_table_round_trip(fixture_dir, tmp_path):
    contigs = load_fasta(fixture_dir / "sampleA.fasta")
    pair = PrimerPair("VNTR1", FWD, REV, 70)
    amps = extract_amplicons([Hit("sampleA", "contig1", 301, 370, "+", "VNTR1")], contigs, {"VNTR1": pair})
    path = tmp_path / "hits.tsv"
    write_hits_table(amps, path)
    back = read_hits_table(path)
    assert back == amps
    written = write_marker_fastas(amps, tmp_path / "sequences")
    assert (tmp_path / "sequences" / "VNTR1.fasta").read_text().startswith(">sampleA_1\n")
    assert set(written) == {"VNTR1"}


def test_hits_table_old_format_readable(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text("sample\tmarker\tcopy\tcontig\tstart\tend\tstrand\tlength\thas_N\tsequence\n"
                    "s\tm\t1\tc\t1\t4\t+\t4\t0\tACGT\n")
    amps = read_hits_table(path)
    assert amps[0].hit.fwd_mismatches == -1 and amps[0].hit.n_alternatives == 0
