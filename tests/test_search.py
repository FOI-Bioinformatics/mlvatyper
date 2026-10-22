from mlvatyper.search import MerPcrSearcher
from mlvatyper.sts import read_sts, write_sts


def make_searcher(fixture_dir, margin, mismatches):
    valid, _ = read_sts(fixture_dir / "primers.sts")
    write_sts(valid, fixture_dir / "clean.sts")
    return MerPcrSearcher(fixture_dir / "clean.sts", margin=margin, mismatches=mismatches)


def test_merpcr_finds_planted_amplicons(fixture_dir):
    searcher = make_searcher(fixture_dir, margin=50, mismatches=0)
    hits_a = searcher.search_fasta(fixture_dir / "sampleA.fasta")
    assert len(hits_a) == 1
    assert hits_a[0].strand == "+" and hits_a[0].sample == "sampleA" and hits_a[0].contig == "contig1"
    assert hits_a[0].start == 301 and hits_a[0].length == 70

    hits_b = searcher.search_fasta(fixture_dir / "sampleB.fasta")
    strands = {(h.contig, h.strand) for h in hits_b}
    assert ("c1", "-") in strands  # reverse-strand copy with 7 repeats
    assert len(hits_b) == 2  # the N lies between the primers, so that amplicon is also found


def test_merpcr_size_window_rejects_long_product(fixture_dir):
    searcher = make_searcher(fixture_dir, margin=5, mismatches=0)
    assert searcher.search_fasta(fixture_dir / "sampleB.fasta") == []


def test_hits_match_merpcr_text_output(fixture_dir):
    """Our Hit objects must correspond to the me-PCR formatted lines merPCR writes."""
    searcher = make_searcher(fixture_dir, margin=50, mismatches=0)
    hits = searcher.search_fasta(fixture_dir / "sampleB.fasta")
    records = searcher.engine.load_fasta_file(str(fixture_dir / "sampleB.fasta"))
    lines = [h.format() for h in searcher.engine.find_hits(records)]
    assert lines == [f"{h.contig}\t{h.start}..{h.end}\t{h.primer_id}\t({h.strand})" for h in hits]


DEGENERATE_FWD = "AGAGATAAYAATCAACTCWGTTT"  # longest unambiguous run 9, no hash word at W=11
PLAIN_REV = "GTCATTCGACTCATTCAGTATC"


def write_one_sided_fixture(tmp_path):
    from Bio.Seq import Seq
    fwd_instance = DEGENERATE_FWD.replace("Y", "T").replace("W", "T")
    product = fwd_instance + "ACGT" * 40 + str(Seq(PLAIN_REV).reverse_complement())
    genome = "GGGG" * 30 + str(Seq(product).reverse_complement()) + "CCCC" * 30
    (tmp_path / "s.fasta").write_text(">g\n" + genome + "\n")
    (tmp_path / "p.sts").write_text(f"OneSided\t{DEGENERATE_FWD}\t{PLAIN_REV}\t{len(product)}\n")
    return len(product)


def test_one_sided_option_recovers_marker(tmp_path):
    size = write_one_sided_fixture(tmp_path)
    dropped = MerPcrSearcher(tmp_path / "p.sts", margin=50, mismatches=0, iupac_mode=1)
    assert dropped.search_fasta(tmp_path / "s.fasta") == []
    kept = MerPcrSearcher(tmp_path / "p.sts", margin=50, mismatches=0, iupac_mode=1, one_sided=1)
    hits = kept.search_fasta(tmp_path / "s.fasta")
    assert [(h.primer_id, h.strand, h.length) for h in hits] == [("OneSided", "-", size)]


def test_smaller_wordsize_recovers_marker(tmp_path):
    size = write_one_sided_fixture(tmp_path)
    small_w = MerPcrSearcher(tmp_path / "p.sts", margin=50, mismatches=0, iupac_mode=1, wordsize=8)
    hits = small_w.search_fasta(tmp_path / "s.fasta")
    assert len(hits) == 1 and hits[0].length == size
