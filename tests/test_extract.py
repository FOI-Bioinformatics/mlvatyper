from Bio.Seq import Seq

from mlvatyper.extract import (
    extract_amplicons, load_fasta, read_hits_table, write_hits_table, write_marker_fastas,
)
from mlvatyper.models import Hit


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


def test_hits_table_round_trip(fixture_dir, tmp_path):
    contigs = load_fasta(fixture_dir / "sampleA.fasta")
    amps = extract_amplicons([Hit("sampleA", "contig1", 301, 370, "+", "VNTR1")], contigs)
    path = tmp_path / "hits.tsv"
    write_hits_table(amps, path)
    back = read_hits_table(path)
    assert back == amps
    written = write_marker_fastas(amps, tmp_path / "sequences")
    assert (tmp_path / "sequences" / "VNTR1.fasta").read_text().startswith(">sampleA_1\n")
    assert set(written) == {"VNTR1"}
