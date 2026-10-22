import shutil
from pathlib import Path

import pytest
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mlvatyper.align import align_marker, detect_aligner, fasta_to_clustal, filter_n


def write(path: Path, seqs: dict[str, str]) -> None:
    SeqIO.write([SeqRecord(Seq(s), id=k, description="") for k, s in seqs.items()], str(path), "fasta")


def test_fasta_to_clustal_parses_back(tmp_path):
    fa = tmp_path / "a.fasta"
    write(fa, {"x_1": "ACGT--ACGT", "y_1": "ACGTTTACGT"})
    clw = tmp_path / "a.clw"
    fasta_to_clustal(fa, clw)
    aln = AlignIO.read(str(clw), "clustal")
    assert len(aln) == 2 and aln.get_alignment_length() == 10


def test_single_sequence_is_copied(tmp_path):
    fa = tmp_path / "m.fasta"
    write(fa, {"s_1": "ACGTACGT"})
    assert align_marker(None, fa, tmp_path / "m_aligned.fasta", tmp_path / "m.clw")
    assert (tmp_path / "m_aligned.fasta").read_text() == fa.read_text()
    assert AlignIO.read(str(tmp_path / "m.clw"), "clustal").get_alignment_length() == 8


def test_filter_n(tmp_path):
    fa = tmp_path / "m.fasta"
    write(fa, {"a": "ACGT", "b": "ACNT", "c": "acgt"})
    kept, removed = filter_n(fa, tmp_path / "clean.fasta")
    assert (kept, removed) == (2, 1)


def test_detect_aligner_none():
    assert detect_aligner("none") is None
    with pytest.raises(FileNotFoundError):
        detect_aligner("no_such_aligner_binary")


@pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not on PATH")
def test_mafft_alignment(tmp_path):
    fa = tmp_path / "m.fasta"
    write(fa, {"a_1": "GATCGA" * 5, "b_1": "GATCGA" * 7})
    assert align_marker("mafft", fa, tmp_path / "m_aligned.fasta", tmp_path / "m.clw")
    aln = AlignIO.read(str(tmp_path / "m.clw"), "clustal")
    assert len(aln) == 2 and aln.get_alignment_length() == 42


def test_align_all_with_clwclean_and_no_aligner(tmp_path):
    """Without an aligner, single sequences still get CLUSTAL output and clean files are produced."""
    from mlvatyper.align import align_all
    write(tmp_path / "m1.fasta", {"s_1": "ACGTACGT"})
    write(tmp_path / "m2.fasta", {"a_1": "ACGT", "b_1": "ACNT"})
    align_all(tmp_path, None, clwclean=True)
    assert (tmp_path / "m1_aligned.clw").exists()
    assert (tmp_path / "m1_clean_aligned.clw").exists()
    assert not (tmp_path / "m2_aligned.clw").exists()  # two sequences, no aligner
    assert (tmp_path / "m2_clean.fasta").exists()
    assert (tmp_path / "m2_clean_aligned.clw").exists()  # one clean sequence left
    # derived files are not treated as markers on a second pass
    names_before = sorted(p.name for p in tmp_path.iterdir())
    align_all(tmp_path, None, clwclean=True)
    assert names_before == sorted(p.name for p in tmp_path.iterdir())


def test_muscle_command_selection(tmp_path, monkeypatch):
    """The MUSCLE command line follows the detected major version."""
    import subprocess
    from mlvatyper import align as mod

    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(mod.subprocess, "run", fake_run)
    monkeypatch.setattr(mod, "muscle_major_version", lambda: 5)
    mod.run_aligner("muscle", Path("in.fa"), Path("out.fa"))
    assert calls[-1][:2] == ["muscle", "-align"]
    monkeypatch.setattr(mod, "muscle_major_version", lambda: 3)
    mod.run_aligner("muscle", Path("in.fa"), Path("out.fa"))
    assert calls[-1][:2] == ["muscle", "-in"] and "-quiet" in calls[-1]
    with pytest.raises(ValueError):
        mod.run_aligner("clustalw", Path("in.fa"), Path("out.fa"))


def test_muscle_version_parsing(monkeypatch):
    import subprocess
    from mlvatyper import align as mod

    def fake(text):
        return lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0, text, "")

    monkeypatch.setattr(mod.subprocess, "run", fake("muscle 5.1.osx64"))
    assert mod.muscle_major_version() == 5
    monkeypatch.setattr(mod.subprocess, "run", fake("MUSCLE v3.8.1551 by Robert C. Edgar"))
    assert mod.muscle_major_version() == 3


def test_aligner_failure_is_reported(tmp_path, monkeypatch):
    import subprocess
    from mlvatyper import align as mod

    def failing(cmd, **kwargs):
        raise subprocess.CalledProcessError(1, cmd, stderr="boom")

    monkeypatch.setattr(mod.subprocess, "run", failing)
    fa = tmp_path / "m.fasta"
    write(fa, {"a_1": "ACGT", "b_1": "ACGG"})
    assert not align_marker("mafft", fa, tmp_path / "o.fa", tmp_path / "o.clw")
