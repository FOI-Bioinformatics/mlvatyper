import bz2
import gzip
import lzma
from pathlib import Path

import pandas as pd

from vntrfinder.cli import main
from vntrfinder.extract import load_fasta
from vntrfinder.io import open_text, strip_compression
from vntrfinder.search import sample_name
from vntrfinder.sts import read_sts

DATA = Path(__file__).parent / "data"


def test_open_text_handles_each_compression(tmp_path):
    text = "a\tb\n"
    (tmp_path / "x.txt").write_text(text)
    with gzip.open(tmp_path / "x.txt.gz", "wt") as fh:
        fh.write(text)
    with bz2.open(tmp_path / "x.txt.bz2", "wt") as fh:
        fh.write(text)
    with lzma.open(tmp_path / "x.txt.xz", "wt") as fh:
        fh.write(text)
    for name in ["x.txt", "x.txt.gz", "x.txt.bz2", "x.txt.xz"]:
        with open_text(tmp_path / name) as fh:
            assert fh.read() == text


def test_strip_compression_and_sample_name():
    assert strip_compression("s/mini.fasta.gz") == Path("s/mini.fasta")
    assert strip_compression("mini.fasta") == Path("mini.fasta")
    assert sample_name("mini.fasta.gz") == "mini"
    assert sample_name("mini.fa.bz2") == "mini"
    assert sample_name("mini.fasta") == "mini"


def test_gz_fixtures_match_plain():
    plain, _ = read_sts(DATA / "mini.sts")
    gz, _ = read_sts(DATA / "mini.sts.gz")
    assert [p.id for p in gz] == [p.id for p in plain]
    assert load_fasta(DATA / "mini.fasta.gz") == load_fasta(DATA / "mini.fasta")
    assert len(load_fasta(DATA / "mini.fasta")["mini"]) == 4129


def test_search_compressed_inputs(tmp_path):
    out = tmp_path / "out"
    rc = main(["search", "--sts", str(DATA / "mini.sts.gz"), "--fasta", str(DATA / "mini.fasta.gz"),
               "--outdir", str(out), "--align", "none"])
    assert rc == 0
    hits = pd.read_csv(out / "hits.tsv", sep="\t")
    assert set(hits["sample"]) == {"mini"}
    assert list(hits["marker"]) == ["testID146", "testID147", "testID148", "testID149", "testID150"]
    assert list(hits["length"]) == [366, 390, 335, 356, 364]
