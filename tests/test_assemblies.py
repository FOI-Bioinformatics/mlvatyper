import io
import zipfile

import pytest

from mlvatyper import assemblies
from mlvatyper.assemblies import cached_path, download_assembly, fasta_from_zip, is_assembly_accession
from mlvatyper.cli import main
from tests.conftest import FWD, REV, UNIT, make_amplicon, random_dna


def _zip(accession: str, records: dict[str, str]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("README.md", "not a sequence")
        for name, seq in records.items():
            archive.writestr(f"ncbi_dataset/data/{accession}/{name}_genomic.fna", f">{name}\n{seq}\n")
        archive.writestr("ncbi_dataset/data/assembly_data_report.jsonl", "{}")
    return buffer.getvalue()


def test_is_assembly_accession():
    assert is_assembly_accession("GCF_000006945.2") and is_assembly_accession(" GCA_000008485.1 ")
    assert not is_assembly_accession("AE006468.2") and not is_assembly_accession("GCF_0000069.2")


def test_fasta_from_zip_concatenates_genomic_files():
    text = fasta_from_zip(_zip("GCF_000000001.1", {"chr": "ACGT", "plasmid": "GGCC"}), "GCF_000000001.1")
    assert text == ">chr\nACGT\n>plasmid\nGGCC\n"
    with pytest.raises(RuntimeError):
        fasta_from_zip(_zip("GCF_000000001.1", {}), "GCF_000000001.1")


def test_download_assembly_uses_cache(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(assemblies, "_fetch_zip", lambda acc: calls.append(acc) or _zip(acc, {"c": "ACGT"}))
    path = download_assembly("GCF_000000001.1", tmp_path)
    assert path == cached_path("GCF_000000001.1", tmp_path) and path.read_text() == ">c\nACGT\n"
    download_assembly("GCF_000000001.1", tmp_path)
    assert calls == ["GCF_000000001.1"]  # second call served from the cache
    download_assembly("GCF_000000001.1", tmp_path, force=True)
    assert len(calls) == 2
    with pytest.raises(ValueError):
        download_assembly("AE006468", tmp_path)
    monkeypatch.setenv("MLVATYPER_CACHE", str(tmp_path / "env"))
    assert cached_path("GCF_000000001.1").parent == tmp_path / "env"


def test_check_scheme_downloads_the_scheme_assembly(tmp_path, monkeypatch):
    genome = random_dna(200, 7) + make_amplicon(5) + random_dna(200, 8)
    monkeypatch.setattr(assemblies, "_fetch_zip", lambda acc: _zip(acc, {"chr": genome}))
    d = tmp_path / "scheme"
    d.mkdir()
    (d / "primers.sts").write_text(f"VNTR1\t{FWD}\t{REV}\t70\n")
    (d / "markers.tsv").write_text(f"marker\tpatterns\tunit_length\tref_size\tref_repeats\nVNTR1\t{UNIT}\t6\t70\t5\n")
    (d / "scheme.yml").write_text("name: s\ncitation: c\nurl: u\nreference: r\nassembly: GCF_000000002.1\n")
    cache = tmp_path / "cache"
    out = tmp_path / "out"
    assert main(["check-scheme", "--scheme", str(d), "--outdir", str(out), "--cache-dir", str(cache),
                 "-N", "0", "-M", "50"]) == 0
    assert (cache / "GCF_000000002.1.fasta").exists()
    table = (out / "scheme_check.tsv").read_text().splitlines()
    assert table[1].split("\t")[0] == "VNTR1" and "\t1\t" in table[1]  # one product
    assert main(["check-scheme", "--scheme", str(d), "--outdir", str(out), "--cache-dir", str(cache),
                 "--no-download"]) == 0
    # an accession on the command line overrides the scheme's
    monkeypatch.setattr(assemblies, "_fetch_zip", lambda acc: (_ for _ in ()).throw(RuntimeError("offline")))
    assert main(["check-scheme", "--scheme", str(d), "--outdir", str(out), "--cache-dir", str(cache),
                 "--assembly", "GCF_000000003.1"]) == 1
