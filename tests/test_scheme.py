import pandas as pd
import pytest

from mlvatyper.cli import main
from mlvatyper.markers import read_markers, write_markers
from mlvatyper.mlvabank import bins_from_table, export_profiles
from mlvatyper.models import MarkerDef
from mlvatyper.scheme import SchemeError, is_placeholder, list_schemes, load_scheme, parse_flat_yaml
from tests.conftest import FWD, REV, UNIT


def test_bundled_scheme_and_metadata():
    names = [s.name for s in list_schemes()]
    assert "testM" in names
    scheme = load_scheme("testM")
    assert scheme.sts.name == "primers.sts" and scheme.markers is not None and scheme.bins is None
    assert scheme.meta["organism"] == "test organism" and scheme.calling_defaults()["rounding"] == "intermediate"
    assert scheme.calling_defaults()["allele_tolerance"] == 0.25
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == ["assembly"]  # private genome
    with pytest.raises(SchemeError):
        load_scheme("does-not-exist")


def test_parse_flat_yaml():
    meta = parse_flat_yaml("# comment\nname: 'X'\nrounding: nearest  # trailing\nmargin: 300\n\nnotes: a: b\n")
    assert meta == {"name": "X", "rounding": "nearest", "margin": 300, "notes": "a: b"}
    with pytest.raises(SchemeError):
        parse_flat_yaml("margin: many\n")
    with pytest.raises(SchemeError):
        parse_flat_yaml("no separator\n")


def _scheme_dir(fixture_dir, rounding="nearest"):
    d = fixture_dir / "scheme"
    d.mkdir()
    (d / "primers.sts").write_text(f"VNTR1\t{FWD}\t{REV}\t70\n")
    (d / "markers.tsv").write_text(f"marker\tpatterns\tunit_length\tref_size\tref_repeats\tlabel\n"
                                   f"VNTR1\t{UNIT}\t6\t70\t5\tLocus_1\n")
    (d / "bins.tsv").write_text("marker\tsize\tallele\nVNTR1\t82\t9\n")
    (d / "scheme.yml").write_text(f"name: mini\norganism: synthetic\nrounding: {rounding}\nallele_tolerance: 0.1\n"
                                  "citation: Example et al. J Test 2026;1:1-2\nurl: https://example.test/mini\n"
                                  "reference: strain MINI-1\nretrieved: 2026-10-23\n")
    return d


def test_run_with_scheme_directory_applies_files_and_defaults(fixture_dir):
    d = _scheme_dir(fixture_dir)
    out = fixture_dir / "out"
    assert main(["run", "--scheme", str(d), "--fasta", str(fixture_dir / "*.fasta"), "--outdir", str(out),
                 "--prefix", "t", "-N", "0", "-M", "50", "--align", "none", "--no-report", "--mlvabank"]) == 0
    alleles = pd.read_csv(out / "t_alleles.tsv", sep="\t", index_col=0, dtype=str)
    assert alleles.loc["sampleA", "VNTR1"] == "5"
    status = pd.read_csv(out / "t_status.tsv", sep="\t", index_col=0, dtype=str)
    assert "BINNED" in status.loc["sampleB", "VNTR1"]  # the 82 bp amplicon hits the bin
    info = __import__("json").loads((out / "run_info.json").read_text())
    assert info["summarize"]["parameters"]["rounding"] == "nearest"
    assert info["summarize"]["parameters"]["allele_tolerance"] == 0.1
    bank = (out / "t_mlvabank.csv").read_text().splitlines()
    assert bank[0] == "Strain,Locus_1" and bank[1] == "sampleA,5"
    # an explicit option wins over the scheme default
    assert main(["summarize", "--scheme", str(d), "--outdir", str(out), "--prefix", "u", "--rounding", "up"]) == 0
    info = __import__("json").loads((out / "run_info.json").read_text())
    assert info["summarize"]["parameters"]["rounding"] == "up"


def test_sts_or_scheme_required(fixture_dir):
    assert main(["search", "--fasta", str(fixture_dir / "sampleA.fasta"), "--outdir", str(fixture_dir / "o")]) == 1
    assert main(["search", "--scheme", "nope", "--fasta", str(fixture_dir / "sampleA.fasta"),
                 "--outdir", str(fixture_dir / "o")]) == 1


def test_export_profiles_and_missing_token(tmp_path):
    alleles = pd.DataFrame({"m1": ["5", "2/2"], "m2": ["NOT_DETECTED", "3"]}, index=["a", "b"])
    defs = [MarkerDef("m2", label="Second"), MarkerDef("m1"), MarkerDef("absent")]
    export_profiles(alleles, defs, tmp_path / "p.csv")
    assert (tmp_path / "p.csv").read_text().splitlines() == ["Strain,Second,m1", "a,NOT_DETECTED,5", "b,3,2/2"]
    export_profiles(alleles, defs, tmp_path / "q.csv", missing="-", sep="\t")
    assert (tmp_path / "q.csv").read_text().splitlines()[1:] == ["a\t-\t5", "b\t3\t-"]


def test_bins_from_table_and_labels(tmp_path):
    table = tmp_path / "sizes.csv"
    table.write_text("Locus,Allele,Fragment size\nL1,5,310\nL1,6,313.4\nL2,3,440-460\nL2,4,abc\n")
    n = bins_from_table(table, tmp_path / "bins.tsv", labels={"L1": "VNTR1"})
    lines = (tmp_path / "bins.tsv").read_text().splitlines()
    assert n == 3 and lines[1:] == ["VNTR1\t310-310\t5", "VNTR1\t313-313\t6", "L2\t440-460\t3"]
    (tmp_path / "m.tsv").write_text("marker\tlabel\nVNTR1\tL1\n")
    assert main(["bins-from-table", str(table), "--output", str(tmp_path / "b2.tsv"),
                 "--markers", str(tmp_path / "m.tsv")]) == 0
    assert (tmp_path / "b2.tsv").read_text() == (tmp_path / "bins.tsv").read_text()
    (tmp_path / "bad.csv").write_text("a,b\n1,2\n")
    assert main(["bins-from-table", str(tmp_path / "bad.csv"), "--output", str(tmp_path / "b3.tsv")]) == 1


def test_label_round_trip(tmp_path):
    defs = [MarkerDef("m1", ("GATCGA",), 6, 70, 5, label="Locus_1"), MarkerDef("m2")]
    write_markers(defs, tmp_path / "m.tsv")
    back = read_markers(tmp_path / "m.tsv")
    assert back[0].label == "Locus_1" and back[1].label is None and back[0].display_name == "Locus_1"


def test_is_placeholder_and_missing_provenance(fixture_dir):
    for value in ("", "  ", "none", "None.", "N/A", "unknown", "TBD", "to be added", "-"):
        assert is_placeholder(value)
    for value in ("Example et al. 2026", "https://example.test/x", "strain MINI-1", "0"):
        assert not is_placeholder(value)
    d = _scheme_dir(fixture_dir)
    (d / "scheme.yml").write_text("name: mini\ncitation: none\nurl:\nreference: strain MINI-1\n")
    scheme = load_scheme(d)
    assert scheme.missing_provenance() == ["citation", "url"]
    assert scheme.missing_recommended() == ["organism", "retrieved", "assembly"]
    assert scheme.provenance["reference"] == "strain MINI-1"


def test_check_scheme_requires_provenance(fixture_dir, caplog):
    d = _scheme_dir(fixture_dir)
    out = fixture_dir / "chk"
    (d / "scheme.yml").write_text("name: mini\ncitation: TBD\nreference: strain MINI-1\n")
    assert main(["check-scheme", "--scheme", str(d), "--outdir", str(out)]) == 1
    assert "must state citation, url" in caplog.text
    assert (out / "scheme_check.tsv").exists()  # the other checks still ran and were written
    caplog.clear()
    assert main(["check-scheme", "--scheme", str(d), "--outdir", str(out), "--allow-missing-provenance"]) == 0
    assert "cannot be traced without citation, url" in caplog.text
    assert not [r for r in caplog.records if r.levelname == "ERROR"]


def test_check_scheme(fixture_dir):
    d = _scheme_dir(fixture_dir)
    (d / "primers.sts").write_text(f"VNTR1\t{FWD}\t{REV}\t70\nABSENT\t{'G' * 20}\t{'C' * 20}\t100\n"
                                   f"UNHASH\tACGTNACGTNACGTNACGTN\tCCCCCCCCCCCCCCCCCCCC\t100\n")
    out = fixture_dir / "chk"
    assert main(["check-scheme", "--scheme", str(d), "--outdir", str(out)]) == 0
    table = pd.read_csv(out / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert table.loc["UNHASH", "hash_word_11"] == 0 and table.loc["VNTR1", "in_bins"] == 1
    assert table.loc["ABSENT", "in_markers"] == 0 and table["n_products"].eq("").all()
    assert main(["check-scheme", "--scheme", str(d), "--outdir", str(out), "--reference",
                 str(fixture_dir / "sampleA.fasta"), "-N", "0", "-M", "50"]) == 0
    table = pd.read_csv(out / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    row = table.loc["VNTR1"]
    assert (int(row["n_products"]), int(row["length"]), int(row["size_ok"]), int(row["repeats_ok"])) == (1, 70, 1, 1)
    assert int(table.loc["ABSENT", "n_products"]) == 0 and "not found" in table.loc["ABSENT", "note"]
    assert main(["check-scheme", "--scheme", str(d), "--outdir", str(out), "--fetch", "X"]) == 1  # needs --email
    assert main(["check-scheme", "--scheme", str(d), "--outdir", str(out), "--reference", "missing.fa"]) == 1


def test_new_scheme_scaffold(tmp_path, caplog):
    d = tmp_path / "organism_2015"
    assert main(["new-scheme", str(d)]) == 0
    assert sorted(p.name for p in d.iterdir()) == ["bins.tsv", "markers.tsv", "primers.sts", "scheme.yml"]
    scheme = load_scheme(d)
    assert scheme.name == "organism_2015" and scheme.calling_defaults()["rounding"] == "intermediate"
    assert scheme.missing_provenance() == ["citation", "url", "reference"]
    assert read_markers(d / "markers.tsv") == []  # header and comments only
    from mlvatyper.calling import read_bins
    from mlvatyper.sts import read_sts
    assert read_bins(d / "bins.tsv") == {} and read_sts(d / "primers.sts") == ([], [])
    # the template is usable as soon as the files are filled in
    (d / "primers.sts").write_text(f"L1\t{FWD}\t{REV}\t70\n")
    (d / "markers.tsv").write_text("marker\tunit_length\tref_size\tref_repeats\tpatterns\tlabel\n"
                                   f"L1\t6\t70\t5\t{UNIT}\tL1\n")
    assert read_markers(d / "markers.tsv")[0].ref_repeats == 5
    assert main(["new-scheme", str(d)]) == 1 and "already holds a scheme" in caplog.text
    assert main(["new-scheme", str(tmp_path / "other"), "--name", "named"]) == 0
    assert load_scheme(tmp_path / "other").name == "named"


def test_assembly_accession_in_metadata(tmp_path):
    d = _scheme_dir(tmp_path)
    (d / "scheme.yml").write_text("name: x\nassembly: GCF_000006945.2\n")
    assert load_scheme(d).assembly == "GCF_000006945.2"
    (d / "scheme.yml").write_text("name: x\nassembly: AE006468\n")
    with pytest.raises(SchemeError):
        load_scheme(d)
    (d / "scheme.yml").write_text("name: x\nassembly: none\n")
    assert load_scheme(d).assembly is None
