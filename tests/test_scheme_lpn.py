"""The bundled Legionella pneumophila MLVA-12 scheme against the Philadelphia-1 assembly, downloaded by accession."""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCF_000008485.1"
ORDER = ["Lpms01", "Lpms03", "Lpms13", "Lpms19", "Lpms31", "Lpms33", "Lpms34", "Lpms35", "Lpms38", "Lpms39",
         "Lpms40", "Lpms44"]
# repeat unit (bp), Philadelphia-1 size (bp) and copy number as given in the source table
SOURCE = {
    "Lpms01": (45, 494, 8), "Lpms03": (96, 942, 8), "Lpms13": (24, 793, 11), "Lpms19": (21, 128, 4),
    "Lpms31": (45, 1043, 17), "Lpms33": (125, 317, 1), "Lpms34": (125, 265, 1), "Lpms35": (18, 205, 3),
    "Lpms38": (8, 264, 3), "Lpms39": (6, 79, 6), "Lpms40": (6, 198, 4), "Lpms44": (6, 173, 9),
}


def test_scheme_files_match_the_source_table():
    scheme = load_scheme("Lpn_MLVA12")
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    assert list(defs) == ORDER
    for marker, (unit, size, repeats) in SOURCE.items():
        d = defs[marker]
        assert (d.unit_length, d.ref_size, d.ref_repeats) == (unit, size, repeats), marker
        assert d.can_call and not d.has_repeats  # called from length alone
        assert d.offset == size - repeats * unit and d.offset >= 40  # the flank holds both primers
    sts = {line.split("\t")[0]: int(line.split("\t")[3]) for line in scheme.sts.read_text().splitlines()}
    assert sts == {m: SOURCE[m][1] for m in ORDER}


def test_check_scheme_against_philadelphia_1(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", "Lpn_MLVA12", "--outdir", str(tmp_path), "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert list(table.index) == ORDER
    assert (table["n_products"].astype(int) == 1).all() and (table["size_ok"].astype(int) == 1).all()
    assert (table["length"].astype(int) == pd.Series({m: SOURCE[m][1] for m in ORDER})).all()


def test_philadelphia_profile_equals_the_published_copy_numbers(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", "Lpn_MLVA12", "--fasta", str(fixture), "--outdir", str(tmp_path), "--prefix", "p",
                 "--mlvabank", "--align", "none", "--no-report"]) == 0
    alleles = pd.read_csv(tmp_path / "p_alleles.tsv", sep="\t", index_col=0, dtype=str)
    assert alleles.loc["GCF_000008485.1"].to_dict() == {m: str(SOURCE[m][2]) for m in ORDER}
    status = pd.read_csv(tmp_path / "p_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc["GCF_000008485.1"]) == {"OK"}
    long = pd.read_csv(tmp_path / "p_calls_long.tsv", sep="\t")
    assert (long.fwd_mismatches + long.rev_mismatches == 0).all() and (long.residual_bp == 0).all()
    assert (tmp_path / "p_mlvabank.csv").read_text().splitlines()[0] == "Strain," + ",".join(ORDER)
