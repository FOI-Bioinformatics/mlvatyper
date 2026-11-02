"""The bundled Brucella MLVA-16 scheme against the B. melitensis 16M assembly, downloaded by accession."""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCF_000007125.1"
# scheme order: panel 1 (MLVA-8), panel 2A, panel 2B
ORDER = ["Bruce06", "Bruce08", "Bruce11", "Bruce12", "Bruce42", "Bruce43", "Bruce45", "Bruce55",
         "Bruce18", "Bruce19", "Bruce21", "Bruce04", "Bruce07", "Bruce09", "Bruce16", "Bruce30"]
# repeat unit (bp), B. melitensis 16M product (bp) and copy number, from Le Fleche et al. 2006 Table 1
SOURCE = {
    "Bruce06": (134, 408, 3), "Bruce08": (18, 348, 4), "Bruce11": (63, 257, 2), "Bruce12": (15, 392, 13),
    "Bruce42": (125, 539, 4), "Bruce43": (12, 182, 2), "Bruce45": (18, 151, 3), "Bruce55": (40, 273, 3),
    "Bruce18": (8, 146, 5), "Bruce19": (6, 163, 18), "Bruce21": (8, 148, 6),
    "Bruce04": (8, 152, 2), "Bruce07": (8, 158, 5), "Bruce09": (8, 156, 7), "Bruce16": (8, 152, 3),
    "Bruce30": (8, 151, 6),
}


def test_scheme_files_match_the_source_table():
    scheme = load_scheme("Bru_MLVA16")
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    assert list(defs) == ORDER
    for marker, (unit, size, repeats) in SOURCE.items():
        d = defs[marker]
        assert (d.unit_length, d.ref_size, d.ref_repeats) == (unit, size, repeats), marker
        assert d.can_call and not d.has_repeats  # called from length alone, no unit sequences shipped
        assert d.offset == size - repeats * unit and d.offset >= 0, marker
    sts = {line.split("\t")[0]: int(line.split("\t")[3])
           for line in scheme.sts.read_text().splitlines() if line and not line.startswith("#")}
    assert sts == {m: SOURCE[m][1] for m in ORDER}


def test_check_scheme_against_16m(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", "Bru_MLVA16", "--outdir", str(tmp_path), "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert list(table.index) == ORDER
    assert (table["n_products"].astype(int) == 1).all() and (table["size_ok"].astype(int) == 1).all()
    assert (table["length"].astype(int) == pd.Series({m: SOURCE[m][1] for m in ORDER})).all()


def test_16m_profile_equals_the_published_copy_numbers(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", "Bru_MLVA16", "--fasta", str(fixture), "--outdir", str(tmp_path), "--prefix", "b",
                 "--mlvabank", "--align", "none", "--no-report"]) == 0
    alleles = pd.read_csv(tmp_path / "b_alleles.tsv", sep="\t", index_col=0, dtype=str)
    assert alleles.loc[ASSEMBLY].to_dict() == {m: str(SOURCE[m][2]) for m in ORDER}
    status = pd.read_csv(tmp_path / "b_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "b_calls_long.tsv", sep="\t")
    assert (long.fwd_mismatches + long.rev_mismatches == 0).all() and (long.residual_bp == 0).all()
    assert (tmp_path / "b_mlvabank.csv").read_text().splitlines()[0] == "Strain," + ",".join(ORDER)
