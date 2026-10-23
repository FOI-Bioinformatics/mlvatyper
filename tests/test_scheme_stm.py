"""The bundled Salmonella Typhimurium scheme against the committed LT2 amplicon regions."""

from pathlib import Path

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme

FIXTURE = Path(__file__).parent / "data" / "STm_ECDC_reference.fasta"
EXPECTED = {"STTR9": "4", "STTR5": "13", "STTR6": "13", "STTR10pl": "10", "STTR3": "NO_DEFINITION"}


def test_scheme_files_are_consistent():
    scheme = load_scheme("STm_ECDC")
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.calling_defaults()["unit_mismatches"] == 0
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    assert list(defs) == ["STTR9", "STTR5", "STTR6", "STTR10pl", "STTR3"]  # the profile order of the nomenclature
    for d in defs.values():
        if d.can_call:
            assert d.offset == d.ref_size - d.ref_repeats * d.unit_length
            assert len(d.patterns) == 1 and len(d.patterns[0]) == d.unit_length
    assert defs["STTR10pl"].label == "STTR10" and not defs["STTR3"].can_call and defs["STTR3"].ref_size == 490


def test_check_scheme_against_lt2_regions(tmp_path):
    assert main(["check-scheme", "--scheme", "STm_ECDC", "--outdir", str(tmp_path), "--reference", str(FIXTURE)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert (table["n_products"].astype(int) == 1).all() and (table["size_ok"].astype(int) == 1).all()
    assert list(table.loc[["STTR9", "STTR5", "STTR6", "STTR10pl"], "repeats_ok"].astype(int)) == [1, 1, 1, 1]


def test_lt2_profile_is_4_13_13_10(tmp_path):
    assert main(["run", "--scheme", "STm_ECDC", "--fasta", str(FIXTURE), "--outdir", str(tmp_path), "--prefix", "lt2",
                 "--mlvabank", "--align", "none", "--no-report"]) == 0
    alleles = pd.read_csv(tmp_path / "lt2_alleles.tsv", sep="\t", index_col=0, dtype=str)
    assert alleles.loc["STm_ECDC_reference"].to_dict() == EXPECTED
    status = pd.read_csv(tmp_path / "lt2_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc["STm_ECDC_reference", ["STTR9", "STTR5", "STTR6", "STTR10pl"]]) == {"OK"}
    bank = (tmp_path / "lt2_mlvabank.csv").read_text().splitlines()
    assert bank[0] == "Strain,STTR9,STTR5,STTR6,STTR10,STTR3"
    long = pd.read_csv(tmp_path / "lt2_calls_long.tsv", sep="\t")
    found = long[long.marker != "STTR3"]
    assert (found.residual_bp == 0).all() and (found.tandem_copies == found.allele.astype(int)).all()
    assert (found.expected_size == found.length).all()  # STS sizes are the LT2 products
