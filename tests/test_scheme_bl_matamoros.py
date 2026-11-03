"""The bundled Bifidobacterium longum MLVA of Matamoros et al. 2011 against NCC2705.

Every reference value is published: Table 3 gives the repeat units and Table 1
the NCC2705 profile, 6-4-3-7-16. The primers reproduce each amplicon inside the
published size range, and a run on NCC2705 returns the published profile.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from mlvatyper.sts import read_sts
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000007525.1"
SCHEME = "blongum_mlva5_matamoros_2011"
CHROMOSOME = "AE014295.3"

ORDER = ["VNTR12", "VNTR21", "VNTR23", "VNTR25", "VNTR26"]
UNIT = {"VNTR12": 14, "VNTR21": 9, "VNTR23": 10, "VNTR25": 12, "VNTR26": 10}
NCC2705 = {"VNTR12": 6, "VNTR21": 4, "VNTR23": 3, "VNTR25": 7, "VNTR26": 16}
PROFILE = "6-4-3-7-16"
# Table 2 amplicon size ranges; the measured product must fall inside (allowing one bp of slack)
SIZE_RANGE = {"VNTR12": (436, 500), "VNTR21": (374, 392), "VNTR23": (313, 373),
              "VNTR25": (365, 414), "VNTR26": (252, 332)}


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    pairs, _ = read_sts(scheme.sts)
    assert [p.id for p in pairs] == ORDER
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.can_call, d.primer_id
        assert d.unit_length == UNIT[d.primer_id], d.primer_id
        assert d.ref_repeats == NCC2705[d.primer_id], d.primer_id
        assert d.patterns == (), "imperfect arrays, so no stored pattern"
        low, high = SIZE_RANGE[d.primer_id]
        assert low - 1 <= d.ref_size <= high + 1, d.primer_id


def test_check_scheme_against_ncc2705(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 5
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()


def test_ncc2705_profile_is_the_published_one(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "ncc", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "ncc_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    assert PROFILE.split("-") == [str(NCC2705[m]) for m in ORDER]
    status = pd.read_csv(tmp_path / "ncc_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "ncc_calls_long.tsv", sep="\t").set_index("marker")
    assert (long["residual_bp"] == 0).all()
    assert set(long["contig"]) == {CHROMOSOME}
    assert {m: int(long.loc[m, "allele"]) for m in ORDER} == NCC2705
