"""The bundled Mycoplasma agalactiae MLVA of McAuliffe et al. 2008 against type strain PG2.

Table 1 gives the repeat unit and the PG2 copy number for all four VNTRs. Three
counts are fractional because the arrays are imperfect, so the reference alleles
lie on an offset grid. Magal_VNTR5 has a mismatched reverse primer and is found
by the rescue pass.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from mlvatyper.sts import read_sts
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000063605.1"
SCHEME = "magalactiae_mlva4_mcauliffe_2008"
CHROMOSOME = "CU179680.1"

ORDER = ["Magal_VNTR5", "Magal_VNTR14", "Magal_VNTR17", "Magal_VNTR19"]
UNIT = {"Magal_VNTR5": 24, "Magal_VNTR14": 13, "Magal_VNTR17": 14, "Magal_VNTR19": 22}
PG2 = {"Magal_VNTR5": 2.3, "Magal_VNTR14": 2.0, "Magal_VNTR17": 3.6, "Magal_VNTR19": 2.3}
REF_SIZE = {"Magal_VNTR5": 226, "Magal_VNTR14": 158, "Magal_VNTR17": 209, "Magal_VNTR19": 184}
RESCUED = "Magal_VNTR5"  # reverse primer mismatches PG2, found by the aligned/rescue pass


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
        assert d.ref_repeats == PG2[d.primer_id], d.primer_id
        assert d.ref_size == REF_SIZE[d.primer_id], d.primer_id
        assert d.patterns == ()
    # three loci carry a fractional reference count, one is whole
    assert {m for m in ORDER if PG2[m] != int(PG2[m])} == {"Magal_VNTR5", "Magal_VNTR17", "Magal_VNTR19"}


def test_check_scheme_against_pg2(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 4
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    # the mismatched locus is recovered by a later pass, not the main one
    assert table.loc[RESCUED, "search_pass"] != "main"
    assert (table.loc[[m for m in ORDER if m != RESCUED], "search_pass"] == "main").all()


def test_pg2_profile_is_the_published_copy_numbers(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "pg2", "--rescue-indels", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "pg2_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == "2.3-2-3.6-2.3"
    long = pd.read_csv(tmp_path / "pg2_calls_long.tsv", sep="\t").set_index("marker")
    assert (long["residual_bp"] == 0).all()
    assert set(long["contig"]) == {CHROMOSOME}
    assert {m: float(long.loc[m, "allele"]) for m in ORDER} == PG2
