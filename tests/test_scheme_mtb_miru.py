"""The bundled MIRU-24 and MIRU-15 schemes of Supply et al. against the H37Rv assembly.

MIRU-15 is the discriminatory subset of MIRU-24, so the two must agree at every
shared locus. The copy numbers are derived from the reference rather than
published, which the scheme notes state; these tests pin them so a change in
the search or the tandem counter shows up here.
"""

import pandas as pd
import pytest

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000195955.2"
MIRU24 = "mtb_miru24_supply_2006"
MIRU15 = "mtb_miru15_supply_2006"
# ascending locus position, the order of the published digit code
ORDER_24 = ["154", "424", "577", "580", "802", "960", "1644", "1955", "2059", "2163b", "2165", "2347",
            "2401", "2461", "2531", "2687", "2996", "3007", "3171", "3192", "3690", "4052", "4156", "4348"]
AUXILIARY = ["154", "2059", "2347", "2461", "2531", "2687", "3007", "3171", "4348"]
ORDER_15 = [m for m in ORDER_24 if m not in AUXILIARY]
PROFILE_24 = "2-2-3-3-1-3-2-1-2-5-3-3-2-3-6-1-3-3-2-3-2-5-2-2"
PROFILE_15 = "2-3-3-1-3-2-1-5-3-2-3-3-2-5-2"
UNIT_LENGTH = {"154": 53, "424": 51, "577": 58, "580": 77, "802": 54, "960": 53, "1644": 53, "1955": 57,
               "2059": 77, "2163b": 69, "2165": 75, "2347": 57, "2401": 58, "2461": 57, "2531": 53,
               "2687": 54, "2996": 51, "3007": 53, "3171": 54, "3192": 53, "3690": 58, "4052": 111,
               "4156": 59, "4348": 53}
# loci whose array is a single unit in H37Rv, so there is nothing to cross-check
SINGLE_COPY = {"802", "1955", "2687"}
# loci whose two copies differ by a few bases, counted without a stored pattern
IMPERFECT = {"2401", "3690", "4348"}


@pytest.mark.parametrize("name, order", [(MIRU24, ORDER_24), (MIRU15, ORDER_15)])
def test_scheme_files_are_consistent(name, order):
    scheme = load_scheme(name)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    # the mismatch allowance is deliberately left at its default; raising it miscounts 3171 and 4156
    assert "unit_mismatches" not in scheme.calling_defaults()
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == order
    for d in defs:
        assert d.can_call, d.primer_id
        assert d.unit_length == UNIT_LENGTH[d.primer_id]
        assert d.ref_size - d.ref_repeats * d.unit_length >= 0
        if d.primer_id in SINGLE_COPY:
            assert d.ref_repeats == 1 and d.patterns == ()
        elif d.primer_id in IMPERFECT:
            assert d.patterns == ()  # counted from the adjacent-pair test, not stored
        else:
            assert len(d.patterns) == 1 and len(d.patterns[0]) == d.unit_length


def test_miru15_is_the_discriminatory_subset_of_miru24():
    """The two schemes must carry identical primers and reference values at shared loci."""
    by_name = {}
    for name in (MIRU24, MIRU15):
        scheme = load_scheme(name)
        sts = {line.split("\t")[0]: line.split("\t")[1:3]
               for line in scheme.sts.read_text().splitlines() if line.strip() and not line.startswith("#")}
        defs = {d.primer_id: d for d in read_markers(scheme.markers)}
        by_name[name] = (sts, defs)
    sts24, defs24 = by_name[MIRU24]
    sts15, defs15 = by_name[MIRU15]
    assert set(sts15) == set(ORDER_15) and set(sts24) - set(sts15) == set(AUXILIARY)
    for marker in ORDER_15:
        assert sts15[marker] == sts24[marker], marker
        a, b = defs15[marker], defs24[marker]
        assert (a.unit_length, a.ref_size, a.ref_repeats, a.patterns, a.label) == \
               (b.unit_length, b.ref_size, b.ref_repeats, b.patterns, b.label), marker


@pytest.mark.parametrize("name, order", [(MIRU24, ORDER_24), (MIRU15, ORDER_15)])
def test_check_scheme_against_h37rv(name, order, tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", name, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == len(order)
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    assert (table["note"] == "").all()
    checked = [m for m in order if m not in SINGLE_COPY and m not in IMPERFECT]
    assert (table.loc[checked, "repeats_ok"].astype(int) == 1).all()


@pytest.mark.parametrize("name, order, profile", [(MIRU24, ORDER_24, PROFILE_24),
                                                  (MIRU15, ORDER_15, PROFILE_15)])
def test_h37rv_profile(name, order, profile, tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", name, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "h37rv", "--align", "none", "--no-report"]) == 0
    assert pd.read_csv(tmp_path / "h37rv_profile.tsv", sep="\t", index_col=0,
                       dtype=str).loc[ASSEMBLY, "profile"] == profile
    status = pd.read_csv(tmp_path / "h37rv_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "h37rv_calls_long.tsv", sep="\t")
    assert len(long) == len(order) and (long["residual_bp"] == 0).all()


def test_the_two_profiles_agree_at_the_shared_loci():
    """The 15 digit code is the 24 digit code restricted to the discriminatory loci."""
    full = dict(zip(ORDER_24, PROFILE_24.split("-"), strict=True))
    assert PROFILE_15.split("-") == [full[m] for m in ORDER_15]
