"""The bundled Yersinia pestis MLVA-25 panel of Le Fleche et al. against CO92.

Unit length, repeat number and expected product all come from the source table,
and every one of the 25 products is reproduced in CO92 at the position the full
locus name encodes, so these tests check published values rather than derived
ones. Six loci are shared with yersinia_klevytska_2001, which numbers them
differently; that difference is pinned here so it cannot be lost.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000009065.1"
SCHEME = "yersinia_mlva25_lefleche_2001"
KLEVYTSKA = "yersinia_klevytska_2001"
CHROMOSOME = "AL590842.1"
# the order of the source table, by descending repeat number
ORDER = ["ms01", "ms04", "ms05", "ms06", "ms07", "ms09", "ms15", "ms20", "ms21", "ms35", "ms38",
         "ms40", "ms41", "ms44", "ms45", "ms46", "ms51", "ms54", "ms56", "ms62", "ms69", "ms70",
         "ms71", "ms73", "ms74"]
PROFILE = "8-8-11-9-9-33-10-9-9-8-8-7-7-7-7-7-7-7-7-7-6-6-6-6-6"
# the equivalences the source states, confirmed in the reference
EQUIVALENT = {"ms04": "M58", "ms05": "M59", "ms07": "M37", "ms20": "M51", "ms45": "M42", "ms62": "M34"}
# the seven marker subset Pourcel et al. 2004 use for routine work
POURCEL7 = ["ms01", "ms04", "ms06", "ms07", "ms46", "ms62", "ms70"]
# the one locus whose published size range excludes its own published product
INCONSISTENT_RANGE = {"ms40": (214, "220-240")}


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    assert "unit_mismatches" not in scheme.calling_defaults()
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    assert len(defs) == 25, "the panel name is MLVA-25"
    for d in defs:
        assert d.can_call, d.primer_id
        assert d.ref_size - d.ref_repeats * d.unit_length >= 0, d.primer_id
        # the arrays in CO92 are degenerate and no locus reaches its published repeat number,
        # so no motif is stored and the tandem count is not a cross-check anywhere
        assert d.patterns == (), d.primer_id
        assert d.label.startswith("yp") and d.label.endswith(d.primer_id), d.primer_id
    assert set(POURCEL7) <= set(ORDER)


def test_the_inconsistent_size_range_is_recorded_with_the_published_product():
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    text = load_scheme(SCHEME).markers.read_text(encoding="utf-8")
    for marker, (size, published_range) in INCONSISTENT_RANGE.items():
        assert defs[marker].ref_size == size
        assert published_range in text, f"the size range of {marker} must stay on record"


def test_check_scheme_against_co92(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 25
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all(), "every published product must be reproduced"
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    assert (table["note"] == "").all(), "no primer mismatch anywhere in the reference"
    assert (table["repeats_ok"] == "").all(), "no marker stores a motif, so none is cross-checked"


def test_co92_profile_and_locus_positions(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "co92", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "co92_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    status = pd.read_csv(tmp_path / "co92_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "co92_calls_long.tsv", sep="\t").set_index("marker")
    assert len(long) == 25 and (long["residual_bp"] == 0).all()
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    assert {m: int(long.loc[m, "allele"]) for m in ORDER} == {m: defs[m].ref_repeats for m in ORDER}
    assert set(long["contig"]) == {CHROMOSOME}, "all 25 loci are chromosomal"
    # the full locus name gives the position in the CO92 chromosome in kb
    for marker in ORDER:
        kb = int(defs[marker].label[2:6])
        assert kb * 1000 <= int(long.loc[marker, "start"]) < (kb + 1) * 1000, marker


def test_the_shared_loci_agree_on_the_unit_but_not_on_the_copy_number():
    """The equivalences the source states hold for the unit length; the numbering differs."""
    here = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    there = {d.primer_id: d for d in read_markers(load_scheme(KLEVYTSKA).markers)}
    for mine, theirs in EQUIVALENT.items():
        a, b = here[mine], there[theirs]
        assert a.unit_length == b.unit_length, f"{mine}/{theirs}"
        assert a.ref_repeats != b.ref_repeats, f"{mine}/{theirs} must not look interchangeable"
        assert a.ref_size != b.ref_size, f"{mine}/{theirs} use different primers"
