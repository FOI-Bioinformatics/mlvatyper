"""The bundled Coxiella burnetii MLVA panels of Arricau-Bouvery et al. against Nine Mile RSA 493.

The unit lengths and copy numbers are the published values encoded in the locus
identifiers of the source table. Fourteen of the seventeen published amplicon
sizes are reproduced in the reference; the three that are not are pinned here
with their measured sizes, so a change in the search shows up as a test failure
rather than as a silent change of the allele grid.
"""

import re

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000007765.2"
SCHEME = "coxiella_mlva17_arricau_bouvery_2006"
PANEL1 = ["ms01", "ms03", "ms07", "ms12", "ms20", "ms21", "ms22", "ms26", "ms30", "ms36"]
PANEL2 = ["ms23", "ms24", "ms27", "ms28", "ms31", "ms33", "ms34"]
ORDER = PANEL1 + PANEL2
PROFILE = "4-7-8-8-15-6-6-4-6-4-8-27-4-6-5-9-5"
CHROMOSOME = "AE016828.3"
# the identifiers are written NAME_<unit>bp_<n>U_<size>bp, not the order parse_marker_name reads
LABEL_RE = re.compile(r"^(?P<gene>Cbu\d+)_(?P<ms>ms\d+)_(?P<unit>\d+)bp_(?P<n>\d+)U_(?P<size>\d+)bp$")
# the three loci where in silico PCR with these primers disagrees with the published size:
# marker -> (published size, measured size). markers.tsv sets out the evidence.
SIZE_DISAGREES = {"ms30": (215, 306), "ms31": (182, 285), "ms36": (447, 477)}
# the one locus whose detected repeat unit occurs the published number of times
STORED_PATTERN = {"ms31": "ACAGAGG"}
# ms20 and ms36 share a forward primer in Cbu1941, so the ms36 product contains the ms20 product
NESTED = ("ms20", "ms36")


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    assert "unit_mismatches" not in scheme.calling_defaults()
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    assert len(defs) == 17, "ten panel 1 loci and seven panel 2 loci"
    for d in defs:
        assert d.can_call, d.primer_id
        assert d.ref_size - d.ref_repeats * d.unit_length >= 0, d.primer_id
        if d.primer_id in STORED_PATTERN:
            assert d.patterns == (STORED_PATTERN[d.primer_id],), d.primer_id
            assert len(d.patterns[0]) == d.unit_length, d.primer_id
        else:
            assert d.patterns == (), d.primer_id
    # panel 1 holds the minisatellites, panel 2 the 6 and 7 bp microsatellites. The published
    # panel 1 heading says the unit is above 9 bp, but ms26 and ms36 are exactly 9 bp.
    by_name = {d.primer_id: d for d in defs}
    assert all(by_name[m].unit_length >= 9 for m in PANEL1)
    assert {m for m in PANEL1 if by_name[m].unit_length == 9} == {"ms26", "ms36"}
    assert all(by_name[m].unit_length in (6, 7) for m in PANEL2)


def test_the_published_identifier_supplies_the_unit_and_the_copy_number():
    for d in read_markers(load_scheme(SCHEME).markers):
        m = LABEL_RE.match(d.label)
        assert m is not None, d.label
        assert m["ms"] == d.primer_id
        assert (d.unit_length, d.ref_repeats) == (int(m["unit"]), int(m["n"])), d.primer_id
        published = int(m["size"])
        if d.primer_id in SIZE_DISAGREES:
            assert (published, d.ref_size) == SIZE_DISAGREES[d.primer_id], d.primer_id
        else:
            assert d.ref_size == published, d.primer_id


def test_the_two_disagreeing_sizes_are_off_the_published_repeat_grid():
    """No choice of reference values puts ms30 and ms31 on the grid the published sizes imply."""
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    for marker in ("ms30", "ms31"):
        published, measured = SIZE_DISAGREES[marker]
        assert (measured - published) % defs[marker].unit_length != 0, marker
    # ms36 differs by three bases less than four units, which markers.tsv reads as a transposed size
    assert SIZE_DISAGREES["ms36"] == (447, 477)


def test_check_scheme_against_nine_mile(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 17
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    assert (table["note"] == "").all(), "no primer mismatch anywhere in the reference"
    assert (table.loc[sorted(STORED_PATTERN), "repeats_ok"].astype(int) == 1).all()
    assert (table.loc[[m for m in ORDER if m not in STORED_PATTERN], "repeats_ok"] == "").all()


def test_nine_mile_profile_and_replicon(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "nm", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "nm_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    status = pd.read_csv(tmp_path / "nm_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "nm_calls_long.tsv", sep="\t").set_index("marker")
    assert len(long) == 17 and (long["residual_bp"] == 0).all()
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    assert {m: int(long.loc[m, "allele"]) for m in ORDER} == {m: defs[m].ref_repeats for m in ORDER}
    assert set(long["contig"]) == {CHROMOSOME}, "all seventeen loci are chromosomal"
    # the published digit code is panel 1 followed by panel 2, so the string splits at ms36
    digits = PROFILE.split("-")
    assert len(digits[:len(PANEL1)]) == 10 and len(digits[len(PANEL1):]) == 7


def test_ms36_contains_ms20(tmp_path):
    """The two Cbu1941 loci are not independent: the ms36 product spans the whole ms20 product."""
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "nm", "--align", "none", "--no-report"]) == 0
    hits = pd.read_csv(tmp_path / "hits.tsv", sep="\t").set_index("marker")
    inner, outer = (hits.loc[m] for m in NESTED)
    assert inner["contig"] == outer["contig"] == CHROMOSOME
    assert outer["start"] - inner["start"] == 2, "the forward primers are two bases apart"
    assert outer["end"] > inner["end"]
    assert int(outer["length"]) - int(inner["length"]) == 75
