"""The bundled Yersinia pestis Klevytska scheme against the CO92 assembly.

This is a 42 locus marker survey rather than a standardised panel. The unit
lengths and copy numbers are the publication's CO92 array column; the reference
sizes are measured here, because the paper gives no product sizes. Three loci
are numbered with a partial unit as n.5, so their alleles lie on a half-integer
grid, and these tests hold that behaviour in place.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.models import format_allele
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000009065.1"
N_MARKERS = 42
# published as unit by copies, with a partial unit written as n.5
HALF_COPY = {"M23": 6.5, "M58": 7.5, "M65": 2.5}
UNIT_LENGTH = {"M02": 1, "M06": 2, "M09": 3, "M12": 4, "M15": 5, "M18": 6, "M2": 7, "M23": 7, "M33": 9,
               "M36": 10, "M73": 30, "M74": 36, "M76": 41, "M77": 45, "M79": 8}
REF_SIZE = {"M02": 210, "M23": 169, "M58": 395, "M65": 316, "M77": 591, "M79": 252, "M09": 109}
# M02 and M2 are two different loci in the publication, not a padding mistake
DISTINCT_IDS = ("M02", "M2")


def test_scheme_files_are_consistent():
    scheme = load_scheme("yersinia_klevytska_2001")
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    assert scheme.calling_defaults()["rounding"] == "intermediate"  # needed by the half-copy loci
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    assert len(defs) == N_MARKERS
    assert all(i in defs for i in DISTINCT_IDS)
    assert defs["M02"].unit_length == 1 and defs["M2"].unit_length == 7
    for marker, unit in UNIT_LENGTH.items():
        assert defs[marker].unit_length == unit, marker
    for marker, size in REF_SIZE.items():
        assert defs[marker].ref_size == size, marker
    for d in defs.values():
        assert d.can_call and d.patterns == ()  # no consensus units are published
        assert d.ref_size - d.ref_repeats * d.unit_length >= 0


def test_half_copy_loci_keep_their_fractional_reference():
    defs = {d.primer_id: d for d in read_markers(load_scheme("yersinia_klevytska_2001").markers)}
    for marker, repeats in HALF_COPY.items():
        assert defs[marker].ref_repeats == repeats and defs[marker].grid_offset == 0.5, marker
    whole = [d for d in defs.values() if d.primer_id not in HALF_COPY]
    assert all(d.grid_offset == 0.0 for d in whole)
    assert len(whole) == N_MARKERS - len(HALF_COPY)


def test_check_scheme_against_co92(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", "yersinia_klevytska_2001", "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == N_MARKERS
    # the generous 60-1500 search window still gives one product per marker
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    assert (table["repeats_ok"] == "").all()  # no patterns, so no tandem cross-check
    for marker, size in REF_SIZE.items():
        assert int(table.loc[marker, "length"]) == size


def test_co92_calls_every_marker_and_numbers_the_half_copy_loci_on_their_grid(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", "yersinia_klevytska_2001", "--fasta", str(fixture), "--outdir",
                 str(tmp_path), "--prefix", "co92", "--align", "none", "--no-report"]) == 0
    status = pd.read_csv(tmp_path / "co92_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}, "the reference must call cleanly on every marker"
    alleles = pd.read_csv(tmp_path / "co92_alleles.tsv", sep="\t", index_col=0, dtype=str)
    defs = {d.primer_id: d for d in read_markers(load_scheme("yersinia_klevytska_2001").markers)}
    # reference sizes are measured here, so each call returns the published copy number
    for marker, d in defs.items():
        assert alleles.loc[ASSEMBLY, marker] == format_allele(d.ref_repeats), marker
    # the half-copy loci report n.5 and are not flagged as intermediate
    for marker, repeats in HALF_COPY.items():
        assert alleles.loc[ASSEMBLY, marker] == str(repeats)
        assert status.loc[ASSEMBLY, marker] == "OK"
    long = pd.read_csv(tmp_path / "co92_calls_long.tsv", sep="\t")
    assert len(long) == N_MARKERS and (long["residual_bp"] == 0).all()
