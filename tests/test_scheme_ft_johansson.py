"""The bundled Francisella tularensis Johansson scheme against the SCHU S4 assembly.

The reference sizes come from the publication's genomic coordinates and are
confirmed by in silico PCR. The repeat units are not published in the supplied
table and are filled from the reference only where the array is strong, so the
scheme calls 16 of its 25 markers and sizes the rest. These tests hold that
split in place, together with the two primer corrections the audit established.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000008985.1"
N_MARKERS = 25
CALLABLE = ["Ft-M2", "Ft-M3", "Ft-M4", "Ft-M5", "Ft-M6", "Ft-M7", "Ft-M8", "Ft-M9", "Ft-M10", "Ft-M14",
            "Ft-M15", "Ft-M16", "Ft-M17", "Ft-M18", "Ft-M20", "Ft-M21"]
SIZED_ONLY = ["Ft-M1", "Ft-M11", "Ft-M12", "Ft-M13", "Ft-M19", "Ft-M22", "Ft-M23", "Ft-M24", "Ft-M25"]
DUPLICATED = "Ft-M14"  # two copies in SCHU S4, both 442 bp
# the sizes the audit corrected, and the shared repeat family
CORRECTED = {"Ft-M22": 241, "Ft-M25": 151, "Ft-M10": 617}
SHARED_UNIT = "AAAAAGATATTTGTAG"
SHARED_UNIT_MARKERS = ["Ft-M5", "Ft-M7", "Ft-M8", "Ft-M9", "Ft-M10"]
# Genomic location column of the publication, with the abbreviated end coordinate expanded. The span is
# the amplicon, outer primer edge to outer primer edge, so it is what ref_size must equal.
PUBLISHED_SPAN = {"Ft-M1": 310, "Ft-M2": 463, "Ft-M3": 432, "Ft-M4": 434, "Ft-M5": 208, "Ft-M6": 311,
                  "Ft-M7": 207, "Ft-M8": 286, "Ft-M9": 248, "Ft-M10": 646, "Ft-M11": 392, "Ft-M12": 345,
                  "Ft-M13": 325, "Ft-M14": 442, "Ft-M15": 222, "Ft-M16": 221, "Ft-M17": 351, "Ft-M18": 372,
                  "Ft-M19": 316, "Ft-M20": 255, "Ft-M21": 403, "Ft-M22": 241, "Ft-M23": 435, "Ft-M24": 459,
                  "Ft-M25": 151}
# the only marker whose published span the assembly cannot produce; 617 bp is the single possible product
SPAN_DISAGREES = "Ft-M10"


def test_scheme_files_are_consistent():
    scheme = load_scheme("francisella_johansson_2004")
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    assert scheme.known_multiple == (DUPLICATED,)
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    assert len(defs) == N_MARKERS
    assert sorted(defs) == sorted(CALLABLE + SIZED_ONLY)
    for marker in CALLABLE:
        d = defs[marker]
        assert d.can_call, f"{marker} should carry a unit and a repeat count"
        assert d.ref_size - d.ref_repeats * d.unit_length >= 0
        assert len(d.patterns) == 1 and len(d.patterns[0]) == d.unit_length >= 5
        assert d.ref_repeats >= 3
    for marker in SIZED_ONLY:
        d = defs[marker]
        assert not d.can_call and d.patterns == () and d.ref_repeats is None
        assert d.ref_size is not None, f"{marker} must still carry a reference size"
    for marker, size in CORRECTED.items():
        assert defs[marker].ref_size == size
    assert {defs[m].patterns[0] for m in SHARED_UNIT_MARKERS} == {SHARED_UNIT}


def test_reference_sizes_are_the_published_spans():
    """ref_size equals the publication's genomic span for every marker but Ft-M10."""
    defs = {d.primer_id: d for d in read_markers(load_scheme("francisella_johansson_2004").markers)}
    assert set(defs) == set(PUBLISHED_SPAN)
    agree = {m: s for m, s in PUBLISHED_SPAN.items() if m != SPAN_DISAGREES}
    assert {m: defs[m].ref_size for m in agree} == agree
    # Ft-M10 keeps the measured product, which is shorter than the published span
    assert defs[SPAN_DISAGREES].ref_size == CORRECTED[SPAN_DISAGREES] < PUBLISHED_SPAN[SPAN_DISAGREES]
    # the difference is not a whole number of repeat units, so it is not an allele difference
    unit = defs[SPAN_DISAGREES].unit_length
    assert (PUBLISHED_SPAN[SPAN_DISAGREES] - CORRECTED[SPAN_DISAGREES]) % unit != 0


def test_the_two_corrected_primers_carry_no_extra_bases():
    """Ft-M22 lost a 5' G and Ft-M25 lost an 8 base 5' addition; both are genomic as written now."""
    scheme = load_scheme("francisella_johansson_2004")
    rows = {line.split("\t")[0]: line.split("\t")[1:3]
            for line in scheme.sts.read_text().splitlines() if line.strip() and not line.startswith("#")}
    assert rows["Ft-M22"][0] == "TCAAAATCTCAAGATGAGCAAATATTTGAATGGT"
    assert rows["Ft-M25"][1] == "ATCCCATATGTAAGTACAAATGTAGC"
    # the five markers that share one reverse primer, which is not locus-specific
    shared = {rows[m][1] for m in SHARED_UNIT_MARKERS}
    assert shared == {"CAGCTCGAACTCCGTCATAC"}


def test_check_scheme_against_schu_s4(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", "francisella_johansson_2004", "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == N_MARKERS
    assert (table["size_ok"].astype(int) == 1).all(), "every marker must reproduce its reference size"
    assert (table.loc[CALLABLE, "repeats_ok"].astype(int) == 1).all()
    assert (table.loc[SIZED_ONLY, "repeats_ok"] == "").all()  # no unit, so no tandem cross-check
    # the duplicated locus is the only marker with more than one product, and it is declared
    multiple = table.index[table["n_products"].astype(int) != 1].tolist()
    assert multiple == [DUPLICATED]
    assert multiple == list(load_scheme("francisella_johansson_2004").known_multiple)


def test_schu_s4_run_calls_the_strong_markers_and_flags_the_duplicate(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", "francisella_johansson_2004", "--fasta", str(fixture), "--outdir",
                 str(tmp_path), "--prefix", "schu", "--align", "none", "--no-report"]) == 0
    status = pd.read_csv(tmp_path / "schu_status.tsv", sep="\t", index_col=0, dtype=str)
    row = status.loc[ASSEMBLY]
    assert set(row[SIZED_ONLY]) == {"NO_DEFINITION"}
    assert "MULTIPLE" in row[DUPLICATED]
    unambiguous = [m for m in CALLABLE if m != DUPLICATED]
    assert set(row[unambiguous]) == {"OK"}
    alleles = pd.read_csv(tmp_path / "schu_alleles.tsv", sep="\t", index_col=0, dtype=str)
    defs = {d.primer_id: d for d in read_markers(load_scheme("francisella_johansson_2004").markers)}
    # the reference values were measured here, so each call must return its own repeat count
    assert {m: alleles.loc[ASSEMBLY, m] for m in unambiguous} == \
           {m: str(defs[m].ref_repeats) for m in unambiguous}
    long = pd.read_csv(tmp_path / "schu_calls_long.tsv", sep="\t")
    called = long[long.marker.isin(unambiguous)]
    assert (called["residual_bp"] == 0).all()
    assert (called["tandem_copies"] == called["allele"].astype(int)).all()
