"""The bundled Francisella tularensis Johansson scheme against the SCHU S4 assembly.

Every reference value is the published one: Table 2 of Johansson et al. 2004
gives the repeat motif, its size and the SCHU S4 copy number for all 25
markers, so the scheme calls every locus on the published scale. A stored
repeat pattern is kept only where counting the motif in the reference returns
the published number. Ft-M14 is duplicated in SCHU S4 and reports MULTIPLE.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000008985.1"
N_MARKERS = 25
DUPLICATED = "Ft-M14"  # two copies in SCHU S4, both 442 bp

# Table 2 of Johansson et al. 2004: unit length and SCHU S4 copy number
UNIT = {"Ft-M1": 3, "Ft-M2": 6, "Ft-M3": 9, "Ft-M4": 5, "Ft-M5": 16, "Ft-M6": 21, "Ft-M7": 16,
        "Ft-M8": 16, "Ft-M9": 16, "Ft-M10": 16, "Ft-M11": 10, "Ft-M12": 10, "Ft-M13": 12,
        "Ft-M14": 6, "Ft-M15": 6, "Ft-M16": 10, "Ft-M17": 6, "Ft-M18": 6, "Ft-M19": 13,
        "Ft-M20": 12, "Ft-M21": 7, "Ft-M22": 6, "Ft-M23": 23, "Ft-M24": 21, "Ft-M25": 2}
SCHU_COUNT = {"Ft-M1": 3, "Ft-M2": 4, "Ft-M3": 25, "Ft-M4": 3, "Ft-M5": 3, "Ft-M6": 4, "Ft-M7": 4,
              "Ft-M8": 4, "Ft-M9": 4, "Ft-M10": 18, "Ft-M11": 5, "Ft-M12": 2, "Ft-M13": 2,
              "Ft-M14": 3, "Ft-M15": 2, "Ft-M16": 2, "Ft-M17": 3, "Ft-M18": 4, "Ft-M19": 2,
              "Ft-M20": 3, "Ft-M21": 3, "Ft-M22": 2, "Ft-M23": 2, "Ft-M24": 1, "Ft-M25": 5}
ALL_MARKERS = list(UNIT)

# a stored pattern is kept only where the motif counts back to the published number
PATTERNED = ["Ft-M4", "Ft-M5", "Ft-M7", "Ft-M8", "Ft-M9", "Ft-M10", "Ft-M20", "Ft-M21"]
SHARED_UNIT = "AAAAAGATATTTGTAG"
SHARED_UNIT_MARKERS = ["Ft-M5", "Ft-M7", "Ft-M8", "Ft-M9", "Ft-M10"]

# eight counts the published table corrects from the earlier derived values
CORRECTED_COUNT = {"Ft-M2", "Ft-M3", "Ft-M6", "Ft-M14", "Ft-M15", "Ft-M16", "Ft-M17", "Ft-M18"}
# sizes the earlier audit established, kept here
CORRECTED_SIZE = {"Ft-M22": 241, "Ft-M25": 151, "Ft-M10": 617}


def test_scheme_files_are_consistent():
    scheme = load_scheme("francisella_johansson_2004")
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    assert scheme.known_multiple == (DUPLICATED,)
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    assert len(defs) == N_MARKERS and sorted(defs) == sorted(ALL_MARKERS)
    for marker in ALL_MARKERS:
        d = defs[marker]
        assert d.can_call, f"{marker} should call on the published scale"
        assert d.unit_length == UNIT[marker], marker
        assert d.ref_repeats == SCHU_COUNT[marker], marker
        assert d.ref_size - d.ref_repeats * d.unit_length >= 0, marker
    patterned = {m for m in ALL_MARKERS if defs[m].patterns}
    assert patterned == set(PATTERNED)
    assert {defs[m].patterns[0] for m in SHARED_UNIT_MARKERS} == {SHARED_UNIT}
    for marker, size in CORRECTED_SIZE.items():
        assert defs[marker].ref_size == size


def test_check_scheme_against_schu_s4(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", "francisella_johansson_2004", "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == N_MARKERS
    assert (table["size_ok"].astype(int) == 1).all(), "every published reference size reproduces"
    # Ft-M14 gives two products and is declared; the rest give one
    assert int(table.loc[DUPLICATED, "n_products"]) > 1
    assert (table.loc[[m for m in ALL_MARKERS if m != DUPLICATED], "n_products"].astype(int) == 1).all()
    # the eight patterned loci carry a tandem cross-check and it passes
    assert (table.loc[PATTERNED, "repeats_ok"].astype(int) == 1).all()


def test_schu_s4_run_returns_the_published_profile(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", "francisella_johansson_2004", "--fasta", str(fixture),
                 "--outdir", str(tmp_path), "--prefix", "schu", "--align", "none", "--no-report"]) == 0
    long = pd.read_csv(tmp_path / "schu_calls_long.tsv", sep="\t")
    assert (long["residual_bp"].dropna() == 0).all()
    # every locus calls its published SCHU S4 count; Ft-M14's two products both do
    by_marker = long.groupby("marker")["allele"].apply(lambda s: {int(float(a)) for a in s if pd.notna(a)})
    for marker in ALL_MARKERS:
        assert SCHU_COUNT[marker] in by_marker[marker], marker
    status = pd.read_csv(tmp_path / "schu_status.tsv", sep="\t", index_col=0, dtype=str)
    assert status.loc[ASSEMBLY, DUPLICATED] == "MULTIPLE"


def test_eight_counts_are_corrected_from_the_earlier_derived_values():
    """The published table moved these counts off the derived ones; the units were already right."""
    defs = {d.primer_id: d for d in read_markers(load_scheme("francisella_johansson_2004").markers)}
    # the corrected loci now hold the published count and have no stored pattern
    for marker in CORRECTED_COUNT:
        assert defs[marker].ref_repeats == SCHU_COUNT[marker], marker
        assert defs[marker].patterns == (), f"{marker} drops its pattern, which would disagree"
