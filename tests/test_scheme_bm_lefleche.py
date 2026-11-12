"""The bundled Brucella MLVA-15 scheme of Le Fleche et al. against the 16M assembly.

The reference numbers are read from the MLVAbank identifier in the label column
and agree with the publication's B. melitensis size column, so 16M is the anchor.
These tests hold the panel order, the reference sizes and the 16M profile.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import parse_marker_name, read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000007125.1"
PANEL_1 = ["Bruce06", "Bruce08", "Bruce11", "Bruce12", "Bruce42", "Bruce43", "Bruce45", "Bruce55"]
PANEL_2 = ["Bruce04", "Bruce07", "Bruce09", "Bruce16", "Bruce18", "Bruce21", "Bruce30"]
# the published profile order is panel 1 then panel 2
ORDER = PANEL_1 + PANEL_2
PROFILE = "3-4-2-13-4-2-3-3-2-5-7-3-5-6-6"
REF_SIZE = {"Bruce06": 408, "Bruce08": 348, "Bruce11": 257, "Bruce12": 392, "Bruce42": 539, "Bruce43": 182,
            "Bruce45": 151, "Bruce55": 273, "Bruce04": 152, "Bruce07": 158, "Bruce09": 156, "Bruce16": 152,
            "Bruce18": 146, "Bruce21": 148, "Bruce30": 151}
# loci the table places on chromosome II; the rest are on chromosome I
CHROMOSOME_II = {"Bruce12", "Bruce16", "Bruce18", "Bruce21"}
CHR_I, CHR_II = "AE008917.1", "AE008918.1"
# large units whose reference size leaves less room than the two primers need
PRIMERS_INSIDE_ARRAY = {"Bruce06", "Bruce42"}


def test_scheme_files_are_consistent():
    scheme = load_scheme("brucella_mlva15_lefleche_2006")
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER, "panel 1 must precede panel 2"
    assert len(defs) == 15
    for d in defs:
        assert d.can_call and d.patterns == ()  # no consensus units are published
        assert d.ref_size == REF_SIZE[d.primer_id]


def test_reference_numbers_match_the_mlvabank_identifier():
    """The label encodes unit length, reference size and copy number; the columns must agree with it."""
    defs = read_markers(load_scheme("brucella_mlva15_lefleche_2006").markers)
    for d in defs:
        assert d.label and d.label.startswith("BRU"), d.primer_id
        encoded = parse_marker_name(d.label)
        assert encoded is not None, d.label
        assert (d.unit_length, d.ref_size, d.ref_repeats) == \
               (encoded.unit_length, encoded.ref_size, encoded.ref_repeats), d.primer_id
        slack = d.ref_size - d.ref_repeats * d.unit_length
        assert slack >= 0
        if d.primer_id in PRIMERS_INSIDE_ARRAY:
            assert slack < 44  # the primers must sit inside the array for these two
        else:
            assert slack >= 44


def test_check_scheme_against_16m(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", "brucella_mlva15_lefleche_2006", "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 15
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    assert (table["note"] == "").all(), "no marker should report a primer mismatch in the reference"
    assert (table["repeats_ok"] == "").all()  # no patterns, so no tandem cross-check


def test_16m_profile_and_chromosomes(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", "brucella_mlva15_lefleche_2006", "--fasta", str(fixture), "--outdir",
                 str(tmp_path), "--prefix", "16m", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "16m_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    status = pd.read_csv(tmp_path / "16m_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"BINNED"}  # every 16M size is in bins.tsv
    long = pd.read_csv(tmp_path / "16m_calls_long.tsv", sep="\t").set_index("marker")
    assert (long["residual_bp"] == 0).all()
    # every locus sits on the chromosome the publication states
    for marker in ORDER:
        expected = CHR_II if marker in CHROMOSOME_II else CHR_I
        assert long.loc[marker, "contig"] == expected, marker
    # the reference call is the published copy number for each locus
    defs = {d.primer_id: d for d in read_markers(load_scheme("brucella_mlva15_lefleche_2006").markers)}
    assert {m: int(long.loc[m, "allele"]) for m in ORDER} == {m: defs[m].ref_repeats for m in ORDER}
