"""The bundled Bacillus anthracis MLVA-8 scheme of Keim et al. against Ames Ancestor.

The reference numbers are read from the MLVAbank identifier in the label column,
so they are published values; the primer table of Keim et al. itself carries
none. These tests pin them, the replicon each locus sits on and the reference
profile, so a change in the search or the tandem counter shows up here.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import parse_marker_name, read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000008445.1"
SCHEME = "bacillus_mlva8_keim_2000"
# the published marker order: six chromosomal loci, then the two plasmid loci
ORDER = ["vrrA", "vrrB1", "vrrB2", "vrrC1", "vrrC2", "CG3", "pXO1-aat", "pXO2-at"]
PROFILE = "4-20-13-53-17-2-7-9"
CHROMOSOME, PXO1, PXO2 = "AE017334.2", "AE017336.2", "AE017335.3"
# the loci whose published copy number a count of the published motif does not reproduce, so no
# pattern is stored and the cross-check is skipped rather than made to contradict the panel
NO_PATTERN = {"vrrB1", "vrrB2", "vrrC1", "vrrC2", "pXO2-at"}
# published copy numbers that leave less room than the two primers need, so the primers sit
# inside the repeat array
PRIMERS_INSIDE_ARRAY = {"vrrB1", "vrrB2"}
# the two adjacent paralogue pairs, whose amplicons overlap in the reference
OVERLAP = [("vrrB2", "vrrB1", 27), ("vrrC2", "vrrC1", 36)]


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    # the allowance stays at its default; counting the vrrC arrays would need one that is too loose
    assert "unit_mismatches" not in scheme.calling_defaults()
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.can_call, d.primer_id
        slack = d.ref_size - d.ref_repeats * d.unit_length
        assert slack >= 0, d.primer_id
        if d.primer_id in PRIMERS_INSIDE_ARRAY:
            assert slack < 52  # less than the two primers occupy
        else:
            assert slack >= 44
        if d.primer_id in NO_PATTERN:
            assert d.patterns == ()
        else:
            assert len(d.patterns) == 1 and len(d.patterns[0]) == d.unit_length, d.primer_id


def test_reference_numbers_match_the_mlvabank_identifier():
    """The label encodes unit length, reference size and copy number; the columns must agree."""
    for d in read_markers(load_scheme(SCHEME).markers):
        encoded = parse_marker_name(d.label)
        assert encoded is not None, d.label
        assert (d.unit_length, d.ref_size, d.ref_repeats) == \
               (encoded.unit_length, encoded.ref_size, encoded.ref_repeats), d.primer_id


def test_the_stored_patterns_are_published_motifs_that_the_count_reproduces():
    """A pattern is stored only where counting the published motif gives the published number."""
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    stored = {"vrrA": "CAATATCAACAA", "CG3": "TAATA", "pXO1-aat": "AAT"}
    for marker, motif in stored.items():
        assert defs[marker].patterns == (motif,), marker
        assert defs[marker].unit_length == len(motif), marker
    # pXO2-at has a published motif too, but counting AT gives 10 where the panel numbers 9
    assert defs["pXO2-at"].patterns == () and defs["pXO2-at"].ref_repeats == 9
    assert set(stored) | NO_PATTERN == set(ORDER)


def test_check_scheme_against_ames_ancestor(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == len(ORDER)
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    assert (table["note"] == "").all(), "no marker should report a primer mismatch in the reference"
    checked = [m for m in ORDER if m not in NO_PATTERN]
    assert (table.loc[checked, "repeats_ok"].astype(int) == 1).all()
    assert (table.loc[sorted(NO_PATTERN), "repeats_ok"] == "").all()


def test_ames_ancestor_profile_and_replicons(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "ames", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "ames_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    status = pd.read_csv(tmp_path / "ames_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "ames_calls_long.tsv", sep="\t").set_index("marker")
    assert len(long) == len(ORDER) and (long["residual_bp"] == 0).all()
    # each locus lies on the replicon its name implies
    assert long.loc["pXO1-aat", "contig"] == PXO1
    assert long.loc["pXO2-at", "contig"] == PXO2
    for marker in ORDER:
        if not marker.startswith("pXO"):
            assert long.loc[marker, "contig"] == CHROMOSOME, marker
    # the call is the derived reference copy number at every locus
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    assert {m: int(long.loc[m, "allele"]) for m in ORDER} == {m: defs[m].ref_repeats for m in ORDER}


def test_the_paralogue_pairs_overlap_in_the_reference(tmp_path):
    """vrrB1/vrrB2 and vrrC1/vrrC2 are adjacent, so a rearrangement there can affect both."""
    fixture = reference_assembly(ASSEMBLY)
    assert main(["search", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path)]) == 0
    hits = pd.read_csv(tmp_path / "hits.tsv", sep="\t").set_index("marker")
    for lower, upper, overlap in OVERLAP:
        assert hits.loc[lower, "contig"] == hits.loc[upper, "contig"] == CHROMOSOME
        assert int(hits.loc[lower, "end"]) - int(hits.loc[upper, "start"]) + 1 == overlap, lower
