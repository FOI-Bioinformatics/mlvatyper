"""The bundled Bacillus anthracis MLVA-31 panel of MLVAbank against Ames Ancestor.

Every reference number in this scheme is read from the MLVAbank identifier in
the label column, and all 31 are reproduced in the reference assembly, so these
tests are a complete check of the published values rather than of derived ones.
The scheme also supplies the published values that bacillus_mlva8_keim_2000 and
bacillus_mlva15_vanert_2007 use, so the three must agree at their shared loci.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import parse_marker_name, read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000008445.1"
SCHEME = "bacillus_mlva31_mlvabank"
MLVA8, MLVA15 = "bacillus_mlva8_keim_2000", "bacillus_mlva15_vanert_2007"
KEIM = ["vrrA", "vrrB1", "vrrB2", "vrrC1", "vrrC2", "CG3", "pXO1-aat", "pXO2-at"]
BAMS = ["bams01", "bams03", "bams05", "bams13", "bams15", "bams21", "bams22", "bams23", "bams24",
        "bams25", "bams28", "bams30", "bams31", "bams34", "bams44", "bams51", "bams53"]
VNTR = ["vntr12", "vntr16", "vntr17", "vntr19", "vntr23", "vntr35"]
ORDER = KEIM + BAMS + VNTR
PROFILE = "4-20-13-53-17-2-7-9-16-26-5-70-24-10-16-11-11-13-14-57-64-11-8-9-8-6-20-4-4-4-5"
# the locus in the source table that the panel name leaves out
HELD_OUT = "bams07"
CHROMOSOME, PXO1, PXO2 = "AE017334.2", "AE017336.2", "AE017335.3"
ON_PXO1, ON_PXO2 = {"pXO1-aat"}, {"pXO2-at", "vntr16", "vntr17"}
# a pattern is stored only where counting the published motif gives the published copy number
STORED_MOTIF = {"vrrA": "CAATATCAACAA", "CG3": "TAATA", "pXO1-aat": "AAT", "vntr12": "TA",
                "vntr17": "ATGTGTAA", "vntr19": "TAT", "vntr23": "ATGAAAAAGAAA", "vntr35": "GATTGT"}
# the five vntr reverse primers whose non-genomic 5' G mismatches; vntr19 happens to match
ONE_MISMATCH = {"vntr12", "vntr16", "vntr17", "vntr23", "vntr35"}
PRIMERS_INSIDE_ARRAY = {"vrrB1", "vrrB2"}
# the loci the three Bacillus schemes share, as named in each
SHARED_WITH_MLVA15 = {"vntr12": "BaVNTR12", "vntr16": "BaVNTR16", "vntr17": "BaVNTR17",
                      "vntr19": "BaVNTR19", "vntr23": "BaVNTR23", "vntr35": "BaVNTR35",
                      **{m: m for m in KEIM}}


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    assert "unit_mismatches" not in scheme.calling_defaults()
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    assert len(defs) == 31, "the panel name is MLVA-31"
    for d in defs:
        assert d.can_call, d.primer_id
        slack = d.ref_size - d.ref_repeats * d.unit_length
        assert slack >= 0, d.primer_id
        if d.primer_id in PRIMERS_INSIDE_ARRAY:
            assert slack < 52  # less than the two primers occupy
        if d.primer_id in STORED_MOTIF:
            assert d.patterns == (STORED_MOTIF[d.primer_id],), d.primer_id
            assert d.unit_length == len(d.patterns[0]), d.primer_id
        else:
            assert d.patterns == (), d.primer_id


def test_every_reference_number_comes_from_the_mlvabank_identifier():
    for d in read_markers(load_scheme(SCHEME).markers):
        encoded = parse_marker_name(d.label)
        assert encoded is not None, d.label
        assert d.label.startswith(d.primer_id + "_"), d.primer_id
        assert (d.unit_length, d.ref_size, d.ref_repeats) == \
               (encoded.unit_length, encoded.ref_size, encoded.ref_repeats), d.primer_id


def test_the_held_out_locus_is_recorded_but_not_in_the_panel():
    """bams07 is in the source table and commented out in both files, so it can be restored."""
    scheme = load_scheme(SCHEME)
    assert HELD_OUT not in {d.primer_id for d in read_markers(scheme.markers)}
    for path in (scheme.sts, scheme.markers):
        commented = [line for line in path.read_text().splitlines()
                     if line.startswith("#") and line.lstrip("#").startswith(HELD_OUT + "\t")]
        assert len(commented) == 1, path.name


def test_check_scheme_against_ames_ancestor(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 31
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all(), "every published reference size must be reproduced"
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    assert (table.loc[sorted(STORED_MOTIF), "repeats_ok"].astype(int) == 1).all()
    assert (table.loc[[m for m in ORDER if m not in STORED_MOTIF], "repeats_ok"] == "").all()
    assert {m for m in ORDER if table.loc[m, "note"]} == ONE_MISMATCH


def test_ames_ancestor_profile_and_replicons(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "ames", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "ames_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    status = pd.read_csv(tmp_path / "ames_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "ames_calls_long.tsv", sep="\t").set_index("marker")
    assert len(long) == 31 and (long["residual_bp"] == 0).all()
    # the call is the published copy number at every locus
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    assert {m: int(long.loc[m, "allele"]) for m in ORDER} == {m: defs[m].ref_repeats for m in ORDER}
    for marker in ORDER:
        expected = PXO1 if marker in ON_PXO1 else PXO2 if marker in ON_PXO2 else CHROMOSOME
        assert long.loc[marker, "contig"] == expected, marker


def test_the_three_bacillus_schemes_agree_at_their_shared_loci():
    """MLVA-8 and MLVA-15 take their published values from this panel, so they must match."""
    here = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    for name in (MLVA8, MLVA15):
        there = {d.primer_id: d for d in read_markers(load_scheme(name).markers)}
        for mine, theirs in SHARED_WITH_MLVA15.items():
            if theirs not in there:
                continue
            a, b = here[mine], there[theirs]
            assert (a.unit_length, a.ref_size, a.ref_repeats, a.patterns, a.label) == \
                   (b.unit_length, b.ref_size, b.ref_repeats, b.patterns, b.label), f"{name}:{theirs}"
    # BaVNTR32 is this panel's bams01 under other primers, so only its size differs
    other = {d.primer_id: d for d in read_markers(load_scheme(MLVA15).markers)}["BaVNTR32"]
    assert other.label == here["bams01"].label
    assert (other.unit_length, other.ref_repeats) == (here["bams01"].unit_length, here["bams01"].ref_repeats)
    assert other.ref_size != here["bams01"].ref_size
