"""The bundled Bacillus anthracis MLVA-15 scheme of Van Ert et al. against Ames Ancestor.

MLVA-15 carries the eight MLVA-8 loci of Keim et al. forward with the same
primers, so the two bundles must agree at those loci even though their orders
differ. The source table publishes a repeat motif and a location per locus but
no allele definitions, so the copy numbers are derived and the scheme notes say
so; these tests pin them, and they pin the comparison with the published spans.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000008445.1"
SCHEME = "bacillus_mlva15_vanert_2007"
MLVA8 = "bacillus_mlva8_keim_2000"
# the order of the published primer table, which is by multiplex reaction
ORDER = ["BaVNTR12", "BaVNTR19", "BaVNTR35", "BaVNTR16", "BaVNTR23", "BaVNTR17", "BaVNTR32",
         "vrrA", "vrrB1", "CG3", "pXO1-aat", "pXO2-at", "vrrB2", "vrrC1", "vrrC2"]
SHARED = ["vrrA", "vrrB1", "vrrB2", "vrrC1", "vrrC2", "CG3", "pXO1-aat", "pXO2-at"]
PROFILE = "6-4-5-4-4-4-11-4-3-2-7-10-2-7-8"
REF_SIZE = {"BaVNTR12": 115, "BaVNTR19": 96, "BaVNTR35": 115, "BaVNTR16": 273, "BaVNTR23": 197,
            "BaVNTR17": 386, "BaVNTR32": 577, "vrrA": 314, "vrrB1": 229, "CG3": 158,
            "pXO1-aat": 126, "pXO2-at": 141, "vrrB2": 153, "vrrC1": 580, "vrrC2": 532}
# the repeat motifs the publication prints, with BaVNTR32 given only as Complex
PUBLISHED_MOTIF = {"BaVNTR12": "TA", "BaVNTR19": "TAT", "BaVNTR35": "GATTGT", "BaVNTR16": "ACATTTAC",
                   "BaVNTR23": "ATGAAAAAGAAA", "BaVNTR17": "ATGTGTAA", "vrrA": "CAATATCAACAA",
                   "CG3": "TAATA", "pXO1-aat": "AAT", "pXO2-at": "AT"}
# units derived from the reference: BaVNTR32 is printed as Complex, vrrB1 as a 9 bp degenerate
# repeat and vrrB2 as Complex, so none of the three has a motif to take from the table
DERIVED_MOTIF = {"BaVNTR32": "GAAGAAACAGAAAAGCCAGGT", "vrrB1": "CAAGGTCAC", "vrrB2": "CAACAGCAC"}
# arrays too imperfect to count from a stored pattern, so the cross-check is skipped
NO_PATTERN = {"vrrC1", "vrrC2"}
# span of the published location, which is in the A2012 numbering
PUBLISHED_SPAN = {"BaVNTR12": 114, "BaVNTR19": 124, "BaVNTR35": 114, "BaVNTR16": 272,
                  "BaVNTR23": 197, "BaVNTR17": 385, "BaVNTR32": 577, "vrrA": 314, "vrrB1": 229,
                  "CG3": 158, "pXO1-aat": 126, "pXO2-at": 141, "vrrB2": 153, "vrrC1": 580,
                  "vrrC2": 532}
# the new reverse primers carry a 5' G the genome lacks; in BaVNTR19 it happens to match
ONE_MISMATCH = {"BaVNTR12", "BaVNTR35", "BaVNTR16", "BaVNTR23", "BaVNTR17", "BaVNTR32"}
CHROMOSOME, PXO1, PXO2 = "AE017334.2", "AE017336.2", "AE017335.3"
ON_PXO1, ON_PXO2 = {"pXO1-aat"}, {"BaVNTR16", "BaVNTR17", "pXO2-at"}


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    assert "unit_mismatches" not in scheme.calling_defaults()
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.can_call, d.primer_id
        assert d.ref_size == REF_SIZE[d.primer_id], d.primer_id
        assert d.ref_size - d.ref_repeats * d.unit_length >= 41, d.primer_id  # room for both primers
        if d.primer_id in NO_PATTERN:
            assert d.patterns == ()
        else:
            motif = {**PUBLISHED_MOTIF, **DERIVED_MOTIF}[d.primer_id]
            assert d.patterns == (motif,), d.primer_id
            assert d.unit_length == len(motif), d.primer_id


def test_the_stored_units_are_the_published_motifs():
    """Only BaVNTR32, printed as Complex, has a unit derived from the reference."""
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    for marker, motif in PUBLISHED_MOTIF.items():
        assert defs[marker].patterns == (motif,), marker
    assert set(DERIVED_MOTIF) | set(PUBLISHED_MOTIF) | NO_PATTERN == set(ORDER)


def test_mlva8_is_an_exact_subset():
    """The eight older loci keep their primers and reference values, in a different order."""
    by_name = {}
    for name in (SCHEME, MLVA8):
        scheme = load_scheme(name)
        sts = {line.split("\t")[0]: line.split("\t")[1:3]
               for line in scheme.sts.read_text().splitlines() if line.strip() and not line.startswith("#")}
        by_name[name] = (sts, {d.primer_id: d for d in read_markers(scheme.markers)})
    sts15, defs15 = by_name[SCHEME]
    sts8, defs8 = by_name[MLVA8]
    assert set(sts8) == set(SHARED) and set(sts15) - set(sts8) == set(ORDER[:7])
    for marker in SHARED:
        assert sts15[marker] == sts8[marker], marker
        a, b = defs15[marker], defs8[marker]
        assert (a.unit_length, a.ref_size, a.ref_repeats, a.patterns) == \
               (b.unit_length, b.ref_size, b.ref_repeats, b.patterns), marker
    # the two orders differ, so the profile strings are not aligned
    assert [m for m in ORDER if m in SHARED] != SHARED


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
    checked = [m for m in ORDER if m not in NO_PATTERN]
    assert (table.loc[checked, "repeats_ok"].astype(int) == 1).all()
    assert (table.loc[sorted(NO_PATTERN), "repeats_ok"] == "").all()
    # the non-genomic 5' G of the new reverse primers, and nowhere else
    flagged = {m for m in ORDER if table.loc[m, "note"]}
    assert flagged == ONE_MISMATCH


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
    # each locus lies on the replicon the publication states
    for marker in ORDER:
        expected = PXO1 if marker in ON_PXO1 else PXO2 if marker in ON_PXO2 else CHROMOSOME
        assert long.loc[marker, "contig"] == expected, marker
    # only the new reverse primers mismatch, by one base each
    assert (long["fwd_mismatches"] == 0).all()
    assert {m for m in ORDER if long.loc[m, "rev_mismatches"] == 1} == ONE_MISMATCH


def test_sizes_against_the_published_locations():
    """The published spans are A2012 coordinates, so they are compared by length only.

    The eight older loci match exactly. Of the new ones, the product is one base longer than
    the genomic span because of the non-genomic 5' G, and the publication includes that base
    in some spans and not in others, so a difference of one base either way is expected.
    BaVNTR19 is the known exception.
    """
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    for marker in SHARED:
        assert defs[marker].ref_size == PUBLISHED_SPAN[marker], marker
    for marker in ORDER[:7]:
        delta = defs[marker].ref_size - PUBLISHED_SPAN[marker]
        if marker == "BaVNTR19":
            assert delta == -28
        else:
            assert delta in (0, 1), marker
