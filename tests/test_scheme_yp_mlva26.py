"""The bundled Yersinia pestis 14+12 scheme of Li et al. 2013 against CO92.

Primers and reference values come from Table S2 of the paper, which encodes the
repeat unit length, the CO92 amplicon size and the CO92 repeat number in every
locus name. Five of the encoded sizes are not the product their own primers
give in CO92, and in each case the measured length matches what an
independently bundled scheme records for the same locus with the same primers,
so the measured size is used and the published identifier is kept in the label.
These tests pin that correction, the membership of the two typing steps, and
the places where this panel disagrees with the other two Yersinia schemes.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import parse_marker_name, read_markers
from mlvatyper.scheme import load_scheme
from mlvatyper.sts import read_sts
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000009065.1"
SCHEME = "yersinia_mlva26_li_2013"
KLEVYTSKA, LEFLECHE = "yersinia_klevytska_2001", "yersinia_mlva25_lefleche_2001"
CHROMOSOME = "AL590842.1"

FIRST14 = ["yp0718ms41", "yp3057ms09", "yp3060ms56", "yp3236ms73", "yp4425ms38",
           "M15", "M21", "M58", "M61", "N0865", "N2486", "N2896", "N2976", "N3779"]
SECOND12 = ["M22", "M23", "M25", "M28", "M29", "M33", "M34", "M43",
            "N1606", "N2117", "N2577", "N3773"]
ORDER = FIRST14 + SECOND12

PROFILE = "7-33-7-6-8-5-5-8-5-2-3-3-3-4-13-7-22-7-8-22-6-5-4-3-7-5"

# Table 1 gives a motif length for each locus independently of the Table S2 locus name
MOTIF_LENGTH = {"yp0718ms41": 17, "yp3057ms09": 18, "yp3060ms56": 16, "yp3236ms73": 18,
                "yp4425ms38": 16, "M15": 5, "M21": 7, "M58": 17, "M61": 18, "N0865": 18,
                "N2486": 17, "N2896": 16, "N2976": 18, "N3779": 17, "M22": 7, "M23": 7,
                "M25": 7, "M28": 8, "M29": 8, "M33": 9, "M34": 9, "M43": 12, "N1606": 9,
                "N2117": 9, "N2577": 7, "N3773": 9}

# locus -> (size measured in CO92, size encoded in the Table S2 name, corroborating scheme)
CORRECTED = {"yp0718ms41": (217, 218, LEFLECHE), "yp4425ms38": (233, 251, LEFLECHE),
             "M61": (356, 357, KLEVYTSKA), "M29": (231, 233, KLEVYTSKA),
             "M34": (225, 255, KLEVYTSKA)}

# the Pourcel name of the loci this panel shares with the Le Fleche panel
POURCEL_ALIAS = {"yp0718ms41": "ms41", "yp3057ms09": "ms09", "yp3060ms56": "ms56",
                 "yp3236ms73": "ms73", "yp4425ms38": "ms38", "M58": "ms04", "M34": "ms62"}
# the two that keep a Pourcel name but are amplified with the Klevytska primers
DIFFERENT_AMPLICON = {"M58", "M34"}

# shared M loci this panel numbers differently from the 42 locus survey: locus -> (here, there)
RENUMBERED = {"M15": (5, 2), "M22": (13, 12), "M23": (7, 6.5), "M25": (22, 20), "M58": (8, 7.5)}

ONE_MISMATCH = {"M29"}


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    pairs, _ = read_sts(scheme.sts)
    assert [p.id for p in pairs] == ORDER
    assert len(FIRST14) == 14 and len(SECOND12) == 12 and len(pairs) == 26
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.can_call, d.primer_id
        assert d.patterns == (), "no repeat unit sequences are published for these loci"
        assert d.ref_size - d.unit_length * d.ref_repeats >= 0, d.primer_id


def test_the_first_fourteen_are_the_first_typing_step():
    """The panel order puts the 14 first, so a profile truncates to the first step result."""
    pairs, _ = read_sts(load_scheme(SCHEME).sts)
    assert [p.id for p in pairs][:14] == FIRST14
    assert [p.id for p in pairs][14:] == SECOND12
    assert not set(FIRST14) & set(SECOND12)
    assert PROFILE.count("-") == 25
    assert len(PROFILE.split("-")[:14]) == 14


def test_unit_lengths_agree_with_the_motif_lengths_of_table_1():
    for d in read_markers(load_scheme(SCHEME).markers):
        assert d.unit_length == MOTIF_LENGTH[d.primer_id], d.primer_id


def test_labels_keep_the_published_identifier_and_its_repeat_number():
    """The label is the Table S2 name; only its size is overruled, and only at five loci."""
    differ = {}
    for d in read_markers(load_scheme(SCHEME).markers):
        encoded = parse_marker_name(d.label)
        assert encoded is not None, d.label
        assert d.label.startswith(d.primer_id + "_"), d.primer_id
        assert encoded.unit_length == d.unit_length, d.primer_id
        assert encoded.ref_repeats == d.ref_repeats, d.primer_id
        if encoded.ref_size != d.ref_size:
            differ[d.primer_id] = (d.ref_size, encoded.ref_size)
    assert differ == {m: (measured, published) for m, (measured, published, _) in CORRECTED.items()}


def test_each_corrected_size_is_corroborated_by_another_bundled_scheme():
    """The measured length must equal what another scheme records with the same primer pair."""
    here_sts = {p.id: p for p in read_sts(load_scheme(SCHEME).sts)[0]}
    here = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    for marker, (measured, published, other) in CORRECTED.items():
        assert measured != published
        name = POURCEL_ALIAS.get(marker, marker) if other == LEFLECHE else marker
        there_sts = {p.id: p for p in read_sts(load_scheme(other).sts)[0]}
        there = {d.primer_id: d for d in read_markers(load_scheme(other).markers)}
        assert (here_sts[marker].forward, here_sts[marker].reverse) == \
               (there_sts[name].forward, there_sts[name].reverse), f"{marker}: primers must match"
        assert there[name].ref_size == measured, marker
        # the repeat number agrees too, so correcting the size shifts no allele
        assert there[name].ref_repeats == here[marker].ref_repeats, marker


def test_where_the_pourcel_panel_uses_other_primers_the_amplicon_differs():
    here_sts = {p.id: p for p in read_sts(load_scheme(SCHEME).sts)[0]}
    here = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    there_sts = {p.id: p for p in read_sts(load_scheme(LEFLECHE).sts)[0]}
    there = {d.primer_id: d for d in read_markers(load_scheme(LEFLECHE).markers)}
    for marker, alias in POURCEL_ALIAS.items():
        same_primers = (here_sts[marker].forward, here_sts[marker].reverse) == \
                       (there_sts[alias].forward, there_sts[alias].reverse)
        if marker in DIFFERENT_AMPLICON:
            assert not same_primers, marker
            assert here[marker].ref_size != there[alias].ref_size, marker
        else:
            assert same_primers, marker
            assert here[marker].ref_size == there[alias].ref_size, marker


def test_the_repeat_numbering_differs_from_the_klevytska_survey():
    here = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    there = {d.primer_id: d for d in read_markers(load_scheme(KLEVYTSKA).markers)}
    shared = [m for m in ORDER if m in there]
    assert len(shared) == 11, "M21 is not in the 42 locus survey"
    differ = {m: (here[m].ref_repeats, there[m].ref_repeats)
              for m in shared if here[m].ref_repeats != there[m].ref_repeats}
    assert differ == RENUMBERED
    # the sizes agree at every shared locus, so only the counting convention differs
    for marker in shared:
        assert here[marker].ref_size == there[marker].ref_size, marker


def test_check_scheme_against_co92(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 26
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all(), "every reference size must be reproduced"
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    assert {m for m in ORDER if table.loc[m, "note"]} == ONE_MISMATCH


def test_co92_profile_is_the_published_repeat_numbers(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "co92", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "co92_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    status = pd.read_csv(tmp_path / "co92_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "co92_calls_long.tsv", sep="\t").set_index("marker")
    assert len(long) == 26 and (long["residual_bp"] == 0).all()
    assert set(long["contig"]) == {CHROMOSOME}, "all 26 loci are chromosomal"
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    assert {m: float(long.loc[m, "allele"]) for m in ORDER} == \
           {m: float(defs[m].ref_repeats) for m in ORDER}
