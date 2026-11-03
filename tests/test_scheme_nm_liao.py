"""The bundled Neisseria meningitidis 12 locus MLVA of Liao et al. 2006 against Z2491.

This is the best specified of the bundled schemes. Table 1 gives the consensus
repeat unit, its length and the repeat number in each of three sequenced
genomes, and Table 5 gives the primers and the predicted amplicon size in the
same three. The two tables agree exactly: at every locus the size minus the
unit length times the repeat number is one flanking length shared by all three
genomes. These tests pin that arithmetic, the reproduction of the reference,
and the fact that the panel also reproduces a genome it is not anchored on.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from mlvatyper.sts import read_sts
from tests.conftest import reference_assembly

SCHEME = "nmeningitidis_mlva12_liao_2006"
Z2491, MC58 = "GCA_000009105.1", "GCA_000008805.1"
CHROMOSOME = "AL157959.1"

# the order the header of Table 2 gives for the MLVA allelic profile
ORDER = ["NMTR1", "NMTR2", "NMTR6", "NMTR7", "NMTR9", "NMTR9a",
         "NMTR9b", "NMTR9c", "NMTR10", "NMTR12", "NMTR18", "NMTR19"]
COMPOUND = "NMTR12"
CALLABLE = [m for m in ORDER if m != COMPOUND]

UNIT = {"NMTR1": 7, "NMTR2": 7, "NMTR6": 21, "NMTR7": 4, "NMTR9": 9, "NMTR9a": 30,
        "NMTR9b": 18, "NMTR9c": 15, "NMTR10": 6, "NMTR18": 9, "NMTR19": 13}

# Table 1 and Table 5, per genome: repeat number and predicted amplicon size
REPEATS = {
    Z2491: {"NMTR1": 25, "NMTR2": 17, "NMTR6": 2, "NMTR7": 12, "NMTR9": 9, "NMTR9a": 2,
            "NMTR9b": 2, "NMTR9c": 3, "NMTR10": 4, "NMTR18": 2, "NMTR19": 2},
    MC58: {"NMTR1": 35, "NMTR2": 4, "NMTR6": 2, "NMTR7": 11, "NMTR9": 7, "NMTR9a": 1,
           "NMTR9b": 2, "NMTR9c": 2, "NMTR10": 9, "NMTR18": 2, "NMTR19": 2},
}
SIZES = {
    Z2491: {"NMTR1": 351, "NMTR2": 327, "NMTR6": 165, "NMTR7": 231, "NMTR9": 243,
            "NMTR9a": 218, "NMTR9b": 200, "NMTR9c": 202, "NMTR10": 221, "NMTR12": 218,
            "NMTR18": 182, "NMTR19": 186},
    MC58: {"NMTR1": 421, "NMTR2": 236, "NMTR6": 165, "NMTR7": 227, "NMTR9": 225,
           "NMTR9a": 188, "NMTR9b": 200, "NMTR9c": 187, "NMTR10": 251, "NMTR12": 548,
           "NMTR18": 182, "NMTR19": 186},
}
# FAM18, the third genome, which is not bundled as an assembly but checks the arithmetic
FAM18_REPEATS = {"NMTR1": 29, "NMTR2": 19, "NMTR6": 4, "NMTR7": 8, "NMTR9": 26, "NMTR9a": 3,
                 "NMTR9b": 1, "NMTR9c": 3, "NMTR10": 4, "NMTR18": 5, "NMTR19": 5}
FAM18_SIZES = {"NMTR1": 379, "NMTR2": 341, "NMTR6": 207, "NMTR7": 215, "NMTR9": 396,
               "NMTR9a": 248, "NMTR9b": 182, "NMTR9c": 202, "NMTR10": 221, "NMTR18": 209,
               "NMTR19": 225}

PROFILE = {Z2491: "25-17-2-12-9-2-2-3-4-NO_DEFINITION-2-2",
           MC58: "35-4-2-11-7-1-2-2-9-NO_DEFINITION-2-2"}

# the motif is stored only where counting it in Z2491 returns the published number
NO_PATTERN = {"NMTR1", "NMTR6", COMPOUND}
DYES = ("6-FAM", "VIC", "NED", "PET")


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == Z2491
    pairs, _ = read_sts(scheme.sts)
    assert [p.id for p in pairs] == ORDER
    assert len(pairs) == 12
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.ref_size == SIZES[Z2491][d.primer_id], d.primer_id
        if d.primer_id == COMPOUND:
            continue
        assert d.can_call, d.primer_id
        assert d.unit_length == UNIT[d.primer_id], d.primer_id
        assert d.ref_repeats == REPEATS[Z2491][d.primer_id], d.primer_id


def test_the_compound_locus_is_sized_only():
    """NMTR12 mixes a 13 bp and a 12 bp unit, which the single unit formula cannot express."""
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    compound = defs[COMPOUND]
    assert not compound.can_call
    assert compound.unit_length == 0 and compound.ref_repeats is None
    assert compound.ref_size == SIZES[Z2491][COMPOUND]
    assert all(defs[m].can_call for m in CALLABLE)


def test_the_two_source_tables_agree_in_all_three_genomes():
    """size - unit * repeats is one flanking length per locus, shared by Z2491, MC58 and FAM18."""
    for marker in CALLABLE:
        unit = UNIT[marker]
        flanks = {g: SIZES[g][marker] - unit * REPEATS[g][marker] for g in (Z2491, MC58)}
        flanks["FAM18"] = FAM18_SIZES[marker] - unit * FAM18_REPEATS[marker]
        assert len(set(flanks.values())) == 1, f"{marker}: {flanks}"
        assert next(iter(flanks.values())) > 0, marker
    # the compound locus too, counted as a 13 bp and a 12 bp unit together
    compound = {Z2491: (2, 1), MC58: (20, 9), "FAM18": (13, 11)}
    sizes = {Z2491: SIZES[Z2491][COMPOUND], MC58: SIZES[MC58][COMPOUND], "FAM18": 481}
    flanks = {g: sizes[g] - (a * 13 + b * 12) for g, (a, b) in compound.items()}
    assert set(flanks.values()) == {180}


def test_patterns_are_stored_only_where_they_count_back():
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    without = {m for m, d in defs.items() if d.patterns == ()}
    assert without == NO_PATTERN
    for marker, d in defs.items():
        if marker not in NO_PATTERN:
            assert len(d.patterns) == 1 and len(d.patterns[0]) == d.unit_length, marker


def test_no_primer_keeps_its_dye_label():
    for pair in read_sts(load_scheme(SCHEME).sts)[0]:
        for seq in (pair.forward, pair.reverse):
            assert set(seq) <= set("ACGTRYSWKMBDHVN"), pair.id
            assert not any(seq.startswith(dye) for dye in DYES), pair.id


def test_check_scheme_against_z2491(tmp_path):
    fixture = reference_assembly(Z2491)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 12
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all(), "all twelve published sizes must reproduce"
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    counted = [m for m in ORDER if m not in NO_PATTERN]
    assert (table.loc[counted, "repeats_ok"].astype(int) == 1).all()
    assert (table.loc[sorted(NO_PATTERN), "repeats_ok"] == "").all()


def test_z2491_profile_is_the_published_repeat_numbers(tmp_path):
    fixture = reference_assembly(Z2491)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "z", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "z_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[Z2491, "profile"] == PROFILE[Z2491]
    long = pd.read_csv(tmp_path / "z_calls_long.tsv", sep="\t").set_index("marker")
    assert set(long["contig"]) == {CHROMOSOME}
    assert (long.loc[CALLABLE, "residual_bp"] == 0).all()
    assert {m: float(long.loc[m, "allele"]) for m in CALLABLE} == \
           {m: float(REPEATS[Z2491][m]) for m in CALLABLE}


def test_mc58_reproduces_its_published_numbers_too(tmp_path):
    """MC58 is not the reference, so typing it checks the panel rather than the anchor."""
    fixture = reference_assembly(MC58)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "mc", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "mc_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[MC58, "profile"] == PROFILE[MC58]
    long = pd.read_csv(tmp_path / "mc_calls_long.tsv", sep="\t").set_index("marker")
    assert {m: int(long.loc[m, "length"]) for m in ORDER} == SIZES[MC58]
    assert {m: float(long.loc[m, "allele"]) for m in CALLABLE} == \
           {m: float(REPEATS[MC58][m]) for m in CALLABLE}
