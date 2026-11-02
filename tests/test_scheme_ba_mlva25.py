"""The bundled Bacillus anthracis MLVA-25 panel of Lista et al. 2006 against Ames Ancestor.

The publication gives no primer sequences, so the panel reuses the MLVAbank
primers already bundled as bacillus_mlva31_mlvabank, of which these 25 loci are
the first 25. What it does give is the panel composition, the repeat unit
lengths, the fragment size ranges of Table 2 and an explicit table of sequenced
amplicon size against repeat number, Table 3, which ships here as bins.tsv.

The repeat numbering is the paper's own and differs from the MLVAbank one at
five loci, so these tests pin both the agreement in size and the disagreement in
copy number, to keep the two numbering systems from being merged by accident.
"""

import pandas as pd

from mlvatyper.calling import read_bins
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from mlvatyper.sts import read_sts
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000008445.1"
SCHEME = "bacillus_mlva25_lista_2006"
MLVA31 = "bacillus_mlva31_mlvabank"
MLVA8 = "bacillus_mlva8_keim_2000"

KEIM = ["vrrA", "vrrB1", "vrrB2", "vrrC1", "vrrC2", "CG3", "pXO1-aat", "pXO2-at"]
BAMS = ["bams01", "bams03", "bams05", "bams13", "bams15", "bams21", "bams22", "bams23", "bams24",
        "bams25", "bams28", "bams30", "bams31", "bams34", "bams44", "bams51", "bams53"]
ORDER = KEIM + BAMS

PROFILE = "10-16-6-53-17-2-7-10-16-26-5-70-24-10-16-11-11-13-14-57-64-9-8-9-8"
# the Methods state the Ames MLVA-8 code in both units, which is the first eight digits above
AMES_MLVA8_BP = [314, 229, 153, 580, 532, 158, 126, 141]
AMES_MLVA8_REPEATS = [10, 16, 6, 53, 17, 2, 7, 10]

# a pattern is stored only where the published repeat number equals the number of motif copies
STORED_MOTIF = {"CG3": "TAATA", "pXO1-aat": "AAT"}

# the loci this panel numbers differently from the MLVAbank identifiers: name -> (here, MLVAbank)
RENUMBERED = {"vrrA": (10, 4), "vrrB1": (16, 20), "vrrB2": (6, 13),
              "pXO2-at": (10, 9), "bams34": (9, 11)}

# Table 2, the fragment size range observed over the 160 isolates
SIZE_RANGE = {
    "CG3": (153, 158), "bams44": (339, 417), "bams03": (474, 609), "vrrB2": (135, 171),
    "bams05": (229, 385), "bams15": (409, 616), "bams01": (380, 485), "vrrC1": (517, 616),
    "bams13": (337, 868), "vrrB1": (193, 256), "bams28": (373, 493), "vrrC2": (532, 604),
    "bams53": (212, 236), "bams31": (331, 952), "vrrA": (290, 326), "bams25": (376, 391),
    "bams21": (541, 676), "bams34": (238, 581), "bams24": (469, 637), "bams51": (358, 538),
    "bams22": (555, 735), "bams23": (399, 651), "bams30": (268, 916), "pXO1-aat": (123, 135),
    "pXO2-at": (132, 145),
}

# Table 3 is not on a single grid at these loci; see bins.tsv for what each one is
OFF_GRID = {"bams05": {332}, "bams13": {391}, "bams23": {639}, "bams24": {582},
            "bams34": {238, 503, 581}}

CHROMOSOME, PXO1, PXO2 = "AE017334.2", "AE017336.2", "AE017335.3"


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    assert scheme.bins is not None, "the panel is numbered by a published size table"
    pairs, _ = read_sts(scheme.sts)
    assert [p.id for p in pairs] == ORDER
    assert len(pairs) == 25, "the panel name is MLVA-25"
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.can_call, d.primer_id
        assert d.ref_size - d.ref_repeats * d.unit_length >= 0, d.primer_id
        if d.primer_id in STORED_MOTIF:
            assert d.patterns == (STORED_MOTIF[d.primer_id],), d.primer_id
            assert d.unit_length == len(d.patterns[0]), d.primer_id
        else:
            assert d.patterns == (), d.primer_id


def test_the_methods_mlva8_code_is_reproduced():
    """The Methods print the Ames code in base pairs and in repeats; both must match."""
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    assert [defs[m].ref_size for m in KEIM] == AMES_MLVA8_BP
    assert [defs[m].ref_repeats for m in KEIM] == AMES_MLVA8_REPEATS
    assert PROFILE.split("-")[:8] == [str(n) for n in AMES_MLVA8_REPEATS]


def test_every_reference_size_is_inside_the_published_range():
    for d in read_markers(load_scheme(SCHEME).markers):
        low, high = SIZE_RANGE[d.primer_id]
        assert low <= d.ref_size <= high, d.primer_id


def test_primers_and_sizes_are_the_mlvabank_ones_but_the_numbering_is_not():
    """Same loci and same products as the 31 locus panel; five differ in repeat number."""
    here_sts = {p.id: p for p in read_sts(load_scheme(SCHEME).sts)[0]}
    there_sts = {p.id: p for p in read_sts(load_scheme(MLVA31).sts)[0]}
    assert set(here_sts) < set(there_sts), "this panel is a subset of the 31 locus one"
    for marker, mine in here_sts.items():
        theirs = there_sts[marker]
        assert (mine.forward, mine.reverse) == (theirs.forward, theirs.reverse), marker

    here = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    there = {d.primer_id: d for d in read_markers(load_scheme(MLVA31).markers)}
    for marker, mine in here.items():
        assert mine.ref_size == there[marker].ref_size, marker
        assert mine.unit_length == there[marker].unit_length, marker
    differ = {m for m in here if here[m].ref_repeats != there[m].ref_repeats}
    assert differ == set(RENUMBERED)
    for marker, (mine, theirs) in RENUMBERED.items():
        assert (here[marker].ref_repeats, there[marker].ref_repeats) == (mine, theirs), marker


def test_the_mlva8_loci_are_shared_with_the_keim_scheme_in_the_same_order():
    """The first eight digits line up with bacillus_mlva8_keim_2000, with a different numbering."""
    there = [d.primer_id for d in read_markers(load_scheme(MLVA8).markers)]
    assert there == KEIM
    keim = {d.primer_id: d for d in read_markers(load_scheme(MLVA8).markers)}
    here = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    for marker in KEIM:
        assert here[marker].ref_size == keim[marker].ref_size, marker
    assert {m for m in KEIM if here[m].ref_repeats != keim[m].ref_repeats} == \
           {m for m in RENUMBERED if m in KEIM}


def test_bins_come_from_table_3_and_agree_with_the_marker_definitions():
    scheme = load_scheme(SCHEME)
    bins = read_bins(scheme.bins)
    assert set(bins) == set(ORDER), "every locus has published allele definitions"
    assert sum(len(v) for v in bins.values()) == 132
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    for marker, rows in bins.items():
        assert all(b.low == b.high for b in rows), f"{marker}: sizes are exact, not ranges"
        sizes = [b.low for b in rows]
        assert sizes == sorted(set(sizes)), f"{marker}: sizes are unique and ascending"
        # the reference product is one of the published rows, and the two agree on its allele
        d = defs[marker]
        at_reference = [b for b in rows if b.low == d.ref_size]
        assert len(at_reference) == 1, marker
        assert at_reference[0].allele == str(d.ref_repeats), marker


def test_bins_lie_on_the_repeat_grid_except_where_table_3_does_not():
    """size - unit_length * repeats is the flanking length and is constant within a locus."""
    scheme = load_scheme(SCHEME)
    bins, defs = read_bins(scheme.bins), {d.primer_id: d for d in read_markers(scheme.markers)}
    for marker, rows in bins.items():
        d = defs[marker]
        flanks = {b.low: b.low - d.unit_length * float(b.allele) for b in rows}
        on_grid = {size for size, flank in flanks.items() if size not in OFF_GRID.get(marker, set())}
        assert len({flanks[size] for size in on_grid}) == 1, marker
        for size in OFF_GRID.get(marker, set()):
            assert size in flanks, f"{marker}: {size} is published and kept"


def test_check_scheme_against_ames_ancestor(tmp_path):
    from mlvatyper.cli import main
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 25
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all(), "every reference size must be reproduced"
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    assert (table.loc[sorted(STORED_MOTIF), "repeats_ok"].astype(int) == 1).all()


def test_ames_ancestor_profile_is_called_from_the_published_table(tmp_path):
    from mlvatyper.cli import main
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "ames", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "ames_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    long = pd.read_csv(tmp_path / "ames_calls_long.tsv", sep="\t").set_index("marker")
    assert len(long) == 25
    # every locus is called from bins.tsv, and the formula agrees with it exactly
    assert (long["status"] == "BINNED").all()
    assert (long["residual_bp"] == 0).all()
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    assert {m: int(long.loc[m, "allele"]) for m in ORDER} == {m: defs[m].ref_repeats for m in ORDER}
    for marker in ORDER:
        expected = PXO1 if marker == "pXO1-aat" else PXO2 if marker == "pXO2-at" else CHROMOSOME
        assert long.loc[marker, "contig"] == expected, marker
