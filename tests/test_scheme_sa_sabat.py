"""The bundled Staphylococcus aureus five locus MLVA of Sabat et al. 2003 against NCTC 8325.

The paper is a gel fingerprint. It publishes primers and numbers the band
patterns, but gives no repeat counts, so no published scale exists. The scheme
therefore defines its own: the allele is the number of whole repeat units that
fit between the two primer sites of the reference product. The primers lie
inside the repeat arrays, so there is no flanking landmark to anchor a count on,
and this convention is the one rule that needs no outside value, is identical
for every strain and gives a negative result only for a product shorter than its
own primers. These tests pin that convention and the properties it relies on.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from mlvatyper.sts import read_sts
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000013425.1"
SCHEME = "saureus_mlva5_sabat_2003"
SCHOULS = "saureus_mlva8_schouls_2009"
CHROMOSOME = "CP000253.1"

ORDER = ["clfA", "clfB", "sdr", "spa", "ssp"]
MULTIPLE = "sdr"
SINGLE = [m for m in ORDER if m != MULTIPLE]

# the unit each locus steps by, checked against eight S. aureus genomes
UNIT = {"clfA": 6, "clfB": 6, "sdr": 6, "spa": 24, "ssp": 9}
# the Introduction calls the Ser-Asp repeat 18 bp, which is three dipeptides
PUBLISHED_UNIT = {"clfA": 18, "clfB": 18, "sdr": 18, "spa": 24, "ssp": 9}
SER_ASP = ("clfA", "clfB", "sdr")

# the reference products in NCTC 8325, and the repeat count the convention gives them
REF_SIZE = {"clfA": 997, "clfB": 828, "sdr": 796, "spa": 308, "ssp": 155}
REF_REPEATS = {"clfA": 159, "clfB": 130, "sdr": 125, "spa": 11, "ssp": 12}
SDR_PRODUCTS = {640, 796}  # NCTC 8325 gives two sdr products; the paper expects three

# every product size seen at each locus across N315, Mu50, MW2, COL, MRSA252, MSSA476, NCTC 8325, USA300
OBSERVED = {
    "clfA": [997, 1003, 1015, 1021, 1057, 1183, 1303],
    "clfB": [816, 828, 894, 912, 918, 936],
    "sdr": [556, 568, 580, 604, 634, 640, 652, 670, 676, 682, 688, 736, 748, 796],
    "spa": [236, 284, 308],
    "ssp": [128, 146, 155, 173, 218],
}

# the sspA primers carry IUPAC codes, so that marker is not searchable at word size 11
AMBIGUOUS = "ssp"


def _primer_length(pair):
    return len(pair.forward) + len(pair.reverse)


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    assert scheme.bins is None
    pairs, _ = read_sts(scheme.sts)
    assert [p.id for p in pairs] == ORDER
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.unit_length == UNIT[d.primer_id], d.primer_id
        assert d.ref_size == REF_SIZE[d.primer_id], d.primer_id
        assert d.ref_repeats == REF_REPEATS[d.primer_id], d.primer_id
        assert d.patterns == (), "the paper publishes no repeat unit sequences"
        assert d.label == d.primer_id


def test_all_five_loci_are_callable():
    defs = read_markers(load_scheme(SCHEME).markers)
    assert all(d.can_call for d in defs)


def test_ref_repeats_are_the_whole_units_between_the_primer_sites():
    """The convention: floor((reference product - both primers) / unit)."""
    scheme = load_scheme(SCHEME)
    pairs = {p.id: p for p in read_sts(scheme.sts)[0]}
    for d in read_markers(scheme.markers):
        interior = d.ref_size - _primer_length(pairs[d.primer_id])
        assert interior > 0, d.primer_id
        assert d.ref_repeats == interior // d.unit_length, d.primer_id


def test_one_remainder_per_locus_so_every_strain_gets_an_exact_integer():
    """Every product of a locus leaves the same remainder, so the formula never gives a fraction.

    This is what makes the convention exact across strains, and it includes sdr, whose three genes
    give 14 distinct sizes that still share one remainder.
    """
    scheme = load_scheme(SCHEME)
    pairs = {p.id: p for p in read_sts(scheme.sts)[0]}
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    for marker, sizes in OBSERVED.items():
        d = defs[marker]
        assert {(s - d.ref_size) % d.unit_length for s in sizes} == {0}, marker
        interior = [s - _primer_length(pairs[marker]) for s in sizes]
        assert len({i % d.unit_length for i in interior}) == 1, marker
        alleles = [d.ref_repeats + (s - d.ref_size) // d.unit_length for s in sizes]
        assert alleles == [i // d.unit_length for i in interior], marker
        assert min(alleles) > 0, f"{marker}: no observed product is shorter than its primers"


def test_the_convention_values_at_the_extremes():
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    def allele(marker, size):
        d = defs[marker]
        return d.ref_repeats + (size - d.ref_size) // d.unit_length
    assert allele("clfA", 1303) == 210 and allele("clfA", 997) == 159
    assert allele("spa", 236) == 8 and allele("spa", 308) == 11
    assert allele("ssp", 128) == 9 and allele("ssp", 218) == 19
    assert allele("sdr", 556) == 85 and allele("sdr", 796) == 125


def test_sdr_is_declared_as_a_locus_with_several_products():
    """One primer pair amplifies sdrC, sdrD and sdrE together, so it gives several products."""
    scheme = load_scheme(SCHEME)
    assert scheme.known_multiple == (MULTIPLE,)
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    # the reference size is the larger of the two NCTC 8325 products; the same rule fits both
    assert defs[MULTIPLE].ref_size == max(SDR_PRODUCTS)
    d = defs[MULTIPLE]
    assert d.ref_repeats + (min(SDR_PRODUCTS) - d.ref_size) // d.unit_length == 99


def test_the_ssp_primers_keep_their_ambiguity_codes():
    pairs = {p.id: p for p in read_sts(load_scheme(SCHEME).sts)[0]}
    for pair in pairs.values():
        assert set(pair.forward) <= set("ACGTRYSWKMBDHVN"), pair.id
        assert set(pair.reverse) <= set("ACGTRYSWKMBDHVN"), pair.id
    assert set(pairs[AMBIGUOUS].forward) - set("ACGT"), "the published sspA primer is degenerate"
    for marker in ORDER:
        if marker != AMBIGUOUS:
            assert not set(pairs[marker].forward) - set("ACGT"), marker
            assert not set(pairs[marker].reverse) - set("ACGT"), marker


def test_check_scheme_against_nctc_8325(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 5
    # every locus, sdr included, reproduces its reference size
    assert (table["size_ok"].astype(int) == 1).all()
    for marker in SINGLE:
        assert int(table.loc[marker, "n_products"]) == 1, marker
    assert int(table.loc[MULTIPLE, "n_products"]) > 1
    # no stored pattern, so there is no tandem cross-check
    assert (table["repeats_ok"] == "").all()
    # the degenerate primer pair is the one that needs the rescue pass
    assert int(table.loc[AMBIGUOUS, "hash_word_11"]) == 0
    assert table.loc[AMBIGUOUS, "search_pass"] == "rescue"
    assert all(table.loc[m, "search_pass"] == "main" for m in ORDER if m != AMBIGUOUS)


def test_a_run_on_nctc_8325_calls_the_convention_alleles(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "n8325", "--align", "none", "--no-report"]) == 0
    long = pd.read_csv(tmp_path / "n8325_calls_long.tsv", sep="\t")
    assert set(long["contig"]) == {CHROMOSOME}
    assert (long["residual_bp"].fillna(0) == 0).all()
    assert (long["allele"].astype(float) > 0).all()
    by_marker = long.groupby("marker")["allele"].apply(lambda s: sorted(int(float(a)) for a in s))
    for marker in SINGLE:
        assert by_marker[marker] == [REF_REPEATS[marker]], marker
    assert by_marker[MULTIPLE] == [99, 125], "both sdr products call, each by the same rule"
    status = long.groupby("marker")["status"].apply(set)
    assert status[MULTIPLE] == {"MULTIPLE"}
    assert all(status[m] == {"OK"} for m in SINGLE)
    alleles = pd.read_csv(tmp_path / "n8325_alleles.tsv", sep="\t", index_col=0, dtype=str)
    assert alleles.loc[ASSEMBLY, SINGLE].tolist() == [str(REF_REPEATS[m]) for m in SINGLE]


def test_the_two_saureus_schemes_overlap_at_spa_with_different_primers():
    here = {p.id: p for p in read_sts(load_scheme(SCHEME).sts)[0]}
    there = {p.id: p for p in read_sts(load_scheme(SCHOULS).sts)[0]}
    assert "VNTR24_01" in there, "the Schouls panel reads spa as VNTR24_01"
    assert (here["spa"].forward, here["spa"].reverse) != \
           (there["VNTR24_01"].forward, there["VNTR24_01"].reverse)
    # no locus id is shared, so the two schemes cannot be confused by name either
    assert not set(here) & set(there)
    # they read the same 24 bp array but on different scales and from different products, so a count
    # from one is not a count from the other: the reference products differ in size, and the same
    # strain, N315, whose spa product here is 284 bp, gets 10 under this convention and 11 in Schouls
    d_here = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}["spa"]
    d_there = {d.primer_id: d for d in read_markers(load_scheme(SCHOULS).markers)}["VNTR24_01"]
    assert d_here.unit_length == d_there.unit_length == 24
    assert d_here.ref_size != d_there.ref_size
    n315_here = d_here.ref_repeats + (284 - d_here.ref_size) // d_here.unit_length
    n315_there = d_there.ref_repeats  # N315 is the Schouls reference
    assert (n315_here, n315_there) == (10, 11)


def test_the_ser_asp_loci_step_by_one_dipeptide():
    """The array gains and loses Ser-Asp dipeptides singly, so the unit is 6 bp, not 18."""
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    differ = {m: (d.unit_length, PUBLISHED_UNIT[m])
              for m, d in defs.items() if d.unit_length != PUBLISHED_UNIT[m]}
    assert set(differ) == set(SER_ASP)
    for marker in SER_ASP:
        assert differ[marker] == (6, 18), marker
        assert PUBLISHED_UNIT[marker] % defs[marker].unit_length == 0, marker
    # the two loci outside that family keep the published unit
    for marker in ("spa", "ssp"):
        assert defs[marker].unit_length == PUBLISHED_UNIT[marker], marker
