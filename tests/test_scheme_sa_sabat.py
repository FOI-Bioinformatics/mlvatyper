"""The bundled Staphylococcus aureus five locus MLVA of Sabat et al. 2003 against NCTC 8325.

This is the thinnest of the bundled schemes and the tests say so. The paper
publishes primers and nothing else quantitative: no amplicon sizes, no repeat
numbers, no allele definitions. Its MLVA pattern numbers are arbitrary labels
for band patterns on a gel, not per locus alleles. So the scheme sizes its
markers without calling them, and it is the first bundled scheme with a locus
that gives more than one product on purpose.
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

# the repeat unit the Introduction states for each gene family
UNIT = {"clfA": 18, "clfB": 18, "sdr": 18, "spa": 24, "ssp": 9}
# the products measured in NCTC 8325; none of these is a published value
REF_SIZES = {"clfA": 997, "clfB": 828, "spa": 308, "ssp": 155}
SDR_PRODUCTS = {640, 796}

# the sspA primers carry IUPAC codes, so that marker is not searchable at word size 11
AMBIGUOUS = "ssp"


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
        assert d.patterns == (), "the paper publishes no repeat unit sequences"
        assert d.label == d.primer_id


def test_nothing_can_be_called():
    """No reference value of any kind is published, so no marker may be callable."""
    defs = read_markers(load_scheme(SCHEME).markers)
    assert all(not d.can_call for d in defs)
    assert all(d.ref_repeats is None for d in defs)


def test_sdr_is_declared_as_a_locus_with_several_products():
    """One primer pair amplifies sdrC, sdrD and sdrE together, so it has no single size."""
    scheme = load_scheme(SCHEME)
    assert scheme.known_multiple == (MULTIPLE,)
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    assert defs[MULTIPLE].ref_size is None, "a locus with several products has no reference size"
    for marker in SINGLE:
        assert defs[marker].ref_size == REF_SIZES[marker], marker


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
    # the four single product loci reproduce their measured size
    for marker in SINGLE:
        assert int(table.loc[marker, "n_products"]) == 1, marker
        assert int(table.loc[marker, "size_ok"]) == 1, marker
    assert int(table.loc[MULTIPLE, "n_products"]) > 1
    assert table.loc[MULTIPLE, "size_ok"] == "", "a declared multiple locus has no size to check"
    # nothing is callable anywhere
    assert (table["repeats_ok"] == "").all()
    # the degenerate primer pair is the one that needs the rescue pass
    assert int(table.loc[AMBIGUOUS, "hash_word_11"]) == 0
    assert table.loc[AMBIGUOUS, "search_pass"] == "rescue"
    assert all(table.loc[m, "search_pass"] == "main" for m in ORDER if m != AMBIGUOUS)


def test_a_run_on_nctc_8325_sizes_the_single_product_loci(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "n8325", "--align", "none", "--no-report"]) == 0
    hits = pd.read_csv(tmp_path / "hits.tsv", sep="\t", keep_default_na=False)
    assert set(hits["contig"]) == {CHROMOSOME}
    by_marker = hits.groupby("marker")["length"].apply(lambda s: set(s.astype(int)))
    for marker in SINGLE:
        assert by_marker[marker] == {REF_SIZES[marker]}, marker
    assert by_marker[MULTIPLE] == SDR_PRODUCTS, "both sdr genes amplify"
    assert not (tmp_path / "n8325_alleles.tsv").exists()
    assert not (tmp_path / "n8325_status.tsv").exists()


def test_the_two_saureus_schemes_overlap_at_spa_with_different_primers():
    here = {p.id: p for p in read_sts(load_scheme(SCHEME).sts)[0]}
    there = {p.id: p for p in read_sts(load_scheme(SCHOULS).sts)[0]}
    assert "VNTR24_01" in there, "the Schouls panel reads spa as VNTR24_01"
    assert (here["spa"].forward, here["spa"].reverse) != \
           (there["VNTR24_01"].forward, there["VNTR24_01"].reverse)
    # no locus id is shared, so the two schemes cannot be confused by name either
    assert not set(here) & set(there)
