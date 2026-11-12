"""The bundled HOOF-Prints scheme (Bricker et al. 2003) against the three genomes of its Figure 1.

Alleles are the number of complete octameric repeat units in the amplicon, called from length against
B. abortus 9-941. Figure 1 of the paper draws the units of all three sequenced genomes, which is the
only published allele table the assay has.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

SCHEME = "brucella_hoofprints_bricker_2003"
ABORTUS = "GCF_000008145.1"  # B. abortus bv. 1 str. 9-941, the reference
MELITENSIS = "GCF_000007125.1"  # 16M
SUIS = "GCF_000007505.1"  # 1330
ORDER = [f"Locus-{i}" for i in range(1, 9)]
# complete repeat units drawn in Figure 1 for each of the three sequenced genomes
FIGURE = {"Locus-1": (2, 5, 4), "Locus-2": (5, 5, 2), "Locus-3": (4, 1, 4), "Locus-4": (4, 7, 1),
          "Locus-5": (2, 7, 7), "Locus-6": (2, 1, 5), "Locus-7": (14, 8, 11), "Locus-8": (2, 6, 4)}
# product of 9-941 measured in the reference assembly
REF_SIZE = {"Locus-1": 140, "Locus-2": 125, "Locus-3": 151, "Locus-4": 155,
            "Locus-5": 138, "Locus-6": 158, "Locus-7": 183, "Locus-8": 91}


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ABORTUS
    assert set(scheme.interrupted_arrays) == {"Locus-3", "Locus-6", "Locus-8"}
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    assert list(defs) == ORDER
    for marker, d in defs.items():
        assert d.unit_length == 8 and d.ref_size == REF_SIZE[marker]
        assert d.ref_repeats == FIGURE[marker][0], marker
        assert d.offset == d.ref_size - 8 * d.ref_repeats and d.offset > 0, marker
        assert d.has_repeats  # one octamer with two variants, so the count is a cross-check


def test_check_scheme_against_9941(tmp_path):
    fixture = reference_assembly(ABORTUS)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path), "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert list(table.index) == ORDER
    assert (table["n_products"].astype(int) == 1).all() and (table["size_ok"].astype(int) == 1).all()
    assert (table["length"].astype(int) == pd.Series(REF_SIZE)).all()


def test_the_three_sequenced_genomes_give_the_figure_1_counts(tmp_path):
    """Seven of eight loci reproduce Figure 1 in all three genomes; the two exceptions are published."""
    fastas = [str(reference_assembly(a)) for a in (ABORTUS, MELITENSIS, SUIS)]
    assert main(["run", "--scheme", SCHEME, "--fasta", *fastas, "--outdir", str(tmp_path),
                 "--prefix", "h", "--align", "none", "--no-report"]) == 0
    alleles = pd.read_csv(tmp_path / "h_alleles.tsv", sep="\t", index_col=0, dtype=str)
    assert list(alleles.columns) == ORDER and len(alleles) == 3
    expected = {a: {m: str(FIGURE[m][i]) for m in ORDER} for i, a in enumerate((ABORTUS, MELITENSIS, SUIS))}
    assert alleles.loc[ABORTUS].to_dict() == expected[ABORTUS]
    assert alleles.loc[MELITENSIS].to_dict() == expected[MELITENSIS]
    # B. suis 1330 departs from the figure at two loci, both for reasons the paper states
    differs = {m for m in ORDER if alleles.loc[SUIS, m] != expected[SUIS][m]}
    assert differs == {"Locus-1", "Locus-4"}


def test_a_flank_deletion_is_caught_by_the_repeat_count(tmp_path):
    """B. suis 1330 has a 16 bp deletion in the Locus-1 flank, so its 4 unit array sizes as 2 units."""
    fixture = reference_assembly(SUIS)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "s", "--align", "none", "--no-report"]) == 0
    long = pd.read_csv(tmp_path / "s_calls_long.tsv", sep="\t").set_index("marker")
    assert long.loc["Locus-1", "length"] == 140  # the same product size as 2 units in 9-941
    assert long.loc["Locus-1", "allele"] == 2 and long.loc["Locus-1", "tandem_copies"] == 4
    assert "COUNT_DISAGREES" in long.loc["Locus-1", "status"]
    # the 11 bp offset the paper reports for B. suis biovar 1 leaves Locus-4 between two alleles
    assert long.loc["Locus-4", "length"] == 142 and "INTERMEDIATE" in long.loc["Locus-4", "status"]
    # the three loci whose arrays open with a partial unit are declared, so their counts are not compared
    for marker in ("Locus-3", "Locus-6", "Locus-8"):
        assert "COUNT_DISAGREES" not in str(long.loc[marker, "status"]), marker
