"""The bundled Brucella MLVA-16 scheme against the B. melitensis 16M assembly, downloaded by accession.

Alleles come from the published size to copy number table of Al Dahouk et al. 2007, shipped as bins.tsv,
with the repeat formula behind it for sizes the table does not list.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCF_000007125.1"
# scheme order: panel 1 (MLVA-8), panel 2A, panel 2B
ORDER = ["Bruce06", "Bruce08", "Bruce11", "Bruce12", "Bruce42", "Bruce43", "Bruce45", "Bruce55",
         "Bruce18", "Bruce19", "Bruce21", "Bruce04", "Bruce07", "Bruce09", "Bruce16", "Bruce30"]
# repeat unit (bp), B. melitensis 16M product (bp) and copy number, from Le Fleche et al. 2006 Table 1;
# Bruce19 is numbered in 3 bp units as in MLVAbank allele numbering version 3.6 (16M = 36, twice the 6 bp value)
SOURCE = {
    "Bruce06": (134, 408, 3), "Bruce08": (18, 348, 4), "Bruce11": (63, 257, 2), "Bruce12": (15, 392, 13),
    "Bruce42": (125, 539, 4), "Bruce43": (12, 182, 2), "Bruce45": (18, 151, 3), "Bruce55": (40, 273, 3),
    "Bruce18": (8, 146, 5), "Bruce19": (3, 163, 36), "Bruce21": (8, 148, 6),
    "Bruce04": (8, 152, 2), "Bruce07": (8, 158, 5), "Bruce09": (8, 156, 7), "Bruce16": (8, 152, 3),
    "Bruce30": (8, 151, 6),
}


def test_scheme_files_match_the_source_table():
    scheme = load_scheme("Bru_MLVA16")
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    assert list(defs) == ORDER
    for marker, (unit, size, repeats) in SOURCE.items():
        d = defs[marker]
        assert (d.unit_length, d.ref_size, d.ref_repeats) == (unit, size, repeats), marker
        assert d.can_call and not d.has_repeats  # called from length alone, no unit sequences shipped
        assert d.offset == size - repeats * unit and d.offset >= 0, marker
    sts = {line.split("\t")[0]: int(line.split("\t")[3])
           for line in scheme.sts.read_text().splitlines() if line and not line.startswith("#")}
    assert sts == {m: SOURCE[m][1] for m in ORDER}


def test_check_scheme_against_16m(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", "Bru_MLVA16", "--outdir", str(tmp_path), "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert list(table.index) == ORDER
    assert (table["n_products"].astype(int) == 1).all() and (table["size_ok"].astype(int) == 1).all()
    assert (table["length"].astype(int) == pd.Series({m: SOURCE[m][1] for m in ORDER})).all()


def test_16m_profile_equals_the_published_copy_numbers(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", "Bru_MLVA16", "--fasta", str(fixture), "--outdir", str(tmp_path), "--prefix", "b",
                 "--mlvabank", "--align", "none", "--no-report"]) == 0
    alleles = pd.read_csv(tmp_path / "b_alleles.tsv", sep="\t", index_col=0, dtype=str)
    assert alleles.loc[ASSEMBLY].to_dict() == {m: str(SOURCE[m][2]) for m in ORDER}
    status = pd.read_csv(tmp_path / "b_status.tsv", sep="\t", index_col=0, dtype=str)
    # every reference size is in the published allele table, so each locus is called from the bins
    assert set(status.loc[ASSEMBLY]) == {"BINNED"}
    long = pd.read_csv(tmp_path / "b_calls_long.tsv", sep="\t")
    assert (long.fwd_mismatches + long.rev_mismatches == 0).all() and (long.residual_bp == 0).all()
    assert (tmp_path / "b_mlvabank.csv").read_text().splitlines()[0] == "Strain," + ",".join(ORDER)


# rows of bins.tsv taken from the MLVAbank version 3.6 binning file that are not points on the ladder
V36_RANGES = {"Bruce12": {(327, 332), (342, 347), (392, 397), (407, 411)}, "Bruce21": {(172, 175)}}
V36_OFF_LADDER = {"Bruce16": {254, 270}, "Bruce55": {533}}


def test_the_allele_table_agrees_with_the_reference_values():
    """Apart from the listed version 3.6 rows, each ladder steps by the repeat unit and meets 16M."""
    from mlvatyper.calling import read_bins
    scheme = load_scheme("Bru_MLVA16")
    bins = read_bins(scheme.bins)
    assert set(bins) == set(ORDER), "the table covers all 16 loci"
    for marker, (unit, ref_size, ref_repeats) in SOURCE.items():
        ranges = {(b.low, b.high) for b in bins[marker] if b.low != b.high}
        assert ranges == V36_RANGES.get(marker, set()), marker
        for b in bins[marker]:
            if b.low in V36_OFF_LADDER.get(marker, set()):
                continue
            # a range must contain its ladder size; a point must be on the ladder
            on_ladder = [s for s in range(b.low, b.high + 1) if (s - ref_size) % unit == 0]
            assert len(on_ladder) == 1, f"{marker} {b.low}-{b.high}"
            assert int(b.allele) == ref_repeats + (on_ladder[0] - ref_size) // unit, f"{marker} {b.low}"
        assert any(b.low == ref_size and int(b.allele) == ref_repeats for b in bins[marker]), marker


def test_bruce19_3bp_numbering_gives_integer_alleles_off_the_6bp_ladder(tmp_path):
    """B. abortus 82 gives 184 bp at Bruce19, half a 6 bp unit off the B. melitensis ladder; it is 43."""
    fixture = reference_assembly("GCF_000473805.1")  # B. abortus 82
    assert main(["run", "--scheme", "Bru_MLVA16", "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "a", "--align", "none", "--no-report"]) == 0
    long = pd.read_csv(tmp_path / "a_calls_long.tsv", sep="\t").set_index("marker")
    assert long.loc["Bruce19", "length"] == 184
    assert long.loc["Bruce19", "allele"] == 43 and long.loc["Bruce19", "status"] == "BINNED"
    # 3 loci have no product in this draft assembly; of the 13 that do, 12 are on the table
    typed = long[long["length"].notna()]
    assert len(typed) == 13 and (typed["status"] == "BINNED").sum() == 12
    off = typed[typed["status"].str.contains("OFF_TABLE", na=False)]
    assert list(off.index) == ["Bruce12"] and off.loc["Bruce12", "allele"] == 12  # 375 bp, 2 bp under 377


def test_primer_site_substitution_near_the_3_prime_end_is_rescued(tmp_path):
    """B. canis carries a substitution 5 nt from the 3' end of the Bruce11 forward primer.

    The main search misses it; the scheme turns on rescue_indels, which finds the 698 bp product,
    the published B. canis allele 9 (Al Dahouk et al. 2007, Table 2). Bruce21 at 175 bp is 9 by
    the version 3.6 range 172-175 rather than 9.5.
    """
    assert load_scheme("Bru_MLVA16").meta.get("rescue_indels") is True
    fixture = reference_assembly("GCF_000018525.1")  # B. canis ATCC 23365
    assert main(["run", "--scheme", "Bru_MLVA16", "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "c", "--align", "none", "--no-report"]) == 0
    long = pd.read_csv(tmp_path / "c_calls_long.tsv", sep="\t").set_index("marker")
    assert long.loc["Bruce11", "length"] == 698 and long.loc["Bruce11", "allele"] == 9
    assert long.loc["Bruce11", "search_pass"] != "main" and long.loc["Bruce11", "fwd_mismatches"] == 1
    assert long.loc["Bruce21", "length"] == 175 and long.loc["Bruce21", "allele"] == 9
    assert long["allele"].notna().all()
