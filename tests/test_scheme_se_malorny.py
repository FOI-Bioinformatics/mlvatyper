"""The bundled Salmonella Enteritidis nine locus MLVA of Malorny et al. 2008 against NCTC 13349.

The panel is unusual in two ways that these tests pin. Seven of its nine repeat
numbers are fractional, because the array ends in a partial unit, so the
alleles of those loci lie on a grid offset from the integers and a call on that
grid must come out as n.m rather than be flagged INTERMEDIATE. And four of the
sequence lengths Table 2 publishes are not the product their own primers give
in the assembly of the strain they are attributed to, so the measured length is
used instead.
"""

import pandas as pd
import pytest

from mlvatyper.calling import round_repeats
from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from mlvatyper.sts import read_sts
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000009505.1"
SCHEME = "salmonella_enteritidis_mlva9_malorny_2008"
CHROMOSOME = "AM933172.1"

# the Methods give this order for the MLVA type
ORDER = ["SENTR7", "SENTR5", "SENTR6", "SENTR2", "SENTR3", "SENTR4", "SE-3", "SENTR1", "SE-7"]
PROFILE = "3.1-10-5.1-9.3-4.4-4.9-3-7.8-8.5"

# Table 1: repeat size and the repeat number in NCTC 13349
UNIT = {"SENTR1": 45, "SENTR2": 39, "SENTR3": 93, "SENTR4": 7, "SENTR5": 6,
        "SENTR6": 7, "SENTR7": 9, "SE-3": 12, "SE-7": 61}
REPEATS = {"SENTR1": 7.8, "SENTR2": 9.3, "SENTR3": 4.4, "SENTR4": 4.9, "SENTR5": 10.0,
           "SENTR6": 5.1, "SENTR7": 3.1, "SE-3": 3.0, "SE-7": 8.5}
WHOLE = {"SENTR5", "SE-3"}

# locus -> (measured in the reference assembly, DNA sequence length of Table 2)
CORRECTED = {"SENTR2": (518, 519), "SENTR3": (493, 494), "SENTR4": (119, 120), "SE-7": (719, 717)}
PUBLISHED_SIZE = {"SENTR1": 440, "SENTR5": 271, "SENTR6": 180, "SENTR7": 135, "SE-3": 308,
                  **{m: pub for m, (_, pub) in CORRECTED.items()}}

# Table 1 gives the Boxrud et al. name for five of the SENTR loci
BOXRUD = {"SENTR1": "SE-10", "SENTR4": "SE-1", "SENTR5": "SE-5", "SENTR6": "SE-2",
          "SENTR7": "SE-9", "SE-3": "SE-3", "SE-7": "SE-7"}

DYES = ("FAM", "HEX", "TET")


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    pairs, _ = read_sts(scheme.sts)
    assert [p.id for p in pairs] == ORDER, "the order the Methods give for the MLVA type"
    assert len(pairs) == 9
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.can_call, d.primer_id
        assert d.patterns == (), "no repeat unit sequences are published"
        assert d.unit_length == UNIT[d.primer_id], d.primer_id
        assert d.ref_repeats == REPEATS[d.primer_id], d.primer_id
        assert d.ref_size - d.unit_length * d.ref_repeats >= 0, d.primer_id


def test_no_primer_carries_its_dye_label():
    """The 5' dye is an assay property, not part of the oligonucleotide."""
    for pair in read_sts(load_scheme(SCHEME).sts)[0]:
        for seq in (pair.forward, pair.reverse):
            assert set(seq) <= set("ACGTRYSWKMBDHVN"), pair.id
            assert not any(seq.startswith(dye) for dye in DYES), pair.id


def test_seven_loci_sit_on_a_fractional_grid():
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    fractional = {m for m, d in defs.items() if d.grid_offset != 0.0}
    assert fractional == set(ORDER) - WHOLE
    assert len(fractional) == 7
    for marker in fractional:
        d = defs[marker]
        assert round(d.grid_offset, 3) == round(d.ref_repeats - int(d.ref_repeats), 3), marker


def test_a_call_on_the_fractional_grid_is_exact_not_intermediate():
    """One repeat either side of the reference must come out as n.m with no INTERMEDIATE flag."""
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    for marker in ORDER:
        d = defs[marker]
        for step in (-1, 0, 1, 4):
            raw = d.ref_repeats + step
            allele, intermediate = round_repeats(raw, 0.25, "intermediate", d.grid_offset)
            assert not intermediate, f"{marker} at {raw}"
            assert float(allele) == pytest.approx(raw), f"{marker} at {raw}"


def test_four_published_sequence_lengths_are_overruled_by_the_measurement():
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    differ = {m: (d.ref_size, PUBLISHED_SIZE[m]) for m, d in defs.items()
              if d.ref_size != PUBLISHED_SIZE[m]}
    assert differ == CORRECTED
    # each difference must be far inside the allele tolerance, so no published size changes grid point
    for marker, (measured, published) in CORRECTED.items():
        shift = abs(published - measured) / defs[marker].unit_length
        assert shift < 0.25, f"{marker} moves alleles by {shift:.3f} repeats"


def test_check_scheme_against_nctc_13349(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 9
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    assert not any(table.loc[m, "note"] for m in ORDER), "every primer matches the reference exactly"


def test_reference_profile_is_the_published_repeat_numbers(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "nctc", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "nctc_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    status = pd.read_csv(tmp_path / "nctc_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}, "no locus is flagged, including the fractional ones"
    long = pd.read_csv(tmp_path / "nctc_calls_long.tsv", sep="\t").set_index("marker")
    assert len(long) == 9 and (long["residual_bp"] == 0).all()
    assert set(long["contig"]) == {CHROMOSOME}
    assert {m: float(long.loc[m, "allele"]) for m in ORDER} == REPEATS


def test_labels_keep_the_boxrud_names():
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    assert {m: d.label for m, d in defs.items() if d.label != m} == \
           {m: name for m, name in BOXRUD.items() if name != m}
    assert defs["SENTR2"].label == "SENTR2", "SENTR2 has no Boxrud name"
    assert defs["SENTR3"].label == "SENTR3"


def test_it_shares_no_locus_with_the_typhimurium_scheme():
    here = {p.id for p in read_sts(load_scheme(SCHEME).sts)[0]}
    there = {p.id for p in read_sts(load_scheme("STm_ECDC").sts)[0]}
    assert not here & there
