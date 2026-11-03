"""The bundled Staphylococcus aureus eight locus MLVA of Schouls et al. 2009 against N315.

All eight loci call on the published scale. Figure 7 of the paper gives the MLVA
profile of 14 sequenced strains, and the N315 profile anchors every locus. The
anchor is validated against the other sequenced strains that have an assembly:
seven of eight loci reproduce each published profile exactly, the exception
being VNTR61_01, whose flank the paper reports as heterogeneous.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from mlvatyper.sts import read_sts
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000009645.1"
SCHEME = "saureus_mlva8_schouls_2009"
CHROMOSOME = "BA000018.3"

ORDER = ["VNTR09_01", "VNTR61_01", "VNTR61_02", "VNTR67_01",
         "VNTR21_01", "VNTR24_01", "VNTR63_01", "VNTR81_01"]

# Table 2: start and end of each product in NC_002745, this assembly's chromosome
COORDINATES = {"VNTR09_01": (1021501, 1021864), "VNTR61_01": (2440537, 2440899),
               "VNTR61_02": (659826, 660093), "VNTR67_01": (680255, 680600),
               "VNTR21_01": (526755, 526904), "VNTR24_01": (122873, 123289),
               "VNTR63_01": (1056187, 1056833), "VNTR81_01": (266137, 266632)}

# the unit each locus steps by, checked across eight genomes (VNTR63_01 is 63, not the table's 64)
UNIT = {"VNTR09_01": 9, "VNTR61_01": 60, "VNTR61_02": 66, "VNTR67_01": 67,
        "VNTR21_01": 21, "VNTR24_01": 24, "VNTR63_01": 63, "VNTR81_01": 81}
PUBLISHED_UNIT = {**UNIT, "VNTR63_01": 64}

# Figure 7: the N315 MLVA profile, the reference count for each locus in ORDER
N315_PROFILE = {"VNTR09_01": 14, "VNTR61_01": 2, "VNTR61_02": 1, "VNTR67_01": 3,
                "VNTR21_01": 2, "VNTR24_01": 11, "VNTR63_01": 7, "VNTR81_01": 5}
PROFILE_STRING = "14-2-1-3-2-11-7-5"
# the locus the paper flags as having a heterogeneous flank; it alone may deviate off-reference
HETEROGENEOUS = "VNTR61_01"

DYE_PREFIX = ("F", "N", "V", "P")


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    assert scheme.bins is None
    pairs, _ = read_sts(scheme.sts)
    assert [p.id for p in pairs] == ORDER and len(pairs) == 8
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.can_call, d.primer_id
        assert d.unit_length == UNIT[d.primer_id], d.primer_id
        assert d.ref_repeats == N315_PROFILE[d.primer_id], d.primer_id
        assert d.patterns == (), "no repeat unit sequences are published"


def test_all_eight_loci_call_on_the_figure_7_scale():
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    assert all(d.can_call for d in defs.values())
    assert {m: d.ref_repeats for m, d in defs.items()} == N315_PROFILE
    assert PROFILE_STRING.split("-") == [str(N315_PROFILE[m]) for m in ORDER]


def test_the_corrected_unit_length_is_kept():
    """VNTR63_01 steps by 63 bp across eight genomes, not the 64 bp Table 1 averages to."""
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    differ = {m: (d.unit_length, PUBLISHED_UNIT[m])
              for m, d in defs.items() if d.unit_length != PUBLISHED_UNIT[m]}
    assert differ == {"VNTR63_01": (63, 64)}


def test_reference_sizes_equal_the_span_of_the_published_coordinates():
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    for marker, (start, end) in COORDINATES.items():
        assert defs[marker].ref_size == end - start + 1, marker


def test_no_forward_primer_keeps_its_dye_character():
    for pair in read_sts(load_scheme(SCHEME).sts)[0]:
        for seq in (pair.forward, pair.reverse):
            assert set(seq) <= set("ACGTRYSWKMBDHVN"), pair.id
            assert "-" not in seq, pair.id
        assert not any(pair.forward.startswith(d + "-") for d in DYE_PREFIX), pair.id


def test_check_scheme_against_n315(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 8
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    assert {m for m in ORDER if table.loc[m, "ref_repeats"]} == set(ORDER), "every locus has a count"


def test_n315_run_returns_the_published_figure_7_profile(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "n315", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "n315_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE_STRING
    long = pd.read_csv(tmp_path / "n315_calls_long.tsv", sep="\t").set_index("marker")
    assert (long["residual_bp"] == 0).all()
    assert set(long["contig"]) == {CHROMOSOME}
    assert {m: int(long.loc[m, "allele"]) for m in ORDER} == N315_PROFILE


def test_it_shares_two_arrays_with_the_pourcel_panel():
    """VNTR24_01 and VNTR63_01 are the Pourcel Sa0122 and Sa1132 arrays, numbered differently."""
    here = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    pourcel = {d.primer_id: d for d in read_markers(load_scheme("saureus_mlva14_pourcel_2009").markers)}
    assert here["VNTR24_01"].unit_length == pourcel["Sa0122"].unit_length == 24
    assert here["VNTR63_01"].unit_length == pourcel["Sa1132"].unit_length == 63
    # the two panels number the same arrays on their own scales, so the reference counts differ
    assert here["VNTR24_01"].ref_repeats != pourcel["Sa0122"].ref_repeats
