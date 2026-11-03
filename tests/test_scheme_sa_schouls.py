"""The bundled Staphylococcus aureus eight locus MLVA of Schouls et al. 2009 against N315.

This scheme sizes its markers without calling them. The paper derives an allele
as the product length minus the flanking length over the repeat size, but
publishes neither the flanking lengths nor the repeat unit sequences and does
not state the repeat number of the sequenced reference, and the arrays are too
imperfect to count back. These tests pin that the scheme stays sizing only, and
check the primers against the coordinates the paper gives for the same pairs.
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

# Table 2: the start and end of each product in NC_002745, which is this assembly's chromosome
COORDINATES = {"VNTR09_01": (1021501, 1021864), "VNTR61_01": (2440537, 2440899),
               "VNTR61_02": (659826, 660093), "VNTR67_01": (680255, 680600),
               "VNTR21_01": (526755, 526904), "VNTR24_01": (122873, 123289),
               "VNTR63_01": (1056187, 1056833), "VNTR81_01": (266137, 266632)}

# Table 1: average repeat size and the repeat range over the 1681 isolates
UNIT = {"VNTR09_01": 9, "VNTR61_01": 60, "VNTR61_02": 66, "VNTR67_01": 67,
        "VNTR21_01": 21, "VNTR24_01": 24, "VNTR63_01": 64, "VNTR81_01": 81}
REPEAT_RANGE = {"VNTR09_01": (6, 25), "VNTR61_01": (0, 6), "VNTR61_02": (0, 5),
                "VNTR67_01": (0, 10), "VNTR21_01": (0, 13), "VNTR24_01": (1, 18),
                "VNTR63_01": (0, 11), "VNTR81_01": (1, 9)}

DYE_PREFIX = ("F", "N", "V", "P")


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    assert scheme.bins is None
    pairs, _ = read_sts(scheme.sts)
    assert [p.id for p in pairs] == ORDER
    assert len(pairs) == 8
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.unit_length == UNIT[d.primer_id], d.primer_id
        assert d.patterns == (), "no repeat unit sequences are published"
        assert d.label == d.primer_id


def test_the_scheme_sizes_its_markers_without_calling_them():
    """No reference repeat number is published or derivable, so no marker may be callable."""
    defs = read_markers(load_scheme(SCHEME).markers)
    assert all(not d.can_call for d in defs)
    assert all(d.ref_repeats is None for d in defs)
    # the sizes and unit lengths are still there, so only the one missing number blocks calling
    assert all(d.ref_size is not None and d.unit_length > 0 for d in defs)


def test_reference_sizes_equal_the_span_of_the_published_coordinates():
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    for marker, (start, end) in COORDINATES.items():
        assert defs[marker].ref_size == end - start + 1, marker


def test_no_forward_primer_keeps_its_dye_character():
    """Table 2 prints a dye letter before each forward primer, separated by a hyphen."""
    for pair in read_sts(load_scheme(SCHEME).sts)[0]:
        for seq in (pair.forward, pair.reverse):
            assert set(seq) <= set("ACGTRYSWKMBDHVN"), pair.id
            assert "-" not in seq, pair.id
        assert not any(pair.forward.startswith(d + "-") for d in DYE_PREFIX), pair.id


def test_the_window_covers_every_published_allele():
    """N315 may sit anywhere in the published repeat range, so the window spans it both ways."""
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    windows = {}
    for line in load_scheme(SCHEME).sts.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        windows[parts[0]] = tuple(int(x) for x in parts[3].split("-"))
    assert set(windows) == set(ORDER)
    for marker in ORDER:
        low, high = windows[marker]
        reach = (REPEAT_RANGE[marker][1] - REPEAT_RANGE[marker][0]) * UNIT[marker]
        size = defs[marker].ref_size
        assert low <= size <= high, marker
        assert high - size >= reach, marker
        assert size - low >= min(reach, size - 60), marker


def test_check_scheme_against_n315(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 8
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()
    assert not any(table.loc[m, "note"] for m in ORDER), "every primer matches N315 exactly"
    # nothing is callable, so the repeat column stays empty for every marker
    assert (table["repeats_ok"] == "").all()


def test_a_run_on_n315_sizes_every_amplicon_but_calls_no_allele(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "n315", "--align", "none", "--no-report"]) == 0
    sizes = {m: d.ref_size for m, d in
             {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}.items()}
    lengths = pd.read_csv(tmp_path / "n315_lengths_summary.tsv", sep="\t", index_col=0)
    assert lengths.loc[ASSEMBLY].astype(int).to_dict() == sizes
    # nothing can be called, so no allele table is written at all
    assert not (tmp_path / "n315_alleles.tsv").exists()
    assert not (tmp_path / "n315_status.tsv").exists()
    assert not (tmp_path / "n315_profile.tsv").exists()
    missing = pd.read_csv(tmp_path / "missing_markers.tsv", sep="\t", keep_default_na=False)
    assert len(missing) == 0, "every marker amplified"
    hits = pd.read_csv(tmp_path / "hits.tsv", sep="\t", keep_default_na=False)
    assert hits.set_index("marker")["length"].astype(int).to_dict() == sizes
    assert set(hits["contig"]) == {CHROMOSOME}
