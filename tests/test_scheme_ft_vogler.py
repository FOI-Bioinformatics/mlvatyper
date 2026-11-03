"""The bundled Francisella tularensis Vogler scheme against the LVS assembly, downloaded by accession.

The scheme carries primers and reference amplicon sizes but no repeat units or
copy numbers, because the publication's table does not give them. These tests
fix what is known: the primers amplify one product of the recorded size in LVS,
and no allele is reported until the missing definitions are supplied.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000009245.1"
# Table S1 order, which is the PCR mix order: 1A, 1B, 2, 3A, 3B
MARKERS = ["Ft-M22", "Ft-M23", "Ft-M24", "Ft-M03", "Ft-M20A", "Ft-M10", "Ft-M05", "Ft-M04", "Ft-M20B",
           "Ft-M06", "Ft-M02"]
REF_SIZES = {"Ft-M22": 178, "Ft-M23": 327, "Ft-M24": 417, "Ft-M03": 353, "Ft-M20A": 319, "Ft-M10": 198,
             "Ft-M05": 300, "Ft-M04": 222, "Ft-M20B": 153, "Ft-M06": 272, "Ft-M02": 339}
SHARED_REVERSE = ("Ft-M05", "Ft-M10")  # the same oligo, separated by product size alone


def test_scheme_files_are_consistent():
    scheme = load_scheme("francisella_vogler_2008")
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    assert list(defs) == MARKERS
    for marker, size in REF_SIZES.items():
        d = defs[marker]
        assert d.ref_size == size
        # the publication gives no copy number, so no allele can be computed yet
        assert not d.can_call and d.patterns == () and d.ref_repeats is None
        assert d.label == marker


# the unit each locus steps by, where a cross strain test and the Johansson panel agree
DERIVED_UNIT = {"Ft-M02": 6, "Ft-M03": 9, "Ft-M04": 5, "Ft-M05": 16,
                "Ft-M06": 21, "Ft-M10": 16, "Ft-M20A": 12}
# the Johansson spelling of the same locus, which pads the number differently
JOHANSSON = {"Ft-M02": "Ft-M2", "Ft-M03": "Ft-M3", "Ft-M04": "Ft-M4", "Ft-M05": "Ft-M5",
             "Ft-M06": "Ft-M6", "Ft-M10": "Ft-M10", "Ft-M20A": "Ft-M20"}


def test_unit_lengths_are_filled_only_where_two_sources_agree():
    defs = {d.primer_id: d for d in read_markers(load_scheme("francisella_vogler_2008").markers)}
    filled = {m: d.unit_length for m, d in defs.items() if d.unit_length}
    assert filled == DERIVED_UNIT
    # the other four stay empty, so nothing is invented for them
    assert {m for m, d in defs.items() if not d.unit_length} == set(MARKERS) - set(DERIVED_UNIT)


def test_the_two_francisella_schemes_agree_on_the_shared_units():
    """The locus names differ only in padding, and the units were derived independently."""
    here = {d.primer_id: d for d in read_markers(load_scheme("francisella_vogler_2008").markers)}
    there = {d.primer_id: d for d in read_markers(load_scheme("francisella_johansson_2004").markers)}
    for mine, theirs in JOHANSSON.items():
        assert there[theirs].unit_length == here[mine].unit_length, mine
        # different primers and different references, so the products are not the same fragment
        assert there[theirs].ref_size != here[mine].ref_size, mine


def test_every_primer_carries_the_published_5_prime_g():
    """The 5' G is not genomic; ref_size is only valid for the primers as written."""
    scheme = load_scheme("francisella_vogler_2008")
    rows = [line.split("\t") for line in scheme.sts.read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    assert len(rows) == len(MARKERS)
    assert all(row[1].startswith("G") and row[2].startswith("G") for row in rows)
    shared = {row[0]: row[2] for row in rows if row[0] in SHARED_REVERSE}
    assert len(set(shared.values())) == 1  # Ft-M05-2R and Ft-M10-2R are one oligo


def test_check_scheme_against_lvs(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", "francisella_vogler_2008", "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert (table["n_products"].astype(int) == 1).all(), "each primer pair must give one product in LVS"
    assert (table["size_ok"].astype(int) == 1).all()
    assert table["length"].astype(int).to_dict() == REF_SIZES
    assert (table["repeats_ok"] == "").all()  # no units, so the tandem cross-check is skipped


def test_lvs_run_sizes_every_amplicon_but_calls_no_allele(tmp_path):
    """With no units or copy numbers, calling is skipped and only the lengths are reported."""
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", "francisella_vogler_2008", "--fasta", str(fixture), "--outdir",
                 str(tmp_path), "--prefix", "lvs", "--align", "none", "--no-report"]) == 0
    lengths = pd.read_csv(tmp_path / "lvs_lengths_summary.tsv", sep="\t", index_col=0)
    assert lengths.loc[ASSEMBLY].astype(int).to_dict() == REF_SIZES
    assert not (tmp_path / "lvs_alleles.tsv").exists()
    assert not (tmp_path / "lvs_status.tsv").exists()
    # every marker amplified, so nothing is reported as missing
    missing = pd.read_csv(tmp_path / "missing_markers.tsv", sep="\t", keep_default_na=False)
    assert len(missing) == 0
    hits = pd.read_csv(tmp_path / "hits.tsv", sep="\t", keep_default_na=False)
    assert hits.set_index("marker")["length"].astype(int).to_dict() == REF_SIZES
