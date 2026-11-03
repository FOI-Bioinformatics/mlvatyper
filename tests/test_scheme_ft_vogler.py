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
from mlvatyper.sts import read_sts
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
        assert d.patterns == (), "the publication gives no repeat unit sequences"
        assert d.label == marker


# the unit each locus steps by, where a cross strain test and the Johansson panel agree
DERIVED_UNIT = {"Ft-M02": 6, "Ft-M03": 9, "Ft-M04": 5, "Ft-M05": 16,
                "Ft-M06": 21, "Ft-M10": 16, "Ft-M20A": 12}
# the Johansson spelling of the same locus, which pads the number differently
JOHANSSON = {"Ft-M02": "Ft-M2", "Ft-M03": "Ft-M3", "Ft-M04": "Ft-M4", "Ft-M05": "Ft-M5",
             "Ft-M06": "Ft-M6", "Ft-M10": "Ft-M10", "Ft-M20A": "Ft-M20"}


# the LVS copy number calibrated through SCHU S4, where the Johansson panel counts the same array
CALIBRATED = {"Ft-M02": 3, "Ft-M03": 12, "Ft-M04": 4, "Ft-M05": 2,
              "Ft-M06": 3, "Ft-M10": 2, "Ft-M20A": 4}
SIZED_ONLY = ["Ft-M20B", "Ft-M22", "Ft-M23", "Ft-M24"]


def test_seven_loci_call_and_four_are_sized_only():
    defs = {d.primer_id: d for d in read_markers(load_scheme("francisella_vogler_2008").markers)}
    assert {m for m, d in defs.items() if d.can_call} == set(CALIBRATED)
    assert {m for m, d in defs.items() if not d.can_call} == set(SIZED_ONLY)
    for marker, count in CALIBRATED.items():
        assert defs[marker].ref_repeats == count, marker
    for marker in SIZED_ONLY:
        assert defs[marker].ref_repeats is None, marker


def test_the_calibrated_flanks_are_whole_repeats_and_fit_the_primers():
    """A flank is the product minus the array, so it must be positive and a whole unit apart."""
    scheme = load_scheme("francisella_vogler_2008")
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    pairs = {p.id: p for p in read_sts(scheme.sts)[0]}
    for marker in CALIBRATED:
        d = defs[marker]
        flank = d.ref_size - d.unit_length * d.ref_repeats
        assert flank > 0, marker
        assert flank >= len(pairs[marker].forward) + len(pairs[marker].reverse), marker


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
        # but both now count the same array, so the reference alleles sit on one scale
        assert there[theirs].can_call and here[mine].can_call, mine


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


def test_lvs_run_sizes_every_amplicon_and_calls_the_seven(tmp_path):
    """Seven loci return their calibrated count; the other four report NO_DEFINITION."""
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", "francisella_vogler_2008", "--fasta", str(fixture), "--outdir",
                 str(tmp_path), "--prefix", "lvs", "--align", "none", "--no-report"]) == 0
    # every marker amplified, so nothing is reported as missing
    missing = pd.read_csv(tmp_path / "missing_markers.tsv", sep="\t", keep_default_na=False)
    assert len(missing) == 0
    hits = pd.read_csv(tmp_path / "hits.tsv", sep="\t", keep_default_na=False)
    assert hits.set_index("marker")["length"].astype(int).to_dict() == REF_SIZES

    long = pd.read_csv(tmp_path / "lvs_calls_long.tsv", sep="\t").set_index("marker")
    assert len(long) == len(MARKERS)
    # the reference calls back as the count recorded for it, with nothing left over
    assert {m: int(long.loc[m, "allele"]) for m in CALIBRATED} == CALIBRATED
    assert (long.loc[list(CALIBRATED), "residual_bp"] == 0).all()
    status = pd.read_csv(tmp_path / "lvs_status.tsv", sep="\t", index_col=0, dtype=str)
    assert {status.loc[ASSEMBLY, m] for m in CALIBRATED} == {"OK"}
    assert {status.loc[ASSEMBLY, m] for m in SIZED_ONLY} == {"NO_DEFINITION"}
