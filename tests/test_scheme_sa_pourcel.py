"""The bundled Staphylococcus aureus 14 locus MLVA of Pourcel et al. 2009 against Mu50.

Every reference value is published: Pourcel et al. Table 1 gives the repeat unit
and the Mu50 copy number, and Bogut et al. 2024 Table 2 gives the Mu50 amplicon
size. The primers reproduce those amplicons, the Mu50 profile matches the
published one exactly, and the scheme is cross-checked against the 31 clinical
profiles of Bogut et al. Table 3. This is the first fully callable S. aureus
scheme in the set, and it shares two arrays with the Schouls panel.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from mlvatyper.sts import read_sts
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000009665.1"
SCHEME = "saureus_mlva14_pourcel_2009"
CHROMOSOME = "BA000017.4"

PANEL1 = ["Sa0122", "Sa0266", "Sa0311", "Sa0704", "Sa1132", "Sa1194", "Sa1291",
          "Sa1729", "Sa1866", "Sa2039"]
PANEL2 = ["Sa0906", "Sa1213", "Sa1425", "Sa1756"]
ORDER = PANEL1 + PANEL2

# Pourcel Table 1: unit length and Mu50 copy number; Bogut Table 2: Mu50 amplicon
UNIT = {"Sa0122": 24, "Sa0266": 81, "Sa0311": 55, "Sa0704": 67, "Sa1132": 63, "Sa1194": 67,
        "Sa1291": 64, "Sa1729": 56, "Sa1866": 159, "Sa2039": 56, "Sa0906": 56, "Sa1213": 56,
        "Sa1425": 58, "Sa1756": 131}
REPEATS = {"Sa0122": 10, "Sa0266": 6, "Sa0311": 3, "Sa0704": 4, "Sa1132": 6, "Sa1194": 7,
           "Sa1291": 4, "Sa1729": 5, "Sa1866": 3, "Sa2039": 3, "Sa0906": 3, "Sa1213": 5,
           "Sa1425": 4, "Sa1756": 2}
AMPLICON = {"Sa0122": 392, "Sa0266": 630, "Sa0311": 272, "Sa0704": 380, "Sa1132": 532,
            "Sa1194": 524, "Sa1291": 369, "Sa1729": 499, "Sa1866": 607, "Sa2039": 282,
            "Sa0906": 864, "Sa1213": 868, "Sa1425": 630, "Sa1756": 365}
PROFILE = "10-6-3-4-6-7-4-5-3-3-3-5-4-2"


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    pairs, _ = read_sts(scheme.sts)
    assert [p.id for p in pairs] == ORDER
    assert len(PANEL1) == 10 and len(PANEL2) == 4
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.can_call, d.primer_id
        assert d.unit_length == UNIT[d.primer_id], d.primer_id
        assert d.ref_repeats == REPEATS[d.primer_id], d.primer_id
        assert d.ref_size == AMPLICON[d.primer_id], d.primer_id
        assert d.patterns == (), "no repeat unit sequences are published"
        assert d.ref_size - d.unit_length * d.ref_repeats >= 0, d.primer_id


def test_check_scheme_reproduces_every_mu50_amplicon(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 14
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all(), "every published Mu50 amplicon must reproduce"
    assert (table["hash_word_11"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()


def test_mu50_profile_is_the_published_one(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "mu50", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "mu50_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    assert PROFILE.split("-") == [str(REPEATS[m]) for m in ORDER]
    status = pd.read_csv(tmp_path / "mu50_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "mu50_calls_long.tsv", sep="\t").set_index("marker")
    assert (long["residual_bp"] == 0).all()
    assert set(long["contig"]) == {CHROMOSOME}
    assert {m: int(long.loc[m, "allele"]) for m in ORDER} == REPEATS


def test_the_calling_model_matches_the_published_clinical_profiles():
    """Bogut Table 3 pairs an amplicon with an allele; the formula must reproduce the allele.

    The source sized amplicons on agarose gels, so a cell is accepted when the computed allele is
    within half a repeat of the published one; most agree exactly.
    """
    defs = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    # a sample of Bogut Table 3 panel 1 cells: locus -> (amplicon, published allele)
    cells = {
        "Sa0122": [(296, 6), (248, 4), (392, 10), (416, 11), (176, 1), (512, 15), (320, 7)],
        "Sa0311": [(162, 1), (272, 3), (327, 4), (217, 2), (437, 6), (382, 5)],
        "Sa0704": [(447, 5), (380, 4), (246, 2), (514, 6), (313, 3)],
        "Sa1132": [(217, 1), (532, 6), (343, 3), (784, 10), (658, 8)],
        "Sa1291": [(305, 3), (369, 4), (177, 1), (241, 2)],
        "Sa1729": [(499, 5), (387, 3), (331, 2), (275, 1), (555, 6)],
    }
    exact = total = 0
    for marker, rows in cells.items():
        d = defs[marker]
        for amplicon, published in rows:
            predicted = d.ref_repeats + (amplicon - d.ref_size) / d.unit_length
            total += 1
            exact += abs(predicted - published) < 1e-9
            assert abs(predicted - published) <= 0.5, f"{marker} {amplicon}: {predicted} vs {published}"
    assert exact >= total - 2, "the great majority of clinical cells match exactly"


def test_it_shares_two_arrays_with_the_schouls_panel():
    """Sa0122 and Sa1132 are the arrays the Schouls panel reads as VNTR24_01 and VNTR63_01."""
    here = {d.primer_id: d for d in read_markers(load_scheme(SCHEME).markers)}
    schouls = {d.primer_id: d for d in read_markers(load_scheme("saureus_mlva8_schouls_2009").markers)}
    # same repeat unit, the same physical arrays
    assert here["Sa0122"].unit_length == schouls["VNTR24_01"].unit_length == 24
    assert here["Sa1132"].unit_length == schouls["VNTR63_01"].unit_length == 63
    # each panel numbers the shared array on its own published scale, so the counts differ by a constant
    assert here["Sa0122"].ref_repeats == 10 and schouls["VNTR24_01"].ref_repeats == 11
    assert here["Sa1132"].ref_repeats == 6 and schouls["VNTR63_01"].ref_repeats == 7
    # different primers, so the Mu50 products are not the same fragment
    hp = {p.id: p for p in read_sts(load_scheme(SCHEME).sts)[0]}
    sp = {p.id: p for p in read_sts(load_scheme("saureus_mlva8_schouls_2009").sts)[0]}
    assert (hp["Sa0122"].forward, hp["Sa0122"].reverse) != (sp["VNTR24_01"].forward, sp["VNTR24_01"].reverse)
