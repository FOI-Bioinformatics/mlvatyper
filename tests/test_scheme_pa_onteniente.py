"""The bundled Pseudomonas aeruginosa MLVA of Onteniente et al. 2003 against PAO1.

Table 1 gives the repeat unit and the number of units in PAO1 for all eight
loci. The primers reproduce the published PAO1 products, and a run returns the
published unit profile 11-12-5-8-7-12-14-45.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from mlvatyper.sts import read_sts
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000006765.1"
SCHEME = "paeruginosa_mlva8_onteniente_2003"
CHROMOSOME = "AE004091.2"

ORDER = ["ms010", "ms061", "ms077", "ms127", "ms142", "ms172", "ms173", "ms194"]
UNIT = {"ms010": 6, "ms061": 6, "ms077": 39, "ms127": 15, "ms142": 115,
        "ms172": 54, "ms173": 243, "ms194": 12}
PAO1 = {"ms010": 11, "ms061": 12, "ms077": 5, "ms127": 8, "ms142": 7,
        "ms172": 12, "ms173": 14, "ms194": 45}
REF_SIZE = {"ms010": 167, "ms061": 126, "ms077": 442, "ms127": 210, "ms142": 890,
            "ms172": 789, "ms173": 3503, "ms194": 690}
PROFILE = "11-12-5-8-7-12-14-45"


def test_scheme_files_are_consistent():
    scheme = load_scheme(SCHEME)
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    pairs, _ = read_sts(scheme.sts)
    assert [p.id for p in pairs] == ORDER
    defs = read_markers(scheme.markers)
    assert [d.primer_id for d in defs] == ORDER
    for d in defs:
        assert d.can_call, d.primer_id
        assert d.unit_length == UNIT[d.primer_id], d.primer_id
        assert d.ref_repeats == PAO1[d.primer_id], d.primer_id
        assert d.ref_size == REF_SIZE[d.primer_id], d.primer_id
        assert d.patterns == ()


def test_check_scheme_against_pao1(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 8
    assert (table["n_products"].astype(int) == 1).all(), "every locus gives one product, ms173 included"
    assert (table["size_ok"].astype(int) == 1).all()


def test_pao1_profile_is_the_published_unit_counts(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "pao1", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "pao1_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    assert PROFILE.split("-") == [str(PAO1[m]) for m in ORDER]
    status = pd.read_csv(tmp_path / "pao1_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "pao1_calls_long.tsv", sep="\t").set_index("marker")
    assert (long["residual_bp"] == 0).all()
    assert set(long["contig"]) == {CHROMOSOME}
    assert {m: int(long.loc[m, "allele"]) for m in ORDER} == PAO1
