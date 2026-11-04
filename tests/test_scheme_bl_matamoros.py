"""The bundled Bifidobacterium longum MLVA of Matamoros et al. 2011 against NCC2705.

Table 3 gives the repeat units and Table 1 the NCC2705 profile, 6-4-3-7-16. Four
of the five loci use the published copy number and a run on NCC2705 returns it.
VNTR23 does not: the published 3 implies a flank longer than the shortest product
seen in other strains, which is impossible, so the scheme uses the 7 units its
amplicon holds and that locus is offset by +4 from the publication.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from mlvatyper.sts import read_sts
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000007525.1"
SCHEME = "blongum_mlva5_matamoros_2011"
CHROMOSOME = "AE014295.3"

ORDER = ["VNTR12", "VNTR21", "VNTR23", "VNTR25", "VNTR26"]
UNIT = {"VNTR12": 14, "VNTR21": 9, "VNTR23": 10, "VNTR25": 12, "VNTR26": 10}
PUBLISHED = {"VNTR12": 6, "VNTR21": 4, "VNTR23": 3, "VNTR25": 7, "VNTR26": 16}
PUBLISHED_PROFILE = "6-4-3-7-16"
# VNTR23 is numbered on the units its amplicon holds, which is 4 more than the published count
OFFSET = {"VNTR23": 4}
NCC2705 = {m: v + OFFSET.get(m, 0) for m, v in PUBLISHED.items()}
PROFILE = "-".join(str(NCC2705[m]) for m in ("VNTR12", "VNTR21", "VNTR23", "VNTR25", "VNTR26"))
# Table 2 amplicon size ranges; the measured product must fall inside (allowing one bp of slack)
SIZE_RANGE = {"VNTR12": (436, 500), "VNTR21": (374, 392), "VNTR23": (313, 373),
              "VNTR25": (365, 414), "VNTR26": (252, 332)}


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
        assert d.ref_repeats == NCC2705[d.primer_id], d.primer_id
        assert d.patterns == (), "imperfect arrays, so no stored pattern"
        low, high = SIZE_RANGE[d.primer_id]
        assert low - 1 <= d.ref_size <= high + 1, d.primer_id


def test_check_scheme_against_ncc2705(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", SCHEME, "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False).set_index("marker")
    assert len(table) == 5
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    assert (table["search_pass"] == "main").all()


def test_ncc2705_profile_is_the_published_one(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", SCHEME, "--fasta", str(fixture), "--outdir", str(tmp_path),
                 "--prefix", "ncc", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "ncc_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    assert PROFILE.split("-") == [str(NCC2705[m]) for m in ORDER] == ["6", "4", "7", "7", "16"]
    assert PUBLISHED_PROFILE == "6-4-3-7-16", "the published profile is what the offset is measured from"
    status = pd.read_csv(tmp_path / "ncc_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "ncc_calls_long.tsv", sep="\t").set_index("marker")
    assert (long["residual_bp"] == 0).all()
    assert set(long["contig"]) == {CHROMOSOME}
    assert {m: int(long.loc[m, "allele"]) for m in ORDER} == NCC2705


def test_vntr23_is_numbered_on_the_units_its_amplicon_holds():
    """The published 3 implies a flank longer than products seen in other strains, so it cannot be right."""
    scheme = load_scheme(SCHEME)
    d = {x.primer_id: x for x in read_markers(scheme.markers)}["VNTR23"]
    assert d.ref_repeats == PUBLISHED["VNTR23"] + OFFSET["VNTR23"] == 7
    pair = {p.id: p for p in read_sts(scheme.sts)[0]}["VNTR23"]
    published_flank = d.ref_size - PUBLISHED["VNTR23"] * d.unit_length
    flank = d.ref_size - d.ref_repeats * d.unit_length
    assert published_flank > 314, "a 314 bp product occurs, so the published flank is impossible"
    assert len(pair.forward) + len(pair.reverse) <= flank <= 314
    assert "shifted by -4" in scheme.meta["notes"], "the offset from the publication is recorded"
