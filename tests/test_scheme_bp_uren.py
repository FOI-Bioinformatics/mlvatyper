"""The bundled Burkholderia pseudomallei U'Ren scheme against the K96243 assembly, downloaded by accession.

ref_size and ref_repeats of this scheme were measured in K96243 rather than
taken from the publication, whose copy numbers could not be reproduced. These
tests hold that measurement in place: every marker must give one product of the
recorded size with the recorded repeat count, so a change in the search or the
tandem counter shows up here.
"""

import pandas as pd

from mlvatyper.cli import main
from mlvatyper.markers import read_markers
from mlvatyper.scheme import load_scheme
from tests.conftest import reference_assembly

ASSEMBLY = "GCA_000011545.1"
PROFILE = ("9-8-8-7-6-18-15-7-6-12-17-8-30-8-7-9-19-23-13-7-24-6-8-15-12-12-7-15-13-16-15-7")
N_MARKERS = 32
# the two loci the metadata records as limitations
COMPOUND = "2050k"  # two adjacent arrays with different units, reported as one call
# markers that use the published copy number: the measured count gave alleles of 0 or below in other strains
PUBLISHED_COUNT = {"2115k": 15, "1500k": 6, "3564k": 30, "1788k": 8, "1217k": 15}
HIGHER_ORDER = "3564k"  # TRF reports a 75 bp period, five times the published 15 bp unit


def test_scheme_files_are_consistent():
    scheme = load_scheme("burkholderia_uren_2007")
    assert scheme.missing_provenance() == [] and scheme.missing_recommended() == []
    assert scheme.assembly == ASSEMBLY
    defs = read_markers(scheme.markers)
    assert len(defs) == N_MARKERS
    for d in defs:
        assert d.can_call, f"{d.primer_id} must carry a reference size and repeat count"
        assert d.ref_size - d.ref_repeats * d.unit_length >= 0
        assert all(len(p) == d.unit_length for p in d.patterns)
        assert d.label == d.primer_id.replace("k", " k")  # the published name keeps the space
    by_id = {d.primer_id: d for d in defs}
    assert len(by_id[COMPOUND].patterns) == 2  # both consensus units are listed
    assert by_id[HIGHER_ORDER].ref_repeats == 30  # the published copy number, see PUBLISHED_COUNT


def test_check_scheme_against_k96243(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["check-scheme", "--scheme", "burkholderia_uren_2007", "--outdir", str(tmp_path),
                 "--reference", str(fixture)]) == 0
    table = pd.read_csv(tmp_path / "scheme_check.tsv", sep="\t", keep_default_na=False)
    assert len(table) == N_MARKERS
    assert (table["n_products"].astype(int) == 1).all()
    assert (table["size_ok"].astype(int) == 1).all()
    # the count is compared only for markers whose arrays the patterns can follow
    skip = set(load_scheme("burkholderia_uren_2007").interrupted_arrays)
    table = table.set_index("marker")
    assert (table.loc[sorted(skip), "repeats_ok"] == "").all()
    assert (table.drop(index=sorted(skip))["repeats_ok"].astype(int) == 1).all()
    table = table.reset_index()
    assert (table["hash_word_11"].astype(int) == 1).all()  # no marker needs the rescue pass
    assert (table["search_pass"] == "main").all()


def test_k96243_profile_is_reproduced(tmp_path):
    fixture = reference_assembly(ASSEMBLY)
    assert main(["run", "--scheme", "burkholderia_uren_2007", "--fasta", str(fixture), "--outdir",
                 str(tmp_path), "--prefix", "k96243", "--align", "none", "--no-report"]) == 0
    profile = pd.read_csv(tmp_path / "k96243_profile.tsv", sep="\t", index_col=0, dtype=str)
    assert profile.loc[ASSEMBLY, "profile"] == PROFILE
    status = pd.read_csv(tmp_path / "k96243_status.tsv", sep="\t", index_col=0, dtype=str)
    assert set(status.loc[ASSEMBLY]) == {"OK"}
    long = pd.read_csv(tmp_path / "k96243_calls_long.tsv", sep="\t")
    assert len(long) == N_MARKERS
    assert (long["residual_bp"] == 0).all()
    # the length based allele and the independent tandem count agree wherever the count is compared
    skip = set(load_scheme("burkholderia_uren_2007").interrupted_arrays)
    kept = long[~long["marker"].isin(skip)]
    assert len(kept) == N_MARKERS - len(skip)
    assert (kept["tandem_copies"] == kept["allele"].astype(int)).all()


def test_published_copy_numbers_are_recorded_and_used_at_five_markers():
    """The publication's copy numbers are kept as comments; five of them replace the measured count."""
    scheme = load_scheme("burkholderia_uren_2007")
    text = scheme.markers.read_text()
    assert "differ for 20 of 32" in text
    comments = [line for line in text.splitlines() if line.startswith("#")]
    assert len(comments) > 8 and any("3564k 30" in line for line in comments)
    assert len(read_markers(scheme.markers)) == N_MARKERS


def test_the_scheme_turns_on_the_indel_rescue():
    scheme = load_scheme("burkholderia_uren_2007")
    assert scheme.calling_defaults()["rescue_indels"] is True


def test_the_indel_rescue_recovers_loci_in_a_second_strain(tmp_path):
    """In B. pseudomallei 1106a two loci are missed by the main search and found by the aligned rescue.

    The default search finds 28 of the 32 loci in this strain; the rescue gives 1764k and 2445k, the
    remaining four have primer sites that differ too much or a product outside the size window.
    """
    fixture = reference_assembly("GCF_000756145.1")
    assert main(["run", "--scheme", "burkholderia_uren_2007", "--fasta", str(fixture), "--outdir",
                 str(tmp_path), "--prefix", "s", "--align", "none", "--no-report"]) == 0
    hits = pd.read_csv(tmp_path / "hits.tsv", sep="\t").set_index("marker")
    assert set(hits.index[hits["search_pass"] != "main"]) == {"1764k", "2445k"}
    assert len(hits) == 28


def test_five_markers_use_the_published_count_and_sixteen_skip_the_count_check():
    scheme = load_scheme("burkholderia_uren_2007")
    defs = {d.primer_id: d for d in read_markers(scheme.markers)}
    assert {m: int(defs[m].ref_repeats) for m in PUBLISHED_COUNT} == PUBLISHED_COUNT
    skip = scheme.interrupted_arrays
    assert len(skip) == 16 and set(PUBLISHED_COUNT) <= set(skip) and set(skip) <= set(defs)
    marked = {d.primer_id for d in scheme.marker_defs() if not d.compare_count}
    assert marked == set(skip)


def test_a_second_strain_gives_positive_alleles_and_no_flag_at_interrupted_arrays(tmp_path):
    fixture = reference_assembly("GCF_000756145.1")
    assert main(["run", "--scheme", "burkholderia_uren_2007", "--fasta", str(fixture), "--outdir",
                 str(tmp_path), "--prefix", "s", "--align", "none", "--no-report"]) == 0
    long = pd.read_csv(tmp_path / "s_calls_long.tsv", sep="\t")
    alleles = pd.to_numeric(long["allele"], errors="coerce")
    assert (alleles.dropna() > 0).all()
    assert not long["status"].str.contains("NEGATIVE", na=False).any()
    skip = set(load_scheme("burkholderia_uren_2007").interrupted_arrays)
    flagged = long[long["marker"].isin(skip)]["status"].str.contains("COUNT_DISAGREES", na=False)
    assert not flagged.any()
