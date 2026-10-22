"""End-to-end check of status codes on a synthetic assembly."""

import pandas as pd
from Bio.Seq import Seq

from mlvatyper.cli import main
from tests.conftest import FWD, REV, UNIT, make_amplicon, random_dna

FWD2 = "GGCATTACGTTAGCCATGCA"
REV2 = "CTAGGCATTCGGATCACGTA"
UNIT2 = "ATCGGAT"
FLANK2 = "GGCCTAA"  # 7 bp between the forward primer and the array, equal to one unit length


def amplicon2(repeats: int, flank: str = FLANK2) -> str:
    return FWD2 + flank + UNIT2 * repeats + str(Seq(REV2).reverse_complement())


def test_status_codes(tmp_path):
    # reference: VNTR1 with 5 units, VNTR2 with 4 units
    (tmp_path / "ref.fasta").write_text(
        ">c1\n" + random_dna(200, 11) + make_amplicon(5) + random_dna(200, 12)
        + amplicon2(4) + random_dna(200, 13) + "\n"
    )
    # split: VNTR1 forward primer at the very end of a contig, reverse primer at the start of the next;
    # deletion: VNTR2 has 6 units but the 7 bp flank between primer and array is missing
    amp1 = make_amplicon(6)
    cut = len(FWD) + 2 * len(UNIT)
    deleted = amplicon2(6, flank="")
    (tmp_path / "split.fasta").write_text(
        ">left\n" + random_dna(300, 14) + amp1[:cut] + "\n"
        ">right\n" + amp1[cut:] + random_dna(300, 15) + "\n"
        ">other\n" + random_dna(300, 16) + deleted + random_dna(300, 17) + "\n"
    )
    (tmp_path / "p.sts").write_text(
        f"VNTR1\t{FWD}\t{REV}\t{len(make_amplicon(5))}\nVNTR2\t{FWD2}\t{REV2}\t{len(amplicon2(4))}\n"
    )
    (tmp_path / "m.tsv").write_text(f"marker\tpatterns\tunit_length\nVNTR1\t{UNIT}\t6\nVNTR2\t{UNIT2}\t7\n")
    out = tmp_path / "out"
    rc = main(["run", "--sts", str(tmp_path / "p.sts"), "--fasta", str(tmp_path / "*.fasta"), "--outdir", str(out),
               "--prefix", "t", "-N", "2", "-M", "100", "--markers", str(tmp_path / "m.tsv"),
               "--reference-sample", "ref", "--align", "none"])
    assert rc == 0
    status = pd.read_csv(out / "t_status.tsv", sep="\t", index_col=0, dtype=str)
    alleles = pd.read_csv(out / "t_alleles.tsv", sep="\t", index_col=0, dtype=str)
    assert status.loc["ref", "VNTR1"] == "OK" and alleles.loc["ref", "VNTR1"] == "5"
    assert status.loc["ref", "VNTR2"] == "OK" and alleles.loc["ref", "VNTR2"] == "4"
    assert status.loc["split", "VNTR1"] == "NEAR_CONTIG_END"
    assert alleles.loc["split", "VNTR1"] == "NEAR_CONTIG_END"
    # the 7 bp flank deletion equals one unit, so the length says 5 while the array holds 6
    assert alleles.loc["split", "VNTR2"] == "5"
    assert "COUNT_DISAGREES" in status.loc["split", "VNTR2"]
    long = pd.read_csv(out / "t_calls_long.tsv", sep="\t")
    row = long[(long["sample"] == "split") & (long["marker"] == "VNTR2")].iloc[0]
    assert row["tandem_copies"] == 6 and row["fwd_mismatches"] == 0 and row["residual_bp"] == 0
    near = (out / "near_contig_end.tsv").read_text().splitlines()
    assert near[1:] == ["split\tVNTR1"]
