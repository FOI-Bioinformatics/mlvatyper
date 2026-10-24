import os
import random
from pathlib import Path

import pytest
from Bio.Seq import Seq

FWD = "ACGTACGGTACCAGTTGACA"
REV = "TTGCAGGCTAGCATCCTTCA"  # as written 5'->3' on the opposite strand; its reverse complement must not contain UNIT
UNIT = "GATCGA"


def make_amplicon(repeats: int) -> str:
    """Forward primer, repeat units, then reverse complement of the reverse primer."""
    return FWD + UNIT * repeats + str(Seq(REV).reverse_complement())


def random_dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture
def fixture_dir(tmp_path: Path) -> Path:
    """Two samples. sampleA has 5 repeats on + strand; sampleB has 7 repeats on - strand and an amplicon with N."""
    amp_a = make_amplicon(5)
    amp_b = str(Seq(make_amplicon(7)).reverse_complement())
    amp_n = make_amplicon(3).replace(UNIT, "GATNGA", 1)  # one unit contains N
    (tmp_path / "sampleA.fasta").write_text(
        ">contig1 some description\n" + random_dna(300, 1) + amp_a + random_dna(300, 2) + "\n"
    )
    (tmp_path / "sampleB.fasta").write_text(
        ">c1\n" + random_dna(200, 3) + amp_b + random_dna(200, 4) + "\n"
        ">c2\n" + random_dna(150, 5) + amp_n + random_dna(150, 6) + "\n"
    )
    expected = len(make_amplicon(5))
    (tmp_path / "primers.sts").write_text(
        f"VNTR1\t{FWD}\t{REV}\t{expected}\n"
        f"BAD1\t1{FWD}\t{REV}\t{expected}\n"
        f"ABSENT\tGGGGGGGGGGGGGGGGGGGG\tCCCCCCCCCCCCCCCCCCCC\t100\n"
    )
    (tmp_path / "patterns.tsv").write_text(
        "marker\tpattern\tlength\n"
        f"VNTR1\t{UNIT}\t{len(UNIT)}\n"
        "ABSENT\tN/A\n"
    )
    return tmp_path


def reference_assembly(accession: str) -> Path:
    """Path of a cached NCBI assembly for tests that need a real reference.

    The file is never committed. It is used when already in the cache, downloaded when MLVATYPER_ONLINE is
    set (as in CI), and the test is skipped otherwise.
    """
    from mlvatyper.assemblies import cached_path, download_assembly

    path = cached_path(accession)
    if path.exists():
        return path
    if os.environ.get("MLVATYPER_ONLINE"):
        return download_assembly(accession)
    pytest.skip(f"reference assembly {accession} not in the cache; set MLVATYPER_ONLINE=1 to download it")
