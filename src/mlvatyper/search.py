"""Primer pair search using merPCR.

merPCR (https://github.com/FOI-Bioinformatics/merpcr) is a Python
reimplementation of me-PCR. Its find_hits method (version 1.0.2 or newer)
returns SearchHit objects with 1-based inclusive coordinates. The reverse
primer is searched only within the STS expected size plus or minus the
margin, so for VNTR markers the expected size and the margin together define
the allowed amplicon size window.
"""

from pathlib import Path

from .io import strip_compression
from .log import get_logger
from .models import Hit

log = get_logger()


class SearchError(RuntimeError):
    pass


def sample_name(fasta_path: str | Path) -> str:
    """Sample name is the FASTA file name without its last suffix and any compression suffix."""
    return strip_compression(fasta_path).stem


class MerPcrSearcher:
    """Thin wrapper around merpcr.MerPCR for repeated searches with one STS file."""

    def __init__(self, sts_path: str | Path, margin: int, mismatches: int, threads: int = 1, **kwargs):
        try:
            from merpcr import MerPCR
        except ImportError as exc:  # pragma: no cover - depends on environment
            raise SearchError(
                "merpcr is not installed. Install with: pip install git+https://github.com/FOI-Bioinformatics/merpcr"
            ) from exc
        self.engine = MerPCR(margin=margin, mismatches=mismatches, threads=threads, **kwargs)
        if not self.engine.load_sts_file(str(sts_path)):
            raise SearchError(f"merPCR could not load STS file {sts_path}")

    def search_fasta(self, fasta_path: str | Path) -> list[Hit]:
        sample = sample_name(fasta_path)
        records = self.engine.load_fasta_file(str(fasta_path))
        if not records:
            log.warning("%s: no FASTA records found", fasta_path)
            return []
        hits = [
            Hit(sample=sample, contig=h.seq_label, start=h.start, end=h.end, strand=h.strand, primer_id=h.sts_id,
                expected_size=h.sts.pcr_size)
            for h in self.engine.find_hits(records)
        ]
        log.info("%s: %d hits", sample, len(hits))
        return hits
