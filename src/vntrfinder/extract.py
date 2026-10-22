"""Amplicon extraction and the hits table."""

import csv
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .io import open_text
from .log import get_logger
from .models import Amplicon, Hit

log = get_logger()

HITS_COLUMNS = [
    "sample", "marker", "copy", "contig", "start", "end", "strand", "length", "has_N", "sequence",
]


def load_fasta(path: str | Path) -> dict[str, str]:
    """Load a FASTA file once, keyed by record id (text up to the first whitespace).

    Files ending in .gz, .bz2 or .xz are decompressed on the fly.
    """
    with open_text(path) as handle:
        return {rec.id: str(rec.seq) for rec in SeqIO.parse(handle, "fasta")}


def extract_amplicons(hits: list[Hit], contigs: dict[str, str]) -> list[Amplicon]:
    """Slice amplicons for one sample. Copy numbers follow hit order per marker."""
    counter: dict[str, int] = defaultdict(int)
    amplicons = []
    for hit in hits:
        seq = contigs.get(hit.contig)
        if seq is None:
            log.warning("%s: contig %s not found in FASTA, hit for %s skipped", hit.sample, hit.contig, hit.primer_id)
            continue
        fragment = seq[hit.start - 1 : hit.end]
        if hit.strand == "-":
            fragment = str(Seq(fragment).reverse_complement())
        counter[hit.primer_id] += 1
        amplicons.append(Amplicon(hit=hit, copy=counter[hit.primer_id], sequence=fragment))
    return amplicons


def write_hits_table(amplicons: list[Amplicon], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HITS_COLUMNS)
        for amp in amplicons:
            h = amp.hit
            writer.writerow([
                h.sample, h.primer_id, amp.copy, h.contig, h.start, h.end, h.strand,
                len(amp.sequence), int(amp.has_n), amp.sequence,
            ])


def read_hits_table(path: str | Path) -> list[Amplicon]:
    amplicons = []
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            hit = Hit(
                sample=row["sample"], contig=row["contig"], start=int(row["start"]),
                end=int(row["end"]), strand=row["strand"], primer_id=row["marker"],
            )
            amplicons.append(Amplicon(hit=hit, copy=int(row["copy"]), sequence=row["sequence"]))
    return amplicons


def write_samples(samples: list[str], path: str | Path) -> None:
    Path(path).write_text("\n".join(samples) + "\n", encoding="utf-8")


def read_samples(path: str | Path) -> list[str]:
    return [s for s in Path(path).read_text(encoding="utf-8").splitlines() if s.strip()]


def group_by_marker(amplicons: list[Amplicon]) -> dict[str, list[Amplicon]]:
    groups: dict[str, list[Amplicon]] = defaultdict(list)
    for amp in amplicons:
        groups[amp.hit.primer_id].append(amp)
    return dict(groups)


def write_marker_fastas(amplicons: list[Amplicon], seq_dir: str | Path) -> dict[str, Path]:
    """Write sequences/<marker>.fasta with record ids <sample>_<copy>."""
    seq_dir = Path(seq_dir)
    seq_dir.mkdir(parents=True, exist_ok=True)
    written = {}
    for marker, group in group_by_marker(amplicons).items():
        records = [SeqRecord(Seq(a.sequence), id=a.record_id, description="") for a in group]
        path = seq_dir / f"{marker}.fasta"
        SeqIO.write(records, str(path), "fasta")
        written[marker] = path
    return written
