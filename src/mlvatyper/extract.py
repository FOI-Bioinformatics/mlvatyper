"""Amplicon extraction, hit deduplication and the hits table."""

import csv
import re
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .io import open_text
from .log import get_logger
from .models import Amplicon, Hit, PrimerPair

log = get_logger()

HITS_COLUMNS = [
    "sample", "marker", "copy", "contig", "start", "end", "strand", "length", "has_N",
    "expected_size", "fwd_mismatches", "rev_mismatches", "n_alternatives", "dist_to_contig_end", "sequence",
]


def load_fasta(path: str | Path) -> dict[str, str]:
    """Load a FASTA file once, keyed by record id (text up to the first whitespace).

    Files ending in .gz, .bz2 or .xz are decompressed on the fly.
    """
    with open_text(path) as handle:
        return {rec.id: str(rec.seq) for rec in SeqIO.parse(handle, "fasta")}


def primer_mismatches(observed: str, primer: str) -> int:
    """Positions where the observed base is not among the bases a primer IUPAC code stands for."""
    if len(observed) != len(primer):
        return max(len(observed), len(primer))
    count = 0
    for o, p in zip(observed.upper(), primer.upper(), strict=True):
        if o not in ambiguous_dna_values.get(p, p):
            count += 1
    return count


def annotate_hit(hit: Hit, fragment: str, contig_length: int, pair: PrimerPair | None) -> Hit:
    """Fill mismatch counts, expected size and distance to the nearest contig end."""
    dist = min(hit.start - 1, contig_length - hit.end)
    if pair is None:
        return replace(hit, dist_to_contig_end=dist)
    rev_rc = str(Seq(pair.reverse).reverse_complement())
    return replace(
        hit,
        expected_size=pair.expected_size if hit.expected_size < 0 else hit.expected_size,
        fwd_mismatches=primer_mismatches(fragment[: len(pair.forward)], pair.forward),
        rev_mismatches=primer_mismatches(fragment[len(fragment) - len(rev_rc) :], rev_rc),
        dist_to_contig_end=dist,
    )


def _forward_primer_position(hit: Hit) -> int:
    return hit.start if hit.strand == "+" else hit.end


def _overlap(a: Hit, b: Hit) -> int:
    return max(0, min(a.end, b.end) - max(a.start, b.start) + 1)


def _rank(hit: Hit) -> tuple[int, int, int]:
    mismatches = hit.total_mismatches if hit.total_mismatches >= 0 else 0
    deviation = abs(hit.size_deviation) if hit.size_deviation is not None else 0
    return (mismatches, deviation, hit.start)


def deduplicate_hits(hits: list[Hit], primer_length: int = 20) -> list[Hit]:
    """Collapse hits that describe the same locus.

    Hits for one marker on one contig and strand are the same locus when they
    share the forward primer position or overlap by more than one primer
    length. The hit with fewest primer mismatches, then smallest deviation
    from the expected size, is kept; the number of discarded hits is recorded
    in n_alternatives.
    """
    groups: dict[tuple[str, str, str, str], list[Hit]] = defaultdict(list)
    for h in hits:
        groups[(h.sample, h.primer_id, h.contig, h.strand)].append(h)
    kept: list[Hit] = []
    for group in groups.values():
        clusters: list[list[Hit]] = []
        for h in sorted(group, key=lambda x: x.start):
            for cluster in clusters:
                if any(_forward_primer_position(c) == _forward_primer_position(h) or _overlap(c, h) > primer_length
                       for c in cluster):
                    cluster.append(h)
                    break
            else:
                clusters.append([h])
        for cluster in clusters:
            best = min(cluster, key=_rank)
            kept.append(replace(best, n_alternatives=len(cluster) - 1))
    return sorted(kept, key=lambda h: (h.sample, h.primer_id, h.contig, h.start))


def extract_amplicons(hits: list[Hit], contigs: dict[str, str], primers: dict[str, PrimerPair] | None = None,
                      deduplicate: bool = True) -> list[Amplicon]:
    """Slice amplicons for one sample. Copy numbers follow hit order per marker."""
    annotated: list[tuple[Hit, str]] = []
    for hit in hits:
        seq = contigs.get(hit.contig)
        if seq is None:
            log.warning("%s: contig %s not found in FASTA, hit for %s skipped", hit.sample, hit.contig, hit.primer_id)
            continue
        fragment = seq[hit.start - 1 : hit.end]
        if hit.strand == "-":
            fragment = str(Seq(fragment).reverse_complement())
        pair = primers.get(hit.primer_id) if primers else None
        annotated.append((annotate_hit(hit, fragment, len(seq), pair), fragment))
    def key(h: Hit) -> tuple[str, str, int, int, str]:
        return (h.primer_id, h.contig, h.start, h.end, h.strand)

    fragments = {key(h): f for h, f in annotated}
    final = [h for h, _ in annotated]
    if deduplicate:
        plen = max((len(p.forward) for p in primers.values()), default=20) if primers else 20
        before = len(final)
        final = deduplicate_hits(final, plen)
        if len(final) < before:
            log.info("%s: %d overlapping hits collapsed", final[0].sample if final else "", before - len(final))
    counter: dict[str, int] = defaultdict(int)
    amplicons = []
    for hit in final:
        counter[hit.primer_id] += 1
        amplicons.append(Amplicon(hit=hit, copy=counter[hit.primer_id], sequence=fragments[key(hit)]))
    return amplicons


def _iupac_regex(primer: str) -> str:
    return "".join(f"[{ambiguous_dna_values.get(b, b)}]" if len(ambiguous_dna_values.get(b, b)) > 1 else b
                   for b in primer.upper())


def primers_near_contig_ends(contigs: dict[str, str], pairs: list[PrimerPair], window: int) -> set[str]:
    """Markers for which either primer (any orientation) lies within window bp of a contig end.

    A missing marker with a primer near a contig end is more likely split by
    the assembly than deleted. Exact IUPAC-aware matching is used.
    """
    found: set[str] = set()
    ends = []
    for seq in contigs.values():
        seq = seq.upper()
        if len(seq) <= 2 * window:
            ends.append(seq)
        else:
            ends.append(seq[:window])
            ends.append(seq[-window:])
    text = "\n".join(ends)
    for pair in pairs:
        for primer in (pair.forward, pair.reverse):
            rc = str(Seq(primer).reverse_complement())
            if re.search(_iupac_regex(primer), text) or re.search(_iupac_regex(rc), text):
                found.add(pair.id)
                break
    return found


def write_hits_table(amplicons: list[Amplicon], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HITS_COLUMNS)
        for amp in amplicons:
            h = amp.hit
            writer.writerow([
                h.sample, h.primer_id, amp.copy, h.contig, h.start, h.end, h.strand,
                len(amp.sequence), int(amp.has_n),
                h.expected_size, h.fwd_mismatches, h.rev_mismatches, h.n_alternatives, h.dist_to_contig_end,
                amp.sequence,
            ])


def _int(row: dict[str, str], key: str, default: int) -> int:
    value = row.get(key)
    return int(value) if value not in (None, "") else default


def read_hits_table(path: str | Path) -> list[Amplicon]:
    amplicons = []
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            hit = Hit(
                sample=row["sample"], contig=row["contig"], start=int(row["start"]),
                end=int(row["end"]), strand=row["strand"], primer_id=row["marker"],
                expected_size=_int(row, "expected_size", -1),
                fwd_mismatches=_int(row, "fwd_mismatches", -1),
                rev_mismatches=_int(row, "rev_mismatches", -1),
                dist_to_contig_end=_int(row, "dist_to_contig_end", -1),
                n_alternatives=_int(row, "n_alternatives", 0),
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
