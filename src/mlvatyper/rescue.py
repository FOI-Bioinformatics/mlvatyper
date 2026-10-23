"""Rescue searches for markers without a hit, and reasons for markers that stay missing.

merPCR needs a run of wordsize unambiguous bases in a primer to hash it and
drops the record otherwise; it also allows no indels in a primer. Two rescue
steps run only for the markers that the main search did not find in a
sample, so that the main scan keeps its speed and specificity:

1. A second merPCR pass with a smaller word size, one-sided search (a record
   whose other primer cannot be hashed is still searched) and a wider size
   window. Hits carry search_pass "rescue".
2. An alignment of one primer with up to one indel, anchored on an exact
   IUPAC-aware match of the other primer, within the size window. Hits
   carry search_pass "aligned" and the indel counts.

Markers that remain missing are given one reason:

  NOT_SEARCHED     no primer holds a hash word of the word size used
  NEAR_CONTIG_END  a primer lies within the size window of a contig end
  ONE_PRIMER_ONLY  one primer matches exactly somewhere in the assembly
  NOT_FOUND        neither primer matches exactly
"""

import re
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path

from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq

from .log import get_logger
from .models import Hit, PrimerPair
from .search import MerPcrSearcher, SearchError, sample_name
from .sts import has_hash_word, write_sts

log = get_logger()

NOT_SEARCHED = "NOT_SEARCHED"
NEAR_CONTIG_END = "NEAR_CONTIG_END"
ONE_PRIMER_ONLY = "ONE_PRIMER_ONLY"
NOT_FOUND = "NOT_FOUND"


@dataclass(frozen=True)
class RescueOptions:
    wordsize: int = 8  # 0 disables the second merPCR pass
    margin: int = 1000
    indels: bool = False
    mismatches: int = 2
    iupac: bool = True
    three_prime_match: int = 1
    threads: int = 1


def _iupac_regex(primer: str) -> str:
    return "".join(f"[{ambiguous_dna_values.get(b, b)}]" if len(ambiguous_dna_values.get(b, b)) > 1 else b
                   for b in primer.upper())


def _matches(base: str, code: str) -> bool:
    return base in ambiguous_dna_values.get(code, code)


def rescue_pass(fasta: Path, missing: list[PrimerPair], options: RescueOptions) -> list[Hit]:
    """Second merPCR pass for the missing markers of one sample."""
    if options.wordsize <= 0 or not missing:
        return []
    with tempfile.TemporaryDirectory() as tmp:
        sts = Path(tmp) / "rescue.sts"
        write_sts(missing, sts)
        try:
            searcher = MerPcrSearcher(sts, margin=options.margin, mismatches=options.mismatches,
                                      threads=options.threads, iupac_mode=int(options.iupac),
                                      wordsize=options.wordsize, one_sided=1,
                                      three_prime_match=options.three_prime_match)
            hits = searcher.search_fasta(fasta)
        except SearchError as exc:
            log.warning("%s: rescue pass failed: %s", fasta, exc)
            return []
    if hits:
        log.info("%s: rescue pass (word size %d, one sided) found %d hits for %d markers", sample_name(fasta),
                 options.wordsize, len(hits), len({h.primer_id for h in hits}))
    return [replace(h, search_pass="rescue") for h in hits]


def align_primer(primer: str, window: str) -> tuple[int, int, int, int] | None:
    """Best placement of the primer in the window with at most one indel.

    Returns (mismatches, indels, start, end) of the aligned segment in the
    window (0-based, half open), or None. Every length of len(primer) - 1,
    len(primer) and len(primer) + 1 is tried at every position; a segment one
    base shorter or longer than the primer is aligned with one gap placed
    where it costs the fewest mismatches.
    """
    n = len(primer)
    best: tuple[int, int, int, int] | None = None
    for start in range(len(window)):
        for length in (n, n - 1, n + 1):
            end = start + length
            if end > len(window) or length < 1:
                continue
            segment = window[start:end]
            if length == n:
                result = (sum(not _matches(b, c) for b, c in zip(segment, primer, strict=True)), 0)
            else:
                result = _gapped(primer, segment)
            if best is None or (result[0] + 2 * result[1]) < (best[0] + 2 * best[1]):
                best = (result[0], result[1], start, end)
    return best


def _gapped(primer: str, segment: str) -> tuple[int, int]:
    """Mismatches with one gap in the primer (segment longer) or in the segment (segment shorter)."""
    longer, shorter = (segment, primer) if len(segment) > len(primer) else (primer, segment)
    primer_is_longer = longer is primer
    best = len(longer)
    for gap in range(len(longer)):
        skipped = longer[:gap] + longer[gap + 1 :]
        if primer_is_longer:
            mismatches = sum(not _matches(b, c) for b, c in zip(shorter, skipped, strict=True))
        else:
            mismatches = sum(not _matches(b, c) for b, c in zip(skipped, shorter, strict=True))
        best = min(best, mismatches)
    return best, 1


def _products(text: str, fwd: str, rev_rc: str, pair: PrimerPair, margin: int, max_mismatch: int):
    """Products on one strand as (start, end, fwd_mm, rev_mm, fwd_indels, rev_indels), 0-based half open.

    The product size window is the STS size range, or the expected size,
    widened by the margin on both sides.
    """
    low = pair.expected_size if pair.size_low is None else pair.size_low
    high = pair.expected_size if pair.size_high is None else pair.size_high
    found = []
    for m in re.finditer(f"(?={_iupac_regex(fwd)})", text):
        p = m.start()
        lo, hi = max(p + len(fwd), p + low - margin - len(rev_rc)), min(len(text), p + high + margin)
        if hi <= lo:
            continue
        best = align_primer(rev_rc, text[lo:hi])
        if best and best[0] <= max_mismatch and best[1] <= 1:
            found.append((p, lo + best[3], 0, best[0], 0, best[1]))
    for m in re.finditer(f"(?={_iupac_regex(rev_rc)})", text):
        e = m.start() + len(rev_rc)
        lo, hi = max(0, e - high - margin), min(e - len(rev_rc), e - low + margin + len(fwd))
        if hi <= lo:
            continue
        best = align_primer(fwd, text[lo:hi])
        if best and best[0] <= max_mismatch and best[1] <= 1:
            found.append((lo + best[2], e, best[0], 0, best[1], 0))
    return found


def find_with_indels(contigs: dict[str, str], pair: PrimerPair, sample: str, margin: int,
                     max_mismatch: int = 2) -> list[Hit]:
    """Hits for one marker where one primer matches exactly and the other with up to one indel.

    A product without indels is also reported, since the marker reached this
    step only because merPCR did not find it.
    """
    fwd = pair.forward.upper()
    rev_rc = str(Seq(pair.reverse.upper()).reverse_complement())
    hits: list[Hit] = []
    for contig, seq in contigs.items():
        seq = seq.upper()
        for strand, text in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
            for s, e, fmm, rmm, find, rind in _products(text, fwd, rev_rc, pair, margin, max_mismatch):
                start, end = (s + 1, e) if strand == "+" else (len(seq) - e + 1, len(seq) - s)
                hits.append(Hit(sample, contig, start, end, strand, pair.id, pair.expected_size, fmm, rmm,
                                fwd_indels=find, rev_indels=rind, search_pass="aligned"))
    return sorted(hits, key=lambda h: (h.contig, h.start))


def _primer_present(text: str, primer: str) -> bool:
    rc = str(Seq(primer).reverse_complement())
    return bool(re.search(_iupac_regex(primer), text) or re.search(_iupac_regex(rc), text))


def missing_reasons(contigs: dict[str, str], missing: list[PrimerPair], near_end: set[str], wordsize: int,
                    one_sided: bool) -> dict[str, str]:
    """Reason per missing marker; wordsize and one_sided describe the most permissive merPCR pass run."""
    text = "\n".join(s.upper() for s in contigs.values())
    reasons: dict[str, str] = {}
    for pair in missing:
        hashable = [has_hash_word(pair.forward, wordsize), has_hash_word(pair.reverse, wordsize)]
        searched = any(hashable) if one_sided else all(hashable)
        if not searched:
            reasons[pair.id] = NOT_SEARCHED
        elif pair.id in near_end:
            reasons[pair.id] = NEAR_CONTIG_END
        elif _primer_present(text, pair.forward) or _primer_present(text, pair.reverse):
            reasons[pair.id] = ONE_PRIMER_ONLY
        else:
            reasons[pair.id] = NOT_FOUND
    return reasons


MISSING_COLUMNS = ["sample", "marker", "reason"]


def write_missing(reasons: dict[tuple[str, str], str], path: Path) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(MISSING_COLUMNS) + "\n")
        for (sample, marker), reason in sorted(reasons.items()):
            handle.write(f"{sample}\t{marker}\t{reason}\n")


def read_missing(path: Path) -> dict[tuple[str, str], str]:
    if not path.exists():
        return {}
    reasons = {}
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        fields = line.split("\t")
        if len(fields) >= 3:
            reasons[(fields[0], fields[1])] = fields[2]
    return reasons
