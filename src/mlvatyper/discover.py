"""Discovery of repeat units from the amplicons of one sample.

For each amplicon the dominant tandem array is located by comparing the
sequence with itself shifted by every candidate period. A run of positions
where seq[i] == seq[i + p], with at most max_mismatch differences in any
window of p positions, marks an array of period p. The smallest period whose
array is nearly as long as the longest one is preferred, so that a unit is
not reported as its own multiple. The repeat unit is the per-column majority
over the windows of the array, and the final count is taken from
count_tandem_array so that it is consistent with allele calling.

This is a simple detector, not a replacement for Tandem Repeats Finder. It
is meant to seed a marker definition file that is then reviewed.
"""

import math
from collections import Counter
from dataclasses import dataclass

from Bio.Seq import Seq

from .log import get_logger
from .models import Amplicon, MarkerDef
from .tandem import TandemArray, count_tandem_array, effective_mismatch

MATCH_SCORE = 2
MISMATCH_PENALTY = 7

log = get_logger()


@dataclass(frozen=True)
class Discovery:
    marker: str
    sample: str
    amplicon_length: int
    unit: str
    array: TandemArray
    left_flank: int
    right_flank: int
    score: int = 0
    engine: str = "builtin"
    engine_copies: float | None = None  # copy number as reported by the engine, fractional for TRF
    pct_match: int | None = None
    pct_indel: int | None = None

    @property
    def unit_length(self) -> int:
        return len(self.unit)

    @property
    def entropy(self) -> float:
        return unit_entropy(self.unit)

    @property
    def copies(self) -> int:
        return self.array.copies


def _longest_self_match_run(seq: str, period: int, max_mismatch: int) -> tuple[int, int]:
    """Longest run of positions i with seq[i] == seq[i + period], tolerating mismatches.

    Returns (start, length) of the run in match positions; the array spans
    start to start + length + period. A mismatch is tolerated while the last
    period positions hold at most max_mismatch mismatches; two mismatches
    in a row end the run.
    """
    n = len(seq) - period
    best_start = best_len = 0
    i = 0
    while i < n:
        if seq[i] != seq[i + period]:
            i += 1
            continue
        start = i
        recent: list[int] = []  # positions of mismatches within the run
        j = i
        while j < n:
            if seq[j] != seq[j + period]:
                if recent and j - recent[-1] == 1:
                    break
                recent = [r for r in recent if j - r < period] + [j]
                if len(recent) > max_mismatch:
                    break
            j += 1
        # trim trailing mismatches
        end = j
        while end > start and seq[end - 1] != seq[end - 1 + period]:
            end -= 1
        if end - start > best_len:
            best_start, best_len = start, end - start
        i = max(end, start + 1)
    return best_start, best_len


def unit_entropy(unit: str) -> float:
    counts = Counter(unit.upper())
    n = len(unit)
    return -sum(c / n * math.log2(c / n) for c in counts.values()) if n else 0.0


def array_score(seq: str, array: TandemArray, period: int) -> int:
    """Tandem Repeats Finder style score of the self-match over the array."""
    matches = sum(1 for i in range(array.start, array.end - period) if seq[i] == seq[i + period])
    compared = max(0, array.span - period)
    return MATCH_SCORE * matches - MISMATCH_PENALTY * (compared - matches)


def _consensus(seq: str, start: int, period: int, copies: int) -> str:
    windows = [seq[start + k * period : start + (k + 1) * period] for k in range(copies)]
    return "".join(Counter(col).most_common(1)[0][0] for col in zip(*windows, strict=True))


def _candidate(seq: str, period: int, start: int, span: int, max_mismatch: int) -> tuple[str, TandemArray, int] | None:
    unit = _consensus(seq, start, period, span // period)
    array = count_tandem_array(seq, (unit,), max_mismatch)
    if array.copies < 2:
        return None
    if array.orientation == "-":
        unit = str(Seq(unit).reverse_complement())
        array = count_tandem_array(seq, (unit,), max_mismatch)
    unit = seq[array.start : array.start + period]  # phase of the array as found
    array = count_tandem_array(seq, (unit,), max_mismatch)
    return unit, array, array_score(seq, array, period)


def find_repeat_unit(
    sequence: str, min_unit: int = 3, max_unit: int = 100, min_copies: float = 2.0,
    max_mismatch: int = 1, min_score: int = 30, period_slack: float = 0.85,
) -> tuple[str, TandemArray, int] | None:
    """Dominant tandem repeat unit in a sequence as (unit, array, score), or None.

    Candidates need at least min_copies copies and a score of min_score. The
    smallest period whose score is within period_slack of the best score
    wins, so that a unit is not reported as its own multiple.
    """
    seq = sequence.upper()
    candidates: list[tuple[int, int, str, TandemArray]] = []  # (score, period, unit, array)
    for period in range(min_unit, min(max_unit, len(seq) // 2) + 1):
        allowed = effective_mismatch(period, max_mismatch)
        start, length = _longest_self_match_run(seq, period, allowed)
        span = length + period
        if not length or span / period < min_copies:
            continue
        found = _candidate(seq, period, start, span, max_mismatch)
        if found is None:
            continue
        unit, array, score = found
        if array.copies >= min_copies and score >= min_score:
            candidates.append((score, period, unit, array))
    if not candidates:
        return None
    best = max(c[0] for c in candidates)
    score, period, unit, array = min((c for c in candidates if c[0] >= period_slack * best), key=lambda c: c[1])
    return unit, array, score


def discover_units(
    amplicons: list[Amplicon], sample: str, min_unit: int = 3, max_unit: int = 100,
    min_copies: float = 2.0, max_mismatch: int = 1, min_score: int = 30,
) -> tuple[list[Discovery], list[str]]:
    """Discoveries for the amplicons of one sample, and markers where no array was found.

    Markers with several amplicons in the sample use copy 1.
    """
    found: list[Discovery] = []
    missing: list[str] = []
    seen: set[str] = set()
    for amp in sorted((a for a in amplicons if a.hit.sample == sample), key=lambda a: (a.hit.primer_id, a.copy)):
        marker = amp.hit.primer_id
        if marker in seen:
            continue
        seen.add(marker)
        result = find_repeat_unit(amp.sequence, min_unit, max_unit, min_copies, max_mismatch, min_score)
        if result is None:
            missing.append(marker)
            continue
        unit, array, score = result
        found.append(Discovery(marker, sample, len(amp.sequence), unit, array,
                               array.start, len(amp.sequence) - array.end, score))
    log.info("%s: repeat units found for %d markers, none for %d", sample, len(found), len(missing))
    return found, missing


def discoveries_to_definitions(found: list[Discovery], missing: list[str], amplicons: list[Amplicon],
                               sample: str) -> list[MarkerDef]:
    """Marker definitions seeded from discoveries; markers without an array keep only ref_size."""
    lengths = {a.hit.primer_id: len(a.sequence) for a in amplicons if a.hit.sample == sample and a.copy == 1}
    defs = [MarkerDef(d.marker, (d.unit,), d.unit_length, d.amplicon_length, d.copies) for d in found]
    defs.extend(MarkerDef(m, (), 0, lengths.get(m), None) for m in missing)
    return defs


REPORT_COLUMNS = [
    "marker", "sample", "amplicon_length", "unit_length", "unit", "copies", "exact_copies", "variant_copies",
    "score", "entropy", "array_start", "array_end", "array_span", "left_flank", "right_flank",
    "engine", "engine_copies", "pct_match", "pct_indel",
]


def write_report(found: list[Discovery], missing: list[str], path) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(REPORT_COLUMNS) + "\n")
        for d in found:
            a = d.array
            handle.write("\t".join(str(x) for x in [
                d.marker, d.sample, d.amplicon_length, d.unit_length, d.unit, a.copies, a.exact, a.variant,
                d.score, f"{d.entropy:.2f}", a.start, a.end, a.span, d.left_flank, d.right_flank,
                d.engine, "" if d.engine_copies is None else d.engine_copies,
                "" if d.pct_match is None else d.pct_match, "" if d.pct_indel is None else d.pct_indel,
            ]) + "\n")
        for m in missing:
            handle.write(m + "\t" * 5 + "0" + "\t" * (len(REPORT_COLUMNS) - 6) + "\n")
