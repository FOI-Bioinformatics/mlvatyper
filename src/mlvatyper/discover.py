"""Discovery of repeat units from the amplicons of one sample.

Every tandem array that an engine reports in an amplicon becomes a
candidate. The builtin engine compares the sequence with itself shifted by
every candidate period and keeps, per period, the longest run of positions
where seq[i] == seq[i + p] with at most max_mismatch differences in any
window of p positions. Tandem Repeats Finder (trf.py) is the alternative
engine. Candidates that describe the same array at several periods are
merged, keeping the smallest period whose score is within period_slack of
the best one, so that a unit is not reported as its own multiple.

Each candidate is counted with count_tandem_array so that its copy number is
the integer that allele calling would see, and is then placed in a category:

  A  unit of at least strong_unit bp, at least strong_copies integer copies,
     score at least min_score, entropy at least min_entropy, and no overlap
     with a primer
  B  unit of at least strong_unit bp, at least two integer copies and score
     at least min_score
  C  everything else above the weak score floor: short units, weak scores
     and arrays reaching into a primer

The candidates of a marker are ranked by category, cross-sample support
(discover_multi.py), integer copies and score; the first one is the primary
and seeds the marker definition when its category is accepted. A marker with
two or more non-overlapping A or B candidates is reported as mixed (M).

A single reference genome typically holds two or three copies of each unit,
where chance tandems score the same as the true VNTR. Categories only order
the evidence; the comparison of several assemblies in discover_multi.py is
what identifies the array that varies.
"""

import math
from collections import Counter
from dataclasses import dataclass, replace

from Bio.Seq import Seq

from .log import get_logger
from .models import Amplicon, MarkerDef
from .tandem import TandemArray, count_tandem_array, effective_mismatch, rotations

MATCH_SCORE = 2
MISMATCH_PENALTY = 7
CATEGORY_RANK = {"A": 0, "B": 1, "C": 2}

log = get_logger()


@dataclass(frozen=True)
class Thresholds:
    """Category thresholds; see the module docstring."""

    min_score: int = 30
    weak_score: int = 20
    min_copies: float = 1.8  # engine copies, fractional
    strong_unit: int = 5
    strong_copies: int = 3
    min_entropy: float = 1.0
    accept: str = "B"  # least category that seeds a definition


@dataclass(frozen=True)
class RawArray:
    """A tandem array as reported by an engine, 0-based half-open coordinates."""

    start: int
    end: int
    period: int
    consensus: str
    score: int
    copies: float  # fractional copy number from the engine
    pct_match: int | None = None
    pct_indel: int | None = None

    @property
    def span(self) -> int:
        return self.end - self.start


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
    primer_overlap: bool = False
    category: str = ""
    rank: int = 1
    primary: bool = True
    n_samples: int | None = None  # filled by discover_multi
    copy_min: int | None = None
    copy_max: int | None = None
    explains: float | None = None

    @property
    def unit_length(self) -> int:
        return len(self.unit)

    @property
    def entropy(self) -> float:
        return unit_entropy(self.unit)

    @property
    def copies(self) -> int:
        return self.array.copies

    @property
    def unit_key(self) -> str:
        return canonical_unit(self.unit)


def canonical_unit(unit: str) -> str:
    """Lexicographically smallest rotation of the unit or of its reverse complement."""
    unit = unit.upper()
    rc = str(Seq(unit).reverse_complement())
    return min(rotations(unit) | rotations(rc))


def canonical_rotation(unit: str) -> str:
    """Lexicographically smallest rotation of the unit in its own orientation."""
    return min(rotations(unit.upper()))


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


def builtin_arrays(
    sequence: str, min_unit: int = 3, max_unit: int = 100, max_mismatch: int = 1, min_score: int = 20,
) -> list[RawArray]:
    """One array per period from the self-comparison, before harmonic merging."""
    seq = sequence.upper()
    arrays: list[RawArray] = []
    for period in range(min_unit, min(max_unit, len(seq) // 2) + 1):
        allowed = effective_mismatch(period, max_mismatch)
        start, length = _longest_self_match_run(seq, period, allowed)
        span = length + period
        if not length or span < 2 * period:
            continue
        copies = span // period
        unit = _consensus(seq, start, period, copies)
        array = count_tandem_array(seq, (unit,), max_mismatch)
        if array.copies < 2:
            continue
        score = array_score(seq, array, period)
        if score < min_score:
            continue
        arrays.append(RawArray(array.start, array.end, period, unit, score, round(array.span / period, 1)))
    return arrays


def _overlap(a: RawArray, b: RawArray) -> int:
    return max(0, min(a.end, b.end) - max(a.start, b.start))


def merge_harmonics(arrays: list[RawArray], period_slack: float = 0.85) -> list[RawArray]:
    """One array per locus: cluster arrays overlapping by more than half the shorter span.

    Within a cluster the smallest period whose score is within period_slack of
    the best score is kept, as in select_repeat for TRF.
    """
    clusters: list[list[RawArray]] = []
    for a in sorted(arrays, key=lambda r: (r.start, r.end)):
        for cluster in clusters:
            if any(_overlap(a, c) > 0.5 * min(a.span, c.span) for c in cluster):
                cluster.append(a)
                break
        else:
            clusters.append([a])
    merged = []
    for cluster in clusters:
        best = max(r.score for r in cluster)
        merged.append(min((r for r in cluster if r.score >= period_slack * best), key=lambda r: (r.period, -r.score)))
    return merged


def _oriented(seq: str, unit: str, max_mismatch: int) -> tuple[str, TandemArray]:
    """Count the unit in the amplicon and return it in the orientation of the amplicon."""
    array = count_tandem_array(seq, (unit,), max_mismatch)
    if array.orientation == "-":
        unit = str(Seq(unit).reverse_complement())
        array = count_tandem_array(seq, (unit,), max_mismatch)
    return canonical_rotation(unit), array


def make_discovery(
    marker: str, sample: str, seq: str, raw: RawArray, engine: str, max_mismatch: int = 1,
    primer_lengths: tuple[int, int] | None = None, min_copies: int = 2,
) -> Discovery | None:
    """Discovery for an engine array, or None when the tandem counter finds fewer than min_copies copies."""
    unit, array = _oriented(seq, raw.consensus, max_mismatch)
    if array.copies < min_copies:
        return None
    overlap = False
    if primer_lengths:  # at least one whole unit inside a primer; a one base shift of phase is not an overlap
        overlap = (primer_lengths[0] - array.start >= len(unit)
                   or array.end - (len(seq) - primer_lengths[1]) >= len(unit))
    return Discovery(
        marker, sample, len(seq), unit, array, array.start, len(seq) - array.end, raw.score, engine,
        raw.copies, raw.pct_match, raw.pct_indel, overlap,
    )


def categorize(d: Discovery, t: Thresholds) -> str:
    strong_unit = d.unit_length >= t.strong_unit
    if (strong_unit and d.copies >= t.strong_copies and d.score >= t.min_score and d.entropy >= t.min_entropy
            and not d.primer_overlap):
        return "A"
    if strong_unit and d.copies >= 2 and d.score >= t.min_score and not d.primer_overlap:
        return "B"
    return "C"


def support(d: Discovery) -> int:
    """Cross-sample support: 2 when the unit varies between samples and its period explains the
    length differences, 1 when the unit is seen in several samples, 0 otherwise."""
    if d.n_samples is None or d.n_samples < 2:
        return 0
    varies = d.copy_min is not None and d.copy_max is not None and d.copy_max > d.copy_min
    if varies and d.explains is not None and d.explains >= 0.8:
        return 2
    return 1


def rank_candidates(cands: list[Discovery]) -> list[Discovery]:
    """Rank the candidates of one marker and mark the primary; returns them in rank order."""
    ordered = sorted(cands, key=lambda d: (CATEGORY_RANK[d.category], -support(d), -d.copies, -d.score, d.unit))
    return [replace(d, rank=n + 1, primary=(n == 0)) for n, d in enumerate(ordered)]


def _distinct_strong(cands: list[Discovery]) -> int:
    """Number of non-overlapping A or B candidates."""
    strong = sorted((d for d in cands if d.category in "AB"), key=lambda d: d.array.start)
    n = 0
    last_end = -1
    for d in strong:
        if d.array.start >= last_end:
            n += 1
            last_end = d.array.end
    return n


def is_mixed(cands: list[Discovery]) -> bool:
    return _distinct_strong(cands) >= 2


def finish_candidates(cands: list[Discovery], t: Thresholds) -> list[Discovery]:
    """Categorise and rank the candidates of one marker."""
    return rank_candidates([replace(d, category=categorize(d, t)) for d in cands])


def discover_units(
    amplicons: list[Amplicon], sample: str, min_unit: int = 3, max_unit: int = 100,
    min_copies: float = 1.8, max_mismatch: int = 1, min_score: int = 30,
    thresholds: Thresholds | None = None, primer_lengths: dict[str, tuple[int, int]] | None = None,
) -> tuple[dict[str, list[Discovery]], list[str]]:
    """Candidates per marker for the amplicons of one sample with the builtin engine.

    Returns ({marker: ranked candidates}, markers without any candidate).
    Markers with several amplicons in the sample use copy 1.
    """
    t = thresholds or Thresholds(min_score=min_score, min_copies=min_copies)
    found: dict[str, list[Discovery]] = {}
    missing: list[str] = []
    for marker, amp in chosen_amplicons(amplicons, sample).items():
        seq = amp.sequence.upper()
        raw = merge_harmonics(builtin_arrays(seq, min_unit, max_unit, max_mismatch, t.weak_score))
        cands = [make_discovery(marker, sample, seq, r, "builtin", max_mismatch,
                                (primer_lengths or {}).get(marker)) for r in raw if r.copies >= t.min_copies]
        cands = finish_candidates([c for c in cands if c is not None], t)
        if cands:
            found[marker] = cands
        else:
            missing.append(marker)
    log.info("%s: repeat units found for %d markers, none for %d", sample, len(found), len(missing))
    return found, missing


def chosen_amplicons(amplicons: list[Amplicon], sample: str) -> dict[str, Amplicon]:
    """Copy 1 of every marker in the sample, in marker order."""
    chosen: dict[str, Amplicon] = {}
    for amp in sorted((a for a in amplicons if a.hit.sample == sample), key=lambda a: (a.hit.primer_id, a.copy)):
        chosen.setdefault(amp.hit.primer_id, amp)
    return chosen


def primaries(found: dict[str, list[Discovery]]) -> list[Discovery]:
    return [cands[0] for cands in found.values() if cands]


def accepted(d: Discovery, t: Thresholds) -> bool:
    return CATEGORY_RANK[d.category] <= CATEGORY_RANK[t.accept]


def discoveries_to_definitions(
    found: dict[str, list[Discovery]] | list[Discovery], missing: list[str], amplicons: list[Amplicon],
    sample: str, thresholds: Thresholds | None = None,
) -> list[MarkerDef]:
    """Marker definitions seeded from the primary candidates.

    Markers without an accepted candidate keep only ref_size. A plain list of
    discoveries is accepted for backwards compatibility and treated as one
    primary per marker.
    """
    t = thresholds or Thresholds()
    groups = found if isinstance(found, dict) else {d.marker: [d] for d in found}
    lengths = {a.hit.primer_id: len(a.sequence) for a in amplicons if a.hit.sample == sample and a.copy == 1}
    defs = []
    rest = list(missing)
    for marker, cands in groups.items():
        if not cands:
            continue
        d = cands[0]
        if d.category and not accepted(d, t):
            rest.append(marker)
            continue
        defs.append(MarkerDef(marker, companion_units(cands, t), d.unit_length, d.amplicon_length, d.copies))
    defs.extend(MarkerDef(m, (), 0, lengths.get(m), None) for m in rest)
    return defs


def companion_units(cands: list[Discovery], t: Thresholds) -> tuple[str, ...]:
    """The primary unit followed by accepted, non-overlapping units of the same length.

    In a mixed marker whose arrays share the period, the tandem counter then
    accepts windows of either unit, so the count covers the whole region.
    Units of another period cannot be combined and are left to the report.
    """
    d = cands[0]
    units = [d.unit]
    for other in cands[1:]:
        same_period = other.unit_length == d.unit_length and other.unit_key != d.unit_key
        disjoint = other.array.end <= d.array.start or other.array.start >= d.array.end
        if same_period and disjoint and accepted(other, t) and other.unit not in units:
            units.append(other.unit)
    return tuple(units)


CANDIDATE_COLUMNS = [
    "marker", "sample", "amplicon_length", "unit_length", "unit", "copies", "exact_copies", "variant_copies",
    "score", "entropy", "array_start", "array_end", "array_span", "left_flank", "right_flank",
    "engine", "engine_copies", "pct_match", "pct_indel",
    "primer_overlap", "category", "rank", "primary", "n_samples", "copy_range", "explains",
]
REPORT_COLUMNS = CANDIDATE_COLUMNS[:-6] + ["category", "n_candidates", "mixed", "n_samples", "copy_range", "explains"]
# the report keeps the candidate columns up to primer_overlap, then describes the primary and the marker


def _cells(d: Discovery) -> list[str]:
    a = d.array
    return [str(x) for x in [
        d.marker, d.sample, d.amplicon_length, d.unit_length, d.unit, a.copies, a.exact, a.variant,
        d.score, f"{d.entropy:.2f}", a.start, a.end, a.span, d.left_flank, d.right_flank,
        d.engine, "" if d.engine_copies is None else d.engine_copies,
        "" if d.pct_match is None else d.pct_match, "" if d.pct_indel is None else d.pct_indel,
    ]]


def _multi_cells(d: Discovery) -> list[str]:
    return [
        "" if d.n_samples is None else str(d.n_samples),
        "" if d.copy_min is None else f"{d.copy_min}-{d.copy_max}",
        "" if d.explains is None else f"{d.explains:.2f}",
    ]


def write_candidates(found: dict[str, list[Discovery]], path) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(CANDIDATE_COLUMNS) + "\n")
        for cands in found.values():
            for d in cands:
                row = _cells(d) + [str(int(d.primer_overlap)), d.category, str(d.rank), str(int(d.primary))]
                handle.write("\t".join(row + _multi_cells(d)) + "\n")


def write_report(found: dict[str, list[Discovery]] | list[Discovery], missing: list[str], path) -> None:
    """One row per marker: the primary candidate, or an empty row for markers without one."""
    groups = found if isinstance(found, dict) else {d.marker: [d] for d in found}
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(REPORT_COLUMNS) + "\n")
        for cands in groups.values():
            if not cands:
                continue
            d = cands[0]
            row = _cells(d) + [str(int(d.primer_overlap)), d.category, str(len(cands)), str(int(is_mixed(cands)))]
            handle.write("\t".join(row + _multi_cells(d)) + "\n")
        for m in missing:
            handle.write(m + "\t" * 5 + "0" + "\t" * (len(REPORT_COLUMNS) - 6) + "\n")
