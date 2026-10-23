"""Repeat unit discovery across several samples.

Candidates are found per sample (discover.py or trf.py) and pooled per
marker by unit, where units equal up to rotation and reverse complement are
the same unit. For each pooled unit two quantities summarise the evidence
from the other assemblies:

  copy range  the smallest and largest integer copy number over the samples
              in which the unit was found
  explains    among sample pairs whose amplicons differ in length, the
              fraction whose length difference is a multiple of the unit
              length within the allele tolerance; None when all amplicons
              have the same length

A unit seen in two or more samples with a copy range of at least one and
explains of at least 0.8 is placed in category A whatever its copy number
in the reference sample. When another unit of the marker reaches that
level and this one explains less than half of the differences, it is
demoted by one category. The reference sample provides ref_size and
ref_repeats; when the best supported unit was not found as an array in the
reference amplicon, it is counted there directly so that the definition
can still be anchored.
"""

from collections import defaultdict
from dataclasses import dataclass, replace
from itertools import combinations

from .discover import (
    CATEGORY_RANK, Discovery, RawArray, chosen_amplicons, is_mixed, make_discovery, rank_candidates, support,
)
from .log import get_logger
from .models import Amplicon

log = get_logger()

Found = dict[str, list[Discovery]]  # marker -> ranked candidates


@dataclass(frozen=True)
class PooledUnit:
    key: str
    period: int
    samples: tuple[str, ...]
    copy_min: int
    copy_max: int
    explains: float | None


def explained_fraction(lengths: dict[str, int], period: int, tolerance: float) -> float | None:
    """Fraction of sample pairs with different amplicon lengths whose difference is a multiple of period."""
    pairs = [(a, b) for a, b in combinations(sorted(lengths), 2) if lengths[a] != lengths[b]]
    if not pairs:
        return None
    hits = 0
    for a, b in pairs:
        ratio = abs(lengths[a] - lengths[b]) / period
        if round(ratio) >= 1 and abs(ratio - round(ratio)) < tolerance:
            hits += 1
    return hits / len(pairs)


def pool_units(
    by_sample: dict[str, Found], lengths: dict[str, dict[str, int]], tolerance: float = 0.25,
) -> dict[str, dict[str, PooledUnit]]:
    """Pooled units per marker. lengths maps marker -> sample -> copy 1 amplicon length."""
    seen: dict[str, dict[str, list[Discovery]]] = defaultdict(lambda: defaultdict(list))
    for found in by_sample.values():
        for marker, cands in found.items():
            for d in cands:
                seen[marker][d.unit_key].append(d)
    pooled: dict[str, dict[str, PooledUnit]] = {}
    for marker, units in seen.items():
        pooled[marker] = {}
        for key, cands in units.items():
            period = len(key)
            pooled[marker][key] = PooledUnit(
                key, period, tuple(sorted({d.sample for d in cands})),
                min(d.copies for d in cands), max(d.copies for d in cands),
                explained_fraction(lengths.get(marker, {}), period, tolerance),
            )
    return pooled


def _adjust_category(d: Discovery, best_explains: float) -> str:
    if support(d) == 2:
        return "A"
    if best_explains >= 0.8 and d.explains is not None and d.explains < 0.5 and d.category != "C":
        return "B" if d.category == "A" else "C"
    return d.category


def apply_support(by_sample: dict[str, Found], pooled: dict[str, dict[str, PooledUnit]]) -> dict[str, Found]:
    """Fill the cross-sample columns, promote or demote, and re-rank every sample's candidates."""
    result: dict[str, Found] = {}
    for sample, found in by_sample.items():
        result[sample] = {}
        for marker, cands in found.items():
            units = pooled.get(marker, {})
            best = max((u.explains or 0.0 for u in units.values()), default=0.0)
            updated = []
            for d in cands:
                u = units[d.unit_key]
                d = replace(d, n_samples=len(u.samples), copy_min=u.copy_min, copy_max=u.copy_max, explains=u.explains)
                updated.append(replace(d, category=_adjust_category(d, best)))
            result[sample][marker] = rank_candidates(updated)
    return result


def best_pooled_unit(by_sample: dict[str, Found], marker: str) -> Discovery | None:
    """The top ranked candidate of the marker over all samples."""
    cands = [c[0] for found in by_sample.values() for m, c in found.items() if m == marker and c]
    if not cands:
        return None
    return min(cands, key=lambda d: (CATEGORY_RANK[d.category], -support(d), -d.copies, -d.score, d.unit))


def reference_candidates(
    by_sample: dict[str, Found], reference: str, amplicons: list[Amplicon], max_mismatch: int = 1,
    primer_lengths: dict[str, tuple[int, int]] | None = None,
) -> tuple[Found, list[str]]:
    """Candidates of the reference sample, with the best pooled unit added when it is missing there.

    A unit counted directly in the reference amplicon keeps the category and
    cross-sample columns of its best candidate elsewhere and is marked with
    engine "pooled"; its copy number may be below two.
    """
    found = dict(by_sample.get(reference, {}))
    chosen = chosen_amplicons(amplicons, reference)
    missing: list[str] = []
    markers = {m for f in by_sample.values() for m in f} | set(chosen)
    for marker in sorted(markers):
        amp = chosen.get(marker)
        if amp is None:
            continue
        best = best_pooled_unit(by_sample, marker)
        if best is None:
            missing.append(marker)
            continue
        if any(d.unit_key == best.unit_key for d in found.get(marker, [])):
            continue
        seq = amp.sequence.upper()
        raw = RawArray(0, 0, best.unit_length, best.unit, best.score, 0.0)
        d = make_discovery(marker, reference, seq, raw, "pooled", max_mismatch, (primer_lengths or {}).get(marker),
                           min_copies=0)
        if d is None:
            missing.append(marker)
            continue
        d = replace(d, category=best.category, n_samples=best.n_samples, copy_min=best.copy_min,
                    copy_max=best.copy_max, explains=best.explains)
        found[marker] = rank_candidates(found.get(marker, []) + [d])
    log.info("%s: %d markers with a candidate in the reference, %d without", reference, len(found), len(missing))
    return found, [m for m in missing if m not in found]


def mixed_markers(found: Found) -> list[str]:
    return [m for m, cands in found.items() if is_mixed(cands)]
