"""Mismatch tolerant counting of tandem repeat arrays.

An array is the longest run of consecutive, non-overlapping windows of the
unit length in which every window is within a Hamming distance of
max_mismatch from some rotation of one of the given unit variants. The
allowance is scaled down for short units (one mismatch per six bases at
most), since one difference in a trinucleotide would match a quarter of all
windows. Both
orientations of the units are tried. The count is an independent estimate
of the repeat number that does not depend on the amplicon length and is used
as a cross-check of the length based allele.
"""

from dataclasses import dataclass

from Bio.Seq import Seq


@dataclass(frozen=True)
class TandemArray:
    copies: int
    start: int  # 0-based start of the array in the sequence
    end: int  # 0-based exclusive end
    exact: int  # windows identical to a variant
    variant: int  # windows with 1..max_mismatch differences
    orientation: str  # "+" if the units matched as given, "-" if reverse complemented

    @property
    def span(self) -> int:
        return self.end - self.start


EMPTY = TandemArray(0, 0, 0, 0, 0, "+")


def rotations(unit: str) -> set[str]:
    return {unit[i:] + unit[:i] for i in range(len(unit))}


def _hamming(a: str, b: str, limit: int) -> int:
    d = 0
    for x, y in zip(a, b, strict=True):
        if x != y:
            d += 1
            if d > limit:
                return d
    return d


def _best_distance(window: str, variants: set[str], limit: int) -> int:
    best = limit + 1
    for v in variants:
        d = _hamming(window, v, limit)
        if d < best:
            best = d
            if best == 0:
                break
    return best


def _longest_run(seq: str, variants: set[str], length: int, max_mismatch: int, orientation: str) -> TandemArray:
    n = len(seq)
    distances = [_best_distance(seq[i : i + length], variants, max_mismatch) for i in range(n - length + 1)]
    best = EMPTY
    for start in range(len(distances)):
        if distances[start] > max_mismatch:
            continue
        copies = exact = 0
        pos = start
        while pos < len(distances) and distances[pos] <= max_mismatch:
            copies += 1
            exact += distances[pos] == 0
            pos += length
        if copies > best.copies or (copies == best.copies and exact > best.exact):
            best = TandemArray(copies, start, start + copies * length, exact, copies - exact, orientation)
    return best


def effective_mismatch(unit_length: int, max_mismatch: int) -> int:
    """Mismatches allowed per window: at most one per six bases of unit length."""
    return max(0, min(max_mismatch, unit_length // 6))


def count_tandem_array(sequence: str, units: tuple[str, ...], max_mismatch: int = 1) -> TandemArray:
    """Longest tandem array of the given units in either orientation.

    Ties are resolved towards the array with more exact windows, then towards
    the orientation in which the units were given.
    """
    units = tuple(u.upper() for u in units if u)
    if not units:
        return EMPTY
    length = len(units[0])
    if any(len(u) != length for u in units) or length == 0 or len(sequence) < length:
        return EMPTY
    seq = sequence.upper()
    max_mismatch = effective_mismatch(length, max_mismatch)
    fwd = set().union(*(rotations(u) for u in units))
    rev = {str(Seq(u).reverse_complement()) for u in fwd}
    best = _longest_run(seq, fwd, length, max_mismatch, "+")
    other = _longest_run(seq, rev, length, max_mismatch, "-")
    if other.copies > best.copies or (other.copies == best.copies and other.exact > best.exact):
        best = other
    return best
