"""Reading and validation of STS primer files.

The STS format is tab separated: id, forward primer, reverse primer, expected
product size, and an optional alias. The size may be a range low-high, which
merPCR turns into a per-marker size window. merPCR silently skips primers it
cannot hash, so validation is done here first to make such problems visible.

Primers may contain IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N).
Whether they are honoured during the search is controlled by the iupac option
of the search step.
"""

from pathlib import Path

from .io import open_text
from .log import get_logger
from .models import PrimerPair

log = get_logger()

IUPAC_BASES = set("ACGTRYSWKMBDHVN")
STRICT_BASES = set("ACGT")


class StsError(ValueError):
    """Raised when the STS file contains invalid entries."""


def _clean_primer(primer: str, label: str, primer_id: str, line_number: int, fix: bool) -> str:
    primer = primer.upper()
    bad = [(n + 1, c) for n, c in enumerate(primer) if c not in IUPAC_BASES]
    if bad and fix:
        for position, char in bad:
            log.warning("line %d (%s): removed %r at position %d of the %s primer", line_number, primer_id, char,
                        position, label)
        primer = "".join(c for c in primer if c in IUPAC_BASES)
    elif bad:
        where = ", ".join(f"{c!r} at position {p}" for p, c in bad)
        raise StsError(
            f"line {line_number} ({primer_id}): {label} primer {primer!r} contains non-IUPAC characters: {where}"
        )
    if not primer:
        raise StsError(f"line {line_number} ({primer_id}): empty {label} primer")
    return primer


def parse_size(size: str, primer_id: str, line_number: int) -> tuple[int, int | None, int | None]:
    """Expected size as (size, low, high); a range low-high gives the midpoint as in merPCR."""
    low = high = None
    try:
        if "-" in size:
            low_text, high_text = size.split("-", 1)
            low, high = int(low_text), int(high_text)
            if low < 0 or high < low:
                raise ValueError
            expected = (low + high) // 2
        else:
            expected = int(size)
            if expected < 0:
                raise ValueError
    except ValueError as exc:
        raise StsError(f"line {line_number} ({primer_id}): expected size {size!r} is not a non-negative "
                       "integer or a low-high range") from exc
    return expected, low, high


def parse_sts_line(line: str, line_number: int, fix: bool = False) -> PrimerPair:
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 4:
        raise StsError(f"line {line_number}: expected at least 4 tab-separated fields, got {len(fields)}")
    primer_id, fwd, rev, size = (f.strip() for f in fields[:4])
    alias = fields[4].strip() if len(fields) > 4 else None
    if not primer_id:
        raise StsError(f"line {line_number}: empty id")
    fwd = _clean_primer(fwd, "forward", primer_id, line_number, fix)
    rev = _clean_primer(rev, "reverse", primer_id, line_number, fix)
    expected, low, high = parse_size(size, primer_id, line_number)
    return PrimerPair(primer_id, fwd, rev, expected, alias, line_number, low, high)


def read_sts(path: str | Path, fix: bool = False) -> tuple[list[PrimerPair], list[tuple[int, str, str]]]:
    """Return (valid primer pairs, invalid entries as (line_number, raw_line, reason)).

    With fix, characters that are not IUPAC codes are removed from the primers
    and each removal is logged.
    """
    valid: list[PrimerPair] = []
    invalid: list[tuple[int, str, str]] = []
    seen: set[str] = set()
    with open_text(path) as handle:
        for number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            try:
                pair = parse_sts_line(raw, number, fix)
            except StsError as exc:
                invalid.append((number, raw.rstrip("\r\n"), str(exc)))
                continue
            if pair.id in seen:
                invalid.append((number, raw.rstrip("\r\n"), f"line {number}: duplicate id {pair.id}"))
                continue
            seen.add(pair.id)
            valid.append(pair)
    return valid, invalid


def is_degenerate(pair: PrimerPair) -> bool:
    return bool((set(pair.forward) | set(pair.reverse)) - STRICT_BASES)


def has_hash_word(primer: str, wordsize: int) -> bool:
    """Whether merPCR can hash the primer: it needs wordsize consecutive unambiguous bases."""
    run = 0
    for base in primer.upper():
        run = run + 1 if base in STRICT_BASES else 0
        if run >= wordsize:
            return True
    return False


def write_sts(pairs: list[PrimerPair], path: str | Path) -> None:
    """Write primer pairs in STS format with LF line endings."""
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for pair in pairs:
            fields = [pair.id, pair.forward, pair.reverse, pair.size_text]
            if pair.alias is not None:
                fields.append(pair.alias)
            handle.write("\t".join(fields) + "\n")


def write_invalid(invalid: list[tuple[int, str, str]], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("line\treason\traw\n")
        for number, raw, reason in invalid:
            handle.write(f"{number}\t{reason}\t{raw}\n")
