"""Reading and validation of STS primer files.

The STS format is tab separated: id, forward primer, reverse primer, expected
product size, and an optional alias. merPCR silently skips primers it cannot
hash, so validation is done here first to make such problems visible.

Primers may contain IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N).
Whether they are honoured during the search is controlled by the iupac option
of the search step.
"""

from pathlib import Path

from .io import open_text
from .models import PrimerPair

IUPAC_BASES = set("ACGTRYSWKMBDHVN")
STRICT_BASES = set("ACGT")


class StsError(ValueError):
    """Raised when the STS file contains invalid entries."""


def parse_sts_line(line: str, line_number: int) -> PrimerPair:
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 4:
        raise StsError(f"line {line_number}: expected at least 4 tab-separated fields, got {len(fields)}")
    primer_id, fwd, rev, size = (f.strip() for f in fields[:4])
    alias = fields[4].strip() if len(fields) > 4 else None
    if not primer_id:
        raise StsError(f"line {line_number}: empty id")
    for label, primer in (("forward", fwd), ("reverse", rev)):
        bad = set(primer.upper()) - IUPAC_BASES
        if not primer or bad:
            raise StsError(
                f"line {line_number} ({primer_id}): {label} primer {primer!r} "
                f"contains non-IUPAC characters {sorted(bad)}"
            )
    try:
        expected = int(size)
    except ValueError as exc:
        raise StsError(f"line {line_number} ({primer_id}): expected size {size!r} is not an integer") from exc
    if expected < 0:
        raise StsError(f"line {line_number} ({primer_id}): expected size must not be negative")
    return PrimerPair(primer_id, fwd.upper(), rev.upper(), expected, alias, line_number)


def read_sts(path: str | Path) -> tuple[list[PrimerPair], list[tuple[int, str, str]]]:
    """Return (valid primer pairs, invalid entries as (line_number, raw_line, reason))."""
    valid: list[PrimerPair] = []
    invalid: list[tuple[int, str, str]] = []
    seen: set[str] = set()
    with open_text(path) as handle:
        for number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            try:
                pair = parse_sts_line(raw, number)
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


def write_sts(pairs: list[PrimerPair], path: str | Path) -> None:
    """Write primer pairs in STS format with LF line endings."""
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        for pair in pairs:
            fields = [pair.id, pair.forward, pair.reverse, str(pair.expected_size)]
            if pair.alias is not None:
                fields.append(pair.alias)
            handle.write("\t".join(fields) + "\n")


def write_invalid(invalid: list[tuple[int, str, str]], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("line\treason\traw\n")
        for number, raw, reason in invalid:
            handle.write(f"{number}\t{reason}\t{raw}\n")
