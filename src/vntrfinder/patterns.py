"""Reader for the repeat pattern file.

Tab separated with a header line. Columns: marker id, pattern(s) separated by
semicolons or the literal N/A, and an optional repeat unit length. Further
columns present in older files are ignored.
"""

from pathlib import Path

from .log import get_logger
from .models import PatternDef

log = get_logger()


def parse_pattern_line(line: str, line_number: int) -> PatternDef | None:
    """Parse one data line. Returns None for lines with fewer than two columns."""
    parts = line.rstrip("\r\n").split("\t")
    if len(parts) < 2:
        log.warning("pattern file line %d: fewer than 2 columns, skipped", line_number)
        return None
    primer_id = parts[0].strip()
    pattern_str = parts[1].strip()
    if not pattern_str or pattern_str.upper() == "N/A":
        return PatternDef(primer_id)
    patterns = tuple(p.strip().upper() for p in pattern_str.split(";") if p.strip())
    length = len(patterns[0])
    if len(parts) > 2 and parts[2].strip():
        try:
            length = int(parts[2])
        except ValueError:
            log.warning(
                "pattern file line %d (%s): pattern length %r not an integer, using %d",
                line_number, primer_id, parts[2], length,
            )
    return PatternDef(primer_id, patterns, length)


def read_patterns(path: str | Path) -> list[PatternDef]:
    results: list[PatternDef] = []
    with open(path, encoding="utf-8") as handle:
        handle.readline()  # header
        for number, line in enumerate(handle, start=2):
            if not line.strip() or line.startswith("#"):
                continue
            parsed = parse_pattern_line(line, number)
            if parsed is not None:
                results.append(parsed)
    log.info("Loaded %d pattern definitions from %s", len(results), path)
    return results
