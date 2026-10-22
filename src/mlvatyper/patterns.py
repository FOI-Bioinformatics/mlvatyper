"""Backwards compatible entry points for the older pattern file format.

The pattern file (marker, pattern, length) is a subset of the marker
definition file; see markers.py.
"""

from pathlib import Path

from .markers import parse_header, parse_marker_line, read_markers
from .models import MarkerDef

LEGACY_COLUMNS = parse_header("marker\tpattern\tlength")


def parse_pattern_line(line: str, line_number: int) -> MarkerDef | None:
    """Parse one data line of a header-less three column pattern file."""
    if len(line.rstrip("\r\n").split("\t")) < 2:
        return None
    return parse_marker_line(line, LEGACY_COLUMNS, line_number)


def read_patterns(path: str | Path) -> list[MarkerDef]:
    return read_markers(path)
