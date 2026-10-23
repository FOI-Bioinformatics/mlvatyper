"""Reader for marker definition files and MLVAbank style marker names.

The definition file is tab separated with a header line. Recognised columns
(case insensitive): marker, patterns (or pattern), unit_length (or length),
ref_size, ref_repeats, offset. Only marker is required. patterns holds one or
more repeat unit sequences separated by semicolons, or N/A. Other columns are
ignored, which keeps the older pattern file format (marker, pattern, length)
readable.

Marker ids written in the MLVAbank convention, NAME_<unit>bp_<size>bp_<n>u,
carry the unit length, the amplicon size in the reference strain and the
repeat number in that strain. They are parsed into the same definition.
"""

import re
from pathlib import Path

from .log import get_logger
from .models import MarkerDef

log = get_logger()

NAME_RE = re.compile(r"^(?P<name>.+?)_(?P<unit>\d+)bp_(?P<size>\d+)bp_(?P<n>\d+)[uU]$")

COLUMN_ALIASES = {
    "marker": "marker", "id": "marker", "primer_id": "marker", "locus": "marker",
    "pattern": "patterns", "patterns": "patterns", "unit": "patterns",
    "length": "unit_length", "unit_length": "unit_length", "pattern_length": "unit_length",
    "ref_size": "ref_size", "ref_repeats": "ref_repeats", "offset": "offset",
    "label": "label", "name": "label", "mlvabank": "label", "locus_name": "label",
}


class MarkerError(ValueError):
    """Raised for unusable marker definitions."""


def parse_marker_name(marker_id: str) -> MarkerDef | None:
    """Return a definition from an MLVAbank style id, or None if the id has no such suffix."""
    m = NAME_RE.match(marker_id)
    if not m:
        return None
    return MarkerDef(marker_id, (), int(m["unit"]), int(m["size"]), int(m["n"]))


def _parse_int(value: str, field: str, marker: str, line_number: int) -> int | None:
    value = value.strip()
    if not value or value.upper() in {"NA", "N/A"}:
        return None
    try:
        return int(value)
    except ValueError as exc:
        raise MarkerError(f"line {line_number} ({marker}): {field} {value!r} is not an integer") from exc


def _parse_patterns(value: str) -> tuple[str, ...]:
    value = value.strip()
    if not value or value.upper() == "N/A":
        return ()
    return tuple(p.strip().upper() for p in value.split(";") if p.strip())


def parse_header(line: str) -> dict[str, int]:
    """Map canonical column names to column indices."""
    columns: dict[str, int] = {}
    for index, raw in enumerate(line.rstrip("\r\n").split("\t")):
        canonical = COLUMN_ALIASES.get(raw.strip().lower().lstrip("#"))
        if canonical and canonical not in columns:
            columns[canonical] = index
    if "marker" not in columns:
        raise MarkerError("marker definition file needs a header with a 'marker' column")
    return columns


def parse_marker_line(line: str, columns: dict[str, int], line_number: int) -> MarkerDef | None:
    fields = line.rstrip("\r\n").split("\t")

    def get(name: str) -> str:
        index = columns.get(name)
        return fields[index] if index is not None and index < len(fields) else ""

    marker = get("marker").strip()
    if not marker:
        log.warning("marker file line %d: empty marker id, skipped", line_number)
        return None
    patterns = _parse_patterns(get("patterns"))
    unit_length = _parse_int(get("unit_length"), "unit_length", marker, line_number)
    if unit_length is None:
        unit_length = len(patterns[0]) if patterns else 0
    ref_size = _parse_int(get("ref_size"), "ref_size", marker, line_number)
    ref_repeats = _parse_int(get("ref_repeats"), "ref_repeats", marker, line_number)
    offset = _parse_int(get("offset"), "offset", marker, line_number)
    from_name = parse_marker_name(marker)
    if from_name is not None:
        unit_length = unit_length or from_name.unit_length
        ref_size = from_name.ref_size if ref_size is None else ref_size
        ref_repeats = from_name.ref_repeats if ref_repeats is None else ref_repeats
    if unit_length < 0:
        raise MarkerError(f"line {line_number} ({marker}): unit_length must not be negative")
    bad = [p for p in patterns if len(p) != unit_length]
    if bad:
        raise MarkerError(
            f"line {line_number} ({marker}): patterns {bad} do not have unit_length {unit_length}"
        )
    if (ref_size is not None and ref_repeats is not None and offset is None and unit_length > 0
            and ref_size - ref_repeats * unit_length < 0):
        raise MarkerError(f"line {line_number} ({marker}): ref_repeats * unit_length exceeds ref_size")
    return MarkerDef(marker, patterns, unit_length, ref_size, ref_repeats, offset, get("label").strip() or None)


LEGACY_COLUMNS = {"marker": 0, "patterns": 1}  # the length column is checked, later columns are ignored
SEQUENCE_RE = re.compile(r"^(N/A|[ACGTRYSWKMBDHVN;]+)$")  # upper case only, so that a header word is not a unit


def looks_like_data(line: str) -> bool:
    """A headerless legacy pattern line: an id followed by upper case unit sequences or N/A."""
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 2 or not fields[0].strip() or COLUMN_ALIASES.get(fields[0].strip().lower().lstrip("#")):
        return False
    return bool(SEQUENCE_RE.match(fields[1].strip()))


def read_markers(path: str | Path) -> list[MarkerDef]:
    """Read a marker definition file.

    A file without a header is read in the legacy layout (marker, patterns,
    unit length, then ignored columns), with its first line as data; the
    legacy script skipped that line. The unit length is then taken from the
    patterns, and a differing stated length is only reported.
    """
    results: list[MarkerDef] = []
    with open(path, encoding="utf-8") as handle:
        header = handle.readline()
        legacy = looks_like_data(header)
        if legacy:
            log.warning("%s has no header; read as marker, patterns, unit length", path)
            columns = dict(LEGACY_COLUMNS)
            lines = [(1, header)]
        else:
            columns = parse_header(header)
            lines = []
        lines.extend(enumerate(handle, start=2))
        for number, line in lines:
            if not line.strip() or line.startswith("#"):
                continue
            parsed = parse_marker_line(line, columns, number)
            if parsed is not None:
                if legacy:
                    _check_legacy_length(line, parsed, number)
                results.append(parsed)
    log.info("Loaded %d marker definitions from %s", len(results), path)
    return results


def _check_legacy_length(line: str, d: MarkerDef, number: int) -> None:
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) > 2 and fields[2].strip().isdigit() and d.has_repeats and int(fields[2]) != d.unit_length:
        log.warning("line %d (%s): stated unit length %s differs from the pattern length %d; the pattern is used",
                    number, d.primer_id, fields[2].strip(), d.unit_length)


def markers_from_ids(marker_ids: list[str]) -> list[MarkerDef]:
    """Definitions parsed from MLVAbank style ids; ids without the suffix are skipped."""
    defs = [d for d in (parse_marker_name(i) for i in marker_ids) if d is not None]
    if defs:
        log.info("%d marker ids carry unit, reference size and repeats in their name", len(defs))
    return defs


def merge_definitions(primary: list[MarkerDef], secondary: list[MarkerDef]) -> list[MarkerDef]:
    """Combine definitions; fields missing in primary are taken from secondary for the same marker."""
    by_id = {d.primer_id: d for d in secondary}
    merged = []
    seen = set()
    for d in primary:
        s = by_id.get(d.primer_id)
        if s is not None:
            d = MarkerDef(
                d.primer_id, d.patterns or s.patterns, d.unit_length or s.unit_length,
                s.ref_size if d.ref_size is None else d.ref_size,
                s.ref_repeats if d.ref_repeats is None else d.ref_repeats,
                s.offset if d.offset is None else d.offset,
                d.label or s.label,
            )
        merged.append(d)
        seen.add(d.primer_id)
    merged.extend(d for d in secondary if d.primer_id not in seen)
    return merged


def check_against_sts(defs: list[MarkerDef], sts_ids: set[str]) -> list[MarkerDef]:
    """Warn about definitions for markers absent from the STS file and drop them."""
    kept = []
    for d in defs:
        if d.primer_id in sts_ids:
            kept.append(d)
        else:
            log.warning("marker %s is defined but not in the STS file; ignored", d.primer_id)
    return kept


def write_markers(defs: list[MarkerDef], path: str | Path) -> None:
    labels = any(d.label for d in defs)
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("marker\tunit_length\tref_size\tref_repeats\toffset\tpatterns" + ("\tlabel" if labels else "")
                     + "\n")
        for d in defs:
            handle.write("\t".join([
                d.primer_id, str(d.unit_length),
                "" if d.ref_size is None else str(d.ref_size),
                "" if d.ref_repeats is None else str(d.ref_repeats),
                "" if d.offset is None else str(d.offset),
                ";".join(d.patterns) if d.patterns else "N/A",
            ] + ([d.label or ""] if labels else [])) + "\n")
