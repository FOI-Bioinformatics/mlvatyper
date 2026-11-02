"""Allele calling from amplicon lengths.

The repeat number is anchored to a reference: raw = ref_repeats +
(length - ref_size) / unit_length. Values within a tolerance of an integer
are reported as that integer. Other values are rounded according to the
chosen policy and flagged as INTERMEDIATE; the policy "intermediate" reports
them as n.5 in the manner of MLVA_finder. An optional bin table maps sizes or
size ranges directly to alleles and takes precedence over the formula.

Every call carries a list of status codes. A call without problems has the
single code OK. Codes are text so that they cannot be confused with allele
numbers.
"""

import math
from collections import defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path

import pandas as pd

from .log import get_logger
from .models import Amplicon, MarkerDef
from .tandem import EMPTY, TandemArray, count_tandem_array

log = get_logger()

ROUNDING_POLICIES = ("nearest", "up", "down", "intermediate")

OK = "OK"
NOT_DETECTED = "NOT_DETECTED"
NEAR_CONTIG_END = "NEAR_CONTIG_END"
NOT_SEARCHED = "NOT_SEARCHED"
ONE_PRIMER_ONLY = "ONE_PRIMER_ONLY"
PRIMER_INDEL = "PRIMER_INDEL"
MISSING_CODES = {NEAR_CONTIG_END, NOT_SEARCHED, ONE_PRIMER_ONLY}
INTERMEDIATE = "INTERMEDIATE"
NEGATIVE = "NEGATIVE"
OFF_TABLE = "OFF_TABLE"
BINNED = "BINNED"
MULTIPLE = "MULTIPLE"
AMBIGUOUS = "AMBIGUOUS"
CENSORED = "CENSORED"
HAS_N = "HAS_N"
COUNT_DISAGREES = "COUNT_DISAGREES"
NO_DEFINITION = "NO_DEFINITION"


@dataclass
class Bin:
    marker: str
    low: int
    high: int
    allele: str


@dataclass
class Call:
    """Allele call for one amplicon (or for an absent marker when amplicon is None)."""

    sample: str
    marker: str
    amplicon: Amplicon | None
    raw: float | None = None
    allele: str | None = None
    residual: float | None = None
    array: TandemArray = EMPTY
    status: list[str] = field(default_factory=list)

    @property
    def length(self) -> int | None:
        return None if self.amplicon is None else len(self.amplicon.sequence)

    @property
    def status_text(self) -> str:
        return ";".join(self.status) if self.status else OK

    @property
    def cell(self) -> str:
        """Value for the sample by marker allele table."""
        if self.allele is not None:
            return self.allele
        return self.status[0] if self.status else OK


# hyphen, non-breaking hyphen, figure dash, en dash, em dash and minus sign, as pasted from publications
DASHES = str.maketrans({c: "-" for c in "\u2010\u2011\u2012\u2013\u2014\u2212"})


def parse_size_range(text: str, decimals: bool = False) -> tuple[int, int]:
    """(low, high) from a size or a low-high range; any dash character separates the bounds.

    With decimals, sizes are rounded to whole base pairs, halves up.
    Raises ValueError for anything else, including a range whose high
    bound is below its low bound.
    """
    text = str(text).translate(DASHES).replace(" ", "")
    parts = text.split("-") if "-" in text else [text, text]
    if len(parts) != 2:
        raise ValueError(text)
    low, high = (math.floor(float(p) + 0.5) if decimals else int(p) for p in parts)
    if high < low:
        raise ValueError(text)
    return low, high


def read_bins(path: str | Path) -> dict[str, list[Bin]]:
    """Bin file: marker, size or low-high range, allele. Tab separated, header optional.

    The header, when present, is the first line that is not a comment. A
    line without an allele is skipped with a warning, as is a size that is
    not understood.
    """
    bins: dict[str, list[Bin]] = defaultdict(list)
    first = True
    with open(path, encoding="utf-8-sig") as handle:  # -sig: drop a byte-order mark
        for number, line in enumerate(handle, start=1):
            parts = [p.strip() for p in line.rstrip("\r\n").split("\t")]
            if not parts[0] or parts[0].startswith("#"):
                continue
            header, first = first, False
            if len(parts) < 3 or not parts[2]:
                if not header:
                    log.warning("bin file line %d: no allele, skipped", number)
                continue
            try:
                low, high = parse_size_range(parts[1])
            except ValueError:
                if not header:
                    log.warning("bin file line %d: size %r not understood, skipped", number, parts[1])
                continue
            bins[parts[0]].append(Bin(parts[0], low, high, parts[2]))
    return dict(bins)


def round_repeats(raw: float, tolerance: float, policy: str) -> tuple[str, bool]:
    """Return (allele text, is_intermediate).

    Halves round up (2.5 to 3), not to the even neighbour as Python's round
    does. The intermediate value is floor(raw) + 0.5, so -0.4 gives -0.5.
    """
    nearest = math.floor(raw + 0.5)
    # strict comparison as in MLVA_finder: a value exactly tolerance away is intermediate
    if abs(raw - nearest) < tolerance:
        return str(nearest), False
    if policy == "nearest":
        return str(nearest), True
    if policy == "up":
        return str(math.ceil(raw)), True
    if policy == "down":
        return str(math.floor(raw)), True
    if policy == "intermediate":
        return str(math.floor(raw) + 0.5), True
    raise ValueError(f"unknown rounding policy {policy}")


def _allele_value(text: str) -> float | None:
    """Numeric value of an allele, or None for a non-numeric bin label such as '3a'."""
    try:
        return float(text.lstrip(">="))
    except ValueError:
        return None


def _assign_allele(call: Call, d: MarkerDef, length: int, bins: list[Bin] | None,
                   tolerance: float, policy: str) -> None:
    call.raw = d.ref_repeats + (length - d.ref_size) / d.unit_length
    binned = None
    if bins:
        binned = next((b.allele for b in bins if b.low <= length <= b.high), None)
        if binned is None:
            call.status.append(OFF_TABLE)
    if binned is not None:
        call.allele = binned
        call.status.append(BINNED)
    else:
        call.allele, intermediate = round_repeats(call.raw, tolerance, policy)
        if intermediate:
            call.status.append(INTERMEDIATE)
    if call.raw < 0:
        call.status.append(NEGATIVE)
    value = _allele_value(call.allele)
    expected = None if value is None else d.expected_size(value)
    call.residual = None if expected is None else length - expected


def _check_array(call: Call, d: MarkerDef, max_array_length: int | None) -> None:
    allele = _allele_value(call.allele)
    if d.has_repeats and call.array.copies and allele is not None and abs(call.array.copies - allele) >= 1:
        call.status.append(COUNT_DISAGREES)
    if call.array.copies:
        array_length = call.array.span
    elif allele is not None:
        array_length = allele * d.unit_length
    else:
        array_length = None
    if max_array_length is not None and array_length is not None and array_length >= max_array_length:
        call.allele = f">={call.allele}"
        call.status.append(CENSORED)


def call_amplicon(
    amp: Amplicon, d: MarkerDef, tolerance: float, policy: str,
    bins: list[Bin] | None = None, max_mismatch: int = 1, max_array_length: int | None = None,
) -> Call:
    call = Call(amp.hit.sample, amp.hit.primer_id, amp)
    if amp.has_n:
        call.status.append(HAS_N)
    if amp.hit.n_alternatives:
        call.status.append(AMBIGUOUS)
    if amp.hit.fwd_indels or amp.hit.rev_indels:
        call.status.append(PRIMER_INDEL)
    if d.has_repeats:
        call.array = count_tandem_array(amp.sequence, d.patterns, max_mismatch)
    if not d.can_call:
        call.status.insert(0, NO_DEFINITION)  # first, so that the allele table shows it rather than e.g. HAS_N
        return call
    _assign_allele(call, d, len(amp.sequence), bins, tolerance, policy)
    _check_array(call, d, max_array_length)
    return call


def call_all(
    amplicons: list[Amplicon], samples: list[str], defs: list[MarkerDef], tolerance: float, policy: str,
    bins: dict[str, list[Bin]] | None = None, max_mismatch: int = 1, max_array_length: int | None = None,
    near_contig_end: set[tuple[str, str]] | None = None,
    missing_reasons: dict[tuple[str, str], str] | None = None,
) -> list[Call]:
    """One Call per amplicon, plus one Call per sample and marker without amplicon.

    missing_reasons maps (sample, marker) to one of MISSING_CODES and takes
    precedence over the near_contig_end set, which is kept for older runs.
    """
    by_marker = {d.primer_id: d for d in defs}
    present: dict[tuple[str, str], list[Call]] = defaultdict(list)
    calls: list[Call] = []
    for amp in sorted(amplicons, key=lambda a: (a.hit.sample, a.hit.primer_id, a.copy)):
        d = by_marker.get(amp.hit.primer_id, MarkerDef(amp.hit.primer_id))
        c = call_amplicon(amp, d, tolerance, policy, (bins or {}).get(d.primer_id), max_mismatch, max_array_length)
        present[(amp.hit.sample, amp.hit.primer_id)].append(c)
        calls.append(c)
    for group in present.values():
        if len(group) > 1:
            for c in group:
                c.status.append(MULTIPLE)
    for sample in samples:
        for d in defs:
            if (sample, d.primer_id) not in present:
                code = (missing_reasons or {}).get((sample, d.primer_id))
                if code not in MISSING_CODES:
                    near = near_contig_end and (sample, d.primer_id) in near_contig_end
                    code = NEAR_CONTIG_END if near else NOT_DETECTED
                calls.append(Call(sample, d.primer_id, None, status=[code]))
    return calls


def _frame(cells: dict[tuple[str, str], str], samples: list[str], markers: list[str]) -> pd.DataFrame:
    data: dict[str, dict[str, str]] = defaultdict(dict)
    for (sample, marker), value in cells.items():
        data[sample][marker] = value
    df = pd.DataFrame.from_dict(data, orient="index")
    return df.reindex(index=list(samples), columns=markers)


def allele_table(calls: list[Call], samples: list[str], markers: list[str]) -> pd.DataFrame:
    """Sample by marker allele values; multiple amplicons joined with '/'."""
    grouped: dict[tuple[str, str], list[str]] = defaultdict(list)
    for c in calls:
        grouped[(c.sample, c.marker)].append(c.cell)
    return _frame({k: "/".join(v) for k, v in grouped.items()}, samples, markers)


def status_table(calls: list[Call], samples: list[str], markers: list[str]) -> pd.DataFrame:
    grouped: dict[tuple[str, str], list[str]] = defaultdict(list)
    for c in calls:
        for s in c.status or [OK]:
            if s not in grouped[(c.sample, c.marker)]:
                grouped[(c.sample, c.marker)].append(s)
    return _frame({k: ";".join(v) for k, v in grouped.items()}, samples, markers)


def tandem_count_table(calls: list[Call], samples: list[str], markers: list[str]) -> pd.DataFrame:
    grouped: dict[tuple[str, str], list[str]] = defaultdict(list)
    for c in calls:
        if c.amplicon is not None:
            grouped[(c.sample, c.marker)].append(str(c.array.copies))
    return _frame({k: ";".join(v) for k, v in grouped.items()}, samples, markers)


def profile_table(calls: list[Call], samples: list[str], markers: list[str]) -> pd.DataFrame:
    """Allele per marker in the given order plus a hyphen joined profile string."""
    df = allele_table(calls, samples, markers).fillna(NOT_DETECTED)
    df.insert(0, "profile", df[markers].astype(str).agg("-".join, axis=1) if markers else "")
    return df


LONG_COLUMNS = [
    "sample", "marker", "copy", "contig", "start", "end", "strand", "length", "expected_size",
    "fwd_mismatches", "rev_mismatches", "n_alternatives", "dist_to_contig_end", "fwd_indels", "rev_indels",
    "search_pass", "has_N", "raw_repeats", "allele", "residual_bp", "tandem_copies", "tandem_exact", "tandem_variant",
    "tandem_orientation", "status",
]
INT_COLUMNS = [
    "copy", "start", "end", "length", "expected_size", "fwd_mismatches", "rev_mismatches", "n_alternatives",
    "dist_to_contig_end", "fwd_indels", "rev_indels", "has_N", "tandem_copies", "tandem_exact", "tandem_variant",
]


def long_table(calls: list[Call]) -> pd.DataFrame:
    rows = []
    for c in calls:
        h = c.amplicon.hit if c.amplicon is not None else None
        rows.append({
            "sample": c.sample, "marker": c.marker,
            "copy": c.amplicon.copy if c.amplicon else None,
            "contig": h.contig if h else None, "start": h.start if h else None, "end": h.end if h else None,
            "strand": h.strand if h else None, "length": c.length,
            "expected_size": h.expected_size if h and h.expected_size >= 0 else None,
            "fwd_mismatches": h.fwd_mismatches if h and h.fwd_mismatches >= 0 else None,
            "rev_mismatches": h.rev_mismatches if h and h.rev_mismatches >= 0 else None,
            "n_alternatives": h.n_alternatives if h else None,
            "dist_to_contig_end": h.dist_to_contig_end if h and h.dist_to_contig_end >= 0 else None,
            "fwd_indels": h.fwd_indels if h else None, "rev_indels": h.rev_indels if h else None,
            "search_pass": h.search_pass if h else None,
            "has_N": int(c.amplicon.has_n) if c.amplicon else None,
            "raw_repeats": None if c.raw is None else round(c.raw, 3),
            "allele": c.allele,
            "residual_bp": None if c.residual is None else round(c.residual, 1),
            "tandem_copies": c.array.copies if c.amplicon else None,
            "tandem_exact": c.array.exact if c.amplicon else None,
            "tandem_variant": c.array.variant if c.amplicon else None,
            "tandem_orientation": c.array.orientation if c.amplicon and c.array.copies else None,
            "status": c.status_text,
        })
    df = pd.DataFrame(rows, columns=LONG_COLUMNS)
    for column in INT_COLUMNS:
        df[column] = df[column].astype("Int64")
    return df


def derive_definitions(
    amplicons: list[Amplicon], defs: list[MarkerDef], reference_sample: str, max_mismatch: int = 1,
) -> list[MarkerDef]:
    """Fill ref_size and ref_repeats from the amplicons of one sample.

    ref_size is the amplicon length in that sample; ref_repeats is the tandem
    array count when patterns are known, otherwise 0, in which case alleles
    are offsets from the reference rather than repeat numbers. Markers with
    several amplicons in the reference sample use the first copy and are
    logged. Existing values in defs are kept.
    """
    by_marker: dict[str, list[Amplicon]] = defaultdict(list)
    for a in amplicons:
        if a.hit.sample == reference_sample:
            by_marker[a.hit.primer_id].append(a)
    if not by_marker:
        raise ValueError(f"reference sample {reference_sample!r} has no amplicons")
    derived = []
    for d in defs:
        amps = by_marker.get(d.primer_id)
        if d.can_call or not amps or d.unit_length == 0:
            derived.append(d)
            continue
        amps.sort(key=lambda a: a.copy)
        if len(amps) > 1:
            log.warning("%s: %d amplicons in reference sample %s, copy 1 used",
                        d.primer_id, len(amps), reference_sample)
        amp = amps[0]
        ref_size = d.ref_size if d.ref_size is not None else len(amp.sequence)
        if d.ref_repeats is not None:
            ref_repeats = d.ref_repeats
        elif d.has_repeats:
            ref_repeats = count_tandem_array(amp.sequence, d.patterns, max_mismatch).copies
        else:
            ref_repeats = 0
        derived.append(replace(d, ref_size=ref_size, ref_repeats=ref_repeats))  # keeps offset and label
    n = sum(1 for d in derived if d.can_call)
    log.info("Reference sample %s: %d of %d markers have reference size and repeats", reference_sample, n, len(defs))
    return derived
