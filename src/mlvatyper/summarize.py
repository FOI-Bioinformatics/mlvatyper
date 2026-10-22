"""Summary tables built from the hits table.

All tables have one row per sample and one column per marker. Samples without
a hit for a marker are reported as NA. When a sample has several amplicons
for a marker, the per-copy values are joined with "/" (lengths) or ";"
(counts and length differences), in copy order.

The pattern count and length difference tables are the legacy estimators of
the original vntr_searcher script. They are kept for comparison under the
--legacy-tables option; allele calling lives in calling.py.
"""

from collections import defaultdict
from pathlib import Path

import pandas as pd

from .log import get_logger
from .models import Amplicon, MarkerDef

log = get_logger()


def _frame(data: dict[str, dict[str, str | None]], samples: list[str], markers: list[str]) -> pd.DataFrame:
    df = pd.DataFrame.from_dict(data, orient="index")
    df = df.reindex(index=sorted(samples), columns=sorted(markers))
    return df


def _select(amplicons: list[Amplicon], clean: bool) -> list[Amplicon]:
    return [a for a in amplicons if not (clean and a.has_n)]


def _by_marker_sample(amplicons: list[Amplicon]) -> dict[str, dict[str, list[Amplicon]]]:
    out: dict[str, dict[str, list[Amplicon]]] = defaultdict(lambda: defaultdict(list))
    for a in sorted(amplicons, key=lambda x: (x.hit.primer_id, x.hit.sample, x.copy)):
        out[a.hit.primer_id][a.hit.sample].append(a)
    return out


def lengths_table(amplicons: list[Amplicon], samples: list[str], clean: bool = False) -> pd.DataFrame:
    """Amplicon length per sample and marker; multiple copies joined with '/'."""
    data: dict[str, dict[str, str | None]] = defaultdict(dict)
    grouped = _by_marker_sample(_select(amplicons, clean))
    for marker, per_sample in grouped.items():
        for sample, amps in per_sample.items():
            data[sample][marker] = "/".join(str(len(a.sequence)) for a in amps)
    markers = sorted({a.hit.primer_id for a in amplicons})
    return _frame(data, samples, markers)


def count_repeats(sequence: str, patterns: tuple[str, ...]) -> int:
    """Legacy count: sum of exact, non-overlapping occurrences of each pattern on the amplicon strand."""
    seq = sequence.upper()
    total = 0
    for pattern in patterns:
        pat = pattern.upper()
        pos = 0
        while True:
            pos = seq.find(pat, pos)
            if pos == -1:
                break
            total += 1
            pos += len(pat)
    return total


def pattern_counts_table(
    amplicons: list[Amplicon], samples: list[str], defs: list[MarkerDef], clean: bool = False
) -> pd.DataFrame:
    """Legacy exact pattern counts per copy, joined with ';'."""
    valid = {d.primer_id: d for d in defs if d.has_repeats}
    data: dict[str, dict[str, str | None]] = defaultdict(dict)
    grouped = _by_marker_sample(_select(amplicons, clean))
    for marker, d in valid.items():
        per_sample = grouped.get(marker)
        if not per_sample:
            log.debug("%s: no amplicons for pattern counting", marker)
            continue
        for sample, amps in per_sample.items():
            data[sample][marker] = ";".join(str(count_repeats(a.sequence, d.patterns)) for a in amps)
    return _frame(data, samples, list(valid))


def length_diff_from_shortest_table(
    amplicons: list[Amplicon], samples: list[str], defs: list[MarkerDef], clean: bool = False
) -> pd.DataFrame:
    """Legacy (length - shortest length for the marker) / repeat unit length, per copy, joined with ';'.

    The reference is the shortest amplicon in the current dataset, so the
    values are relative to the run and not comparable between runs.
    """
    valid = {d.primer_id: d for d in defs if d.has_repeats and d.pattern_length > 0}
    data: dict[str, dict[str, str | None]] = defaultdict(dict)
    grouped = _by_marker_sample(_select(amplicons, clean))
    for marker, d in valid.items():
        per_sample = grouped.get(marker)
        if not per_sample:
            continue
        shortest = min(len(a.sequence) for amps in per_sample.values() for a in amps)
        for sample, amps in per_sample.items():
            data[sample][marker] = ";".join(
                f"{(len(a.sequence) - shortest) / d.pattern_length:.3f}" for a in amps
            )
    return _frame(data, samples, list(valid))


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    df.to_csv(path, sep="\t", na_rep="NA", lineterminator="\n")
    log.info("Wrote %s (%d samples, %d markers)", path, len(df.index), len(df.columns))


def write_all_summaries(amplicons: list[Amplicon], samples: list[str], outdir: Path, prefix: str,
                        defs: list[MarkerDef] | None, clwclean: bool, legacy_tables: bool = False) -> list[Path]:
    """Lengths tables, plus the legacy count tables when requested."""
    written = []
    variants = [("", False)] + ([("_clean", True)] if clwclean else [])
    for suffix, clean in variants:
        p = outdir / f"{prefix}_lengths_summary{suffix}.tsv"
        write_table(lengths_table(amplicons, samples, clean), p)
        written.append(p)
        if defs is not None and legacy_tables:
            p = outdir / f"{prefix}_count_summary{suffix}.tsv"
            write_table(pattern_counts_table(amplicons, samples, defs, clean), p)
            written.append(p)
            p = outdir / f"{prefix}_count_lengthdiff{suffix}.tsv"
            write_table(length_diff_from_shortest_table(amplicons, samples, defs, clean), p)
            written.append(p)
    return written
