"""Exchange with MLVA databases: profile export and bin tables from allele size tables.

MLVAbank style profile files are plain tables with one row per strain, a
key column and one column per locus holding the allele (repeat number).
Markers are written in scheme order under their label when one is defined
in the marker file, otherwise under the marker id. Cells without an
integer allele (status codes, intermediate values, several amplicons) are
written as given, so that they are visible rather than silently dropped;
`--mlvabank-missing` replaces them with a fixed token.

A bin table for `--bins` can be built from a long table that lists one
row per locus, allele and amplicon size (or size range), which is how
allele size tables are usually published.
"""

import csv
from pathlib import Path

import pandas as pd

from .log import get_logger
from .models import MarkerDef

log = get_logger()

LOCUS_ALIASES = ("locus", "marker", "vntr", "name")
ALLELE_ALIASES = ("allele", "repeats", "copies", "repeat_number", "code")
SIZE_ALIASES = ("size", "fragment", "fragment_size", "length", "amplicon", "bp", "range")


def export_profiles(alleles: pd.DataFrame, defs: list[MarkerDef], path: str | Path, key: str = "Strain",
                    missing: str | None = None, sep: str = ",") -> int:
    """Write the allele table with one column per marker in defs order; returns the row count.

    alleles has samples as index and marker ids as columns.
    """
    order = [d for d in defs if d.primer_id in alleles.columns]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=sep, lineterminator="\n")
        writer.writerow([key] + [d.display_name for d in order])
        for sample, row in alleles.iterrows():
            cells = []
            for d in order:
                value = "" if pd.isna(row[d.primer_id]) else str(row[d.primer_id])
                if missing is not None and not _is_integer(value):
                    value = missing
                cells.append(value)
            writer.writerow([sample] + cells)
    log.info("Wrote %s (%d strains, %d loci)", path, len(alleles), len(order))
    return len(alleles)


def _is_integer(value: str) -> bool:
    return value.lstrip("-").isdigit()


def _find_column(columns: list[str], aliases: tuple[str, ...], what: str) -> str:
    lowered = {c.lower().strip().replace(" ", "_"): c for c in columns}
    for alias in aliases:
        if alias in lowered:
            return lowered[alias]
    raise ValueError(f"no {what} column found; expected one of {', '.join(aliases)} in {columns}")


def bins_from_table(path: str | Path, out: str | Path, locus_column: str | None = None,
                    allele_column: str | None = None, size_column: str | None = None,
                    labels: dict[str, str] | None = None) -> int:
    """Convert a long allele size table to the bin file format; returns the number of bins.

    The table is tab or comma separated with a header. Sizes may be single
    numbers or low-high ranges; decimals are rounded to whole base pairs.
    labels maps a locus name in the table to a marker id when they differ.
    """
    table = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    cols = list(table.columns)
    locus = locus_column or _find_column(cols, LOCUS_ALIASES, "locus")
    allele = allele_column or _find_column(cols, ALLELE_ALIASES, "allele")
    size = size_column or _find_column(cols, SIZE_ALIASES, "size")
    names = {k.lower(): v for k, v in (labels or {}).items()}
    n = 0
    with open(out, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("marker\tsize\tallele\n")
        for _, row in table.iterrows():
            marker = names.get(str(row[locus]).strip().lower(), str(row[locus]).strip())
            text = str(row[size]).replace(" ", "")
            if not marker or not text:
                continue
            try:
                low, high = (round(float(x)) for x in text.split("-")) if "-" in text else (
                    round(float(text)), round(float(text)))
            except ValueError:
                log.warning("%s: size %r for %s not understood, skipped", path, row[size], marker)
                continue
            handle.write(f"{marker}\t{low}-{high}\t{str(row[allele]).strip()}\n")
            n += 1
    log.info("Wrote %d bins to %s", n, out)
    return n
