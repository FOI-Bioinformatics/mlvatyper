"""Consistency check of a scheme, optionally against its reference assembly.

The scheme must state its provenance: citation, url and reference. A
scheme whose allele numbers cannot be traced to a source is not usable
for comparison with a database, so a missing or placeholder value is an
error rather than a warning.

Without an assembly the check covers the files: STS entries are validated,
marker and bin ids are compared with the STS ids, and primers that merPCR
cannot hash at the default and at the rescue word size are listed. With an
assembly (a FASTA path, or an accession fetched from NCBI) every marker is
searched and the number of products, the product size against ref_size
and the tandem array count against ref_repeats are reported per marker.
"""

from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from .calling import read_bins
from .extract import load_fasta
from .log import get_logger
from .markers import read_markers
from .models import MarkerDef, PrimerPair
from .rescue import RescueOptions, find_with_indels, rescue_pass
from .scheme import PROVENANCE_FIELDS, Scheme
from .search import MerPcrSearcher
from .sts import has_hash_word, read_sts, write_sts
from .tandem import count_tandem_array

log = get_logger()

CHECK_COLUMNS = [
    "marker", "sts_ok", "in_markers", "in_bins", "hash_word_11", "hash_word_6", "n_products", "length",
    "ref_size", "size_ok", "tandem_copies", "ref_repeats", "repeats_ok", "search_pass", "note",
]


def provenance_checks(scheme: Scheme, strict: bool = True) -> list[str]:
    """Log the provenance of the scheme; returns the fields that are missing.

    Missing fields are errors unless the caller accepts them, in which case
    they are warnings like the recommended fields.
    """
    missing = scheme.missing_provenance()
    for field in PROVENANCE_FIELDS:
        if field not in missing:
            log.info("%s: %s: %s", scheme.name, field, scheme.meta[field])
    report = log.error if strict else log.warning
    for field in missing:
        report("%s: scheme.yml gives no %s", scheme.name, field)
    recommended = scheme.missing_recommended()
    if recommended:
        log.warning("%s: scheme.yml gives no %s", scheme.name, " or ".join(recommended))
    return missing


def file_checks(scheme: Scheme) -> tuple[pd.DataFrame, list[PrimerPair], dict[str, MarkerDef]]:
    valid, invalid = read_sts(scheme.sts)
    for _, _raw, reason in invalid:
        log.error("%s: invalid STS entry: %s", scheme.name, reason)
    defs = {d.primer_id: d for d in read_markers(scheme.markers)} if scheme.markers else {}
    bins = read_bins(scheme.bins) if scheme.bins else {}
    sts_ids = {p.id for p in valid}
    for name, ids in (("markers.tsv", set(defs)), ("bins.tsv", set(bins))):
        extra = sorted(ids - sts_ids)
        if extra:
            log.warning("%s: %d ids in %s have no STS entry: %s", scheme.name, len(extra), name, ", ".join(extra[:10]))
    rows = []
    for pair in valid:
        rows.append({
            "marker": pair.id, "sts_ok": 1, "in_markers": int(pair.id in defs), "in_bins": int(pair.id in bins),
            "hash_word_11": int(has_hash_word(pair.forward, 11) and has_hash_word(pair.reverse, 11)),
            "hash_word_6": int(has_hash_word(pair.forward, 6) or has_hash_word(pair.reverse, 6)),
        })
    for number, raw, reason in invalid:
        rows.append({"marker": raw.split("\t")[0] or f"line {number}", "sts_ok": 0, "note": reason})
    return pd.DataFrame(rows, columns=CHECK_COLUMNS), valid, defs


def _search_reference(pairs: list[PrimerPair], fasta: Path, workdir: Path, margin: int, mismatches: int) -> dict:
    """Hits per marker from the main search, the rescue pass and the indel alignment."""
    sts = workdir / "scheme_primers.sts"
    write_sts(pairs, sts)
    searcher = MerPcrSearcher(sts, margin=margin, mismatches=mismatches, iupac_mode=1)
    hits = searcher.search_fasta(fasta)
    contigs = load_fasta(fasta)
    missing = [p for p in pairs if p.id not in {h.primer_id for h in hits}]
    if missing:
        hits += rescue_pass(fasta, missing, RescueOptions(6, 2 * margin, False, mismatches))
        missing = [p for p in pairs if p.id not in {h.primer_id for h in hits}]
    for pair in missing:
        hits += find_with_indels(contigs, pair, "reference", 2 * margin, mismatches)
    by_marker: dict[str, list] = {}
    for h in hits:
        by_marker.setdefault(h.primer_id, []).append(h)
    return {"hits": by_marker, "contigs": contigs}


def _product_row(h, contigs: dict[str, str], d: MarkerDef | None, n_hits: int, unit_mismatches: int) -> dict:
    seq = contigs[h.contig][h.start - 1 : h.end]
    if h.strand == "-":
        seq = str(Seq(seq).reverse_complement())
    row: dict = {"length": h.length, "search_pass": h.search_pass}
    if d is not None and d.ref_size is not None:
        row["size_ok"] = int(h.length == d.ref_size)
    if d is not None and d.has_repeats:
        copies = count_tandem_array(seq, d.patterns, unit_mismatches).copies
        row["tandem_copies"] = copies
        if d.ref_repeats is not None:
            row["repeats_ok"] = int(copies == d.ref_repeats)
    notes = []
    if n_hits > 1:
        notes.append(f"{n_hits} products")
    if h.search_pass != "main":
        notes.append(f"found by the {h.search_pass} pass")
    if h.total_mismatches > 0:
        notes.append(f"{h.total_mismatches} primer mismatches")
    row["note"] = "; ".join(notes)
    return row


def assembly_checks(table: pd.DataFrame, pairs: list[PrimerPair], defs: dict[str, MarkerDef], fasta: Path,
                    workdir: Path, margin: int = 500, mismatches: int = 2, unit_mismatches: int = 1) -> pd.DataFrame:
    """Search the assembly and fill the product columns of the check table."""
    found = _search_reference(pairs, fasta, workdir, margin, mismatches)
    table = table.set_index("marker").astype(object)
    for pair in pairs:
        marker_hits = found["hits"].get(pair.id, [])
        d = defs.get(pair.id)
        row: dict = {"n_products": len(marker_hits)}
        if d is not None:
            row["ref_size"] = "" if d.ref_size is None else d.ref_size
            row["ref_repeats"] = "" if d.ref_repeats is None else d.ref_repeats
        if marker_hits:
            row.update(_product_row(marker_hits[0], found["contigs"], d, len(marker_hits), unit_mismatches))
        else:
            row["note"] = "not found in the reference assembly"
        for key, value in row.items():
            table.loc[pair.id, key] = value
    return table.reset_index()


def fetch_assembly(accession: str, email: str, out: Path) -> Path:
    """Download a nucleotide record from NCBI as FASTA."""
    from Bio import Entrez

    Entrez.email = email
    log.info("Fetching %s from NCBI nucleotide", accession)
    with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
        text = handle.read()
    if not text.startswith(">"):
        raise RuntimeError(f"no FASTA record returned for {accession}")
    out.write_text(text, encoding="utf-8")
    return out


def summary(table: pd.DataFrame, with_assembly: bool) -> str:
    n = len(table)
    parts = [f"{int((table['sts_ok'] == 1).sum())} of {n} STS entries valid",
             f"{int((table['hash_word_11'] == 1).sum())} searchable at word size 11",
             f"{int((table['in_markers'] == 1).sum())} with a marker definition"]
    if with_assembly:
        products = pd.to_numeric(table["n_products"], errors="coerce")
        parts.append(f"{int((products == 1).sum())} with exactly one product in the reference")
        parts.append(f"{int((products == 0).sum())} not found")
        size_ok = pd.to_numeric(table["size_ok"], errors="coerce")
        if size_ok.notna().any():
            parts.append(f"{int((size_ok == 1).sum())} of {int(size_ok.notna().sum())} with the reference size")
        rep_ok = pd.to_numeric(table["repeats_ok"], errors="coerce")
        if rep_ok.notna().any():
            parts.append(f"{int((rep_ok == 1).sum())} of {int(rep_ok.notna().sum())} with the reference repeats")
    return "; ".join(parts)
