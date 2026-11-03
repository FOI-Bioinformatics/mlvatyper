"""Command line interface.

Subcommands:
  run        search, extract, align and call alleles in one go
  search     search and extract only (writes hits.tsv and sequences/)
  summarize  rebuild summary tables and allele calls from an existing hits.tsv
  discover   find repeat units in the amplicons of one sample and write a marker definition file
  report     write a self-contained HTML report from the tables of a run
  schemes    list the bundled schemes, or describe one with its loci and provenance
  new-scheme   create a scheme directory with annotated template files
  check-scheme   check a scheme's files and, with an assembly, its reference values
  bins-from-table   build a bin file from a published allele size table
"""

import argparse
import glob
import hashlib
import json
import logging
import sys
import textwrap
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from datetime import datetime, UTC
from pathlib import Path

from . import __version__
from .align import align_all, detect_aligner
from .discover import Thresholds, discover_units, discoveries_to_definitions, write_candidates, write_report
from .discover_multi import apply_support, pool_units, reference_candidates
from .trf import discover_with_trf, trf_available
from .calling import (
    ROUNDING_POLICIES, allele_table, call_all, derive_definitions, long_table, profile_table, read_bins,
    status_table, tandem_count_table,
)
from .extract import (
    extract_amplicons, load_fasta, primer_end_detail, primers_near_contig_ends, read_hits_table, read_samples,
    write_hits_table, write_marker_fastas, write_samples,
)
from .io import strip_compression
from .log import configure_logging, get_logger
from .markers import MarkerError, check_against_sts, markers_from_ids, merge_definitions, read_markers, write_markers
from .models import MarkerDef, PrimerPair
from .checkscheme import (
    assembly_checks, fetch_assembly, file_checks, provenance_checks, resolve_reference, summary as check_summary,
)
from .mlvabank import bins_from_table, export_profiles
from .report import write_report as write_html_report
from .scheme import SchemeError, apply_scheme, list_schemes, load_scheme, scaffold_scheme
from .rescue import RescueOptions, find_with_indels, missing_reasons, read_missing, rescue_pass, write_missing
from .search import MerPcrSearcher, SearchError, sample_name
from .sts import is_degenerate, read_sts, write_invalid, write_sts
from .summarize import write_all_summaries, write_table

log = get_logger()

HITS_FILE = "hits.tsv"
SAMPLES_FILE = "samples.txt"
SEQ_DIR = "sequences"
STS_USED = "primers_used.sts"
NEAR_END_FILE = "near_contig_end.tsv"
MISSING_FILE = "missing_markers.tsv"
DERIVED_FILE = "marker_defs_derived.tsv"
DISCOVERED_FILE = "marker_defs_discovered.tsv"
DISCOVER_REPORT = "discover_report.tsv"
DISCOVER_CANDIDATES = "discover_candidates.tsv"
RUN_INFO = "run_info.json"
SCHEME_CHECK = "scheme_check.tsv"
# parser defaults that a scheme may override (see scheme.apply_scheme)
CALLING_DEFAULTS = {"rounding": "intermediate", "allele_tolerance": 0.25, "unit_mismatches": 1, "margin": 500,
                    "mismatches": 2, "max_array_length": None}


def _add_search_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--sts", help="STS primer file (id, forward, reverse, expected size[, alias]); "
                                 "not needed with --scheme")
    _add_scheme_arg(p)
    p.add_argument("--fasta", required=True, nargs="+",
                   help="FASTA file(s) or glob pattern(s), one file per sample")
    p.add_argument("-M", "--margin", type=int, default=500,
                   help="size window: the reverse primer is searched within this many bp of the STS "
                        "expected size (default 500)")
    p.add_argument("-N", "--mismatches", type=int, default=2,
                   help="mismatches allowed per primer, outside the protected 3' bases (default 2)")
    p.add_argument("--three-prime-match", type=int, default=1,
                   help="number of 3' primer bases that must match exactly (default 1)")
    p.add_argument("-T", "--threads", type=int, default=1,
                   help="merPCR worker threads within one sample (default 1)")
    p.add_argument("--jobs", type=int, default=1,
                   help="samples searched in parallel, each in its own process (default 1); the total worker "
                        "count is --jobs times --threads")
    p.add_argument("-W", "--wordsize", type=int, default=11,
                   help="merPCR hash word size, 3 to 16; a primer needs this many consecutive unambiguous "
                        "bases to be searched (default 11)")
    p.add_argument("--one-sided", action="store_true",
                   help="search a marker on one strand when only one of its primers has a hash word, "
                        "instead of dropping it as me-PCR does (merPCR -U 1)")
    p.add_argument("--no-iupac", dest="iupac", action="store_false",
                   help="do not interpret IUPAC ambiguity codes in primers; degenerate positions then "
                        "count as mismatches, as in the me-PCR default")
    p.add_argument("--no-dedup", dest="dedup", action="store_false",
                   help="keep overlapping hits for the same locus instead of keeping the best one")
    p.add_argument("--skip-invalid", action="store_true",
                   help="skip invalid STS entries instead of aborting; they are listed in invalid_primers.tsv")
    p.add_argument("--fix-sts", action="store_true",
                   help="remove characters that are not IUPAC codes from primers, logging each removal")
    p.add_argument("--rescue-wordsize", type=int, default=6,
                   help="word size of a second, one-sided merPCR pass for markers without a hit; merPCR needs "
                        "an exact, unambiguous match of the word, so a smaller word finds primers with IUPAC "
                        "codes or a mismatch inside the main word (default %(default)s, 0 disables)")
    p.add_argument("--rescue-margin", type=int,
                   help="size window of the rescue passes (default twice -M)")
    p.add_argument("--max-repeats", type=int, default=50,
                   help="largest repeat number of the per-marker size window written from --markers "
                        "(default 50)")
    p.add_argument("--rescue-indels", action="store_true",
                   help="for markers still missing, anchor on an exact match of one primer and align the other "
                        "with up to one indel")


def _add_scheme_arg(p: argparse.ArgumentParser) -> None:
    if not any(a.dest == "scheme" for a in p._actions):
        p.add_argument("--scheme", help="bundled scheme name or scheme directory (primers.sts, markers.tsv, "
                                        "bins.tsv, scheme.yml); provides --sts, --markers, --bins and the "
                                        "calling defaults of the scheme")


def _add_summary_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--prefix", required=True, help="prefix for summary file names")
    _add_scheme_arg(p)
    p.add_argument("--mlvabank", action="store_true",
                   help="also write <prefix>_mlvabank.csv: one row per sample, one column per marker in scheme "
                        "order under its label")
    p.add_argument("--mlvabank-missing", metavar="TOKEN",
                   help="token written in the MLVAbank file for cells without an integer allele "
                        "(default: the status code or value as reported)")
    p.add_argument("--markers", "--pattern", dest="markers",
                   help="marker definition file: marker, patterns, unit_length, ref_size, ref_repeats, offset")
    p.add_argument("--reference-sample",
                   help="sample whose amplicons provide ref_size and ref_repeats for markers lacking them")
    p.add_argument("--bins", help="bin file: marker, size or low-high range, allele")
    p.add_argument("--rounding", choices=ROUNDING_POLICIES, default="intermediate",
                   help="policy for repeat numbers outside the tolerance (default intermediate: n.5)")
    p.add_argument("--allele-tolerance", type=float, default=0.25,
                   help="repeat numbers within this distance of an integer are that integer (default 0.25)")
    p.add_argument("--unit-mismatches", type=int, default=1,
                   help="mismatches allowed per repeat unit when counting tandem arrays (default 1)")
    p.add_argument("--max-array-length", type=int,
                   help="censor alleles whose repeat array is at least this long (assembly k-mer or read length)")
    p.add_argument("--legacy-tables", action="store_true",
                   help="also write the exact pattern count and shortest-relative length difference tables")
    p.add_argument("--clwclean", action="store_true",
                   help="also produce outputs restricted to sequences without N")


def _add_align_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--align", default="auto", choices=["auto", "none", "mafft", "muscle"],
                   help="aligner for per-marker alignments (default auto)")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mlvatyper", description="Locate VNTR amplicons with merPCR and call MLVA alleles"
    )
    parser.add_argument("--version", action="version", version=f"mlvatyper {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--log-json", metavar="FILE",
                        help="also write one JSON object per log event (JSON Lines) to FILE, for pipeline monitors")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="search, extract, align and call alleles")
    run.add_argument("--outdir", required=True)
    _add_search_args(run)
    _add_align_args(run)
    _add_summary_args(run)
    run.add_argument("--no-report", dest="report", action="store_false", help="do not write the HTML report")
    run.add_argument("--title", help="report title (default: MLVA typing report and the STS file name)")

    search = sub.add_parser("search", help="search and extract amplicons only")
    search.add_argument("--outdir", required=True)
    _add_search_args(search)
    search.add_argument("--markers", "--pattern", dest="markers",
                        help="marker definition file; markers with ref_size and ref_repeats get a per-marker "
                             "size window from 0 to --max-repeats repeats in place of -M")
    _add_align_args(search)
    search.add_argument("--clwclean", action="store_true", help="also align sequences without N")

    summ = sub.add_parser("summarize", help="rebuild tables and allele calls from hits.tsv")
    summ.add_argument("--outdir", required=True, help="directory containing hits.tsv and samples.txt")
    _add_summary_args(summ)

    schemes_p = sub.add_parser("schemes", help="list the bundled schemes, or describe one")
    schemes_p.add_argument("name", nargs="?", help="scheme name, or part of a name or organism; one match is "
                                                   "described in full, several are listed")
    schemes_p.add_argument("--names", action="store_true", help="print the scheme names only, one per line")
    schemes_p.add_argument("--tsv", action="store_true", help="print the table tab separated with the scheme "
                                                             "directory as a last column, for scripts")

    new = sub.add_parser("new-scheme", help="create a scheme directory with annotated template files")
    new.add_argument("directory", help="directory to create")
    new.add_argument("--name", help="scheme name (default: the directory name)")

    chk = sub.add_parser("check-scheme", help="check a scheme's files and reference values")
    chk.add_argument("--scheme", required=True, help="bundled scheme name or scheme directory")
    chk.add_argument("--outdir", required=True, help=f"directory for {SCHEME_CHECK} and fetched sequences")
    chk.add_argument("--reference", help="reference assembly FASTA to search the primers in (default: the "
                                         "assembly accession of the scheme, downloaded from NCBI)")
    chk.add_argument("--assembly", metavar="GCF_or_GCA", help="assembly accession to download instead of the "
                                                             "one named in scheme.yml")
    chk.add_argument("--no-download", dest="download", action="store_false",
                     help="do not download the scheme's assembly; check the files only")
    chk.add_argument("--cache-dir", help="where downloaded assemblies are kept (default $MLVATYPER_CACHE or "
                                         "~/.cache/mlvatyper/assemblies)")
    chk.add_argument("--fetch", metavar="ACCESSION",
                     help="download a nucleotide record (for example a chromosome accession) from NCBI "
                          "instead of an assembly")
    chk.add_argument("--email", help="e-mail address required by NCBI for --fetch")
    chk.add_argument("-M", "--margin", type=int, default=500)
    chk.add_argument("-N", "--mismatches", type=int, default=2)
    chk.add_argument("--unit-mismatches", type=int, default=1)
    chk.add_argument("--allow-missing-provenance", action="store_true",
                     help="report missing citation, url or reference as a warning instead of failing; "
                          "for a scheme still being prepared")

    bft = sub.add_parser("bins-from-table", help="build a bin file from a long allele size table")
    bft.add_argument("table", help="tab or comma separated table with locus, allele and size columns")
    bft.add_argument("--output", required=True, help="bin file to write (marker, size, allele)")
    bft.add_argument("--locus-column")
    bft.add_argument("--allele-column")
    bft.add_argument("--size-column")
    bft.add_argument("--markers", help="marker file whose labels map the table's locus names to marker ids")

    rep = sub.add_parser("report", help="write a self-contained HTML report from the tables of a run")
    rep.add_argument("--outdir", required=True, help="directory with the summarize outputs")
    rep.add_argument("--prefix", required=True, help="prefix used by summarize")
    rep.add_argument("--title", help="report title (default: MLVA typing report and the STS file name)")

    disc = sub.add_parser("discover", help="find repeat units in the amplicons of one or more samples")
    disc.add_argument("--outdir", required=True, help="directory containing hits.tsv")
    disc.add_argument("--sample", action="append", default=None,
                      help="sample to use; repeat for several samples (default: the only sample)")
    disc.add_argument("--all-samples", action="store_true", help="use every sample in hits.tsv")
    disc.add_argument("--reference-sample",
                      help="sample whose amplicons provide ref_size and ref_repeats (default: first sample used)")
    disc.add_argument("--min-unit", type=int, default=3, help="smallest repeat unit length (default 3)")
    disc.add_argument("--max-unit", type=int, default=100, help="largest repeat unit length (default 100)")
    disc.add_argument("--min-copies", type=float, default=1.8,
                      help="smallest array, in engine copies, that becomes a candidate; two integer copies "
                           "are always required (default 1.8)")
    disc.add_argument("--unit-mismatches", type=int, default=1,
                      help="mismatches tolerated per repeat unit, at most one per six bases (default 1)")
    disc.add_argument("--min-score", type=int, default=30,
                      help="array score required for categories A and B, +2 per matching and -7 per differing "
                           "position between consecutive copies (default 30; Tandem Repeats Finder uses 50)")
    disc.add_argument("--weak-score", type=int, default=20,
                      help="score floor for category C candidates (default 20)")
    disc.add_argument("--strong-unit", type=int, default=5,
                      help="shortest unit length for categories A and B (default 5)")
    disc.add_argument("--strong-copies", type=int, default=3,
                      help="integer copies required for category A (default 3)")
    disc.add_argument("--min-entropy", type=float, default=1.0,
                      help="unit base entropy in bits required for category A (default 1.0)")
    disc.add_argument("--accept", choices=["A", "B", "C"], default="B",
                      help="least category whose primary candidate seeds a definition (default B)")
    disc.add_argument("--allele-tolerance", type=float, default=0.25,
                      help="tolerance when testing whether a unit length explains length differences "
                           "between samples (default 0.25)")
    disc.add_argument("--engine", choices=["auto", "builtin", "trf"], default="auto",
                      help="repeat finder: Tandem Repeats Finder (trf) when on PATH, else the built-in "
                           "self-comparison (default auto)")
    disc.add_argument("--output", default=DISCOVERED_FILE,
                      help=f"marker definition file name within outdir (default {DISCOVERED_FILE})")
    return parser


def expand_fasta(patterns: list[str]) -> list[Path]:
    """FASTA files from paths and glob patterns.

    A path without glob characters that does not exist stops the run, since
    a misspelled sample would otherwise be left out without notice; a glob
    pattern that matches nothing is reported as a warning.
    """
    files: list[Path] = []
    for pat in patterns:
        matches = sorted(glob.glob(pat))
        if not matches and Path(pat).exists():
            matches = [pat]
        if not matches:
            if any(c in pat for c in "*?["):
                log.warning("--fasta pattern %s matches no file", pat)
            else:
                raise SystemExit(f"FASTA file {pat} not found")
        files.extend(Path(m) for m in matches)
    return files


def resolve_samples(patterns: list[str]) -> list[Path]:
    """Expand FASTA globs and check that sample names are unique."""
    fasta_files = expand_fasta(patterns)
    if not fasta_files:
        raise SystemExit(f"no FASTA files match {patterns}")
    names = [sample_name(f) for f in fasta_files]
    for f, name in zip(fasta_files, names, strict=True):
        if name != strip_compression(f).stem:
            log.warning("%s: sample name %s, with whitespace replaced by _", f, name)
    if len(set(names)) != len(names):
        raise SystemExit("sample names (FASTA file stems) are not unique")
    return fasta_files


def prepare_sts(sts_path: str, outdir: Path, skip_invalid: bool, iupac: bool,
                fix: bool = False) -> tuple[Path, list[PrimerPair]]:
    """Validate the STS file and write the accepted entries to outdir/primers_used.sts."""
    valid, invalid = read_sts(sts_path, fix)
    for _, _, reason in invalid:
        log.error("invalid STS entry: %s", reason)
    if invalid and not skip_invalid:
        raise SystemExit(f"{len(invalid)} invalid STS entries in {sts_path}; fix them or use --skip-invalid")
    if invalid:
        write_invalid(invalid, outdir / "invalid_primers.tsv")
    if not valid:
        raise SystemExit("no valid STS entries")
    degenerate = sum(is_degenerate(p) for p in valid)
    log.info("%d valid primer pairs loaded from %s (%d with IUPAC ambiguity codes, iupac mode %s)",
             len(valid), sts_path, degenerate, "on" if iupac else "off")
    clean_sts = outdir / STS_USED
    write_sts(valid, clean_sts)
    return clean_sts, valid


def size_windows_from_definitions(pairs: list[PrimerPair], defs: list[MarkerDef],
                                  max_repeats: int) -> list[PrimerPair]:
    """Give every pair with a callable definition the STS size range of 0 to max_repeats repeats.

    merPCR then searches the reverse primer within that range for the marker
    instead of within -M of the single expected size. Hits keep the STS
    expected size rather than the midpoint of the window. A pair whose STS
    size is already a range keeps it.
    """
    by_id = {d.primer_id: d for d in defs if d.can_call}
    result = []
    changed = 0
    for pair in pairs:
        d = by_id.get(pair.id)
        if d is None or pair.size_low is not None:  # an explicit range in the STS file is kept
            result.append(pair)
            continue
        low = max(int(d.expected_size(0)), len(pair.forward) + len(pair.reverse))
        high = int(d.expected_size(max_repeats))
        result.append(replace(pair, size_low=low, size_high=high))  # expected_size stays the STS value
        changed += 1
    log.info("%d markers searched with a size window from their definitions (0 to %d repeats)", changed, max_repeats)
    return result


def _renumber(amplicons: list) -> list:
    """Copy numbers per marker after merging amplicons from several passes."""
    from .models import Amplicon
    counter: dict[str, int] = {}
    result = []
    for amp in sorted(amplicons, key=lambda a: (a.hit.primer_id, a.hit.contig, a.hit.start)):
        counter[amp.hit.primer_id] = counter.get(amp.hit.primer_id, 0) + 1
        result.append(Amplicon(amp.hit, counter[amp.hit.primer_id], amp.sequence))
    return result


def search_sample(searcher: MerPcrSearcher, fasta: Path, pairs: list[PrimerPair], margin: int, dedup: bool,
                  rescue: RescueOptions, wordsize: int, one_sided: bool) -> tuple[list, dict[str, str]]:
    """Main search, rescue passes and reasons for the markers that stay missing, for one sample."""
    primers = {p.id: p for p in pairs}
    sample = sample_name(fasta)
    hits = searcher.search_fasta(fasta)
    contigs = load_fasta(fasta)
    amps = extract_amplicons(hits, contigs, primers, deduplicate=dedup)

    def still_missing() -> list[PrimerPair]:
        found = {a.hit.primer_id for a in amps}
        return [p for p in pairs if p.id not in found]

    missing = still_missing()
    if missing and rescue.wordsize > 0:
        amps += extract_amplicons(rescue_pass(fasta, missing, rescue), contigs, primers, deduplicate=dedup)
        missing = still_missing()
    if missing and rescue.indels:
        extra = [h for p in missing for h in find_with_indels(contigs, p, sample, rescue.margin, rescue.mismatches)]
        if extra:
            log.info("%s: primer alignment found %d hits for %d markers", sample, len(extra),
                     len({h.primer_id for h in extra}))
        amps += extract_amplicons(extra, contigs, primers, deduplicate=dedup)
        missing = still_missing()
    reasons: dict[str, str] = {}
    if missing:
        # per marker: its own expected size plus -M, not the largest one of all missing markers. The top of a
        # definition window (50 repeats by default) would mark almost every primer in a small contig as near an end.
        windows = {p.id: p.expected_size + margin for p in missing}
        near = primers_near_contig_ends(contigs, missing, windows)
        used_wordsize = rescue.wordsize if rescue.wordsize > 0 else wordsize
        details = {p.id: primer_end_detail(contigs, p, windows[p.id]) for p in missing if p.id in near}
        reasons = missing_reasons(contigs, missing, near, used_wordsize, one_sided or rescue.wordsize > 0, details)
    return _renumber(amps), reasons


# fork-safe per-process state for the parallel path, filled by _init_worker
_WORKER: dict = {}

SAMPLE_ERRORS = (SearchError, ValueError, OSError, EOFError)  # not FASTA, or a corrupt/truncated compressed file


def _init_worker(sts_path, searcher_kwargs, params) -> None:
    _WORKER["searcher"] = MerPcrSearcher(sts_path, **searcher_kwargs)
    _WORKER["params"] = params


def _search_one(fasta: Path):
    """Run one sample in a worker; return (fasta, amplicons, reasons) or (fasta, None, error text)."""
    pairs, margin, dedup, rescue, wordsize, one_sided = _WORKER["params"]
    try:
        amps, reasons = search_sample(_WORKER["searcher"], fasta, pairs, margin, dedup, rescue, wordsize, one_sided)
    except SAMPLE_ERRORS as exc:
        return fasta, None, str(exc).splitlines()[0] if str(exc) else type(exc).__name__
    return fasta, amps, reasons


def search_samples(searcher: MerPcrSearcher, fasta_files: list[Path], pairs: list[PrimerPair], margin: int,
                   dedup: bool = True, rescue: RescueOptions | None = None, wordsize: int = 11,
                   one_sided: bool = False, jobs: int = 1, sts_path: str | Path | None = None,
                   searcher_kwargs: dict | None = None) -> tuple[list, dict[tuple[str, str], str], list[Path]]:
    """Search and extract amplicons for every sample.

    With jobs > 1 the samples are searched in parallel, one process each, and
    the results are merged in input order so the outputs are identical to the
    serial run. Each process builds its own searcher from sts_path and
    searcher_kwargs, so those are required when jobs > 1. Returns (amplicons,
    reason per (sample, marker) without a hit, FASTA files that failed).
    """
    rescue = rescue or RescueOptions(wordsize=0)
    params = (pairs, margin, dedup, rescue, wordsize, one_sided)
    workers = max(1, min(jobs, len(fasta_files)))
    if workers == 1:
        results = []
        for fasta in fasta_files:
            _WORKER_SERIAL = searcher  # noqa: N806 - keep the serial path a plain loop
            try:
                amps, sample_reasons = search_sample(_WORKER_SERIAL, fasta, pairs, margin, dedup, rescue,
                                                      wordsize, one_sided)
            except SAMPLE_ERRORS as exc:
                results.append((fasta, None, str(exc).splitlines()[0] if str(exc) else type(exc).__name__))
                continue
            results.append((fasta, amps, sample_reasons))
    else:
        if sts_path is None or searcher_kwargs is None:
            raise ValueError("parallel search needs sts_path and searcher_kwargs")
        log.info("Searching %d samples with %d parallel jobs", len(fasta_files), workers)
        by_fasta: dict[Path, tuple] = {}
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                                 initargs=(str(sts_path), searcher_kwargs, params)) as pool:
            for fasta, amps, reasons in pool.map(_search_one, fasta_files):
                by_fasta[fasta] = (fasta, amps, reasons)
        results = [by_fasta[fasta] for fasta in fasta_files]  # input order, for deterministic output
    amplicons = []
    reasons: dict[tuple[str, str], str] = {}
    failed = []
    for fasta, amps, sample_reasons in results:
        if amps is None:
            log.error("%s: %s", fasta, sample_reasons)
            failed.append(fasta)
            continue
        if workers > 1:
            log.info("%s: %d amplicons", sample_name(fasta), len(amps))
        amplicons.extend(amps)
        reasons.update({(sample_name(fasta), m): r for m, r in sample_reasons.items()})
    return amplicons, reasons, failed


def write_near_end(near_end: set[tuple[str, str]], path: Path) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("sample\tmarker\n")
        for sample, marker in sorted(near_end):
            handle.write(f"{sample}\t{marker}\n")


def read_near_end(path: Path) -> set[tuple[str, str]]:
    if not path.exists():
        return set()
    lines = path.read_text(encoding="utf-8").splitlines()[1:]
    return {tuple(line.split("\t")[:2]) for line in lines if line.strip()}


def file_record(path: Path) -> dict:
    """Name, size and MD5 of an input file for the run record."""
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return {"path": str(path), "size": path.stat().st_size, "md5": digest.hexdigest()}


def write_run_info(outdir: Path, args: argparse.Namespace, stage: str) -> None:
    """Record versions and parameters so a run can be reproduced."""
    path = outdir / RUN_INFO
    info = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    try:
        from importlib.metadata import version
        merpcr_version = version("merpcr")
    except Exception:  # pragma: no cover - metadata missing
        merpcr_version = "unknown"
    info["mlvatyper_version"] = __version__
    info["merpcr_version"] = merpcr_version
    if stage == "search":
        info["inputs"] = [file_record(Path(p)) for p in [args.sts, *resolve_samples(args.fasta)]]
    # marker and bin files are inputs too; they are recorded when the stage that reads them runs
    listed = {record["path"] for record in info.get("inputs", [])}
    for extra in (getattr(args, "markers", None), getattr(args, "bins", None)):
        if extra and str(extra) not in listed and Path(extra).exists():
            info.setdefault("inputs", []).append(file_record(Path(extra)))
            listed.add(str(extra))
    info[stage] = {
        "timestamp": datetime.now(UTC).isoformat(timespec="seconds"),
        "parameters": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()
                       if k not in {"command", "verbose"}},
    }
    if stage == "summarize":  # several prefixes can share an outdir; each report reads its own parameters
        info.setdefault("summarize_by_prefix", {})[args.prefix] = info[stage]
    path.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n", encoding="utf-8")


SCHEME_COLUMNS = ("scheme", "organism", "loci", "called", "assembly")
SCHEME_DETAIL_WIDTH = 96  # the notes of a published scheme are a paragraph, so they are wrapped


def _callable_loci(scheme) -> int:
    """Markers with the reference values an allele call needs; a scheme may size the rest only.

    Reading a marker file logs a line per file, which would bury a listing of
    every bundled scheme, so the logger is quietened for the duration.
    """
    if scheme.markers is None:
        return 0
    binned = set()
    if scheme.bins is not None:
        from .calling import read_bins
        binned = set(read_bins(scheme.bins))
    level = log.level
    log.setLevel(logging.WARNING)
    try:
        # a marker calls when it carries the reference values, or when a bins table defines its alleles
        return sum(1 for d in read_markers(scheme.markers) if d.can_call or d.primer_id in binned)
    finally:
        log.setLevel(level)


def _scheme_row(scheme, with_path: bool = False) -> tuple[str, ...]:
    pairs, _ = read_sts(scheme.sts)
    row = (scheme.name, str(scheme.meta.get("organism", "")).strip() or "-", str(len(pairs)),
           str(_callable_loci(scheme)), scheme.assembly or "-")
    return (*row, str(scheme.path)) if with_path else row


def _print_table(header: tuple[str, ...], rows: list[tuple[str, ...]]) -> None:
    widths = [max(len(row[i]) for row in (header, *rows)) for i in range(len(header))]
    for row in (header, *rows):
        print("  ".join(value.ljust(width) for value, width in zip(row, widths, strict=True)).rstrip())


def _print_scheme(scheme) -> None:
    """Everything one scheme records, including its loci and its notes."""
    pairs, _ = read_sts(scheme.sts)
    print(scheme.name)
    fields = (
        ("organism", str(scheme.meta.get("organism", "")).strip()),
        ("loci", f"{len(pairs)}: " + ", ".join(p.id for p in pairs)),
        ("reference", str(scheme.meta.get("reference", "")).strip()),
        ("assembly", scheme.assembly or ""),
        ("citation", str(scheme.meta.get("citation", "")).strip()),
        ("url", str(scheme.meta.get("url", "")).strip()),
        ("retrieved", str(scheme.meta.get("retrieved", "")).strip()),
        ("defaults", ", ".join(f"{k}={v}" for k, v in scheme.calling_defaults().items())),
        ("duplicated", ", ".join(scheme.known_multiple)),
        ("called", f"{_callable_loci(scheme)} of {len(pairs)} loci carry the reference values an allele "
                   "call needs; the rest are sized only"),
        ("files", ", ".join(p.name for p in (scheme.sts, scheme.markers, scheme.bins) if p is not None)),
        ("directory", str(scheme.path)),
        ("notes", str(scheme.meta.get("notes", "")).strip()),
    )
    for label, value in fields:
        if value:
            print(textwrap.fill(value, width=SCHEME_DETAIL_WIDTH, break_long_words=False,
                                break_on_hyphens=False,
                                initial_indent=f"  {label:<11}", subsequent_indent=" " * 13))


def do_schemes(args: argparse.Namespace) -> int:
    """List the bundled schemes, or describe the one the name selects."""
    schemes = list_schemes()
    if args.name:
        needle = args.name.lower()
        exact = [s for s in schemes if s.name.lower() == needle]
        schemes = exact or [s for s in schemes if needle in s.name.lower()
                            or needle in str(s.meta.get("organism", "")).lower()]
        if not schemes:
            log.error("no bundled scheme matches %r; 'mlvatyper schemes' lists them all", args.name)
            return 1
    if args.names:
        for scheme in schemes:
            print(scheme.name)
        return 0
    if args.tsv:
        for scheme in schemes:
            print("\t".join(_scheme_row(scheme, with_path=True)))
        return 0
    if len(schemes) == 1:
        _print_scheme(schemes[0])
        return 0
    _print_table(SCHEME_COLUMNS, [_scheme_row(s) for s in schemes])
    print(f"\n{len(schemes)} schemes. 'mlvatyper schemes NAME' describes one and lists its loci; "
          "'--scheme NAME' uses it in run, search and check-scheme.")
    return 0


def do_search(args: argparse.Namespace, outdir: Path) -> list[Path]:
    """Search every sample; returns the FASTA files that failed, which are left out of samples.txt."""
    fasta_files = resolve_samples(args.fasta)
    clean_sts, pairs = prepare_sts(args.sts, outdir, args.skip_invalid, args.iupac, args.fix_sts)
    if getattr(args, "markers", None):
        pairs = size_windows_from_definitions(pairs, read_markers(args.markers), args.max_repeats)
        write_sts(pairs, clean_sts)
    searcher_kwargs = dict(margin=args.margin, mismatches=args.mismatches, threads=args.threads,
                           iupac_mode=int(args.iupac), wordsize=args.wordsize, one_sided=int(args.one_sided),
                           three_prime_match=args.three_prime_match)
    jobs = getattr(args, "jobs", 1)
    # the serial path builds the searcher once here; each parallel worker builds its own from clean_sts
    searcher = MerPcrSearcher(clean_sts, **searcher_kwargs) if jobs <= 1 or len(fasta_files) <= 1 else None
    rescue = RescueOptions(args.rescue_wordsize, args.rescue_margin or 2 * args.margin, args.rescue_indels,
                           args.mismatches, args.iupac, args.three_prime_match, args.threads)
    amplicons, reasons, failed = search_samples(searcher, fasta_files, pairs, args.margin, args.dedup, rescue,
                                                args.wordsize, args.one_sided, jobs=jobs, sts_path=clean_sts,
                                                searcher_kwargs=searcher_kwargs)

    # a failed sample is left out rather than listed with every marker NOT_DETECTED, which would read as typed
    write_samples([sample_name(f) for f in fasta_files if f not in failed], outdir / SAMPLES_FILE)
    write_hits_table(amplicons, outdir / HITS_FILE)
    write_missing(reasons, outdir / MISSING_FILE)
    near_end = {key for key, reason in reasons.items() if reason == "NEAR_CONTIG_END"}
    write_near_end(near_end, outdir / NEAR_END_FILE)
    seq_dir = outdir / SEQ_DIR
    written = write_marker_fastas(amplicons, seq_dir)
    log.info("%d amplicons for %d markers written to %s", len(amplicons), len(written), seq_dir)
    rescued = sum(a.hit.search_pass != "main" for a in amplicons)
    if rescued:
        log.info("%d amplicons come from the rescue passes (see search_pass in %s)", rescued, HITS_FILE)
    if reasons:
        counts = {}
        for reason in reasons.values():
            counts[reason] = counts.get(reason, 0) + 1
        log.info("%d sample/marker combinations without a hit: %s", len(reasons),
                 ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    log.info("Processed %d samples, %d failed", len(fasta_files) - len(failed), len(failed))
    write_run_info(outdir, args, "search")

    aligner = detect_aligner(args.align)
    if aligner is not None:
        align_all(seq_dir, aligner, clwclean=getattr(args, "clwclean", False))
    return failed


def _failed_exit(failed: list[Path]) -> None:
    if failed:
        raise SystemExit(f"{len(failed)} samples failed and are not in the results: "
                         f"{', '.join(sample_name(f) for f in failed)}")


def load_definitions(args: argparse.Namespace, outdir: Path, sts_ids: list[str], amplicons: list) -> list[MarkerDef]:
    """Marker definitions from the file, MLVAbank style ids and, if asked, a reference sample."""
    from_file = read_markers(args.markers) if args.markers else []
    defs = merge_definitions(from_file, markers_from_ids(sts_ids))
    if sts_ids:
        defs = check_against_sts(defs, set(sts_ids))
        known = {d.primer_id for d in defs}
        defs.extend(MarkerDef(i) for i in sts_ids if i not in known)
        order = {i: n for n, i in enumerate(sts_ids)}
        defs.sort(key=lambda d: order.get(d.primer_id, len(order)))
    if args.reference_sample:
        try:
            defs = derive_definitions(amplicons, defs, args.reference_sample, args.unit_mismatches)
        except ValueError as exc:  # the reference sample has no amplicons
            raise SystemExit(f"{exc}; choose another --reference-sample") from exc
        write_markers(defs, outdir / DERIVED_FILE)
        log.info("Derived marker definitions written to %s", outdir / DERIVED_FILE)
    return defs


CALL_OUTPUTS = ("_alleles.tsv", "_status.tsv", "_tandem_counts.tsv", "_profile.tsv", "_calls_long.tsv",
                "_mlvabank.csv", "_definitions.tsv", "_report.html")


def _remove_stale_calls(outdir: Path, prefix: str) -> None:
    """Remove the call tables of an earlier summarize with this prefix, which no longer describe the run."""
    for suffix in CALL_OUTPUTS:
        path = outdir / f"{prefix}{suffix}"
        if path.exists():
            path.unlink()
            log.info("Removed %s from an earlier run", path)


def write_calls(calls: list, samples: list[str], markers: list[str], outdir: Path, prefix: str,
                defs: list[MarkerDef] | None = None, mlvabank: bool = False,
                mlvabank_missing: str | None = None) -> None:
    alleles = allele_table(calls, samples, markers)
    write_table(alleles, outdir / f"{prefix}_alleles.tsv")
    if defs:
        write_markers(defs, outdir / f"{prefix}_definitions.tsv")  # what these calls were made with
    if mlvabank:
        export_profiles(alleles, defs or [MarkerDef(m) for m in markers], outdir / f"{prefix}_mlvabank.csv",
                        missing=mlvabank_missing)
    elif (outdir / f"{prefix}_mlvabank.csv").exists():  # from an earlier run with --mlvabank
        (outdir / f"{prefix}_mlvabank.csv").unlink()
        log.info("Removed %s from an earlier run", outdir / f"{prefix}_mlvabank.csv")
    write_table(status_table(calls, samples, markers), outdir / f"{prefix}_status.tsv")
    write_table(tandem_count_table(calls, samples, markers), outdir / f"{prefix}_tandem_counts.tsv")
    write_table(profile_table(calls, samples, markers), outdir / f"{prefix}_profile.tsv")
    long_table(calls).to_csv(outdir / f"{prefix}_calls_long.tsv", sep="\t", index=False, na_rep="NA",
                             lineterminator="\n")
    log.info("Wrote %s_calls_long.tsv (%d rows)", outdir / prefix, len(calls))


def do_summarize(args: argparse.Namespace, outdir: Path) -> None:
    hits_path = outdir / HITS_FILE
    if not hits_path.exists():
        raise SystemExit(f"{hits_path} not found; run 'mlvatyper search' first")
    amplicons = read_hits_table(hits_path)
    samples_path = outdir / SAMPLES_FILE
    samples = read_samples(samples_path) if samples_path.exists() else sorted({a.hit.sample for a in amplicons})
    sts_path = outdir / STS_USED
    sts_ids = [p.id for p in read_sts(sts_path)[0]] if sts_path.exists() else []
    if args.reference_sample and args.reference_sample not in samples:
        raise SystemExit(f"reference sample {args.reference_sample!r} is not among the samples")

    defs = load_definitions(args, outdir, sts_ids, amplicons) if (args.markers or args.reference_sample
                                                                     or sts_ids) else []
    skip = set(getattr(args, "interrupted_arrays", ()))
    defs = [replace(d, compare_count=False) if d.primer_id in skip else d for d in defs]
    legacy_defs = defs if (args.markers or args.reference_sample) else None
    write_all_summaries(amplicons, samples, outdir, args.prefix, legacy_defs, args.clwclean, args.legacy_tables,
                        sts_ids)

    if not any(d.can_call or d.has_repeats for d in defs):
        log.info("No marker has reference size and repeats or repeat patterns; allele calling skipped")
        _remove_stale_calls(outdir, args.prefix)
        write_run_info(outdir, args, "summarize")
        return
    if not any(d.can_call for d in defs):
        log.warning("No marker has reference size and repeats; only tandem counts are reported. "
                    "Give --markers with ref_size and ref_repeats, or --reference-sample")
    bins = read_bins(args.bins) if args.bins else None
    calls = call_all(amplicons, samples, defs, args.allele_tolerance, args.rounding, bins,
                     args.unit_mismatches, args.max_array_length, read_near_end(outdir / NEAR_END_FILE),
                     read_missing(outdir / MISSING_FILE))
    write_calls(calls, samples, [d.primer_id for d in defs], outdir, args.prefix, defs, args.mlvabank,
                args.mlvabank_missing)
    write_run_info(outdir, args, "summarize")


def _discover_samples(args: argparse.Namespace, samples: list[str]) -> list[str]:
    chosen = samples if args.all_samples else args.sample
    if not chosen:
        if len(samples) != 1:
            raise SystemExit(f"--sample or --all-samples is required when hits.tsv holds several samples: {samples}")
        chosen = samples
    unknown = [s for s in chosen if s not in samples]
    if unknown:
        raise SystemExit(f"samples {unknown} have no amplicons in hits.tsv")
    return list(dict.fromkeys(chosen))


def _primer_lengths(outdir: Path) -> dict[str, tuple[int, int]]:
    sts_path = outdir / STS_USED
    if not sts_path.exists():
        return {}
    return {p.id: (len(p.forward), len(p.reverse)) for p in read_sts(sts_path)[0]}


def _discover_one(args: argparse.Namespace, amplicons: list, sample: str, engine: str, thresholds: Thresholds,
                  primer_lengths: dict[str, tuple[int, int]]) -> tuple[dict, list[str]]:
    if engine == "trf":
        if not trf_available():
            raise SystemExit("trf not found on PATH; install it with: conda install -c bioconda trf, "
                             "or use --engine builtin")
        found, missing, command = discover_with_trf(amplicons, sample, args.min_score, args.max_unit,
                                                    args.unit_mismatches, args.min_copies, thresholds, primer_lengths)
        args.trf_command = " ".join(command)
        return found, missing
    return discover_units(amplicons, sample, args.min_unit, args.max_unit, args.min_copies, args.unit_mismatches,
                          args.min_score, thresholds, primer_lengths)


def do_discover(args: argparse.Namespace, outdir: Path) -> None:
    hits_path = outdir / HITS_FILE
    if not hits_path.exists():
        raise SystemExit(f"{hits_path} not found; run 'mlvatyper search' first")
    amplicons = read_hits_table(hits_path)
    samples = _discover_samples(args, sorted({a.hit.sample for a in amplicons}))
    reference = args.reference_sample or samples[0]
    if reference not in samples:
        raise SystemExit(f"reference sample {reference!r} is not among the samples used: {samples}")
    engine = args.engine
    if engine == "auto":
        engine = "trf" if trf_available() else "builtin"
        log.info("Repeat finder: %s", engine)
    args.engine = engine
    thresholds = Thresholds(args.min_score, args.weak_score, args.min_copies, args.strong_unit, args.strong_copies,
                            args.min_entropy, args.accept)
    primer_lengths = _primer_lengths(outdir)
    by_sample = {}
    for sample in samples:
        by_sample[sample] = _discover_one(args, amplicons, sample, engine, thresholds, primer_lengths)[0]
    if len(samples) > 1:
        lengths: dict[str, dict[str, int]] = {}
        for a in amplicons:
            if a.hit.sample in samples and a.copy == 1:
                lengths.setdefault(a.hit.primer_id, {})[a.hit.sample] = len(a.sequence)
        by_sample = apply_support(by_sample, pool_units(by_sample, lengths, args.allele_tolerance))
        found, missing = reference_candidates(by_sample, reference, amplicons, args.unit_mismatches, primer_lengths)
    else:
        found = by_sample[reference]
        missing = [m for m in {a.hit.primer_id for a in amplicons if a.hit.sample == reference} if m not in found]
    sts_path = outdir / STS_USED
    order = {p.id: n for n, p in enumerate(read_sts(sts_path)[0])} if sts_path.exists() else {}
    found = dict(sorted(found.items(), key=lambda kv: order.get(kv[0], len(order))))
    missing.sort(key=lambda m: order.get(m, len(order)))
    write_report(found, missing, outdir / DISCOVER_REPORT)
    all_candidates = {f"{s}:{m}": c for s in samples for m, c in by_sample.get(s, {}).items()}
    if len(samples) > 1:
        all_candidates.update({f"{reference}:{m}": c for m, c in found.items()})
    write_candidates(all_candidates, outdir / DISCOVER_CANDIDATES)
    defs = discoveries_to_definitions(found, missing, amplicons, reference, thresholds)
    defs.sort(key=lambda d: order.get(d.primer_id, len(order)))
    write_markers(defs, outdir / args.output)
    n_def = sum(d.can_call for d in defs)
    log.info("%d marker definitions with a unit written to %s, details in %s and %s", n_def, outdir / args.output,
             outdir / DISCOVER_REPORT, outdir / DISCOVER_CANDIDATES)
    write_run_info(outdir, args, "discover")


def do_check_scheme(args: argparse.Namespace, outdir: Path) -> None:
    scheme = load_scheme(args.scheme)
    log.info("Checking scheme %s (%s)", scheme.name, scheme.path)
    missing = provenance_checks(scheme, strict=not args.allow_missing_provenance)
    table, pairs, defs = file_checks(scheme)
    try:
        if args.fetch:
            if not args.email:
                raise SystemExit("--fetch needs --email, which NCBI requires")
            fasta = fetch_assembly(args.fetch, args.email, outdir / f"{args.fetch}.fasta")
        else:
            fasta = resolve_reference(scheme, args.reference, args.assembly, args.download, args.cache_dir)
    except (FileNotFoundError, ValueError, RuntimeError, OSError) as exc:
        raise SystemExit(f"{scheme.name}: reference assembly: {exc}") from exc
    if fasta is not None:
        table = assembly_checks(table, pairs, defs, fasta, outdir, args.margin, args.mismatches,
                                args.unit_mismatches)
    table.to_csv(outdir / SCHEME_CHECK, sep="\t", index=False, na_rep="")
    log.info("%s: %s; details in %s", scheme.name, check_summary(table, fasta is not None), outdir / SCHEME_CHECK)
    if missing:
        if not args.allow_missing_provenance:
            raise SystemExit(f"{scheme.name}: scheme.yml must state {', '.join(missing)}; the allele numbers of a "
                             f"scheme are only comparable when their source is recorded "
                             f"(--allow-missing-provenance to continue anyway)")
        log.warning("%s: the allele numbers cannot be traced without %s", scheme.name, ", ".join(missing))


def do_bins_from_table(args: argparse.Namespace) -> None:
    labels = {}
    if args.markers:
        labels = {d.label: d.primer_id for d in read_markers(args.markers) if d.label}
    bins_from_table(args.table, args.output, args.locus_column, args.allele_column, args.size_column, labels)


def explicit_options(argv: list[str]) -> set[str]:
    """Destinations of the options given in argv, whatever their value.

    The command line is parsed again by a parser whose defaults are all
    suppressed, so only what was typed ends up in the namespace.
    """
    parser = build_parser()
    parsers = [parser]
    while parsers:
        p = parsers.pop()
        for action in p._actions:
            action.default = argparse.SUPPRESS
            if isinstance(action, argparse._SubParsersAction):
                parsers.extend(action.choices.values())
    return set(vars(parser.parse_args(argv)))


def _resolve_scheme(args: argparse.Namespace, explicit: set[str] | None = None) -> None:
    if getattr(args, "scheme", None):
        apply_scheme(args, load_scheme(args.scheme), CALLING_DEFAULTS, explicit)
    if args.command in ("run", "search") and not args.sts:
        raise SystemExit("--sts or --scheme is required")


def _dispatch(args: argparse.Namespace, outdir: Path, explicit: set[str] | None = None) -> None:
    _resolve_scheme(args, explicit)
    if args.command == "check-scheme":
        do_check_scheme(args, outdir)
    elif args.command == "search":
        _failed_exit(do_search(args, outdir))
    elif args.command == "summarize":
        do_summarize(args, outdir)
    elif args.command == "discover":
        do_discover(args, outdir)
    elif args.command == "report":
        try:
            write_html_report(outdir, args.prefix, args.title)
        except FileNotFoundError as exc:
            raise SystemExit(str(exc)) from exc
    elif args.command == "run":
        # the samples that did not fail are typed and reported; the exit status still reports the failures
        failed = do_search(args, outdir)
        do_summarize(args, outdir)
        if args.report and (outdir / f"{args.prefix}_calls_long.tsv").exists():
            write_html_report(outdir, args.prefix, args.title)
        _failed_exit(failed)


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    args = build_parser().parse_args(argv)
    configure_logging(args.verbose, args.log_json)
    if args.command == "schemes":
        return do_schemes(args)
    if args.command == "new-scheme":
        try:
            scaffold_scheme(args.directory, args.name)
        except SchemeError as exc:
            log.error("%s", exc)
            return 1
        return 0
    if args.command == "bins-from-table":
        try:
            do_bins_from_table(args)
        except (ValueError, FileNotFoundError) as exc:
            log.error("%s", exc)
            return 1
        return 0
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        _dispatch(args, outdir, explicit_options(argv))
    except SystemExit as exc:
        if exc.code not in (None, 0):
            log.error("%s", exc.code)
            return 1
    except (SchemeError, MarkerError) as exc:
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
