"""Command line interface.

Subcommands:
  run        search, extract, align and call alleles in one go
  search     search and extract only (writes hits.tsv and sequences/)
  summarize  rebuild summary tables and allele calls from an existing hits.tsv
  discover   find repeat units in the amplicons of one sample and write a marker definition file
"""

import argparse
import glob
import json
import sys
from datetime import datetime, UTC
from pathlib import Path

from . import __version__
from .align import align_all, detect_aligner
from .discover import discover_units, discoveries_to_definitions, write_report
from .trf import discover_with_trf, trf_available
from .calling import (
    ROUNDING_POLICIES, allele_table, call_all, derive_definitions, long_table, profile_table, read_bins,
    status_table, tandem_count_table,
)
from .extract import (
    extract_amplicons, load_fasta, primers_near_contig_ends, read_hits_table, read_samples, write_hits_table,
    write_marker_fastas, write_samples,
)
from .log import configure_logging, get_logger
from .markers import check_against_sts, markers_from_ids, merge_definitions, read_markers, write_markers
from .models import MarkerDef, PrimerPair
from .search import MerPcrSearcher, SearchError, sample_name
from .sts import is_degenerate, read_sts, write_invalid, write_sts
from .summarize import write_all_summaries, write_table

log = get_logger()

HITS_FILE = "hits.tsv"
SAMPLES_FILE = "samples.txt"
SEQ_DIR = "sequences"
STS_USED = "primers_used.sts"
NEAR_END_FILE = "near_contig_end.tsv"
DERIVED_FILE = "marker_defs_derived.tsv"
DISCOVERED_FILE = "marker_defs_discovered.tsv"
DISCOVER_REPORT = "discover_report.tsv"
RUN_INFO = "run_info.json"


def _add_search_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--sts", required=True,
                   help="STS primer file (id, forward, reverse, expected size[, alias])")
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
                   help="merPCR worker processes (default 1)")
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


def _add_summary_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--prefix", required=True, help="prefix for summary file names")
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
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="search, extract, align and call alleles")
    run.add_argument("--outdir", required=True)
    _add_search_args(run)
    _add_align_args(run)
    _add_summary_args(run)

    search = sub.add_parser("search", help="search and extract amplicons only")
    search.add_argument("--outdir", required=True)
    _add_search_args(search)
    _add_align_args(search)
    search.add_argument("--clwclean", action="store_true", help="also align sequences without N")

    summ = sub.add_parser("summarize", help="rebuild tables and allele calls from hits.tsv")
    summ.add_argument("--outdir", required=True, help="directory containing hits.tsv and samples.txt")
    _add_summary_args(summ)

    disc = sub.add_parser("discover", help="find repeat units in the amplicons of one sample")
    disc.add_argument("--outdir", required=True, help="directory containing hits.tsv")
    disc.add_argument("--sample", help="sample to use (default: the only sample)")
    disc.add_argument("--min-unit", type=int, default=3, help="smallest repeat unit length (default 3)")
    disc.add_argument("--max-unit", type=int, default=100, help="largest repeat unit length (default 100)")
    disc.add_argument("--min-copies", type=float, default=2.0,
                      help="smallest array, in copies, that counts as a repeat (default 2)")
    disc.add_argument("--unit-mismatches", type=int, default=1,
                      help="mismatches tolerated per repeat unit, at most one per six bases (default 1)")
    disc.add_argument("--min-score", type=int, default=30,
                      help="smallest array score, +2 per matching and -7 per differing position between "
                           "consecutive copies (default 30; Tandem Repeats Finder uses 50 for genomes)")
    disc.add_argument("--engine", choices=["auto", "builtin", "trf"], default="auto",
                      help="repeat finder: Tandem Repeats Finder (trf) when on PATH, else the built-in "
                           "self-comparison (default auto)")
    disc.add_argument("--output", default=DISCOVERED_FILE,
                      help=f"marker definition file name within outdir (default {DISCOVERED_FILE})")
    return parser


def expand_fasta(patterns: list[str]) -> list[Path]:
    files: list[Path] = []
    for pat in patterns:
        matches = sorted(glob.glob(pat))
        if not matches and Path(pat).exists():
            matches = [pat]
        files.extend(Path(m) for m in matches)
    return files


def resolve_samples(patterns: list[str]) -> list[Path]:
    """Expand FASTA globs and check that sample names are unique."""
    fasta_files = expand_fasta(patterns)
    if not fasta_files:
        raise SystemExit(f"no FASTA files match {patterns}")
    names = [sample_name(f) for f in fasta_files]
    if len(set(names)) != len(names):
        raise SystemExit("sample names (FASTA file stems) are not unique")
    return fasta_files


def prepare_sts(sts_path: str, outdir: Path, skip_invalid: bool, iupac: bool) -> tuple[Path, list[PrimerPair]]:
    """Validate the STS file and write the accepted entries to outdir/primers_used.sts."""
    valid, invalid = read_sts(sts_path)
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


def search_samples(searcher: MerPcrSearcher, fasta_files: list[Path], pairs: list[PrimerPair], margin: int,
                   dedup: bool = True) -> tuple[list, set[tuple[str, str]], int]:
    """Search and extract amplicons for every sample.

    Returns (amplicons, (sample, marker) pairs with a primer near a contig end
    among the markers without a hit, failed sample count).
    """
    primers = {p.id: p for p in pairs}
    amplicons = []
    near_end: set[tuple[str, str]] = set()
    failed = 0
    for fasta in fasta_files:
        try:
            hits = searcher.search_fasta(fasta)
        except SearchError as exc:
            log.error("%s: %s", fasta, exc)
            failed += 1
            continue
        contigs = load_fasta(fasta)
        amps = extract_amplicons(hits, contigs, primers, deduplicate=dedup)
        amplicons.extend(amps)
        found = {a.hit.primer_id for a in amps}
        missing = [p for p in pairs if p.id not in found]
        if missing:
            window = max(p.expected_size for p in missing) + margin
            sample = sample_name(fasta)
            near_end.update((sample, m) for m in primers_near_contig_ends(contigs, missing, window))
    return amplicons, near_end, failed


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
    info[stage] = {
        "timestamp": datetime.now(UTC).isoformat(timespec="seconds"),
        "parameters": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()
                       if k not in {"command", "verbose"}},
    }
    path.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def do_search(args: argparse.Namespace, outdir: Path) -> None:
    fasta_files = resolve_samples(args.fasta)
    clean_sts, pairs = prepare_sts(args.sts, outdir, args.skip_invalid, args.iupac)
    searcher = MerPcrSearcher(clean_sts, margin=args.margin, mismatches=args.mismatches,
                              threads=args.threads, iupac_mode=int(args.iupac),
                              wordsize=args.wordsize, one_sided=int(args.one_sided),
                              three_prime_match=args.three_prime_match)
    amplicons, near_end, failed = search_samples(searcher, fasta_files, pairs, args.margin, args.dedup)

    write_samples([sample_name(f) for f in fasta_files], outdir / SAMPLES_FILE)
    write_hits_table(amplicons, outdir / HITS_FILE)
    write_near_end(near_end, outdir / NEAR_END_FILE)
    seq_dir = outdir / SEQ_DIR
    written = write_marker_fastas(amplicons, seq_dir)
    log.info("%d amplicons for %d markers written to %s", len(amplicons), len(written), seq_dir)
    if near_end:
        log.info("%d missing sample/marker combinations have a primer near a contig end", len(near_end))
    log.info("Processed %d samples, %d failed", len(fasta_files) - failed, failed)
    write_run_info(outdir, args, "search")

    aligner = detect_aligner(args.align)
    if aligner is not None:
        align_all(seq_dir, aligner, clwclean=getattr(args, "clwclean", False))
    if failed:
        raise SystemExit(f"{failed} samples failed")


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
        defs = derive_definitions(amplicons, defs, args.reference_sample, args.unit_mismatches)
        write_markers(defs, outdir / DERIVED_FILE)
        log.info("Derived marker definitions written to %s", outdir / DERIVED_FILE)
    return defs


def write_calls(calls: list, samples: list[str], markers: list[str], outdir: Path, prefix: str) -> None:
    write_table(allele_table(calls, samples, markers), outdir / f"{prefix}_alleles.tsv")
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
    legacy_defs = defs if (args.markers or args.reference_sample) else None
    write_all_summaries(amplicons, samples, outdir, args.prefix, legacy_defs, args.clwclean, args.legacy_tables)

    if not any(d.can_call or d.has_repeats for d in defs):
        log.info("No marker has reference size and repeats or repeat patterns; allele calling skipped")
        write_run_info(outdir, args, "summarize")
        return
    if not any(d.can_call for d in defs):
        log.warning("No marker has reference size and repeats; only tandem counts are reported. "
                    "Give --markers with ref_size and ref_repeats, or --reference-sample")
    bins = read_bins(args.bins) if args.bins else None
    calls = call_all(amplicons, samples, defs, args.allele_tolerance, args.rounding, bins,
                     args.unit_mismatches, args.max_array_length, read_near_end(outdir / NEAR_END_FILE))
    write_calls(calls, samples, [d.primer_id for d in defs], outdir, args.prefix)
    write_run_info(outdir, args, "summarize")


def do_discover(args: argparse.Namespace, outdir: Path) -> None:
    hits_path = outdir / HITS_FILE
    if not hits_path.exists():
        raise SystemExit(f"{hits_path} not found; run 'mlvatyper search' first")
    amplicons = read_hits_table(hits_path)
    samples = sorted({a.hit.sample for a in amplicons})
    sample = args.sample
    if sample is None:
        if len(samples) != 1:
            raise SystemExit(f"--sample is required when hits.tsv holds several samples: {samples}")
        sample = samples[0]
    elif sample not in samples:
        raise SystemExit(f"sample {sample!r} has no amplicons in {hits_path}")
    engine = args.engine
    if engine == "auto":
        engine = "trf" if trf_available() else "builtin"
        log.info("Repeat finder: %s", engine)
    if engine == "trf":
        found, missing, command = discover_with_trf(amplicons, sample, args.min_score, args.max_unit,
                                                    args.unit_mismatches, args.min_copies)
        args.trf_command = " ".join(command)
    else:
        found, missing = discover_units(amplicons, sample, args.min_unit, args.max_unit, args.min_copies,
                                        args.unit_mismatches, args.min_score)
    args.engine = engine
    write_report(found, missing, outdir / DISCOVER_REPORT)
    defs = discoveries_to_definitions(found, missing, amplicons, sample)
    sts_path = outdir / STS_USED
    if sts_path.exists():
        order = {p.id: n for n, p in enumerate(read_sts(sts_path)[0])}
        defs.sort(key=lambda d: order.get(d.primer_id, len(order)))
    write_markers(defs, outdir / args.output)
    log.info("Marker definitions written to %s, details in %s", outdir / args.output, outdir / DISCOVER_REPORT)
    write_run_info(outdir, args, "discover")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    configure_logging(args.verbose)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        if args.command == "search":
            do_search(args, outdir)
        elif args.command == "summarize":
            do_summarize(args, outdir)
        elif args.command == "discover":
            do_discover(args, outdir)
        elif args.command == "run":
            do_search(args, outdir)
            do_summarize(args, outdir)
    except SystemExit as exc:
        if exc.code not in (None, 0):
            log.error("%s", exc.code)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
