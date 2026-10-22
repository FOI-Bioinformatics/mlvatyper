"""Command line interface.

Subcommands:
  run        search, extract, align and summarise in one go
  search     search and extract only (writes hits.tsv and sequences/)
  summarize  rebuild summary tables from an existing hits.tsv
"""

import argparse
import glob
import sys
from pathlib import Path

from . import __version__
from .align import align_all, detect_aligner
from .extract import (
    extract_amplicons, load_fasta, read_hits_table, read_samples, write_hits_table,
    write_marker_fastas, write_samples,
)
from .log import configure_logging, get_logger
from .patterns import read_patterns
from .search import MerPcrSearcher, SearchError, sample_name
from .sts import is_degenerate, read_sts, write_invalid, write_sts
from .summarize import write_all_summaries

log = get_logger()

HITS_FILE = "hits.tsv"
SAMPLES_FILE = "samples.txt"
SEQ_DIR = "sequences"


def _add_search_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--sts", required=True,
                   help="STS primer file (id, forward, reverse, expected size[, alias])")
    p.add_argument("--fasta", required=True, nargs="+",
                   help="FASTA file(s) or glob pattern(s), one file per sample")
    p.add_argument("-M", "--margin", type=int, default=500,
                   help="allowed deviation from the expected product size in bp (default 500)")
    p.add_argument("-N", "--mismatches", type=int, default=3,
                   help="mismatches allowed per primer (default 3)")
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
    p.add_argument("--skip-invalid", action="store_true",
                   help="skip invalid STS entries instead of aborting; they are listed in invalid_primers.tsv")


def _add_summary_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--prefix", required=True, help="prefix for summary file names")
    p.add_argument("--pattern", help="pattern file with repeat units per marker")
    p.add_argument("--clwclean", action="store_true",
                   help="also produce outputs restricted to sequences without N")


def _add_align_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--align", default="auto", choices=["auto", "none", "mafft", "muscle"],
                   help="aligner for per-marker alignments (default auto)")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mlvatyper", description="Locate VNTR amplicons with merPCR and summarise repeat lengths"
    )
    parser.add_argument("--version", action="version", version=f"mlvatyper {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="search, extract, align and summarise")
    run.add_argument("--outdir", required=True)
    _add_search_args(run)
    _add_align_args(run)
    _add_summary_args(run)

    search = sub.add_parser("search", help="search and extract amplicons only")
    search.add_argument("--outdir", required=True)
    _add_search_args(search)
    _add_align_args(search)
    search.add_argument("--clwclean", action="store_true", help="also align sequences without N")

    summ = sub.add_parser("summarize", help="rebuild summary tables from hits.tsv")
    summ.add_argument("--outdir", required=True, help="directory containing hits.tsv and samples.txt")
    _add_summary_args(summ)
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


def prepare_sts(sts_path: str, outdir: Path, skip_invalid: bool, iupac: bool) -> Path:
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
    clean_sts = outdir / "primers_used.sts"
    write_sts(valid, clean_sts)
    return clean_sts


def search_samples(searcher: MerPcrSearcher, fasta_files: list[Path]) -> tuple[list, int]:
    """Search and extract amplicons for every sample. Returns (amplicons, failed sample count)."""
    amplicons = []
    failed = 0
    for fasta in fasta_files:
        try:
            hits = searcher.search_fasta(fasta)
        except SearchError as exc:
            log.error("%s: %s", fasta, exc)
            failed += 1
            continue
        amplicons.extend(extract_amplicons(hits, load_fasta(fasta)))
    return amplicons, failed


def do_search(args: argparse.Namespace, outdir: Path) -> None:
    fasta_files = resolve_samples(args.fasta)
    clean_sts = prepare_sts(args.sts, outdir, args.skip_invalid, args.iupac)
    searcher = MerPcrSearcher(clean_sts, margin=args.margin, mismatches=args.mismatches,
                              threads=args.threads, iupac_mode=int(args.iupac),
                              wordsize=args.wordsize, one_sided=int(args.one_sided))
    amplicons, failed = search_samples(searcher, fasta_files)

    write_samples([sample_name(f) for f in fasta_files], outdir / SAMPLES_FILE)
    write_hits_table(amplicons, outdir / HITS_FILE)
    seq_dir = outdir / SEQ_DIR
    written = write_marker_fastas(amplicons, seq_dir)
    log.info("%d amplicons for %d markers written to %s", len(amplicons), len(written), seq_dir)
    log.info("Processed %d samples, %d failed", len(fasta_files) - failed, failed)

    aligner = detect_aligner(args.align)
    if aligner is not None:
        align_all(seq_dir, aligner, clwclean=getattr(args, "clwclean", False))
    if failed:
        raise SystemExit(f"{failed} samples failed")


def do_summarize(args: argparse.Namespace, outdir: Path) -> None:
    hits_path = outdir / HITS_FILE
    if not hits_path.exists():
        raise SystemExit(f"{hits_path} not found; run 'mlvatyper search' first")
    amplicons = read_hits_table(hits_path)
    samples_path = outdir / SAMPLES_FILE
    samples = read_samples(samples_path) if samples_path.exists() else sorted({a.hit.sample for a in amplicons})
    defs = read_patterns(args.pattern) if args.pattern else None
    for p in write_all_summaries(amplicons, samples, outdir, args.prefix, defs, args.clwclean):
        log.debug("summary: %s", p)


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
