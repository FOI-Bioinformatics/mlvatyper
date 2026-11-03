"""Searching one sample, running many samples in parallel, and the per-sample results that --resume reads.

A sample is searched with merPCR, then the rescue passes fill in the markers the main pass missed. With
jobs > 1 the samples run one process each and are merged in input order, so the output is identical to a
serial run. Each finished sample is written to the cache directory as its hits, its missing markers and a
done file holding the input file as it was on disk and a hash of the search settings; --resume reads a
sample back when both still match and searches only the rest.
"""

import argparse
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from . import __version__
from .extract import (
    extract_amplicons, load_fasta, primer_end_detail, primers_near_contig_ends, read_hits_table,
    write_hits_table,
)
from .log import get_logger
from .models import PrimerPair
from .rescue import (
    RescueOptions, find_with_indels, missing_reasons, read_missing, rescue_pass, write_missing,
)
from .search import MerPcrSearcher, SearchError, merpcr_version, sample_name

log = get_logger()


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


PER_SAMPLE_DIR = "per_sample"


def _fingerprint(fasta: Path, params_key: str) -> dict:
    """What a per-sample result depends on: the input file as it is on disk and the search settings."""
    stat = Path(fasta).stat()
    return {"fasta": str(Path(fasta).resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "params": params_key}


def _sample_files(cache: Path, sample: str) -> tuple[Path, Path, Path]:
    return cache / f"{sample}.hits.tsv", cache / f"{sample}.missing.tsv", cache / f"{sample}.done.json"


def save_sample(cache: Path, fasta: Path, amps: list, reasons: dict[str, str], params_key: str) -> None:
    """Write the result of one sample; the done file goes last, so a half written sample is never reused."""
    sample = sample_name(fasta)
    hits, missing, done = _sample_files(cache, sample)
    done.unlink(missing_ok=True)
    write_hits_table(amps, hits)
    write_missing({(sample, m): r for m, r in reasons.items()}, missing)
    done.write_text(json.dumps(_fingerprint(fasta, params_key)), encoding="utf-8")


def load_sample(cache: Path, fasta: Path, params_key: str):
    """The saved result of a sample as (amplicons, reasons), or None when it is missing or out of date."""
    sample = sample_name(fasta)
    hits, missing, done = _sample_files(cache, sample)
    if not (hits.exists() and missing.exists() and done.exists()):
        return None
    try:
        if json.loads(done.read_text(encoding="utf-8")) != _fingerprint(fasta, params_key):
            return None
        amps = read_hits_table(hits)
        reasons = {m: r for (_, m), r in read_missing(missing).items()}
    except (OSError, ValueError, KeyError):
        return None
    return amps, reasons


def search_key(sts_path: Path, args: argparse.Namespace) -> str:
    """Hash of everything besides the input file that changes a search result."""
    settings = {k: getattr(args, k, None) for k in (
        "margin", "mismatches", "wordsize", "one_sided", "iupac", "three_prime_match", "dedup", "rescue_wordsize",
        "rescue_margin", "rescue_indels", "max_repeats")}
    digest = hashlib.md5(Path(sts_path).read_bytes())
    digest.update(json.dumps(settings, sort_keys=True, default=str).encode())
    digest.update(f"{__version__} {merpcr_version()}".encode())
    return digest.hexdigest()


def search_samples(searcher: MerPcrSearcher, fasta_files: list[Path], pairs: list[PrimerPair], margin: int,
                   dedup: bool = True, rescue: RescueOptions | None = None, wordsize: int = 11,
                   one_sided: bool = False, jobs: int = 1, sts_path: str | Path | None = None,
                   searcher_kwargs: dict | None = None, cache: Path | None = None, resume: bool = False,
                   params_key: str = "") -> tuple[list, dict[tuple[str, str], str], list[Path]]:
    """Search and extract amplicons for every sample.

    With jobs > 1 the samples are searched in parallel, one process each, and
    the results are merged in input order so the outputs are identical to the
    serial run. Each process builds its own searcher from sts_path and
    searcher_kwargs, so those are required when jobs > 1. Returns (amplicons,
    reason per (sample, marker) without a hit, FASTA files that failed).

    With a cache directory each finished sample is written there as it completes (its hits, its missing
    markers and a done file that records the input file and the settings). With resume, a sample whose done
    file matches is read back instead of searched, so a failed or interrupted run can be repeated for the
    samples that are not finished.
    """
    rescue = rescue or RescueOptions(wordsize=0)
    params = (pairs, margin, dedup, rescue, wordsize, one_sided)
    if cache is not None:
        cache.mkdir(parents=True, exist_ok=True)
    done = _saved_samples(cache, fasta_files, params_key) if resume and cache is not None else {}
    if done:
        log.info("Resuming: %d of %d samples are read from %s", len(done), len(fasta_files), cache)
    todo = [f for f in fasta_files if f not in done]

    def keep(fasta, amps, sample_reasons) -> None:
        if cache is not None and amps is not None:
            save_sample(cache, fasta, amps, sample_reasons, params_key)
        done[fasta] = (fasta, amps, sample_reasons)

    workers = max(1, min(jobs, len(todo))) if todo else 1
    if workers == 1:
        for fasta in todo:
            try:
                amps, sample_reasons = search_sample(searcher, fasta, pairs, margin, dedup, rescue,
                                                     wordsize, one_sided)
            except SAMPLE_ERRORS as exc:
                keep(fasta, None, str(exc).splitlines()[0] if str(exc) else type(exc).__name__)
                continue
            keep(fasta, amps, sample_reasons)
    else:
        _search_parallel(todo, workers, sts_path, searcher_kwargs, params, keep)
    return _merge_results([done[fasta] for fasta in fasta_files], workers > 1)  # input order, deterministic


def _saved_samples(cache: Path, fasta_files: list[Path], params_key: str) -> dict[Path, tuple]:
    done = {}
    for fasta in fasta_files:
        saved = load_sample(cache, fasta, params_key)
        if saved is not None:
            done[fasta] = (fasta, saved[0], saved[1])
    return done


def _search_parallel(todo, workers, sts_path, searcher_kwargs, params, keep) -> None:
    if sts_path is None or searcher_kwargs is None:
        raise ValueError("parallel search needs sts_path and searcher_kwargs")
    log.info("Searching %d samples with %d parallel jobs", len(todo), workers)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(str(sts_path), searcher_kwargs, params)) as pool:
        futures = [pool.submit(_search_one, fasta) for fasta in todo]
        for future in as_completed(futures):  # kept as each sample finishes, so a crash loses only the running ones
            keep(*future.result())


def _merge_results(results: list[tuple], log_counts: bool):
    amplicons = []
    reasons: dict[tuple[str, str], str] = {}
    failed = []
    for fasta, amps, sample_reasons in results:
        if amps is None:
            log.error("%s: %s", fasta, sample_reasons)
            failed.append(fasta)
            continue
        if log_counts:
            log.info("%s: %d amplicons", sample_name(fasta), len(amps))
        amplicons.extend(amps)
        reasons.update({(sample_name(fasta), m): r for m, r in sample_reasons.items()})
    return amplicons, reasons, failed
