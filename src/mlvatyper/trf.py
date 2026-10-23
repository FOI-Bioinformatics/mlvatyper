"""Repeat unit discovery with Tandem Repeats Finder (TRF).

TRF (Benson 1999, https://tandem.bu.edu/trf/trf.html; bioconda package trf)
is run once on a multi-FASTA of the amplicons of one sample with pinned
alignment parameters and the -ngs output format, which is one "@id" line
per sequence followed by one row per repeat: start, end, period, copy
number, consensus size, percent matches, percent indels, score, A, C, G, T
counts, entropy, consensus unit, repeat sequence and the two flanking
sequences. Coordinates are 1-based and inclusive.

TRF reports the same array at several periods and every array in the
sequence. All rows become candidates that discover.merge_harmonics reduces
to one per locus; select_repeat applies the same rule to one list. Each
consensus unit is then counted with the tandem counter so that ref_repeats
is the same integer that allele calling would see.
"""

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .discover import (
    Discovery, RawArray, Thresholds, chosen_amplicons, finish_candidates, make_discovery, merge_harmonics,
)
from .log import get_logger
from .models import Amplicon

log = get_logger()

# match, mismatch, indel, match probability, indel probability: the TRF defaults
TRF_ALIGNMENT = (2, 7, 7, 80, 10)


@dataclass(frozen=True)
class TrfRepeat:
    start: int  # 1-based inclusive
    end: int
    period: int
    copies: float
    consensus_size: int
    pct_match: int
    pct_indel: int
    score: int
    entropy: float
    consensus: str
    sequence: str


def trf_available() -> bool:
    return shutil.which("trf") is not None


def trf_command(fasta: Path, min_score: int, max_period: int) -> list[str]:
    return ["trf", str(fasta), *(str(x) for x in TRF_ALIGNMENT), str(min_score), str(max_period), "-h", "-ngs"]


def parse_trf_ngs(text: str) -> dict[str, list[TrfRepeat]]:
    """Repeats per sequence id from TRF -ngs output."""
    result: dict[str, list[TrfRepeat]] = {}
    current: list[TrfRepeat] | None = None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("@"):
            current = result.setdefault(line[1:].split()[0], [])
            continue
        if current is None:
            raise ValueError("TRF output has a repeat row before any @sequence line")
        f = line.split()
        if len(f) < 15:
            raise ValueError(f"TRF row with {len(f)} fields: {line[:60]}")
        current.append(TrfRepeat(
            int(f[0]), int(f[1]), int(f[2]), float(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]),
            float(f[12]), f[13].upper(), f[14].upper(),
        ))
    return result


def select_repeat(repeats: list[TrfRepeat], period_slack: float = 0.85) -> TrfRepeat | None:
    if not repeats:
        return None
    best = max(r.score for r in repeats)
    return min((r for r in repeats if r.score >= period_slack * best), key=lambda r: (r.period, -r.score))


def run_trf(fasta: Path, min_score: int, max_period: int) -> str:
    """Run TRF and return its -ngs output. TRF exits with the number of sequences processed."""
    cmd = trf_command(fasta, min_score, max_period)
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=fasta.parent)
    if proc.returncode < 0 or (not proc.stdout and proc.stderr):
        raise RuntimeError(f"trf failed ({proc.returncode}): {proc.stderr.strip()[:200]}")
    return proc.stdout


def trf_raw_arrays(repeats: list[TrfRepeat]) -> list[RawArray]:
    """TRF rows as raw arrays with 0-based half-open coordinates."""
    return [RawArray(r.start - 1, r.end, r.period, r.consensus, r.score, r.copies, r.pct_match, r.pct_indel)
            for r in repeats]


def discover_with_trf(
    amplicons: list[Amplicon], sample: str, min_score: int = 30, max_period: int = 100, max_mismatch: int = 1,
    min_copies: float = 1.8, thresholds: Thresholds | None = None,
    primer_lengths: dict[str, tuple[int, int]] | None = None,
) -> tuple[dict[str, list[Discovery]], list[str], list[str]]:
    """Candidates per marker from TRF for the amplicons (copy 1) of one sample.

    TRF is run at the weak score floor so that weak arrays appear as category
    C candidates. Returns ({marker: ranked candidates}, markers without a
    candidate, the command line used).
    """
    if not trf_available():
        raise FileNotFoundError("trf not found on PATH; install with: conda install -c bioconda trf")
    t = thresholds or Thresholds(min_score=min_score, min_copies=min_copies)
    chosen = chosen_amplicons(amplicons, sample)
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "amplicons.fasta"
        with open(fasta, "w", encoding="utf-8") as handle:
            for marker, amp in chosen.items():
                handle.write(f">{marker}\n{amp.sequence}\n")
        floor = min(t.min_score, t.weak_score)
        command = trf_command(fasta, floor, max_period)
        repeats = parse_trf_ngs(run_trf(fasta, floor, max_period))
    found: dict[str, list[Discovery]] = {}
    missing: list[str] = []
    for marker, amp in chosen.items():
        seq = amp.sequence.upper()
        raw = merge_harmonics(trf_raw_arrays(repeats.get(marker, [])))
        cands = [make_discovery(marker, sample, seq, r, "trf", max_mismatch, (primer_lengths or {}).get(marker))
                 for r in raw if r.copies >= t.min_copies]
        cands = finish_candidates([c for c in cands if c is not None], t)
        if cands:
            found[marker] = cands
        else:
            missing.append(marker)
    log.info("%s: TRF found repeat units for %d markers, none for %d", sample, len(found), len(missing))
    command[1] = "amplicons.fasta"
    return found, missing, command
