"""Repeat unit discovery with Tandem Repeats Finder (TRF).

TRF (Benson 1999, https://tandem.bu.edu/trf/trf.html; bioconda package trf)
is run once on a multi-FASTA of the amplicons of one sample with pinned
alignment parameters and the -ngs output format, which is one "@id" line
per sequence followed by one row per repeat: start, end, period, copy
number, consensus size, percent matches, percent indels, score, A, C, G, T
counts, entropy, consensus unit, repeat sequence and the two flanking
sequences. Coordinates are 1-based and inclusive.

TRF reports the same array at several periods. The repeat with the highest
score is chosen, and among repeats scoring within period_slack of it the
smallest period, so that a unit is not reported as its own multiple. The
consensus unit is then counted with the tandem counter so that ref_repeats
is the same integer that allele calling would see.
"""

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .discover import Discovery
from .log import get_logger
from .models import Amplicon
from .tandem import count_tandem_array

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


def discover_with_trf(
    amplicons: list[Amplicon], sample: str, min_score: int = 30, max_period: int = 100, max_mismatch: int = 1,
    min_copies: float = 2.0,
) -> tuple[list[Discovery], list[str], list[str]]:
    """Discoveries from TRF for the amplicons (copy 1) of one sample.

    Returns (found, markers without a repeat, the command line used).
    """
    if not trf_available():
        raise FileNotFoundError("trf not found on PATH; install with: conda install -c bioconda trf")
    chosen: dict[str, Amplicon] = {}
    for amp in sorted((a for a in amplicons if a.hit.sample == sample), key=lambda a: (a.hit.primer_id, a.copy)):
        chosen.setdefault(amp.hit.primer_id, amp)
    with tempfile.TemporaryDirectory() as tmp:
        fasta = Path(tmp) / "amplicons.fasta"
        with open(fasta, "w", encoding="utf-8") as handle:
            for marker, amp in chosen.items():
                handle.write(f">{marker}\n{amp.sequence}\n")
        command = trf_command(fasta, min_score, max_period)
        repeats = parse_trf_ngs(run_trf(fasta, min_score, max_period))
    found: list[Discovery] = []
    missing: list[str] = []
    for marker, amp in chosen.items():
        r = select_repeat(repeats.get(marker, []))
        if r is None or r.copies < min_copies:
            missing.append(marker)
            continue
        array = count_tandem_array(amp.sequence, (r.consensus,), max_mismatch)
        if array.copies == 0:
            missing.append(marker)
            continue
        found.append(Discovery(
            marker, sample, len(amp.sequence), r.consensus, array, array.start, len(amp.sequence) - array.end,
            r.score, "trf", r.copies, r.pct_match, r.pct_indel,
        ))
    log.info("%s: TRF found repeat units for %d markers, none for %d", sample, len(found), len(missing))
    command[1] = "amplicons.fasta"
    return found, missing, command
