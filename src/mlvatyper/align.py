"""Optional multiple sequence alignment of amplicons per marker.

An external aligner (mafft or muscle) is used when available. CLUSTAL output
is written here with the full sequence names, which removes the dependency
on EMBOSS seqret.
"""

import shutil
import subprocess
from pathlib import Path

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

from .log import get_logger

log = get_logger()

ALIGNERS = ("mafft", "muscle")


def detect_aligner(choice: str = "auto") -> str | None:
    """Return the aligner name to use, or None."""
    if choice == "none":
        return None
    if choice == "auto":
        for name in ALIGNERS:
            if shutil.which(name):
                log.info("Using aligner %s", name)
                return name
        log.warning("No aligner (mafft or muscle) found on PATH; alignments are skipped")
        return None
    if shutil.which(choice) is None:
        raise FileNotFoundError(f"aligner {choice!r} not found on PATH")
    return choice


def muscle_major_version() -> int:
    out = subprocess.run(["muscle", "-version"], capture_output=True, text=True)
    text = (out.stdout + out.stderr).lower()
    return 5 if "muscle 5" in text or "muscle v5" in text else 3


def run_aligner(aligner: str, in_fasta: Path, out_fasta: Path) -> None:
    if aligner == "mafft":
        with open(out_fasta, "w", encoding="utf-8") as handle:
            subprocess.run(["mafft", "--quiet", "--auto", str(in_fasta)], stdout=handle, check=True, text=True)
    elif aligner == "muscle":
        if muscle_major_version() >= 5:
            cmd = ["muscle", "-align", str(in_fasta), "-output", str(out_fasta)]
        else:
            cmd = ["muscle", "-in", str(in_fasta), "-out", str(out_fasta), "-quiet"]
        subprocess.run(cmd, capture_output=True, check=True, text=True)
    else:
        raise ValueError(f"unknown aligner {aligner}")


CLUSTAL_BLOCK = 60


def fasta_to_clustal(aligned_fasta: Path, clustal_path: Path) -> None:
    """Write the alignment in CLUSTAL format with the full sequence names.

    Biopython's CLUSTAL writer cuts names to 30 characters, which makes
    sample_1 and sample_2 of a long assembly name the same row; here the
    name column is as wide as the longest name.
    """
    records = list(SeqIO.parse(str(aligned_fasta), "fasta"))
    MultipleSeqAlignment(records)  # raises ValueError when the sequences differ in length
    width = max((len(r.id) for r in records), default=0) + 6
    length = len(records[0].seq) if records else 0
    lines = ["CLUSTAL multiple sequence alignment", "", ""]
    for start in range(0, length, CLUSTAL_BLOCK):
        lines.extend(f"{r.id.ljust(width)}{str(r.seq)[start:start + CLUSTAL_BLOCK]}" for r in records)
        lines.append("")
    Path(clustal_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def align_marker(aligner: str | None, in_fasta: Path, out_fasta: Path, clustal_path: Path) -> bool:
    """Align one marker FASTA. Single sequences are copied unaligned. Returns True on success."""
    records = list(SeqIO.parse(str(in_fasta), "fasta"))
    if not records:
        return False
    if len(records) == 1 or aligner is None:
        if len(records) > 1:
            log.debug("%s: no aligner, writing unaligned CLUSTAL is not possible, skipped", in_fasta.name)
            return False
        shutil.copy(in_fasta, out_fasta)
    else:
        try:
            run_aligner(aligner, in_fasta, out_fasta)
        except subprocess.CalledProcessError as exc:
            log.error("%s: %s failed: %s", in_fasta.name, aligner, exc.stderr)
            return False
    fasta_to_clustal(out_fasta, clustal_path)
    return True


def filter_n(in_fasta: Path, out_fasta: Path) -> tuple[int, int]:
    """Write records without N to out_fasta. Returns (kept, removed)."""
    records = list(SeqIO.parse(str(in_fasta), "fasta"))
    clean = [r for r in records if "N" not in str(r.seq).upper()]
    if clean:
        SeqIO.write(clean, str(out_fasta), "fasta")
    return len(clean), len(records) - len(clean)


def align_all(seq_dir: Path, aligner: str | None, clwclean: bool = False) -> None:
    for in_fasta in sorted(seq_dir.glob("*.fasta")):
        if in_fasta.stem.endswith(("_aligned", "_clean", "_clean_aligned")):
            continue
        marker = in_fasta.stem
        align_marker(aligner, in_fasta, seq_dir / f"{marker}_aligned.fasta", seq_dir / f"{marker}_aligned.clw")
        if clwclean:
            clean_fasta = seq_dir / f"{marker}_clean.fasta"
            kept, removed = filter_n(in_fasta, clean_fasta)
            if kept:
                log.info("%s: %d sequences with N removed", marker, removed)
                align_marker(aligner, clean_fasta,
                             seq_dir / f"{marker}_clean_aligned.fasta", seq_dir / f"{marker}_clean_aligned.clw")
            else:
                log.info("%s: no sequences without N", marker)
