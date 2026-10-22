"""Transparent reading of plain or compressed text files.

Files ending in .gz, .bz2 or .xz are decompressed on the fly. The suffix is
the only signal used, which matches the behaviour of merPCR.
"""

import bz2
import gzip
import lzma
from pathlib import Path
from typing import IO

COMPRESSION_SUFFIXES = {".gz": gzip.open, ".bz2": bz2.open, ".xz": lzma.open}


def open_text(path: str | Path) -> IO[str]:
    """Open a text file for reading, decompressing by suffix if needed."""
    path = Path(path)
    opener = COMPRESSION_SUFFIXES.get(path.suffix.lower())
    if opener is None:
        return open(path, encoding="utf-8")
    return opener(path, "rt", encoding="utf-8")


def strip_compression(path: str | Path) -> Path:
    """Return the path without a trailing compression suffix."""
    path = Path(path)
    if path.suffix.lower() in COMPRESSION_SUFFIXES:
        return path.with_suffix("")
    return path
