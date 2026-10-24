"""Reference assemblies by NCBI accession, downloaded once into a local cache.

A scheme names the assembly its reference values describe by a GenBank or
RefSeq assembly accession (GCA_ or GCF_ followed by nine digits and a
version). The sequence is fetched from the NCBI Datasets API as a zip of
FASTA files and stored as one FASTA under the cache directory, which is
MLVATYPER_CACHE when set, otherwise ~/.cache/mlvatyper/assemblies. The
assembly itself is never kept in a scheme directory, so the repository
holds accessions and not sequences.
"""

import io
import os
import re
import urllib.request
import zipfile
from pathlib import Path

from .log import get_logger

log = get_logger()

ASSEMBLY_RE = re.compile(r"^GC[AF]_\d{9}\.\d+$")
NUCLEOTIDE_RE = re.compile(r"^[A-Z]{1,2}_?\d{5,8}\.\d+$")
DATASETS_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
DEFAULT_CACHE = Path.home() / ".cache" / "mlvatyper" / "assemblies"


def is_assembly_accession(text: str) -> bool:
    return bool(ASSEMBLY_RE.match(text.strip()))


def cache_dir(override: str | Path | None = None) -> Path:
    if override:
        return Path(override)
    return Path(os.environ.get("MLVATYPER_CACHE") or DEFAULT_CACHE)


def cached_path(accession: str, directory: str | Path | None = None) -> Path:
    return cache_dir(directory) / f"{accession}.fasta"


def _fetch_zip(accession: str) -> bytes:
    url = DATASETS_URL.format(accession=accession) + "?include_annotation_type=GENOME_FASTA"
    request = urllib.request.Request(url, headers={"User-Agent": "mlvatyper", "Accept": "application/zip"})
    with urllib.request.urlopen(request, timeout=300) as response:  # noqa: S310 - fixed https host
        return response.read()


def fasta_from_zip(data: bytes, accession: str) -> str:
    """Concatenate the genomic FASTA files of a Datasets zip."""
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        names = sorted(n for n in archive.namelist()
                       if n.endswith((".fna", ".fa", ".fasta")) and f"/data/{accession}/" in n)
        if not names:
            raise RuntimeError(f"the NCBI Datasets archive for {accession} holds no genomic FASTA")
        parts = [archive.read(n).decode("utf-8") for n in names]
    text = "".join(p if p.endswith("\n") else p + "\n" for p in parts)
    if not text.startswith(">"):
        raise RuntimeError(f"the FASTA for {accession} does not start with a record header")
    return text


def download_assembly(accession: str, directory: str | Path | None = None, force: bool = False) -> Path:
    """FASTA path of the assembly, downloaded unless already in the cache."""
    accession = accession.strip()
    if not is_assembly_accession(accession):
        raise ValueError(f"{accession!r} is not an assembly accession (GCA_ or GCF_, nine digits, version)")
    path = cached_path(accession, directory)
    if path.exists() and not force:
        log.info("Using cached assembly %s", path)
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    log.info("Downloading %s from NCBI Datasets", accession)
    text = fasta_from_zip(_fetch_zip(accession), accession)
    tmp = path.with_suffix(".part")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)
    log.info("Stored %s (%d bp in %d records)", path, sum(len(line) for line in text.splitlines()
             if not line.startswith(">")), text.count(">"))
    return path
