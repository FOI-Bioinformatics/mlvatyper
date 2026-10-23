"""Scheme bundles: primers, marker definitions, bins and metadata in one directory.

A scheme directory holds

  primers.sts   the STS file
  markers.tsv   marker definitions (optional; MLVAbank style ids need none)
  bins.tsv      bin table (optional)
  scheme.yml    flat key: value metadata (optional)

Recognised metadata keys: name, organism, reference (strain or accession
the reference values come from), citation, url, retrieved, rounding,
allele_tolerance, unit_mismatches, margin, mismatches, max_array_length,
notes. The calling keys become the defaults of the corresponding options
when the scheme is used; options given on the command line still win.

citation, url and reference are the provenance of the values: who
published them, where the file came from, and which strain or accession
the reference sizes and repeat numbers describe. Allele numbers are only
comparable with a database when these are known, so check-scheme refuses
a scheme that leaves them empty or fills them with a placeholder.

Bundled schemes live under mlvatyper/schemes/<name>. `--scheme` accepts a
bundled name or a directory path.
"""

from dataclasses import dataclass, field
from pathlib import Path

from .log import get_logger

log = get_logger()

SCHEME_DIR = Path(__file__).parent / "schemes"
STS_FILE = "primers.sts"
MARKERS_FILE = "markers.tsv"
BINS_FILE = "bins.tsv"
META_FILE = "scheme.yml"
NUMERIC_KEYS = {"allele_tolerance": float, "unit_mismatches": int, "margin": int, "mismatches": int,
                "max_array_length": int}
PROVENANCE_FIELDS = ("citation", "url", "reference")
RECOMMENDED_FIELDS = ("organism", "retrieved")
# values that state nothing; compared against the whole value, lower case and stripped of trailing punctuation
PLACEHOLDERS = frozenset({"", "-", "?", "na", "n/a", "none", "unknown", "tbd", "todo", "fixme", "xxx",
                          "to be added", "to be determined", "see publication", "see paper"})


class SchemeError(ValueError):
    pass


@dataclass(frozen=True)
class Scheme:
    name: str
    path: Path
    sts: Path
    markers: Path | None
    bins: Path | None
    meta: dict = field(default_factory=dict)

    @property
    def description(self) -> str:
        organism = self.meta.get("organism", "")
        return f"{self.name}: {organism}" if organism else self.name

    @property
    def provenance(self) -> dict:
        return {k: self.meta.get(k, "") for k in PROVENANCE_FIELDS}

    def missing_provenance(self) -> list[str]:
        """Provenance fields that are absent or hold a placeholder."""
        return [k for k, v in self.provenance.items() if is_placeholder(v)]

    def missing_recommended(self) -> list[str]:
        return [k for k in RECOMMENDED_FIELDS if is_placeholder(self.meta.get(k, ""))]

    def calling_defaults(self) -> dict:
        """Metadata values that apply as option defaults."""
        keys = ("rounding", *NUMERIC_KEYS)
        return {k: self.meta[k] for k in keys if k in self.meta}


def is_placeholder(value) -> bool:
    """Whether a metadata value states nothing: empty, or one of the placeholder words."""
    return str(value).strip().rstrip(".!,;:").lower() in PLACEHOLDERS


def parse_flat_yaml(text: str) -> dict:
    """key: value lines; quotes around the value are removed, numbers are typed by key."""
    meta: dict = {}
    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.split("#", 1)[0].rstrip() if not raw.lstrip().startswith("#") else ""
        if not line.strip():
            continue
        if ":" not in line:
            raise SchemeError(f"scheme.yml line {number}: expected 'key: value'")
        key, value = line.split(":", 1)
        key = key.strip().lower()
        value = value.strip().strip("'\"")
        if key in NUMERIC_KEYS and value != "":
            try:
                value = NUMERIC_KEYS[key](value)
            except ValueError as exc:
                raise SchemeError(f"scheme.yml line {number}: {key} must be a number") from exc
        meta[key] = value
    return meta


def list_schemes() -> list[Scheme]:
    if not SCHEME_DIR.exists():
        return []
    return [load_scheme(p) for p in sorted(SCHEME_DIR.iterdir()) if (p / STS_FILE).exists()]


def load_scheme(name_or_path: str | Path) -> Scheme:
    """A bundled scheme by name, or any directory with a primers.sts file."""
    path = Path(name_or_path)
    if not path.is_dir():
        path = SCHEME_DIR / str(name_or_path)
    if not (path / STS_FILE).exists():
        names = ", ".join(p.name for p in list_schemes()) or "none"
        raise SchemeError(f"no scheme at {name_or_path}: a scheme directory needs {STS_FILE} "
                          f"(bundled schemes: {names})")
    meta = parse_flat_yaml((path / META_FILE).read_text(encoding="utf-8")) if (path / META_FILE).exists() else {}
    markers = path / MARKERS_FILE if (path / MARKERS_FILE).exists() else None
    bins = path / BINS_FILE if (path / BINS_FILE).exists() else None
    return Scheme(meta.get("name") or path.name, path, path / STS_FILE, markers, bins, meta)


def apply_scheme(args, scheme: Scheme, parser_defaults: dict) -> None:
    """Fill --sts, --markers and --bins from the scheme and apply its calling defaults.

    An option keeps its command line value when it differs from the parser
    default; otherwise the scheme value is used.
    """
    if getattr(args, "sts", None) is None:
        args.sts = str(scheme.sts)
    if getattr(args, "markers", None) is None and scheme.markers is not None:
        args.markers = str(scheme.markers)
    if getattr(args, "bins", None) is None and scheme.bins is not None:
        args.bins = str(scheme.bins)
    applied = []
    for key, value in scheme.calling_defaults().items():
        if hasattr(args, key) and getattr(args, key) == parser_defaults.get(key):
            setattr(args, key, value)
            applied.append(f"{key}={value}")
    log.info("Scheme %s: %s%s", scheme.name, scheme.description,
             f" (defaults {', '.join(applied)})" if applied else "")
