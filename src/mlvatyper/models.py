"""Data classes shared across the package."""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class PrimerPair:
    """One line of an STS file."""

    id: str
    forward: str
    reverse: str
    expected_size: int
    alias: str | None = None
    line_number: int = 0


@dataclass(frozen=True)
class Hit:
    """One merPCR hit, with 1-based inclusive coordinates on the contig."""

    sample: str
    contig: str
    start: int
    end: int
    strand: str
    primer_id: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class Amplicon:
    """An extracted amplicon. copy is the 1-based index of the hit within a sample for a marker."""

    hit: Hit
    copy: int
    sequence: str

    @property
    def record_id(self) -> str:
        return f"{self.hit.sample}_{self.copy}"

    @property
    def has_n(self) -> bool:
        return "N" in self.sequence.upper()


@dataclass(frozen=True)
class PatternDef:
    """Repeat unit definition for one marker from the pattern file."""

    primer_id: str
    patterns: tuple[str, ...] = field(default_factory=tuple)
    pattern_length: int = 0

    @property
    def has_repeats(self) -> bool:
        return bool(self.patterns)
