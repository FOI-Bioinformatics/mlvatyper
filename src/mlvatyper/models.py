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
    size_low: int | None = None  # set when the STS size was given as a range
    size_high: int | None = None

    @property
    def size_text(self) -> str:
        if self.size_low is not None and self.size_high is not None:
            return f"{self.size_low}-{self.size_high}"
        return str(self.expected_size)


@dataclass(frozen=True)
class Hit:
    """One merPCR hit, with 1-based inclusive coordinates on the contig.

    expected_size is the product size given in the STS file. Mismatch counts
    and the distance to the nearest contig end are filled in during
    extraction; -1 means unknown. n_alternatives is the number of overlapping
    hits for the same locus that were discarded in favour of this one.
    search_pass records which search found the hit.
    """

    sample: str
    contig: str
    start: int
    end: int
    strand: str
    primer_id: str
    expected_size: int = -1
    fwd_mismatches: int = -1
    rev_mismatches: int = -1
    dist_to_contig_end: int = -1
    n_alternatives: int = 0
    fwd_indels: int = 0  # primer indels, only from the indel rescue
    rev_indels: int = 0
    search_pass: str = "main"  # main, rescue (smaller word size, one sided) or aligned (primer indels)

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def total_mismatches(self) -> int:
        if self.fwd_mismatches < 0 or self.rev_mismatches < 0:
            return -1
        return self.fwd_mismatches + self.rev_mismatches

    @property
    def size_deviation(self) -> int | None:
        return None if self.expected_size < 0 else self.length - self.expected_size


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
class MarkerDef:
    """Repeat locus definition for one marker.

    patterns are repeat unit sequences in the orientation of the amplicon
    (forward primer first). unit_length is the repeat unit length in bp.
    ref_size and ref_repeats describe the amplicon in a reference strain
    (outer primer edge to outer primer edge, as merPCR reports it). They are
    None when unknown. offset is the non-repeat part of the amplicon; when not
    given it is derived from ref_size - ref_repeats * unit_length.
    """

    primer_id: str
    patterns: tuple[str, ...] = field(default_factory=tuple)
    unit_length: int = 0
    ref_size: int | None = None
    ref_repeats: int | None = None
    offset: int | None = None
    label: str | None = None  # locus name used in exports, for example the MLVAbank column name

    @property
    def display_name(self) -> str:
        return self.label or self.primer_id

    @property
    def has_repeats(self) -> bool:
        return bool(self.patterns)

    @property
    def pattern_length(self) -> int:
        """Kept for the legacy tables."""
        return self.unit_length

    @property
    def can_call(self) -> bool:
        return self.unit_length > 0 and self.ref_size is not None and self.ref_repeats is not None

    @property
    def flanking(self) -> int | None:
        if self.offset is not None:
            return self.offset
        if self.can_call:
            return self.ref_size - self.ref_repeats * self.unit_length
        return None

    def expected_size(self, repeats: float) -> float | None:
        if not self.can_call:
            return None
        return self.ref_size + (repeats - self.ref_repeats) * self.unit_length


# Backwards compatible name for the pattern-file era.
PatternDef = MarkerDef
