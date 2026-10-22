# mlvatyper

Locate VNTR (variable number tandem repeat) amplicons in genome assemblies with
an STS primer table, and call MLVA (multiple-locus VNTR analysis) alleles that
are anchored to a reference strain. The primer search uses merPCR, a Python
reimplementation of me-PCR. This is a packaged rewrite of an earlier
single-file script (`vntr_searcher_9.py`).

## Method

1. The STS file (id, forward primer, reverse primer, expected product size,
   optional alias) is validated. Primers may contain IUPAC ambiguity codes.
   Entries with other characters, non-integer sizes or duplicate ids are
   reported and, with `--skip-invalid`, written to `invalid_primers.tsv`.
   Both the STS file and the FASTA files may be gzip, bzip2 or xz compressed
   (recognised by the `.gz`, `.bz2` or `.xz` suffix). This requires merPCR 1.0.4
   or newer.
2. Each FASTA file is treated as one sample, named after the file stem with
   any compression suffix removed (`sampleA.fasta.gz` gives `sampleA`).
   Primer pairs are searched with merPCR, allowing `-N` mismatches per primer
   outside the protected 3' bases (`--three-prime-match`). The reverse primer
   is only searched within `-M` bp of the STS expected size, so the expected
   size and the margin together define the allowed amplicon size window.
   IUPAC codes in primers are honoured unless `--no-iupac` is given.
3. Amplicons are extracted from the assemblies (reverse complemented for
   minus-strand hits). Primer mismatches are counted against the amplicon
   ends and the distance to the nearest contig end is recorded. Hits for one
   marker that share the forward primer position or overlap by more than a
   primer length are treated as one locus; the hit with fewest mismatches,
   then smallest deviation from the expected size, is kept and the number of
   discarded alternatives is recorded (`--no-dedup` keeps all). Amplicons are
   written to `hits.tsv` and to one FASTA per marker in `sequences/`. Record
   ids are `<sample>_<copy>`.
4. For markers without a hit in a sample, both primers are searched exactly
   within the expected size plus margin of every contig end. A match is
   recorded in `near_contig_end.tsv` because such a locus is more likely
   split by the assembly than deleted.
5. If MAFFT or MUSCLE is on PATH, amplicons per marker are aligned and written
   as aligned FASTA and CLUSTAL (`.clw`) using Biopython.
6. Summary tables, one row per sample and one column per marker, are written
   with the given prefix. `<prefix>_lengths_summary.tsv` holds amplicon
   lengths, several copies joined by `/`. With `--clwclean` the same table is
   also produced from amplicons without N (`*_clean.tsv`), and clean
   alignments are written.
7. Alleles are called when marker definitions are available (next section).

## Allele calling

The repeat number of an amplicon is anchored to a reference strain:

    raw = ref_repeats + (length - ref_size) / unit_length

where `ref_size` is the amplicon length in the reference strain (outer primer
edge to outer primer edge, as merPCR reports it) and `ref_repeats` the repeat
number in that strain. This is the convention of MLVAbank and of MLVA_finder.
A raw value closer than `--allele-tolerance` (default 0.25) to an integer is
that integer; a value exactly at the tolerance is not, as in MLVA_finder.
Other values depend on `--rounding`: `intermediate` (default)
reports `n.5` as MLVA_finder does, `nearest`, `up` and `down` round as named.
In all cases the call is flagged INTERMEDIATE. Rounding conventions differ
between laboratories, so the policy is recorded in `run_info.json`.

Independently of the length, the longest tandem array of the repeat unit is
counted in each amplicon, allowing `--unit-mismatches` differences per unit
(default 1), in either orientation and at any rotation of the unit. The
count is reported next to the length based allele, and a difference of one
or more repeats is flagged COUNT_DISAGREES. This exposes flanking indels and
wrong unit lengths, which a length alone cannot show.

### Marker definitions

Marker definitions come from three sources, which are merged:

- `--markers FILE`: tab separated with a header. Recognised columns are
  `marker`, `patterns` (repeat unit sequences in the amplicon orientation,
  separated by `;`, or `N/A`), `unit_length`, `ref_size`, `ref_repeats` and
  `offset`. Only `marker` is required; `unit_length` defaults to the length
  of the first pattern. The older pattern file (`marker`, `pattern`,
  `length`) is read by the same code. `--pattern` is accepted as an alias.
- STS ids in the MLVAbank convention `NAME_<unit>bp_<size>bp_<n>u` (for
  example `Bruce06_134bp_408bp_3u`) carry unit length, reference size and
  reference repeats and need no file.
- `--reference-sample NAME`: for markers that still lack `ref_size` and
  `ref_repeats`, the amplicon length of that sample is used as `ref_size`
  and its tandem array count as `ref_repeats` (0 when no pattern is known,
  in which case alleles are offsets from the reference rather than repeat
  numbers). The derived definitions are written to
  `marker_defs_derived.tsv` for review and reuse with `--markers`.

Markers without unit length cannot be called. Without any definition only
the lengths table is produced.

### Discovering repeat units

When neither unit sequences nor unit lengths are known, `mlvatyper discover`
seeds a definition file from the amplicons of one sample:

    mlvatyper search --sts primers.sts --fasta ref.fasta --outdir results
    mlvatyper discover --outdir results --sample ref
    mlvatyper summarize --outdir results --prefix study \
        --markers results/marker_defs_discovered.tsv

For each amplicon the dominant tandem array is located by self-comparison at
every period from `--min-unit` to `--max-unit` (3 to 100 bp), with the
mismatch allowance of the tandem counter. Arrays are scored as in Tandem
Repeats Finder (+2 per matching and -7 per differing position between
consecutive copies) and need `--min-score` (default 30) and `--min-copies`
(default 2). The smallest period whose score is close to the best is
reported so that a unit is not given as its own multiple. The unit is the
per-column majority over the array, written in the phase and orientation of
the amplicon. `discover_report.tsv` lists unit, copies, exact and variant
copies, score, unit entropy and array coordinates for review; markers
without an array keep only their reference size.

The result depends on the sample. A reference strain with one or two copies
at a locus, or with a degenerate array, yields no unit, so discovery is
best run on the sample with the most copies, and the file should be
reviewed and completed by hand where the scheme publication gives the unit.

An optional bin file (`--bins`: marker, size or `low-high` range, allele)
maps amplicon sizes directly to allele names and takes precedence over the
formula, as some schemes (for example Brucella) number alleles by
convention. A size outside every bin for that marker is flagged OFF_TABLE
and falls back to the formula.

`--max-array-length L` censors calls whose repeat array is at least L bp
long, reported as `>=n` with CENSORED. Long tandem arrays collapse in short
read assemblies, so L should be the assembly k-mer or read length.

### Output tables

- `<prefix>_alleles.tsv`: sample by marker allele values. Several amplicons
  are joined with `/`. Missing markers show a status code instead of a
  number.
- `<prefix>_status.tsv`: sample by marker status codes joined by `;`.
- `<prefix>_tandem_counts.tsv`: tandem array counts per amplicon.
- `<prefix>_profile.tsv`: alleles in STS order plus a hyphen joined profile
  string per sample.
- `<prefix>_calls_long.tsv`: one row per amplicon (and per missing marker)
  with coordinates, length, expected size, primer mismatches, discarded
  alternatives, distance to contig end, raw repeat number, allele, residual
  in bp between the length and the expected size of the called allele,
  tandem array count, and status.
- `run_info.json`: versions and all parameters of the search and summary
  steps.

Status codes: OK, NOT_DETECTED (no hit), NEAR_CONTIG_END (no hit but a primer
within the size window of a contig end), INTERMEDIATE, NEGATIVE (fewer bases
than the flanking region), BINNED, OFF_TABLE, MULTIPLE (several loci in the
sample), AMBIGUOUS (overlapping hits were collapsed), CENSORED, HAS_N,
COUNT_DISAGREES and NO_DEFINITION. The codes are text so that they cannot be
confused with the allele 0, which is a valid call.

`--legacy-tables` additionally writes the two tables of the original script:
exact non-overlapping pattern counts (`<prefix>_count_summary.tsv`) and the
length difference to the shortest amplicon in the run divided by the unit
length (`<prefix>_count_lengthdiff.tsv`). The latter is relative to the
dataset and not comparable between runs; both are kept for comparison only.

## Installation

merPCR (1.0.3 or newer) is not on PyPI or conda and requires Python 3.13 or newer.

    git clone https://github.com/FOI-Bioinformatics/mlvatyper
    cd mlvatyper
    conda env create -f environment.yml
    conda activate mlvatyper
    pip install .

or

    pip install git+https://github.com/FOI-Bioinformatics/mlvatyper

An aligner is optional: `conda install -c bioconda mafft` or `muscle`.

## Usage

    mlvatyper run --sts primers.sts --fasta 'assemblies/*.fasta' \
        --outdir results --prefix study -N 2 -M 500 \
        --markers markers.tsv --reference-sample ref_strain \
        --rounding intermediate --clwclean --skip-invalid -T 4

Subcommands: `run` (everything), `search` (hits, sequences, alignments only)
and `summarize` (rebuild tables and calls from an existing `results/hits.tsv`,
for example with other rounding settings or an edited definition file).

## Example

Running on the full example data (one assembly, 208 markers; not included in
the repository, a small slice with five markers is kept under `tests/data`):

    mlvatyper run --sts testID.sts --fasta test.fasta --outdir out_test \
        --prefix test -N 3 -M 500 --skip-invalid --clwclean

Log excerpt:

    ERROR invalid STS entry: line 189 (testID190): forward primer '1ACCATTTGTAGAGCCTGAACT' contains non-IUPAC characters ['1']
    INFO 207 valid primer pairs loaded from testID.sts (24 with IUPAC ambiguity codes, iupac mode on)
    WARNING 3 primers have ambiguities which prevent computation of a hash value: not included in search
    INFO test: 203 hits
    INFO Wrote out_test/test_lengths_summary.tsv (1 samples, 203 markers)

`out_test/` then contains `hits.tsv`, `primers_used.sts`, `invalid_primers.tsv`,
`samples.txt`, `near_contig_end.tsv`, `run_info.json`, `sequences/` with
per-marker FASTA, aligned FASTA and CLUSTAL files, and the lengths table.
Allele tables need marker definitions. First columns of `hits.tsv`:

    sample  marker     copy  contig  start    end      strand  length  has_N  expected_size  fwd_mismatches  rev_mismatches  n_alternatives  dist_to_contig_end
    test    testID001  1     test    1221679  1221902  -       224     0      224            0               0               0               672255

## Markers without hits

A marker can be absent from the output for two reasons.

- The primer pair was not found within the mismatch and size limits. Such a
  marker is still a column in the summary tables, reported as NA in the
  lengths table and NOT_DETECTED or NEAR_CONTIG_END in the allele table.
- merPCR could not build a hash word for a primer because it contains too many
  ambiguity codes close together. merPCR logs
  "primers have ambiguities which prevent computation of a hash value" and
  excludes the affected primer. If this is the forward primer the whole marker
  is dropped; if it is the reverse primer only the reverse-strand search is
  dropped. In the test data testID078 and testID178 are dropped entirely, and
  testID041 and testID043 are searched but not found within the limits.
  Two options help: `--one-sided` keeps such a marker and searches it on the
  strand its other primer allows (merPCR 1.0.3 or newer), and `-W` lowers the
  word size for all markers. In the test data `--one-sided` recovers
  testID078 and testID178, and `-W 8` additionally recovers testID043.

## Limitations

- Alleles are only as good as the reference values. Derived definitions from
  `--reference-sample` assume that the reference amplicon contains a tandem
  array the counter can see; check `marker_defs_derived.tsv`.
- Calls should be validated against isolates typed by conventional MLVA
  before they are compared with published profiles. In silico lengths are
  exact, while capillary sizing offsets are often baked into legacy allele
  tables.
- Short read assemblies collapse long tandem arrays; use `--max-array-length`
  and prefer long reads or large k-mers for such loci.
- Loci sharing a repeat unit are distinguished by their primers only.

## Differences from the legacy script

- me-PCR replaced by merPCR; seqret replaced by Biopython `AlignIO`.
- STS entries are validated up front instead of being silently dropped.
- IUPAC codes in primers are honoured by default (`--no-iupac` restores the
  me-PCR default where they count as mismatches).
- Overlapping hits for one locus are collapsed and primer mismatches are
  reported.
- Repeat numbers are anchored to a reference instead of to the shortest
  amplicon in the run, and non-integer values are flagged. The tandem array
  count tolerates variant units. The legacy tables remain available with
  `--legacy-tables`.
- Summaries are built from `hits.tsv`, not from re-parsed FASTA headers, and
  resuming is explicit through `mlvatyper summarize`.
- The default mismatch allowance per primer is 2 instead of 3.

## Development

    pip install -e ".[dev]"
    ruff check src tests
    pytest --cov=mlvatyper

Tests run in GitHub Actions on every push (`.github/workflows/tests.yml`).

## License

GPL-3.0-or-later, matching merPCR.
