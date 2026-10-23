# mlvatyper

Locate VNTR (variable number tandem repeat) amplicons in genome assemblies with
an STS primer table, and call MLVA (multiple-locus VNTR analysis) alleles that
are anchored to a reference strain. The primer search uses merPCR, a Python
reimplementation of me-PCR. This is a packaged rewrite of an earlier
single-file script (`vntr_searcher_9.py`).

## Method

1. The STS file (id, forward primer, reverse primer, expected product size,
   optional alias) is validated. Primers may contain IUPAC ambiguity codes.
   Entries with other characters, invalid sizes or duplicate ids are
   reported and, with `--skip-invalid`, written to `invalid_primers.tsv`;
   `--fix-sts` instead removes stray characters from primers and logs each
   removal. The size may be a range `low-high`, which merPCR uses as a
   per-marker size window in place of `-M`.
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
   minus-strand hits). Primer mismatches are taken from merPCR, which counts
   them with the rule the search applied, and the distance to the nearest
   contig end is recorded. Hits for one
   marker that share the forward primer position or overlap by more than a
   primer length are treated as one locus; the hit with fewest mismatches,
   then smallest deviation from the expected size, is kept and the number of
   discarded alternatives is recorded (`--no-dedup` keeps all). Amplicons are
   written to `hits.tsv` and to one FASTA per marker in `sequences/`. Record
   ids are `<sample>_<copy>`.
4. Markers without a hit in a sample get two rescue passes. A second merPCR
   pass with a smaller word size (`--rescue-wordsize`, default 8; 0
   disables), one-sided search and a wider window (`--rescue-margin`,
   default twice `-M`) finds primers whose IUPAC codes leave no 11-mer for
   the main pass. With `--rescue-indels`, one primer is then anchored by an
   exact IUPAC-aware match and the other aligned with up to one indel.
   Hits record their origin in the `search_pass` column (`main`, `rescue`,
   `aligned`) and the indel counts. Markers that stay missing are listed in
   `missing_markers.tsv` with a reason: NOT_SEARCHED (no primer holds a
   hash word), NEAR_CONTIG_END (a primer within the size window of a contig
   end, so the locus is more likely split than deleted), ONE_PRIMER_ONLY
   (one primer matches exactly somewhere) or NOT_FOUND.
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
seeds a definition file from the amplicons of one or more samples:

    mlvatyper search --sts primers.sts --fasta 'assemblies/*.fasta' --outdir results
    mlvatyper discover --outdir results --all-samples --reference-sample ref
    mlvatyper summarize --outdir results --prefix study \
        --markers results/marker_defs_discovered.tsv

Two repeat finders are available (`--engine`, default `auto`):

- `trf`: Tandem Repeats Finder (Benson 1999; `conda install -c bioconda trf`,
  included in environment.yml) is run once on the amplicons with the pinned
  alignment parameters 2 7 7 80 10, the weak score floor as its minimum
  score and `--max-unit` as its maximum period. TRF handles indels inside
  arrays and reports fractional copy numbers. The exact command line is
  recorded in `run_info.json`.
- `builtin`: self-comparison at every period from `--min-unit` to
  `--max-unit` (3 to 100 bp), with the mismatch allowance of the tandem
  counter and the same scoring as TRF (+2 per matching and -7 per differing
  position between consecutive copies), without indel handling. Used when
  TRF is not on PATH.

Every array an engine reports becomes a candidate. Arrays that describe the
same locus at several periods are merged, keeping the smallest period whose
score is within 85 percent of the best, so that a unit is not given as its
own multiple. Each candidate is counted again with the tandem counter so
that `ref_repeats` is the integer allele calling will see, and is placed in
a category:

| category | conditions |
|---|---|
| A | unit of at least `--strong-unit` bp (5), at least `--strong-copies` integer copies (3), score at least `--min-score` (30), unit entropy at least `--min-entropy` bits (1.0), no overlap with a primer |
| B | unit of at least `--strong-unit` bp, at least two integer copies (`--min-copies` engine copies, default 1.8), score at least `--min-score`, no primer overlap |
| C | anything else above `--weak-score` (20): short units, weak scores, arrays reaching a whole unit into a primer |

Candidates of a marker are ranked by category, cross-sample support, integer
copies and score. The first is the primary and seeds the definition when its
category is at least `--accept` (default B); markers whose primary is only
category C keep their reference size. A marker with two or more
non-overlapping A or B candidates is flagged as mixed in the report, so that
the choice can be reviewed; when such candidates share the unit length,
their units are all written as patterns of the definition, so the tandem
counter accepts either unit, while `ref_repeats` counts the primary array
alone.

A single reference genome typically holds two or three copies of each unit,
where a chance tandem scores like the true VNTR. The best evidence is
variation: with `--sample` given several times or `--all-samples`,
candidates are pooled per marker by unit (equal up to rotation and reverse
complement) and two columns are added. `copy_range` is the smallest and
largest integer copy number over the samples, and `explains` is the
fraction of sample pairs with different amplicon lengths whose difference
is a whole number of units within `--allele-tolerance`. A unit seen in two
or more samples with a copy range of at least one and `explains` of at
least 0.8 is category A whatever its copy number in the reference; a unit
that explains less than half of the differences while another unit of the
marker explains them is demoted one category. `ref_size` and `ref_repeats`
come from `--reference-sample` (default the first sample). When the best
supported unit was not found as an array in the reference amplicon, it is
counted there directly (engine `pooled`), so the definition is still
anchored to the reference even with one or zero copies.

Outputs: `discover_report.tsv` has one row per marker with the primary
candidate (unit, copies, exact and variant copies, score, entropy, array
coordinates, engine and its own copy number, TRF percent matches and
indels), its category, the number of candidates, the mixed flag and the
cross-sample columns. `discover_candidates.tsv` lists every candidate of
every sample with rank and category. Both are meant for review; the
definition file should be completed by hand where the scheme publication
gives the unit.

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
  alternatives, distance to contig end, primer indels and search pass, raw
  repeat number, allele, residual
  in bp between the length and the expected size of the called allele,
  tandem array count, and status.
- `run_info.json`: versions and all parameters of the search and summary
  steps.

Status codes: OK, NOT_DETECTED (no hit), NEAR_CONTIG_END (no hit but a primer
within the size window of a contig end), NOT_SEARCHED and ONE_PRIMER_ONLY
(from `missing_markers.tsv`), INTERMEDIATE, NEGATIVE (fewer bases than the
flanking region), BINNED, OFF_TABLE, MULTIPLE (several loci in the sample),
AMBIGUOUS (overlapping hits were collapsed), PRIMER_INDEL (hit from the
indel rescue), CENSORED, HAS_N, COUNT_DISAGREES and NO_DEFINITION. The codes are text so that they cannot be
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

`search --markers markers.tsv` gives every marker with `ref_size` and
`ref_repeats` an STS size range from 0 to `--max-repeats` (50) repeats, so
that the size window follows the locus instead of the global `-M`.

Subcommands: `run` (everything), `search` (hits, sequences, alignments only),
`summarize` (rebuild tables and calls from an existing `results/hits.tsv`,
for example with other rounding settings or an edited definition file) and
`discover` (seed a marker definition file from the amplicons).

## Example

Running on the full example data (one assembly, 208 markers; not included in
the repository, a small slice with five markers is kept under `tests/data`):

    mlvatyper run --sts testID.sts --fasta test.fasta --outdir out_test \
        --prefix test -N 2 -M 500 --fix-sts --rescue-indels --clwclean

Log excerpt:

    WARNING line 189 (testID190): removed '1' at position 1 of the forward primer
    INFO 208 valid primer pairs loaded from testID.sts (24 with IUPAC ambiguity codes, iupac mode on)
    WARNING 3 primers have ambiguities which prevent computation of a hash value: not included in search
    INFO test: 203 hits
    INFO test: rescue pass (word size 8, one sided) found 3 hits for 3 markers
    INFO test: primer alignment found 1 hits for 1 markers
    INFO 207 amplicons for 207 markers written to out_test/sequences
    INFO 1 sample/marker combinations without a hit: ONE_PRIMER_ONLY 1

The three markers of the rescue pass (testID043, testID078, testID178) have
IUPAC codes in every 11-mer of one primer; all three are found with no
mismatches at the expected size. testID041 carries one base more in its
reverse primer than the genome and needs `--rescue-indels`. testID062 has
three mismatches in its reverse primer and is found with `-N 3`.

`out_test/` then contains `hits.tsv`, `primers_used.sts`, `samples.txt`,
`missing_markers.tsv`, `near_contig_end.tsv`, `run_info.json`, `sequences/`
with per-marker FASTA, aligned FASTA and CLUSTAL files, and the lengths
table. Allele tables need marker definitions. `mlvatyper discover --outdir
out_test` with TRF gives a candidate in 92 of the 207 amplicons, 7 of
category A and 17 of category B; the report, the candidate table and the
resulting definition file are kept under `examples/testID/`. First columns
of `hits.tsv`:

    sample  marker     copy  contig  start    end      strand  length  has_N  expected_size  fwd_mismatches  rev_mismatches  n_alternatives  dist_to_contig_end  fwd_indels  rev_indels  search_pass
    test    testID001  1     test    1221679  1221902  -       224     0      224            0               0               0               672255              0           0           main

## Markers without hits

A marker without an amplicon in a sample is still a column in the summary
tables, reported as NA in the lengths table and with the reason from
`missing_markers.tsv` in the allele and status tables. The reasons are:

- NOT_SEARCHED: neither primer holds a run of unambiguous bases as long as
  the word size of the most permissive merPCR pass (`-W`, or
  `--rescue-wordsize` when the rescue pass ran). merPCR logs "primers have
  ambiguities which prevent computation of a hash value" for such primers.
- NEAR_CONTIG_END: a primer matches within the size window of a contig end;
  the locus is more likely split by the assembly than deleted.
- ONE_PRIMER_ONLY: one primer matches exactly somewhere in the assembly, so
  the other primer has more than `-N` mismatches or an indel
  (`--rescue-indels` allows one), or the product lies outside the size
  window.
- NOT_FOUND: neither primer matches exactly.

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
