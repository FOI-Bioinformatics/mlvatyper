# vntrfinder

Locate VNTR (variable number tandem repeat) amplicons in genome assemblies with
an STS primer table and summarise amplicon lengths and repeat unit counts across
samples. This is a packaged rewrite of an earlier single-file script
(`vntr_searcher_9.py`), with the external programs me-PCR and EMBOSS seqret
replaced by Python components.

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
   Primer pairs are searched with merPCR, a Python reimplementation of me-PCR,
   allowing `-N` mismatches per primer and a product size within `-M` bp of
   the expected size. IUPAC codes in primers are honoured unless `--no-iupac`
   is given.
3. Amplicons are extracted from the assemblies (reverse complemented for
   minus-strand hits) and written to `hits.tsv` and to one FASTA per marker in
   `sequences/`. Record ids are `<sample>_<copy>`.
4. If MAFFT or MUSCLE is on PATH, amplicons per marker are aligned and written
   as aligned FASTA and CLUSTAL (`.clw`) using Biopython.
5. Summary tables, one row per sample and one column per marker, are written
   with the given prefix:
   - `<prefix>_lengths_summary.tsv`: amplicon lengths, several copies joined by `/`.
   - `<prefix>_count_summary.tsv` (needs `--pattern`): repeat unit counts per
     copy joined by `;`. Counting is non-overlapping, case insensitive and on
     the amplicon strand only.
   - `<prefix>_count_lengthdiff.tsv` (needs `--pattern`): (length - shortest
     length of that marker across all samples) / repeat unit length, per copy.
   With `--clwclean`, the same tables are also produced from amplicons without
   N (`*_clean.tsv`), and clean alignments are written.

The pattern file is tab separated with a header: marker id, repeat unit(s)
separated by `;` or `N/A`, optional repeat unit length. Further columns are
ignored.

## Installation

merPCR (1.0.3 or newer) is not on PyPI or conda and requires Python 3.13 or newer.

    conda env create -f environment.yml
    conda activate vntrfinder

or

    pip install git+https://github.com/FOI-Bioinformatics/merpcr
    pip install .

An aligner is optional: `conda install -c bioconda mafft` or `muscle`.

## Usage

    vntrfinder run --sts primers.sts --fasta 'assemblies/*.fasta' \
        --outdir results --prefix study -N 3 -M 500 \
        --pattern patterns.tsv --clwclean --skip-invalid -T 4

Subcommands: `run` (everything), `search` (hits, sequences, alignments only)
and `summarize` (rebuild tables from an existing `results/hits.tsv`).

## Example

Running on the full example data (one assembly, 208 markers; not included in
the repository, a small slice with five markers is kept under `tests/data`):

    vntrfinder run --sts testID.sts --fasta test.fasta --outdir out_test \
        --prefix test -N 3 -M 500 --skip-invalid --clwclean

Log excerpt:

    ERROR invalid STS entry: line 189 (testID190): forward primer '1ACCATTTGTAGAGCCTGAACT' contains non-IUPAC characters ['1']
    INFO 207 valid primer pairs loaded from testID.sts (24 with IUPAC ambiguity codes, iupac mode on)
    WARNING 3 primers have ambiguities which prevent computation of a hash value: not included in search
    INFO test: 203 hits
    INFO Wrote out_test/test_lengths_summary.tsv (1 samples, 203 markers)

`out_test/` then contains `hits.tsv`, `primers_used.sts`, `invalid_primers.tsv`,
`samples.txt`, `sequences/` with per-marker FASTA, aligned FASTA and CLUSTAL
files, and the summary tables. First lines of `hits.tsv`:

    sample  marker     copy  contig  start  end    strand  length  has_N  sequence
    test    testID165  1     test    15843  16215  -       373     0      TAGAGTTATG...
    test    testID164  1     test    19455  19749  +       295     0      TGACAAAGAA...

## Markers without hits

A marker can be absent from the output for two reasons.

- The primer pair was not found within the mismatch and size limits. Such a
  marker is still a column in the summary tables, reported as NA.
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

## Differences from the legacy script

- me-PCR replaced by merPCR; seqret replaced by Biopython `AlignIO`.
- STS entries are validated up front instead of being silently dropped.
- IUPAC codes in primers are honoured by default (`--no-iupac` restores the
  me-PCR default where they count as mismatches).
- Summaries are built from `hits.tsv`, not from re-parsed FASTA headers, and
  resuming is explicit through `vntrfinder summarize`.
- The length difference table is relative to the shortest amplicon, as the old
  code computed; the old docstring said longest.
- Unused pattern file columns (min and max repeats, score) are ignored.

## Development

    pip install -e ".[dev]"
    ruff check src tests
    pytest --cov=vntrfinder

Tests run in GitHub Actions on every push (`.github/workflows/tests.yml`).

## License

GPL-3.0-or-later, matching merPCR.
