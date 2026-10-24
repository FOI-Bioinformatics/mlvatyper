# mlvatyper

Locate VNTR (variable number tandem repeat) amplicons in genome assemblies with
an STS primer table, and call MLVA (multiple-locus VNTR analysis) alleles that
are anchored to a reference strain. The primer search uses merPCR, a Python
reimplementation of me-PCR. This is a packaged rewrite of an earlier
single-file script (`vntr_searcher_9.py`).

## Documentation

- [Method and output](docs/method.md): the search, allele calling, marker
  definitions, schemes, exchange with MLVA databases, repeat unit discovery,
  every output table, the status codes, the limitations and the differences
  from the legacy script.
- [Run report](docs/report.md): the HTML report written for each run.
- [Building a scheme from a publication](docs/schemes.md): step by step
  from the primer table to a verified scheme directory.
- [Roadmap](docs/roadmap.md): planned work.
- [Audit](docs/audit_2026-10-22.md): comparison with other MLVA tools and the
  evaluation runs behind the current defaults.
- [CHANGELOG](CHANGELOG.md): what changed per release.

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
for example with other rounding settings or an edited definition file),
`discover` (seed a marker definition file from the amplicons), `report`
(write the HTML report from the tables of a run), `schemes`, `check-scheme`
and `bins-from-table` (scheme handling, see [Method and output](docs/method.md)).

## Example

`examples/testM/` holds a 31 marker scheme: `testM.sts` (primers and
expected sizes from the strain the scheme was designed on) and
`testM_pattern.txt`, the legacy pattern file with the repeat unit of every
marker, converted to `testM_markers.tsv` (patterns and unit length, no
reference values). The assembly used below is not included in the
repository; four slices of it around six of the markers are kept under
`tests/data` as `mini.fasta`.

    mlvatyper run --sts examples/testM/testM.sts --fasta test.fasta \
        --outdir out_testM --prefix testM --markers examples/testM/testM_markers.tsv \
        --reference-sample test --rescue-indels
    mlvatyper discover --outdir out_testM
    mlvatyper report --outdir out_testM --prefix testM

Log excerpt:

    INFO 31 valid primer pairs loaded from examples/testM/testM.sts (0 with IUPAC ambiguity codes, iupac mode on)
    INFO test: 31 hits
    INFO test: rescue pass (word size 6, one sided) found 1 hits for 1 markers
    INFO 32 amplicons for 30 markers written to out_testM/sequences
    INFO 1 sample/marker combinations without a hit: ONE_PRIMER_ONLY 1
    WARNING testM05: 2 amplicons in reference sample test, copy 1 used
    INFO Reference sample test: 30 of 31 markers have reference size and repeats

testM18 is found by the rescue pass: its reverse primer has one mismatch
inside the 11 bp word that merPCR must match exactly, so the main pass
misses it whatever `-N` allows. testM25 has six mismatches in its reverse
primer in this assembly and stays ONE_PRIMER_ONLY; testM25_2, a second
primer pair for the same locus, is found. testM05 and testM14 occur twice
in the assembly and are reported as MULTIPLE. The expected sizes in
`testM.sts` are the products of the design strain; the test assembly
differs from it by whole units at testM03 (two 9 bp units fewer),
testM06 (one 21 bp unit fewer) and testM25_2 (four 2 bp units more),
by 40 bp at testM26 and by one base at testM22.

With the assembly as reference every allele is its own tandem count, so
the run mainly checks that the length based and the array based numbers
agree: they do for all 30 markers. `mlvatyper discover` with TRF gives a
candidate in 28 of the 30 markers; for 22 the unit equals the one in the
pattern file, for 2 it has the right length but another sequence (one or
two copies only), 3 are other units, and 3 markers with one copy or a
2 bp unit give nothing. The outputs (allele, status, tandem count and
long tables, derived and discovered definitions, discovery report and
candidate table, missing markers) and the HTML report are kept under
`examples/testM/`.
First columns of `hits.tsv`:

    sample  marker   copy  contig  start   end     strand  length  has_N  expected_size  fwd_mismatches  rev_mismatches  n_alternatives  dist_to_contig_end  fwd_indels  rev_indels  search_pass
    test    testM01  1     test    148691  149000  +       310     0      310            0               0               0               148690              0           0           main

## Development

    pip install -e ".[dev]"
    ruff check src tests
    pytest --cov=mlvatyper

Tests run in GitHub Actions on every push (`.github/workflows/tests.yml`).

## License

GPL-3.0-or-later, matching merPCR.
