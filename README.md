<p align="center">
  <img src="docs/images/logo.jpg" alt="mlvatyper: MLVA typing from sequencing data" width="720">
</p>

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
- [Bundled schemes and their references](docs/scheme_references.md): the
  publication, DOI, reference strain and basis of the values for each of the
  bundled schemes, with a bibliography.
- [Building a scheme from a publication](docs/schemes.md): step by step
  from the primer table to a verified scheme directory.
- [Roadmap](docs/roadmap.md): planned work.
- [Screening check](docs/screening_2026-11-03.md): the U'Ren scheme on 40
  genomes, what it showed and how to repeat it.
- [Audit](docs/audit_2026-10-22.md): comparison with other MLVA tools and the
  evaluation runs behind the current defaults.
- [CHANGELOG](CHANGELOG.md): what changed per release.

## Installation

mlvatyper requires Python 3.13 or newer and merPCR 1.0.5, which is not on
PyPI or conda and is installed from GitHub by both routes below.

    git clone https://github.com/FOI-Bioinformatics/mlvatyper
    cd mlvatyper
    conda env create -f environment.yml
    conda activate mlvatyper

The conda environment installs mlvatyper itself (editable) together with
MAFFT and Tandem Repeats Finder. Without conda:

    pip install git+https://github.com/FOI-Bioinformatics/mlvatyper

Two tools are optional and picked up from PATH: an aligner for the
per-marker alignments (`conda install -c bioconda mafft`, or `muscle`) and
Tandem Repeats Finder for `discover` (`conda install -c bioconda trf`;
a built-in repeat finder is used without it).

## Usage

With a bundled scheme, which supplies the primers, the marker definitions
and the calling defaults. Twenty-five published schemes are bundled, for
*Bacillus anthracis*, *Bifidobacterium longum*, *Brucella*, *Burkholderia
pseudomallei*, *Coxiella burnetii*, *Francisella tularensis*, *Legionella
pneumophila*, *Mycobacterium tuberculosis*, *Mycoplasma agalactiae*,
*Neisseria meningitidis*, *Pseudomonas aeruginosa*, *Salmonella enterica*,
*Staphylococcus aureus* and *Yersinia pestis*. `mlvatyper schemes` lists them,
`mlvatyper schemes NAME` describes one, and the
[scheme references](docs/scheme_references.md) give the publication behind each:

    mlvatyper run --scheme STm_ECDC --fasta 'assemblies/*.fasta' \
        --outdir results --prefix study --mlvabank

With your own primer and marker files:

    mlvatyper run --sts primers.sts --fasta 'assemblies/*.fasta' \
        --outdir results --prefix study -N 2 -M 500 \
        --markers markers.tsv --reference-sample ref_strain \
        --rounding intermediate --clwclean --skip-invalid -T 4

Each FASTA file is one sample, named after the file. `--markers` gives every
marker with `ref_size` and `ref_repeats` a size window from 0 to
`--max-repeats` (50) repeats, so that the window follows the locus instead
of the global `-M`. The main results are `study_alleles.tsv`,
`study_status.tsv` and `study_report.html`; [Method and
output](docs/method.md) lists every file.

Subcommands: `run` (everything), `search` (hits, sequences, alignments only),
`summarize` (rebuild tables and calls from an existing `results/hits.tsv`,
for example with other rounding settings or an edited definition file),
`discover` (seed a marker definition file from the amplicons), `report`
(write the HTML report from the tables of a run), and for schemes
`schemes` (list the bundled MLVA schemes, or describe one with its loci and
provenance), `new-scheme`, `check-scheme` and `bins-from-table` (see [Method
and output](docs/method.md) and [Building a scheme](docs/schemes.md)).
`mlvatyper <subcommand> --help` lists the options.

`search` and `run` take `--jobs N` to search N samples in parallel, one
process each (`--threads` stays merPCR's threads within a sample, so the
total worker count is `--jobs` times `--threads`); the outputs are the same
as a serial run. The global `--log-json FILE`, before the subcommand, also
writes one JSON object per log event to FILE for pipeline monitors.

A run exits with status 1 when an input is invalid or a sample failed; the
other samples are still typed and reported.

## Example

`examples/testM/` holds a 31 marker scheme, also bundled as `--scheme
testM`. It is built on the Johansson *Francisella* panel: 22 of its markers
use a Johansson primer pair and repeat unit, and its expected sizes are those
of a *holarctica* strain. It holds `testM.sts` (primers and expected sizes from
the strain the scheme was designed on) and `testM_pattern.txt`, the legacy pattern file with the
repeat unit of every marker, converted to `testM_markers.tsv` (patterns and
unit length, no reference values). The assembly used below is not included
in the repository; four slices of it around six of the markers are kept
under `tests/data` as `mini.fasta`, which runs from a clone:

    mlvatyper run --scheme testM --fasta tests/data/mini.fasta \
        --outdir out_mini --prefix testM --reference-sample mini --align none

This is the run the CI example job checks: testM01, 02, 03, 06, 3a and 3b
are typed 5, 3, 11, 4, 11 and 11, and the other 25 markers lie outside the
slices. On the full assembly:

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
`examples/testM/`. They were produced with mlvatyper 0.3.0, before the
fixes listed under Unreleased in the [CHANGELOG](CHANGELOG.md), and have
not been regenerated because the assembly is not in the repository; a
current run also writes `testM_definitions.tsv`.
First columns of `hits.tsv`:

    sample  marker   copy  contig  start   end     strand  length  has_N  expected_size  fwd_mismatches  rev_mismatches  n_alternatives  dist_to_contig_end  fwd_indels  rev_indels  search_pass
    test    testM01  1     test    148691  149000  +       310     0      310            0               0               0               148690              0           0           main

## Development

    pip install -e ".[dev]"
    ruff check src tests
    pytest --cov=mlvatyper

Tests that need TRF, an aligner or a reference assembly are skipped when
these are missing. `MLVATYPER_ONLINE=1 pytest` downloads the assemblies of
the bundled schemes from NCBI into the cache, with the few further genomes that
individual scheme tests use, and runs the scheme tests against them. A test
checks that [docs/scheme_references.md](docs/scheme_references.md) lists every
bundled scheme and every DOI its `scheme.yml` cites. Tests run in GitHub Actions on every push
(`.github/workflows/tests.yml`), with all of these enabled.

## License

GPL-3.0-or-later, matching merPCR.
