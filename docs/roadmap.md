# Roadmap

Planned 2026-10-23 after the run report (Part 1 of the plan) was implemented.
Items are ordered by the priorities given by the maintainer; each names the
files it touches so that it can be picked up as its own task.

### A. Scale and pipelines

1. Parallel samples (done): `--jobs N` searches N samples at a time with
   `concurrent.futures.ProcessPoolExecutor`, each process building its own
   searcher; merPCR keeps `-T` for threads within a sample. Results are
   merged in input order so the outputs match a serial run byte for byte.
2. Streaming large sample sets (done): each sample is written under
   `per_sample/` as it finishes and `--resume` skips the samples whose
   input file and search settings are unchanged. The combined tables are
   still built in memory; for very large sets they could be streamed too.
3. Packaging: Bioconda recipe (depends on merpcr, which needs its own
   recipe first, same owner), a Dockerfile and a Singularity definition
   built in CI from the conda environment, versioned with the git tag.
4. nf-core style Nextflow module and a Snakemake wrapper under
   `workflow/`, each running `search` per sample in parallel and
   `summarize`, `discover` and `report` once. Nextflow files without
   Unicode. Tests with nf-test on `tests/data`.
5. Structured log (partly done): `--log-json FILE` writes one JSON object
   per log event (JSON Lines) for pipeline monitors. Distinct exit codes
   (0 ok, 1 input error, 2 some samples failed) are still to be decided;
   the exit codes are now 0, 1 (input error) and 2 (some samples failed).

### B. Scheme and database support (done 2026-10-23)

Implemented: scheme directories with `scheme.yml` defaults and `--scheme`,
the bundled `testM` template, `--mlvabank` export with marker labels,
`bins-from-table`, and `check-scheme` with optional NCBI fetch, which
requires the provenance fields citation, url and reference.

The Salmonella Typhimurium scheme `STm_ECDC` and the Legionella
pneumophila scheme `Lpn_MLVA12` are bundled, verified against LT2 and
Philadelphia-1, and checked in CI against those assemblies, downloaded by
accession. The STTR3 locus of `STm_ECDC` showed a limit of the data model:
compound loci with two unit lengths (27 and 33 bp here) cannot be called
by the single unit formula. Such a locus is now called from a `bins.tsv`
alone, as are NMTR12 of `nmeningitidis_mlva12_liao_2006` and Ft-M20B of
`francisella_vogler_2008`, and a size outside the table gets no allele and
OFF_TABLE. Counting each of the two units separately and naming the allele
from the pair of counts remains a candidate item, since a table covers only
the sizes it lists.

Twenty-five published schemes are bundled and every locus of each calls; the
publication, reference strain and basis of the values of each are in
[scheme_references.md](scheme_references.md). Further published schemes need reference values from their publication,
the provenance fields filled in, and a passing `check-scheme` against the
assembly of the strain the publication used, so that a transcription
error shows up as a marker that misses its stated size or repeat count
([schemes.md](schemes.md) describes the steps). The CI job checks every
bundled scheme against its assembly at each change, so later changes to
the search or calling code cannot shift a published scheme's alleles
unnoticed. The original items follow for reference.

1. Scheme bundles under `src/mlvatyper/schemes/<name>/` with `primers.sts`,
   `markers.tsv`, optional `bins.tsv` and a `scheme.yml` (name, organism,
   reference strain, citation, allele numbering convention). Start with
   the testM scheme as the template; add published Brucella, Francisella
   and Yersinia schemes only with their reference values from the
   publications, each with a note on how the values were obtained.
   `--scheme NAME` replaces `--sts`, `--markers` and `--bins`.
2. MLVAbank export: `profile.tsv` already holds the profile string; add
   `--mlvabank` writing the MLVAbank upload format (one row per sample,
   marker columns in the scheme's order, their locus names) and the
   import of an MLVAbank allele table as a bin file.
3. Allele numbering conventions per scheme: `scheme.yml` carries the
   rounding policy and any locus specific numbering (for example loci
   numbered from the reference allele rather than from 0), applied by
   `calling.py` through the existing `Bin` and `round_repeats` paths.
4. Scheme validation command `mlvatyper check-scheme` that runs the
   primers against the scheme's reference assembly (downloaded by
   accession with Biopython Entrez when `--fetch` is given), confirms
   that every marker gives one product with the reference size and
   repeat count, and reports primers that merPCR cannot hash.

### C. Reads and long repeats

1. Collapsed array detection in assemblies: for every called amplicon
   compare the array length with `--max-array-length` (exists) and, new,
   with the assembly's own evidence when a BAM is given (`--bam`):
   coverage inside the array divided by coverage in the flanks; a ratio
   well above 1 marks a collapsed array (`COLLAPSED` status). Uses pysam,
   optional dependency.
2. Read based calling for long reads: `mlvatyper reads --fastq` extracts
   reads that contain both primers (exact or `-N` mismatches via the
   `rescue.align_primer` code), measures the product length per read,
   clusters the lengths, and calls the allele from the mode with the
   per-read support in the long table (`n_reads`, `fraction_mode`).
   Mixed clusters give `MULTIPLE`; this also separates true duplicates
   from assembly artefacts.
3. k-mer check for short reads: count the repeat unit's k-mers in the
   reads relative to flanking k-mers to estimate copy number for arrays
   longer than the read length; reported as `kmer_estimate` next to the
   assembly call, never as the call itself.

### D. Later: validation against wet-lab MLVA

`mlvatyper concordance --observed table.tsv` compares a conventional
allele table with the calls, reports per-marker agreement, systematic
offsets (capillary sizing) and the list of discordant samples; feeds the
report as a section when present.

### E. Housekeeping before the next release (done 2026-10-23)

Implemented: versions 0.3.0 and 0.4.0 with `CHANGELOG.md` and a release from the tag;
the method description moved to `docs/method.md` with `docs/report.md`
alongside it, leaving the README as usage; and a CI job that checks every
bundled scheme, runs the testM example end to end on the committed
assembly slices, verifies the alleles and the report, and uploads the
outputs as an artifact.

