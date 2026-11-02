# Changelog

Dates are the release dates. Versions follow semantic versioning; the tool
is below 1.0, so the output tables may still change between minor versions.

## Unreleased

### Added

- Bundled scheme `STm_ECDC`, the EU harmonised five locus scheme for
  Salmonella Typhimurium, with reference values verified against strain
  LT2 and a committed fixture of the LT2 amplicon regions that CI checks.
- Bundled scheme `Lpn_MLVA12`, the 12 locus MLVA panel for Legionella
  pneumophila, verified against Philadelphia-1 and called from length
  alone; a fixture of the twelve amplicon regions is checked in CI.
- `mlvatyper new-scheme DIR` scaffolds a scheme directory with annotated
  template files, and `docs/schemes.md` describes how to build a scheme
  from a publication and verify it against the reference strain.

### Changed

- Schemes name their reference assembly by a GenBank or RefSeq accession
  (`assembly: GCF_...` in `scheme.yml`); `check-scheme` downloads it from
  NCBI Datasets into a local cache instead of needing a FASTA file, and the
  reference assemblies that were committed under `tests/data` are removed.
  The two bundled published schemes now carry their accessions.
- `check-scheme` also rejects provenance that opens with a placeholder or
  says it is still to be added, such as `url: not recorded; ... the DOI has
  to be added here`. The Lpn_MLVA12 scheme, which had that `url`, now cites
  Sobral et al. 2011 (doi:10.1128/AEM.05556-11) in full.
- A file that is not FASTA, or a corrupt or truncated compressed file, now
  fails only its own sample instead of stopping the run with a traceback.
  `run` still types and reports the other samples, leaves the failed ones
  out of `samples.txt` and exits with status 1.
- An explicit size range in the STS file is kept when `--markers` derives
  size windows, and hits keep the STS expected size rather than the
  midpoint of a derived window.
- Rounding: `nearest` rounds halves up (2.5 to 3) instead of to the even
  neighbour, and `intermediate` gives -0.5 rather than -1.5 for a raw
  value of -0.4.

### Fixed

- An option given on the command line now always wins over a scheme
  default, also when its value equals the parser default.
- A non-numeric allele in a bin table (for example `3a`) no longer stops
  allele calling.
- `--reference-sample` keeps the marker label and offset, so the MLVAbank
  export uses the locus names.
- The allele table shows `NO_DEFINITION` for a marker without a definition
  instead of another flag such as `HAS_N`.
- A product found by the primer indel alignment, and a product found more
  than once by `check-scheme`, is counted once.
- A `#` inside a `scheme.yml` value, such as a URL fragment, is no longer
  cut off as a comment.
- The command line in the report uses `--no-dedup` and `--no-iupac`
  instead of flags that do not exist, and `NOT_FOUND` has an explanation.

## 0.3.0 - 2026-10-23

### Added

- `mlvatyper report` writes a self-contained HTML report for a run, and
  `run` writes it unless `--no-report` is given. It holds a status summary,
  the allele profiles shaded by status, figures of amplicon lengths,
  residuals and status, the evidence per amplicon, the missing markers with
  their reasons, the discovery report, a methods paragraph generated from
  the recorded parameters, the input checksums and the marker definitions.
- Scheme bundles: a directory with `primers.sts`, `markers.tsv`, optional
  `bins.tsv` and `scheme.yml`, selected with `--scheme`, which also applies
  the scheme's calling defaults. `mlvatyper schemes` lists the bundled
  schemes; `testM` ships as the template.
- `mlvatyper check-scheme` validates a scheme and, with a reference
  assembly (local or fetched from NCBI), reports the products, sizes and
  repeat counts per marker. It requires the provenance fields `citation`,
  `url` and `reference` in `scheme.yml`.
- `--mlvabank` writes an allele table with one column per marker, named by
  the new `label` column of the marker file. `mlvatyper bins-from-table`
  converts a published allele size table into a bin file.
- Rescue passes for markers without a hit: a second one-sided merPCR pass
  with a smaller word size (`--rescue-wordsize`, default 6) and, with
  `--rescue-indels`, an anchored alignment allowing one indel in a primer.
  Hits record their `search_pass` and indel counts.
- `missing_markers.tsv` gives a reason per missing marker (NOT_SEARCHED,
  NEAR_CONTIG_END, ONE_PRIMER_ONLY, NOT_FOUND), reported in the call tables.
- Repeat unit discovery reports every candidate per amplicon with a
  category (A, B, C), a mixed flag and `discover_candidates.tsv`. With
  several samples (`--all-samples`) candidates are pooled per unit and the
  copy range and explained length differences promote the varying unit.
- `--fix-sts` removes characters that are not IUPAC codes from primers.
  STS sizes may be `low-high` ranges, and `search --markers` writes a
  per-marker size window from the definitions.
- `run_info.json` records the size and MD5 sum of every input file.

### Changed

- Primer mismatch counts come from merPCR 1.0.5 instead of being counted
  from the amplicon; the local count remains as a fallback.
- Headerless legacy pattern files keep their first line, which the legacy
  script dropped as a header.
- The README is split: usage stays, the method moves to `docs/method.md`.

### Fixed

- The expected sizes of the bundled testM scheme were one base short of the
  product.

## 0.2.1 - 2026-10-23

- Primer mismatch counts taken from merPCR 1.0.5.
- GitHub Actions moved to the Node 24 releases of checkout and
  setup-miniconda.

## 0.2.0 - 2026-10-22

- Reference-anchored allele calling with a mismatch tolerant tandem array
  count as cross-check, status codes per call, hit deduplication with
  primer mismatch evidence, marker definition files and the `discover`
  subcommand with a built-in detector and Tandem Repeats Finder.
- Renamed from vntrfinder to mlvatyper.
