# Changelog

Dates are the release dates. Versions follow semantic versioning; the tool
is below 1.0, so the output tables may still change between minor versions.

## Unreleased

### Added

- `ref_repeats` may be fractional, for schemes whose reference array ends
  in a partial unit and is numbered `n.5`. The alleles of such a marker
  lie on a half-integer grid, and a call on that grid is reported as `n.5`
  rather than flagged INTERMEDIATE. Whole numbers are unaffected.
- Bundled schemes `mtb_miru24_supply_2006` and `mtb_miru15_supply_2006`,
  the standardised MIRU-VNTR set for Mycobacterium tuberculosis and its
  discriminatory subset, anchored on H37Rv. The source table gives no
  copy numbers, so these are derived from the reference; confirm them
  against the published allele calling table before comparing a profile
  with published MIRU codes.
- Bundled scheme `brucella_mlva15_lefleche_2006`, the Brucella MLVA-15
  panel, anchored on B. melitensis 16M. Every locus reproduces the
  published reference size and copy number on the stated chromosome.
- Bundled scheme `yersinia_klevytska_2001`, the 42 locus VNTR survey of
  Klevytska et al. for Yersinia pestis, anchored on CO92. Three of its
  loci are numbered with a half copy.
- Bundled scheme `burkholderia_uren_2007`, the 32 locus VNTR panel of
  U'Ren et al. for Burkholderia pseudomallei, with reference sizes and
  repeat counts measured in strain K96243 because the publication's copy
  numbers could not be reproduced there.
- Bundled schemes `francisella_johansson_2004` and
  `francisella_vogler_2008` for Francisella tularensis, anchored on SCHU S4
  and LVS. The Johansson panel calls 16 of its 25 markers and sizes the
  rest, because the repeat units of the remaining loci are not in the
  source table and the arrays in the reference are too weak to derive.
- `known_multiple` in `scheme.yml` lists markers whose primers give more
  than one product in the reference because the locus is duplicated there.
  Checks that require one product per marker skip them, so a published
  locus with a known duplication does not have to be dropped.
- Bundled scheme `STm_ECDC`, the EU harmonised five locus scheme for
  Salmonella Typhimurium, with reference values verified against strain
  LT2, whose assembly CI downloads by accession and checks.
- Bundled scheme `Lpn_MLVA12`, the 12 locus MLVA panel for Legionella
  pneumophila (Sobral et al. 2011), verified against Philadelphia-1 and
  called from length alone; CI checks it against that assembly.
- Documentation of every output file, exit status and input rule in
  `docs/method.md`, a runnable example on the committed test slices in the
  README, and the scheme file rules in `docs/schemes.md`.
- Bundled scheme `Bru_MLVA16`, the 16 locus *Brucella* panel (Le Fleche et
  al. 2006; panels of Al Dahouk et al. 2007), with reference values from
  *B. melitensis* 16M; CI checks it against that assembly.
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
- A marker file whose `unit` column holds the unit length (`9`) rather
  than its sequence was read as the pattern `9` with a unit length of 1,
  which made every allele wrong; the length is now used, and a repeat unit
  that is not a DNA sequence is an error.
- A UTF-8 byte-order mark at the start of an STS, marker, bin or scheme
  file, as Excel and Windows editors write, no longer becomes part of the
  first id.
- Bin tables accept ranges written with an en dash or another dash
  character, take the header from the first line that is not a comment,
  and skip rows without an allele instead of writing an empty allele;
  decimal sizes round halves up.
- `report` no longer fails when `hits.tsv` holds markers that the
  definitions used by `summarize` do not list.
- `check-scheme` no longer copies the results of a marker onto the row of
  a rejected duplicate of its id.
- A `--fasta` path that does not exist stops the run instead of being
  left out without notice; a glob pattern that matches nothing is warned
  about.
- CLUSTAL alignments (`*_aligned.clw`) keep the full sequence names;
  names were cut to 30 characters, so copies of one sample or isolates with
  a long shared prefix became indistinguishable.
- `discover --engine trf` finds candidates for markers whose id contains a
  space; TRF cut the name at the space, so those markers had none. Without
  TRF on PATH, `--engine trf` now stops with a message instead of a
  traceback.
- Amplicon sequences are cut from the same text merPCR searched. merPCR
  drops characters such as `X` (masking) before it computes coordinates;
  the sequences were cut from the unfiltered FASTA, so in a masked assembly
  they were shifted while their length and primer mismatches looked right.
  A warning now gives the number of characters dropped.
- A marker defined twice in a marker file is an error instead of a
  duplicate column, and any error in a marker file is reported in one line
  instead of a traceback; so is a `--reference-sample` without amplicons.
- `<prefix>_lengths_summary.tsv` has a column for every STS marker, also
  one without any amplicon, as the other tables do.
- Several summarize runs in one output directory no longer mix in the
  report: parameters are recorded per prefix, the definitions used are
  written to `<prefix>_definitions.tsv`, and a summarize that calls nothing
  removes the call tables an earlier run left under that prefix.
- The report counts only absent markers as missing (not `NO_DEFINITION`)
  and states the size window from the definitions when `--markers` set it.
- `run_info.json` records the size and MD5 sum of the marker and bin files.
- Tandem counts are NA rather than 0 for markers without repeat unit
  sequences, and `COUNT_DISAGREES` is raised when the unit is not found at
  all next to an allele of one or more.
- Whitespace in a sample name (taken from the FASTA file name) becomes `_`,
  since FASTA record ids end at the first space.
- A negative allele is written in parentheses in the profile string,
  `5-(-1)-3`, so that the hyphen joined profile stays unambiguous; the
  allele columns keep the value.
- The built-in `discover` engine reports every separate array of a period,
  as TRF does, so a marker with two arrays of the same period is seen as
  mixed without TRF. Each candidate, from either engine, is counted around
  its own array, so two arrays of one unit are no longer both reported
  with the copy number of the larger one.
- `NEAR_CONTIG_END` uses each missing marker's own expected size plus
  `-M` as its window, instead of the largest expected size of all missing
  markers.
- `--legacy-tables` gives the length difference table for every marker with
  a unit length, also without a unit sequence (all of Lpn_MLVA12).
- A summarize removes the optional tables that an earlier run with the same
  prefix wrote and this run does not (`_mlvabank.csv`, `_clean` and legacy
  tables), so the output directory describes the latest run.

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
