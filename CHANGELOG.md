# Changelog

Dates are the release dates. Versions follow semantic versioning; the tool
is below 1.0, so the output tables may still change between minor versions.

## Unreleased

### Added
- Bundled scheme `nmeningitidis_mlva12_liao_2006`, the 12 locus MLVA of Liao
  et al. for *Neisseria meningitidis*, anchored on Z2491. Every value is
  published: the paper's two tables give the repeat units, unit lengths and
  repeat numbers in three sequenced genomes alongside the primers and the
  predicted amplicon sizes, and they agree at every locus in all three. All
  twelve sizes reproduce in the reference, and typing MC58, which the scheme
  is not anchored on, reproduces its published numbers as well. NMTR12 is a
  compound locus of 13 bp and 12 bp units and is reported by length only.
- Bundled scheme `saureus_mlva5_sabat_2003`, the five locus gel based MLVA
  of Sabat et al. for *Staphylococcus aureus*, anchored on NCTC 8325. The
  paper publishes primers and no other quantitative value, and its pattern
  numbers label gel band patterns rather than per locus alleles, so the
  scheme sizes its markers and calls none. It is the first bundled scheme to
  use `known_multiple`: the `sdr` primer pair amplifies the tandem sdrC,
  sdrD and sdrE genes together by design and gives two products in the
  reference. The degenerate sspA primers are found by the rescue pass.
- Bundled scheme `saureus_mlva8_schouls_2009`, the eight locus MLVA of
  Schouls et al. for *Staphylococcus aureus*, anchored on N315. All eight
  products are reproduced there and each equals the span of the published
  coordinates for the same primer pair. The scheme sizes its markers without
  calling them, because the paper publishes neither the flanking lengths nor
  the repeat unit sequences and does not state the reference repeat number,
  and the arrays are too imperfect to recover it reliably; the scheme notes
  say what would be needed to enable calling.
- Bundled scheme `salmonella_enteritidis_mlva9_malorny_2008`, the nine locus
  multiplex of Malorny et al. for *Salmonella* Enteritidis, anchored on
  NCTC 13349. The loci are in the order the paper gives for the MLVA type.
  Seven of the nine reference repeat numbers are fractional, so those loci
  call on a grid offset from the integers, which is the first bundled use of
  an offset other than a half. Four published sequence lengths are not the
  product their own primers give in the reference assembly and the measured
  length is used; the differences are a fraction of a repeat and move no
  allele. The paper's MLVA type is a string of arbitrary allele indices
  rather than repeat numbers, and that mapping is not published, so a
  profile from this scheme is not directly comparable with a published type.
- Bundled scheme `yersinia_mlva26_li_2013`, the hierarchical 14+12 VNTR
  scheme of Li et al. for *Yersinia pestis*, anchored on CO92. The fourteen
  first step loci come before the twelve that resolve sub-populations, so a
  profile truncates to the first step result. Primers and reference values
  are the paper's supplementary Table S2; five of the amplicon sizes it
  encodes are not the product its own primers give in the reference, and in
  each case the measured size matches what `yersinia_klevytska_2001` or
  `yersinia_mlva25_lefleche_2001` records for the same locus with the same
  primers, so the measured size is used and the published identifier is kept
  in the label. All 26 loci reproduce that size in CO92.
- Bundled scheme `bacillus_mlva25_lista_2006`, the 25 locus MLVA panel of
  Lista et al. for *Bacillus anthracis*, the eight MLVA-8 loci with seventeen
  bams loci, anchored on Ames Ancestor. The paper publishes no primers, so it
  reuses those of `bacillus_mlva31_mlvabank`, of which the panel is the first
  25 loci; all 25 reproduce their reference size. It is the first bundled
  scheme with a `bins.tsv`, holding the paper's table of sequenced amplicon
  size against repeat number for all 25 loci, which is how the assay numbers
  alleles. Its numbering is the paper's own and differs from the MLVAbank one
  at vrrA, vrrB1, vrrB2, pXO2-at and bams34, so profiles from this scheme and
  from the other Bacillus schemes must not be compared.
- Bundled scheme `yersinia_mlva25_lefleche_2001`, the 25 loci of Le Fleche
  et al. polymorphic among *Yersinia pestis* and used as a panel by Pourcel
  et al. 2004, anchored on CO92. Unit length, repeat number and expected
  product are published and all 25 products are reproduced in the reference.
  Six loci are also in `yersinia_klevytska_2001`, which numbers them
  differently, so profiles from the two are not interchangeable.
- `mlvatyper schemes` prints an aligned table of the bundled schemes with
  the organism, the number of loci, how many of them carry the reference
  values a call needs, and the reference assembly accession.
  `mlvatyper schemes NAME` describes one scheme, including its loci,
  provenance, calling defaults and notes; a partial name or an organism
  name filters the list instead. `--names` and `--tsv` give the machine
  readable forms, and the bundled schemes are listed in `docs/method.md`.
- Bundled scheme `coxiella_mlva17_arricau_bouvery_2006`, the Coxiella
  burnetii panels of Arricau-Bouvery et al., ten panel 1 minisatellites
  followed by seven panel 2 microsatellites, anchored on Nine Mile RSA 493.
  Fourteen of the seventeen published amplicon sizes are reproduced there;
  the scheme records the measured size for ms30, ms31 and ms36 and states
  that alleles at ms30 and ms31 should not be compared with published
  allele numbers until their primers are confirmed.
- `ref_repeats` may be fractional, for schemes whose reference array ends
  in a partial unit and is numbered `n.5`. The alleles of such a marker
  lie on a half-integer grid, and a call on that grid is reported as `n.5`
  rather than flagged INTERMEDIATE. Whole numbers are unaffected.
- Bundled scheme `bacillus_mlva31_mlvabank`, the 31 locus composite panel
  for Bacillus anthracis: the eight MLVA-8 loci of Keim et al., seventeen
  bams loci of Lista et al. and six vntr loci of Van Ert et al. All of its
  reference values are read from the MLVAbank locus identifiers, and every
  one of the 31 is reproduced in the Ames Ancestor assembly.
- Bundled scheme `bacillus_mlva15_vanert_2007`, the MLVA-15 panel of Van
  Ert et al. for Bacillus anthracis, anchored on Ames Ancestor. It carries
  the eight MLVA-8 loci forward with the same primers and adds seven
  BaVNTR loci; its published repeat motifs and locations also corroborate
  the units and reference sizes of `bacillus_mlva8_keim_2000`.
- Bundled scheme `bacillus_mlva8_keim_2000`, the eight locus MLVA-8 panel
  of Keim et al. for Bacillus anthracis, anchored on Ames Ancestor. The
  primer table of Keim et al. carries no reference values, so these are
  taken from the MLVAbank locus identifiers and verified in the assembly.
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
- `saureus_mlva8_schouls_2009` calls two of its eight loci, where before it
  called none. Bogut et al. 2024 (doi:10.2478/cipms-2024-0014) tabulate the
  Pourcel MLVA-14 panel for the sequenced reference strain Mu50, and two of
  its VNTRs are the same physical arrays this panel reads under other
  primers, identified by matching genomic position and repeat unit:
  VNTR24_01 is Pourcel Sa0122 (the spa repeat) and VNTR63_01 is Sa1132. The
  shared array carries the same repeat count whichever primers flank it, so
  the published Mu50 counts of 10 and 6, confirmed across that paper's 31
  isolates, anchor these two loci. The scale is the Pourcel convention and
  may differ from Schouls' own by a constant; the other six loci stay sized
  only.
- `francisella_vogler_2008` calls seven of its eleven loci, where before it
  called none. Its source table gives primers only, so both values an allele
  needs come from `francisella_johansson_2004`, which counts the same arrays
  in SCHU S4. The unit is taken where a cross strain test and that scheme's
  independently derived unit agree. The copy number comes from calibrating
  through the shared genome: the array is the same DNA whichever primers
  flank it, so this scheme's product in SCHU S4 minus that array gives its
  flanking length, which fixes its count in LVS. All seven flanks yield a
  whole number of repeats, and across the five genomes where both panels
  amplify, all 35 locus and strain pairs give the same allele. The scale is
  inherited rather than published, so the two schemes now agree with each
  other but may still be offset from published allele names.
- The two *Francisella* schemes were audited across seven genomes, where
  every product of a locus must differ by whole multiples of its repeat
  unit. It confirms 15 of the 16 derived units of
  `francisella_johansson_2004` and fails at Ft-M8, whose steps of 15 and 33
  bases are not multiples of its 16 bp unit because in WY96-3418 its primers
  amplify a different member of the shared repeat family about 1.85 Mb away;
  the scheme now says so, and a fractional allele there is flagged rather
  than rounded. `francisella_vogler_2008` gains a unit length for seven of
  its eleven loci, filled only where that test and the Johansson panel agree
  independently, and still calls nothing because the test cannot supply an
  absolute copy number. Both schemes now record that their locus names are
  the same VNTRs under different spellings, Ft-M02 against Ft-M2, while
  their primers, references and alleles differ.
- The unit lengths of the two *Staphylococcus aureus* schemes were checked
  by amplifying each panel in eight *S. aureus* genomes, where every product
  of a locus must leave the same remainder when divided by its repeat unit,
  and four are corrected. In `saureus_mlva5_sabat_2003` the Ser-Asp arrays
  of clfA, clfB and sdr step by one dipeptide, 6 bp, not the 18 bp the paper
  describes. In `saureus_mlva8_schouls_2009` VNTR63_01 steps by 63 bp, not
  the 64 bp its table averages to. Both schemes now also record which loci
  that test confirms and which fit no single period, including VNTR61_01,
  whose heterogeneous flanking region the paper reports and the data shows.

- The copy numbers of `bacillus_mlva8_keim_2000` and
  `bacillus_mlva15_vanert_2007` are now the published MLVAbank values that
  `bacillus_mlva31_mlvabank` supplies, instead of counts derived from the
  reference. Five of the eight MLVA-8 loci change, including vrrC1 and
  vrrC2, which the panel numbers on a 9 bp and an 18 bp unit rather than
  the 36 bp higher-order period, so profiles from an earlier version of
  these two schemes are not comparable with profiles from this one.
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
