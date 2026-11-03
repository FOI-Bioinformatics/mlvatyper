# Changelog

Dates are the release dates. Versions follow semantic versioning; the tool
is below 1.0, so the output tables may still change between minor versions.

## Unreleased

### Added
- Report figures for five or more samples: a clustered heatmap of the allele profiles with a UPGMA tree, a pairwise distance matrix, allele spectra per locus and call rate bars; the status matrix follows the cluster order. New module `plots.py`; the clustering is written with NumPy and adds no dependency.
- Report layout: the header and the Methods paragraph name the scheme, its publication (DOI), organism and reference strain; table captions come above the table and are numbered in sequence; the allele table uses short codes (ND, CE, 1P, NS, NDef, OFF) with a key and tooltips; header rows and first columns stay in view in wide tables; numbers are right aligned; the status figure has a legend and grid; the per-sample table is collapsed above 20 samples; print styles use landscape A4. The Methods text states the mismatch limit of the aligned rescue and the markers that skip the tandem count comparison.
- The report has an "Across samples" subsection for five or more samples, and two tables, `<prefix>_sample_summary.tsv` and `<prefix>_locus_summary.tsv`: call rate per sample and locus, the most frequent flag per locus, diversity, samples with fewer than 75% of loci OK, and loci to review (one flag on more than half of the typed calls, or an allele of zero or below). New module `screening.py`.
- `missing_markers.tsv` has a `detail` column. For NEAR_CONTIG_END it states where each primer lies relative to the contig ends, so a product split across contigs by a repeat array can be told from a single primer near an end.
- The `interrupted_arrays` key in `scheme.yml` lists markers whose tandem count is reported but not compared with the length based allele, so COUNT_DISAGREES is not raised for arrays that hold variant units the patterns cannot follow. `check-scheme` skips the repeat comparison for them too.
- The `rescue_indels` key in `scheme.yml` (yes or no) sets `--rescue-indels` as a
  scheme default. The `burkholderia_uren_2007` scheme uses it: on five
  *B. pseudomallei* and four *B. mallei* genomes the aligned rescue recovers
  6 of 13 and 6 of 15 missing products at default settings; raising `-N`
  alone changed little. Missing products that remain have primer sites with
  four or more mismatches or a product outside the size window.
- `docs/scheme_references.md`, a page that gives, for every bundled scheme,
  its publication with a link to its DOI, the reference strain and assembly,
  and the basis of its reference values: published, measured, derived,
  calibrated through another scheme, or a convention of the scheme. It also
  groups the schemes that overlap, lists the places where a count is not a
  published value, and ends with a bibliography taken from the DOI registry.
  Two tests keep it in step with the code: every bundled scheme and every DOI
  that a `scheme.yml` cites must appear on it, and its locus counts and
  assemblies must match the scheme files.
- Bundled scheme `paeruginosa_mlva8_onteniente_2003`, the eight locus MLVA
  of Onteniente et al. for *Pseudomonas aeruginosa*, anchored on PAO1. Table
  1 gives the repeat unit and the number of units in PAO1 for each locus;
  the primers reproduce every published PAO1 product and a run returns the
  unit profile 11-12-5-8-7-12-14-45 with zero residual. ms173 is a 243 bp
  unit locus with a product near 3.5 kb. First scheme for this organism, and
  the one the later P. aeruginosa MLVF studies apply.
- Bundled scheme `blongum_mlva5_matamoros_2011`, the five locus MLVA of
  Matamoros et al. for *Bifidobacterium longum* subsp. *longum*, anchored on
  NCC2705. The repeat units and motifs are from the paper's Table 3 and the
  NCC2705 reference profile 6-4-3-7-16 from Table 1; the primers reproduce
  every amplicon in the reference assembly and a run returns the published
  profile with zero residual. First scheme for this organism.
- Bundled scheme `magalactiae_mlva4_mcauliffe_2008`, the four locus MLVA of
  McAuliffe et al. for *Mycoplasma agalactiae*, anchored on type strain PG2.
  Table 1 gives the repeat unit and PG2 copy number for each locus; three
  are fractional, reflecting imperfect arrays, so the reference profile is
  2.3-2-3.6-2.3. Magal_VNTR5 carries a mismatched reverse primer and is
  found by the rescue pass. First scheme for this organism.
- A locus can now be called from a bins table alone, without a repeat unit.
  This makes a compound locus callable, one that interleaves two repeat unit
  sizes and so has no single unit for the allele formula. STTR3 of `STm_ECDC`
  uses the fragment-size-to-allele table of Larsson et al. 2009, and NMTR12
  of `nmeningitidis_mlva12_liao_2006` uses the three sequenced-genome alleles
  of Liao et al. 2006; both reference strains now call their published
  composite allele, and an off-table size is reported OFF_TABLE.
- Bundled scheme `saureus_mlva14_pourcel_2009`, the 14 locus MLVA of Pourcel
  et al. for *Staphylococcus aureus*, in two multiplex panels, anchored on
  Mu50. It is the first fully callable *S. aureus* scheme bundled: the repeat
  units and Mu50 copy numbers are from Pourcel et al. Table 1, and the Mu50
  amplicon sizes and a 31 isolate validation set are from Bogut et al. 2024
  (doi:10.2478/cipms-2024-0014). The primers reproduce every Mu50 amplicon,
  the reference profile matches the published one, and 57 of 66 clinical
  cells reproduce exactly, the rest within one repeat of the gel-sized
  source. Two of its loci are the arrays the Schouls panel reads as
  VNTR24_01 and VNTR63_01.
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
- `search` and `run` take `--jobs N` to search N samples in parallel, one
  process each; the outputs are identical to a serial run, and `--threads`
  stays merPCR's threads within a sample.
- Global `--log-json FILE` streams one JSON object per log event (JSON
  Lines) to FILE for pipeline monitors, alongside the console output.
- `mlvatyper new-scheme DIR` scaffolds a scheme directory with annotated
  template files, and `docs/schemes.md` describes how to build a scheme
  from a publication and verify it against the reference strain.

### Changed
- `francisella_vogler_2008` calls all eleven loci. Ft-M20B, the last, is
  called from a bins table under a convention of the scheme, since no
  publication gives its copy number. Vogler et al. describe it as a roughly
  200 bp insertion with a 15 bp repeat that is present in type A.II and
  *novicida* isolates. In 22 genomes the products are 153 bp without the
  insertion and 351, 366, 381, 396 or 411 bp with it, five sizes in steps of
  exactly 15 bp, which confirms the 15 bp unit and not the 12 bp in one later
  table. The 153 bp state has no array and is allele 0; the others are
  `(size - 306) / 15`, from 3 to 7. The 153 bp state does not lie on the
  15 bp ladder, so the repeat formula cannot link it to the rest, which is
  why the locus is binned, and an unlisted size is reported OFF_TABLE. A
  note that the 153 bp product held one copy of the unit was wrong, since
  the sequence counted lies inside the forward primer.
- `saureus_mlva5_sabat_2003` now calls all five loci, where before it sized
  them only. The paper publishes no repeat counts and its primers lie inside
  the repeat arrays, so no published scale or flanking landmark exists. The
  scheme defines its own allele, the number of whole repeat units between the
  two primer sites, which gives an exact integer for every strain: all 54
  products over eight genomes agree with the repeat formula, with no negative
  alleles. This is a convention of the scheme, documented as such, and not
  Sabat's pattern number or a count comparable with other schemes. The sdr
  locus is called per product and flagged MULTIPLE.
- `saureus_mlva8_schouls_2009` calls all eight loci, up from four, on the
  scheme's own published scale. Figure 7 of Schouls et al. 2009 gives the
  MLVA profile of 14 sequenced strains; the N315 profile, 14-2-1-3-2-11-7-5,
  anchors every locus, and a run on N315 reproduces it with zero residual.
  The four loci previously borrowed from the Pourcel panel are corrected to
  the Schouls convention. Validated against seven more sequenced genomes,
  seven of eight loci reproduce each published profile exactly, including
  the MRSA252 null allele; only VNTR61_01, which the paper reports as having
  a heterogeneous flank, deviates off the reference.
- `francisella_vogler_2008` called ten of its eleven loci, up from seven.
  With the Johansson panel now carrying published counts for Ft-M22, Ft-M23
  and Ft-M24, those three Vogler loci calibrate through the shared genome
  and agree with the Johansson allele in every strain where both amplify.
  Ft-M20B gains its published 15 bp repeat unit from Vogler et al. 2008 but
  stays sized only, because no publication gives its copy number and it is
  monomorphic in the holarctica reference.
- `francisella_johansson_2004` now calls all 25 loci, up from 16, on the
  published allele scale. Table 2 of Johansson et al. 2004 gives the repeat
  motif, its size and the SCHU S4 copy number for every marker, which fills
  the nine loci that were sized only and corrects eight derived counts that
  were off because their arrays are imperfect. A run on SCHU S4 returns the
  published profile with zero residual. The seven calibrated loci of
  `francisella_vogler_2008` are re-derived from these published counts, so
  three of their copy numbers shift onto the published scale; both Vogler
  reference calls stay exact.
- `saureus_mlva8_schouls_2009` now calls four of its eight loci, up from
  two. Rechecking all eight against the full Pourcel panel confirmed two
  more shared arrays: VNTR67_01 is Sa0704 and VNTR81_01 is the coa locus
  Sa0266, anchored at the Mu50 counts 4 and 6. VNTR81_01 gives a whole
  allele in every sequenced genome tested; VNTR67_01 is anchored for the
  reference but its array mixes 67 and 66 bp units in some lineages, which
  are reported INTERMEDIATE rather than rounded.
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
- The two *Francisella* schemes were audited across 22 genomes, after an
  earlier audit reported seven that were really five: two of its accessions
  were not *Francisella* at all and two strain names were wrong. The genomes
  are now identified from their assembly headers, nine subspecies
  *tularensis*, three *holarctica*, one *mediasiatica* and nine *novicida*.
  Where a locus keeps one flank, its products differ by whole multiples of
  the unit; in the 13 genomes outside *novicida* this holds at 19 of the 25
  Johansson loci, and the six exceptions are accounted for in the scheme
  notes, among them the 30 bp *holarctica* deletion at Ft-M19 and the
  insertion at Ft-M20. Ft-M8 is displaced in four of the 13 genomes, at 325,
  1857 and 1871 kb. The Johansson and Vogler panels, which count the same
  arrays with different primers, give the same allele in all 145 comparisons
  at the nine loci other than Ft-M20, and differ at Ft-M20 because the
  insertion that Vogler splits off as Ft-M20B lengthens it.
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
- The aligned rescue (`--rescue-indels`) anchored on an exact match of one primer, so a product whose two primers both carried a mismatch was never found. It now anchors on a gap-free match with up to `-N` mismatches. On 40 *Burkholderia* genomes this recovers 14 products (3671k in 7 strains, 2445k in 7, and others) at plausible sizes, changes no existing hit and raises the run time of that screen from 19 s to 35 s.
- `burkholderia_uren_2007`: on 40 genomes the measured copy number gave alleles of 0 or below at five markers (3564k down to -7) and COUNT_DISAGREES on 39% of calls. 2115k, 1500k, 3564k, 1788k and 1217k now use the copy number of the publication, and sixteen markers are listed under `interrupted_arrays`. In the same genomes no allele is 0 or below and COUNT_DISAGREES falls to 6%. The K96243 profile changes at those five markers.
- Three scheme citations were checked against the DOI registry and corrected.
  The DOI of `blongum_mlva5_matamoros_2011` belonged to a different paper by
  other authors and is now `10.1016/j.mimet.2011.10.005`. The authors of the
  *Mycoplasma agalactiae* paper were given wrongly: they are McAuliffe,
  Churchward, Lawes, Loria, Ayling and Nicholas, and the scheme, which had not
  yet been released, is renamed from `magalactiae_mlva4_nouvel_2008` to
  `magalactiae_mlva4_mcauliffe_2008`. The Vogler paper is cited as 2009, the
  year of its issue, with its full author list; the scheme keeps the id
  `francisella_vogler_2008`.
- Several statements in the documentation no longer matched the code and are
  corrected. `docs/method.md` said that a marker without a unit length cannot
  be called, that `bacillus_mlva25_lista_2006` was the only scheme with a
  bin table, and that the Ft-M8 displacement is seen in one strain, and it
  omitted `known_multiple` from the scheme file keys.
  `docs/schemes.md` gains the third way of numbering alleles, for a locus with
  no single repeat unit, a list of the checks that have caught errors in
  bundled schemes, and checklist items for them. The roadmap and README
  describe the 26 bundled schemes and the way compound loci are now called.

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
