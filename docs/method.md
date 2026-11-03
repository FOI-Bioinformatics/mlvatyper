# Method and output

How mlvatyper locates amplicons, calls alleles and reports them. The
[README](../README.md) covers installation and day to day usage.

## Method

1. The STS file (id, forward primer, reverse primer, expected product size,
   optional alias) is validated. Primers may contain IUPAC ambiguity codes.
   Entries with other characters, invalid sizes or duplicate ids are
   reported and, with `--skip-invalid`, written to `invalid_primers.tsv`;
   `--fix-sts` instead removes stray characters from primers and logs each
   removal. The size may be a range `low-high`, which merPCR uses as a
   per-marker size window in place of `-M`.
   Both the STS file and the FASTA files may be gzip, bzip2 or xz compressed
   (recognised by the `.gz`, `.bz2` or `.xz` suffix). Text inputs (STS,
   marker, bin and scheme files) may start with the byte-order mark that
   Excel and Windows editors write.
2. Each FASTA file is treated as one sample, named after the file stem with
   any compression suffix removed (`sampleA.fasta.gz` gives `sampleA`);
   whitespace in the name becomes `_`, since FASTA record ids end at the
   first space. `--fasta` takes files and glob patterns; a file that does
   not exist stops the run, and a pattern that matches nothing is warned
   about. A sample whose file cannot be read (not FASTA, or a corrupt
   compressed file) fails on its own: the other samples are processed and
   typed, the failed one is left out of the results, and the command exits
   with status 1.
   Primer pairs are searched with merPCR, allowing `-N` mismatches per primer
   outside the protected 3' bases (`--three-prime-match`), with a hash word
   of `-W` bases (default 11) that must match exactly. The reverse primer
   is only searched within `-M` bp of the STS expected size, so the expected
   size and the margin together define the allowed amplicon size window.
   With `--markers`, a marker with `ref_size` and `ref_repeats` is searched
   in a window from the size of 0 to `--max-repeats` (50) repeats instead,
   unless the STS file gives a range. IUPAC codes in primers are honoured
   unless `--no-iupac` is given. merPCR drops a marker when one of its
   primers has no hash word, as me-PCR does; `--one-sided` searches it from
   the primer that has one instead.
3. Amplicons are extracted from the assemblies (reverse complemented for
   minus-strand hits). Like merPCR, the extraction upper-cases the sequence
   and keeps only `GATCNBDHKMRSVWY` and `-`, dropping for example `X`
   masking, before coordinates are counted, so that the coordinates in
   `hits.tsv` refer to the sequence without them; a warning gives the
   number dropped. Primer mismatches are taken from merPCR, which counts
   them with the rule the search applied, and the distance to the nearest
   contig end is recorded. Hits for one
   marker that share the forward primer position or overlap by more than a
   primer length are treated as one locus; the hit with fewest mismatches,
   then smallest deviation from the expected size, is kept and the number of
   discarded alternatives is recorded (`--no-dedup` keeps all). Amplicons are
   written to `hits.tsv` and to one FASTA per marker in `sequences/`. Record
   ids are `<sample>_<copy>`.
4. Markers without a hit in a sample get two rescue passes. A second merPCR
   pass with a smaller word size (`--rescue-wordsize`, default 6; 0
   disables), one-sided search and a wider window (`--rescue-margin`,
   default twice `-M`) finds primers whose IUPAC codes leave no 11-mer for
   the main pass. With `--rescue-indels`, one primer is then anchored by an
   exact IUPAC-aware match and the other aligned with up to one indel.
   Hits record their origin in the `search_pass` column (`main`, `rescue`,
   `aligned`) and the indel counts. Markers that stay missing are listed in
   `missing_markers.tsv` with a reason: NOT_SEARCHED (no primer holds a
   hash word), NEAR_CONTIG_END (a primer within the marker's expected size
   plus `-M` of a contig end, so the locus is more likely split than
   deleted), ONE_PRIMER_ONLY
   (one primer matches exactly somewhere) or NOT_FOUND.
5. If MAFFT or MUSCLE is on PATH, amplicons per marker are aligned and written
   as aligned FASTA and CLUSTAL (`.clw`) with the full sequence names;
   `--align mafft`, `muscle` or `none` overrides the automatic choice.
6. Summary tables, one row per sample and one column per marker, are written
   with the given prefix. `<prefix>_lengths_summary.tsv` holds amplicon
   lengths, several copies joined by `/`, with a column for every STS
   marker. With `--clwclean` the same table is
   also produced from amplicons without N (`*_clean.tsv`), and clean
   alignments are written.
7. Alleles are called when marker definitions are available (see Allele calling).

## Allele calling

The repeat number of an amplicon is anchored to a reference strain:

    raw = ref_repeats + (length - ref_size) / unit_length

where `ref_size` is the amplicon length in the reference strain (outer primer
edge to outer primer edge, as merPCR reports it) and `ref_repeats` the repeat
number in that strain. This is the convention of MLVAbank and of MLVA_finder.
A raw value closer than `--allele-tolerance` (default 0.25) to an integer is
that integer; a value exactly at the tolerance is not, as in MLVA_finder.
Other values depend on `--rounding`: `intermediate` (default)
reports `n.5` as MLVA_finder does (floor plus one half, so -0.4 gives
-0.5), `nearest` rounds halves up (2.5 gives 3), and `up` and `down` round
as named. In all cases the call is flagged INTERMEDIATE. Rounding conventions differ
between laboratories, so the policy is recorded in `run_info.json`.

Independently of the length, the longest tandem array of the repeat unit is
counted in each amplicon, allowing `--unit-mismatches` differences per unit
(default 1), in either orientation and at any rotation of the unit. The
count is reported next to the length based allele, and a difference of one
or more repeats is flagged COUNT_DISAGREES, also when the unit is not
found at all. This exposes flanking indels, wrong unit lengths and wrong
unit sequences, which a length alone cannot show. Markers without a unit
sequence are not counted; their tandem count is NA.

### Marker definitions

Marker definitions come from three sources, which are merged:

- `--markers FILE`: tab separated with a header. Recognised columns are
  `marker`, `patterns` (repeat unit sequences in the amplicon orientation,
  separated by `;`, or `N/A`), `unit_length`, `ref_size`, `ref_repeats`,
  `offset` and `label` (the locus name used in exports). Only `marker` is
  required; `unit_length` defaults to the length of the first pattern, and
  a column named `unit` may hold either the unit sequence or its length.
  Units must be DNA sequences (IUPAC codes allowed), and a marker defined
  twice is an error. The older pattern file (`marker`, `pattern`,
  `length`, then ignored columns) is read by the same code, also without
  a header, in which case the first line counts as data (the legacy script
  skipped it) and a stated length that differs from the pattern is only
  reported. `--pattern` is accepted as an alias.
- STS ids in the MLVAbank convention `NAME_<unit>bp_<size>bp_<n>u` (for
  example `Bruce06_134bp_408bp_3u`) carry unit length, reference size and
  reference repeats and need no file.
- `--reference-sample NAME`: for markers that still lack `ref_size` and
  `ref_repeats`, the amplicon length of that sample is used as `ref_size`
  and its tandem array count as `ref_repeats` (0 when no pattern is known,
  in which case alleles are offsets from the reference rather than repeat
  numbers). The derived definitions are written to
  `marker_defs_derived.tsv` for review and reuse with `--markers`. A
  reference sample without any amplicon is an error.

The definitions a summarize used, from all three sources, are written to
`<prefix>_definitions.tsv`. Markers without unit length cannot be called.
Without any definition only the lengths table is produced, and the call
tables an earlier run wrote under the same prefix are removed.

### Schemes

A scheme directory bundles `primers.sts`, `markers.tsv`, an optional
`bins.tsv` and `scheme.yml`, flat `key: value` metadata (name, organism,
reference, assembly, citation, url, retrieved, notes; `#` starts a comment
at the start of a line or after a space) that may also set the calling defaults
`rounding`, `allele_tolerance`, `unit_mismatches`, `margin`, `mismatches`
and `max_array_length`. `--scheme NAME` (a bundled scheme, listed by
`mlvatyper schemes`) or `--scheme DIR` replaces `--sts`, `--markers` and
`--bins`; `mlvatyper new-scheme DIR` writes the four files as annotated
templates, and [schemes.md](schemes.md) describes how to fill them from a
publication. An option given on the command line wins over the scheme
default, also when it is given at its default value. Schemes are bundled
under `src/mlvatyper/schemes/`. Published schemes should be added only with
reference values taken from their publications, the provenance fields
filled in, and a passing `check-scheme` against the assembly of the strain
the publication used.

Bundled schemes. `mlvatyper schemes` lists them with their organism, locus
count, callable locus count and reference assembly; `mlvatyper schemes NAME`
prints one scheme's loci, provenance, calling defaults and notes, and a
partial name or an organism filters the list. The notes in `scheme.yml` are
the full account of how a scheme's values were obtained, and the summaries
below are shorter than they are.

- `testM`: the example scheme of this repository, without reference
  repeat numbers.
- `STm_ECDC`: the five locus scheme for *Salmonella* Typhimurium
  harmonised in the EU (Lindstedt et al. 2004; nomenclature of Larsson et
  al. 2009; ECDC SOP). Reference sizes are the in silico products in strain
  LT2 (AE006468 and pSLT AE006471) and equal the LT2 sizes of the SOP;
  reference repeat numbers were derived from the SOP size ranges and
  confirmed by exact tandem counts in the LT2 amplicons, which the CI job
  repeats on the LT2 assembly (GCF_000006945.2, downloaded by accession).
  Units are
  counted exactly (`unit_mismatches: 0`) because the nomenclature counts
  complete repeats. STTR3 is a compound locus of 27 bp and 33 bp units
  that the single unit formula cannot express; it is reported by length
  and its allele has to be assigned by hand.

- `salmonella_enteritidis_mlva9_malorny_2008`: the nine locus multiplex for
  *Salmonella* Enteritidis (Malorny et al. 2008), anchored on NCTC 13349,
  which is the sequenced strain P125109. The loci are in the order the paper
  gives for the MLVA type, so a profile reads in the published order. Seven
  of the nine reference repeat numbers are fractional because the array ends
  in a partial unit, so those loci have alleles on a grid offset from the
  integers rather than INTERMEDIATE calls. Four of the published sequence
  lengths are not the product their own primers give in that assembly, and
  the measured length is used; the differences are a fraction of a repeat and
  move no allele. Note that the paper's MLVA type is a string of arbitrary
  allele indices, not repeat numbers, and the index to size table is not
  published, so a profile from this scheme cannot be compared with a
  published type directly.

- `Lpn_MLVA12`: the 12 locus MLVA panel of Orsay for *Legionella
  pneumophila* (Sobral et al. 2011, Appl Environ Microbiol 77:6899-6907,
  doi:10.1128/AEM.05556-11, Table 1; it contains the MLVA-8 markers of
  Pourcel et al. 2007, but the source does not say which eight). Reference
  sizes and copy numbers are those of Philadelphia-1 (AE017354, assembly
  GCF_000008485.1, which the CI job checks); all 12
  products match the source sizes exactly with no primer mismatches. No
  unit sequences are shipped, so alleles are called from amplicon length
  alone and the tandem count is not a cross-check; the arrays in
  Philadelphia-1 are diverged and a counter requiring near perfect units
  gives fewer copies than the scheme counts. In the genomes Paris, Lens and
  Corby the other 11 loci give whole repeat numbers, while Lpms31 gives
  numbers ending in .4 to .6 (INTERMEDIATE) because its length differences
  to Philadelphia-1 are not whole multiples of 45 bp; read that locus
  against the published allele table.

- `Bru_MLVA16`: the 16 locus *Brucella* panel (loci and primers from Le
  Fleche et al. 2006; the panel 1 / 2A / 2B organisation from Al Dahouk et
  al. 2007). Reference sizes and copy numbers are those of *B. melitensis*
  16M (AE008917 and AE008918, assembly GCF_000007125.1, which the CI job
  checks); all 16 products match the Table 1 sizes exactly and typing 16M
  reproduces its published copy numbers. As for `Lpn_MLVA12` no unit
  sequences are shipped, so alleles are called from length alone. The
  panel 2 microsatellites (for example Bruce19, a 6 bp unit) form long
  arrays that collapse in short-read assemblies; use `--max-array-length`
  or long reads and read those loci against the published allele table.

- `bacillus_mlva8_keim_2000`: the eight locus MLVA-8 panel for *Bacillus
  anthracis* (Keim et al. 2000), anchored on Ames Ancestor. The primer table
  carries no reference values, so these are the published MLVAbank values and
  are verified in the assembly.

- `bacillus_mlva15_vanert_2007`: MLVA-15 (Van Ert et al. 2007), the MLVA-8
  loci with the same primers plus seven BaVNTR loci, also on Ames Ancestor.

- `bacillus_mlva25_lista_2006`: MLVA-25 (Lista et al. 2006), the MLVA-8 loci
  with 17 bams loci, run as four multiplex reactions and anchored on Ames
  Ancestor. The paper publishes no primers, so these are the MLVAbank ones of
  `bacillus_mlva31_mlvabank`, of which this panel is the first 25 loci, but
  the repeat numbering is the paper's own and differs from the MLVAbank one at
  vrrA, vrrB1, vrrB2, pXO2-at and bams34; the reference sizes agree at every
  locus. This is the only bundled scheme with a `bins.tsv`: the paper numbers
  alleles from a table of sequenced amplicon size against repeat number, which
  ships here in full, so a profile is comparable with published MLVA-25
  genotypes while one from the other Bacillus schemes is not.

- `bacillus_mlva31_mlvabank`: the 31 locus composite panel of MLVAbank, the
  MLVA-8 loci with 17 bams loci (Lista et al. 2006) and six vntr loci. Every
  reference value is read from the MLVAbank locus identifiers and all 31 are
  reproduced in Ames Ancestor, so this scheme supplies the published values
  the two Bacillus schemes above use.

- `brucella_mlva15_lefleche_2006`: panel 1 and panel 2 of the Brucella
  MLVA-15 assay, anchored on *B. melitensis* 16M. Every locus reproduces the
  published reference size and copy number on the stated chromosome.

- `burkholderia_uren_2007`: the 32 locus panel for *Burkholderia
  pseudomallei*, with reference sizes and copy numbers measured in strain
  K96243 because the published copy numbers could not be reproduced there.

- `coxiella_mlva17_arricau_bouvery_2006`: panel 1 and panel 2 for *Coxiella
  burnetii*, anchored on Nine Mile RSA 493. Fourteen of the 17 published
  amplicon sizes are reproduced; ms30, ms31 and ms36 are recorded with their
  measured sizes, and ms30 and ms31 should not be compared with published
  allele numbers until their primers are confirmed.

- `francisella_johansson_2004` and `francisella_vogler_2008` for *Francisella
  tularensis*, anchored on SCHU S4 and LVS. The Johansson panel calls 16 of
  its 25 markers and sizes the rest, because the repeat units of the others
  are not in the source table; the Vogler table gives amplicon sizes only, so
  that scheme sizes its markers without calling them, though seven of its
  eleven units are now recorded from a cross strain test. The locus names of
  the two are the same VNTRs spelled differently, Ft-M02 against Ft-M2, but
  the primers and references differ so the alleles are not interchangeable.
  Five Johansson loci share a reverse primer and one repeat family and are
  selected by size alone; at Ft-M8 that is shown to pick a different family
  member in one strain, so read its allele with care.

- `mtb_miru24_supply_2006` and `mtb_miru15_supply_2006`: the standardised
  MIRU-VNTR set for *Mycobacterium tuberculosis* and its discriminatory
  subset, anchored on H37Rv. The source table gives no copy numbers, so these
  are derived from the reference; confirm them against the published allele
  calling table before comparing a profile with published MIRU codes.

- `nmeningitidis_mlva12_liao_2006`: the 12 locus MLVA for *Neisseria
  meningitidis* (Liao et al. 2006), anchored on Z2491. It is the best
  specified of the bundled schemes: one table gives the consensus repeat
  unit, its length and the repeat number in each of three sequenced genomes,
  another gives the primers and the predicted amplicon size in the same
  three, and the two agree exactly, so every value is published and none is
  derived. All twelve sizes reproduce in Z2491, and the same primers run
  against MC58 reproduce that genome's twelve sizes and eleven repeat
  numbers as well. NMTR12 is a compound locus of 13 bp and 12 bp units which
  the single unit formula cannot express, so it is reported by length only,
  as STTR3 is in `STm_ECDC`.

- `saureus_mlva5_sabat_2003`: the five locus gel based MLVA for
  *Staphylococcus aureus* (Sabat et al. 2003), anchored on NCTC 8325, the
  control strain the paper types. It is the thinnest bundled scheme: the
  paper publishes primers and nothing else quantitative, and its MLVA
  pattern numbers are arbitrary labels for band patterns rather than per
  locus alleles, so the scheme sizes its markers and calls none. `sdr` is
  one primer pair that amplifies the tandem *sdrC*, *sdrD* and *sdrE* genes
  together, so it is the only bundled locus declared in `known_multiple`;
  it gives three products in six of eight *S. aureus* genomes and two in the
  reference. The *sspA* primers are degenerate and are found by the rescue
  pass rather than at word size 11. Its *spa* locus overlaps
  `saureus_mlva8_schouls_2009` with different primers. The unit lengths were
  checked against eight genomes and three are corrected: the Ser-Asp arrays
  of *clfA*, *clfB* and *sdr* step by one dipeptide, 6 bp, not the 18 bp the
  paper describes.

- `saureus_mlva8_schouls_2009`: the eight locus MLVA for *Staphylococcus
  aureus* (Schouls et al. 2009), anchored on N315. All eight products are
  reproduced there and each equals the span of the coordinates the paper
  publishes for the same primer pair. This scheme sizes its markers but does
  not call them: the paper gives neither the flanking lengths nor the repeat
  unit sequences, and does not state the reference repeat number, so the one
  value that anchors an allele is missing. The arrays are too imperfect to
  count it back reliably, so no copy number is recorded rather than a
  doubtful one. Amplifying the panel in eight *S. aureus* genomes narrows
  each flanking length only to within a whole repeat, so it does not settle
  it either; supplying the flanking lengths or the published size to allele
  table would make it a calling scheme. That comparison did confirm four
  unit lengths, correct VNTR63_01 from 64 bp to 63 bp, and show that three
  loci fit no single period, among them VNTR61_01, whose flanking region the
  paper itself reports as heterogeneous. VNTR24_01 is the *spa* repeat
  region, so this panel is not independent of *spa* typing at that locus.

- `yersinia_klevytska_2001`: the 42 locus VNTR survey for *Yersinia pestis*,
  anchored on CO92. Three loci are numbered with a half copy.

- `yersinia_mlva25_lefleche_2001`: the 25 loci of Le Fleche et al. that are
  polymorphic among *Yersinia pestis*, the panel of Pourcel et al. 2004, also
  anchored on CO92. Unit length, repeat number and expected product are all
  published and every product is reproduced in the reference at the position
  the locus name encodes. Six loci are shared with
  `yersinia_klevytska_2001`, which numbers them differently, so profiles from
  the two schemes are not interchangeable.

- `yersinia_mlva26_li_2013`: the hierarchical 14+12 scheme of Li et al. for
  *Yersinia pestis*, anchored on CO92. The first fourteen loci place an
  isolate in a major population and the twelve that follow resolve it to a
  sub-population, so the panel is bundled in that order and the first
  fourteen digits of a profile are the first step result. Primers and
  reference values are the paper's supplementary Table S2, which encodes unit
  length, amplicon size and repeat number in each locus name; five of the
  encoded sizes are not the product their own primers give in CO92, and the
  measured size is used where another bundled scheme corroborates it with the
  same primer pair. The panel shares loci with both other *Yersinia* schemes
  but numbers several of them differently, so none of the three are
  interchangeable.

`mlvatyper check-scheme --scheme NAME --outdir DIR` first requires the
provenance of the scheme: `citation` (who published the values), `url`
(where the file came from) and `reference` (the strain or accession the
reference sizes and repeat numbers describe). Allele numbers are only
comparable with a database when their source is recorded, so an empty
field or a placeholder such as `none` or `TBD` fails the check, as does
a value that opens with one and explains it (`not recorded; ...`) or says
it is still to be filled in (`the DOI has to be added`); the
remaining checks still run and are written, and
`--allow-missing-provenance` turns the failure into a warning while a
scheme is being prepared. `organism` and `retrieved` are reported as
warnings only. The command then validates the files, compares marker and
bin ids with the STS ids and lists primers that merPCR cannot hash at
word size 11 or 6. When `scheme.yml` names the reference assembly by its `GCA_` or `GCF_`
accession (`assembly`), or `--assembly` or `--reference assembly.fasta`
is given, it also searches every marker and reports the number of
products, the product length against `ref_size` and the tandem count
against `ref_repeats` in `scheme_check.tsv`. Assemblies are downloaded
from NCBI once into `~/.cache/mlvatyper/assemblies` (or `MLVATYPER_CACHE`,
or `--cache-dir`) and are never stored with the scheme; `--no-download`
checks the files only, and `--fetch ACCESSION --email ADDRESS` fetches a
single nucleotide record instead of an assembly. Overlapping hits for one
locus count as one product, as in the search.

### Exchange with MLVA databases

`--mlvabank` (on `run` and `summarize`) writes `<prefix>_mlvabank.csv`:
one row per sample, a `Strain` column and one column per marker in scheme
order, named by the `label` column of the marker file when present.
Cells without an integer allele keep the status code or value as
reported, or are replaced by `--mlvabank-missing TOKEN`. Check the column
names against the target database before uploading.

`mlvatyper bins-from-table TABLE --output bins.tsv` converts a published
allele size table (one row per locus, allele and size or size range; tab
or comma separated; columns found by name, or given with
`--locus-column`, `--allele-column`, `--size-column`) to the bin file
format. Ranges may be written with any dash (`100-110`, `100–110`),
decimal sizes are rounded to whole base pairs (halves up), and rows
without an allele are skipped with a warning. `--markers FILE` maps the
table's locus names to marker ids through the `label` column.

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
  recorded in `run_info.json`. `--engine trf` without TRF on PATH is an
  error; `auto` then falls back to `builtin`.
- `builtin`: self-comparison at every period from `--min-unit` to
  `--max-unit` (3 to 100 bp), with the mismatch allowance of the tandem
  counter and the same scoring as TRF (+2 per matching and -7 per differing
  position between consecutive copies), without indel handling. Every
  separate array of a period is reported, as with TRF. Used when TRF is not
  on PATH.

Every array an engine reports becomes a candidate. Arrays that describe the
same locus at several periods are merged, keeping the smallest period whose
score is within 85 percent of the best, so that a unit is not given as its
own multiple. Each candidate is counted again with the tandem counter,
around its own array, so that `ref_repeats` is the integer allele calling
will see, and is placed in a category:

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

An optional bin file (`--bins`: tab separated marker, size or `low-high`
range, allele; a header and `#` comment lines are optional) maps amplicon
sizes directly to allele names, which need not be numbers, and takes
precedence over the formula, as some schemes (for example Brucella) number alleles by
convention. A size outside every bin for that marker is flagged OFF_TABLE
and falls back to the formula.

`--max-array-length L` censors calls whose repeat array is at least L bp
long, reported as `>=n` with CENSORED. Long tandem arrays collapse in short
read assemblies, so L should be the assembly k-mer or read length.

### Output tables

`search` (and the search step of `run`) writes, in `--outdir`:

- `hits.tsv`: one row per amplicon with its coordinates, length, expected
  size, primer mismatches and indels, discarded alternatives, distance to
  the contig end, search pass and sequence.
- `sequences/<marker>.fasta`: the amplicons of each marker, record ids
  `<sample>_<copy>`, and with an aligner `<marker>_aligned.fasta` and
  `<marker>_aligned.clw` (`_clean` variants with `--clwclean`).
- `missing_markers.tsv`: every sample and marker without a hit, with its
  reason; `near_contig_end.tsv` lists the NEAR_CONTIG_END ones.
- `samples.txt`, `primers_used.sts` (the validated STS entries, with the
  size windows from `--markers`) and `invalid_primers.tsv` (with
  `--skip-invalid`).

`summarize` (and the summary step of `run`) reads those and writes, under
`--prefix`; optional tables that an earlier run with the same prefix wrote
and this one does not are removed:

- `<prefix>_lengths_summary.tsv` (see step 6 above).
- `<prefix>_alleles.tsv`: sample by marker allele values. Several amplicons
  are joined with `/`. Missing markers show a status code instead of a
  number.
- `<prefix>_status.tsv`: sample by marker status codes joined by `;`.
- `<prefix>_tandem_counts.tsv`: tandem array counts per amplicon.
- `<prefix>_profile.tsv`: alleles in STS order plus a hyphen joined profile
  string per sample; a negative allele is written in parentheses there,
  `5-(-1)-3`, so that the string stays unambiguous.
- `<prefix>_calls_long.tsv`: one row per amplicon (and per missing marker)
  with coordinates, length, expected size, primer mismatches, discarded
  alternatives, distance to contig end, primer indels and search pass, raw
  repeat number, allele, residual
  in bp between the length and the expected size of the called allele,
  tandem array count (NA for a marker without repeat unit sequences), and
  status.
- `<prefix>_definitions.tsv`: the marker definitions the calls were made
  with, after merging the marker file, MLVAbank style ids and the reference
  sample; the report shows this table.
- `run_info.json`: versions, all parameters of the search and summary
  steps (the summary parameters also per prefix, so that several prefixes
  can share an output directory), and the size and MD5 sum of every input
  file, including the marker and bin files.
- `<prefix>_report.html`: a self-contained report written by `run` (or by
  `mlvatyper report --outdir DIR --prefix P`) with a status summary, the
  allele profile table shaded by status, figures of amplicon lengths,
  residuals and status per sample and marker, the evidence per amplicon,
  the missing markers with their reasons, the discovery report when
  present, a methods paragraph generated from the recorded parameters,
  the input checksums and the marker definitions. Styles and figures
  (SVG) are embedded, so the file opens offline and prints to PDF from a
  browser. `--no-report` skips it, and `--title` sets its title. See
  [report.md](report.md).

`discover` writes `discover_report.tsv`, `discover_candidates.tsv` and
`marker_defs_discovered.tsv` (see Discovering repeat units), and
`check-scheme` writes `scheme_check.tsv`.

Status codes: OK, NOT_DETECTED (no hit), NEAR_CONTIG_END (no hit but a primer
within the marker's expected size plus `-M` of a contig end), NOT_SEARCHED and ONE_PRIMER_ONLY
(from `missing_markers.tsv`), INTERMEDIATE, NEGATIVE (fewer bases than the
flanking region), BINNED, OFF_TABLE, MULTIPLE (several loci in the sample),
AMBIGUOUS (overlapping hits were collapsed), PRIMER_INDEL (hit from the
indel rescue), CENSORED, HAS_N, COUNT_DISAGREES and NO_DEFINITION (an
amplicon was found, but the marker has no unit length or reference
values). The codes are text so that they cannot be confused with the
allele 0, which is a valid call.

`--legacy-tables` additionally writes the two tables of the original script:
exact non-overlapping pattern counts (`<prefix>_count_summary.tsv`) and the
length difference to the shortest amplicon in the run divided by the unit
length (`<prefix>_count_lengthdiff.tsv`). The latter is relative to the
dataset and not comparable between runs; both are kept for comparison only.

## Markers without hits

A marker without an amplicon in a sample is still a column in the summary
tables, reported as NA in the lengths table and with the reason from
`missing_markers.tsv` in the allele and status tables (NOT_FOUND appears
there as NOT_DETECTED). The reasons are:

- NOT_SEARCHED: neither primer holds a run of unambiguous bases as long as
  the word size of the most permissive merPCR pass (`-W`, or
  `--rescue-wordsize` when the rescue pass ran). merPCR logs "primers have
  ambiguities which prevent computation of a hash value" for such primers.
- NEAR_CONTIG_END: a primer matches within the marker's expected size plus
  `-M` of a contig end; the locus is more likely split by the assembly than
  deleted.
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

- me-PCR replaced by merPCR; seqret replaced by a built-in CLUSTAL writer.
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
- A pattern file without a header keeps its first line; the legacy script
  dropped it as a header.

