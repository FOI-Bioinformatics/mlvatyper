# Building a scheme from a publication

A scheme is a directory with four files that mlvatyper reads through
`--scheme DIR`:

| file | holds | required |
|---|---|---|
| `primers.sts` | marker id, forward primer, reverse primer, expected size | yes |
| `markers.tsv` | repeat unit, unit length, reference size and repeats, export label | recommended |
| `bins.tsv` | amplicon size to allele table | for schemes numbered by a size table, and for loci with no single repeat unit |
| `scheme.yml` | name, organism, provenance, calling defaults | yes |

This page describes how to fill them from an MLVA publication and how to
verify the result before the scheme is used. The order matters: verify
against the reference strain first, then type isolates.

## 1. Create the directory

    mlvatyper new-scheme schemes/organism_2015 --name organism_2015

writes the four files with the column headers and a comment on every
field. Comment lines start with `#` and may stay in place. Name the scheme
after the organism and the publication year, so that two schemes for the
same organism cannot be confused.

## 2. Take the primers from the supplementary table, not the PDF text

Most MLVA papers publish the primers in a supplementary spreadsheet or in
a table in the methods section. Use the spreadsheet when one exists:
copying from a PDF loses characters at line breaks, turns `0` into `O`
and swaps the forward and reverse columns of wide tables. When only the
PDF is available, copy each primer, then compare the primer length with
the length the paper states and read the sequence back once against the
table.

A table printed sideways, or as an image, is the worst case. Read it both ways, from the text layer of
the PDF and from the rendered page, and do not trust either until `check-scheme` has searched the
reference: a primer that does not amplify its own reference strain is a transcription error until
shown otherwise. In one scheme the primers read from the image of a rotated table failed to amplify at
four of eight loci, while the text layer gave primers that amplified at all eight, seven of them at the
published size to the base.

In `primers.sts`, one line per marker, tab separated:

    Locus_1    GAGCTGGTCAAGTTTATTTAAGTA    CAGGATTGCTTGAACATGATA    310

- The id is the name every output uses. Keep the publication's locus
  name without spaces; the `label` column of `markers.tsv` carries the
  exact database name if it differs.
- Primers are written 5' to 3' as published. IUPAC ambiguity codes are
  accepted. A character that is not a base, such as a footnote mark, is
  reported by `check-scheme`; `--fix-sts` removes it during a run.
- The expected size is the product of the reference strain, measured from
  the outer edge of one primer to the outer edge of the other. When the
  paper gives a size range for the locus, write it as `low-high`. A wrong
  size does not change the allele; it narrows the window in which the
  reverse primer is searched, so an amplicon far outside it is missed.

Some papers state the product size without the primers (the repeat region
only) or with one primer excluded. Section 5 shows how to find out.

## 3. Fill markers.tsv

One line per marker with the same id:

    marker    unit_length    ref_size    ref_repeats    patterns    label
    Locus_1   6              310         5              GATCGA      Locus_1

- `patterns`: the repeat unit as published, in the orientation of the
  forward primer. Several variants are separated by `;`. When the paper
  gives only the unit length, write `N/A` and fill `unit_length`; the
  tandem array cross-check is then skipped for that marker.
- `ref_size` and `ref_repeats`: the amplicon size and the repeat number
  of the reference strain, which the paper usually lists for its type
  strain. These two numbers anchor the allele formula, so a mistake in
  either shifts every allele of that marker by a constant. Leave them
  empty when the paper does not give them, and anchor the calls on an
  assembly of the reference strain with `--reference-sample` instead
  (section 5). `ref_repeats` may be fractional when the reference array
  ends in a partial unit and the paper numbers it `n.5`, as several
  Yersinia pestis loci are. The alleles of such a marker then lie on a
  half-integer grid, and a call on that grid is reported as `n.5` without
  being flagged INTERMEDIATE. Do not round these to whole copies: that
  would shift every allele of the marker by half a repeat.
- `label`: the locus name used by the target database, for example the
  MLVAbank column name, when it differs from the id.

Each marker appears once; a second line with the same id is an error, as
is a unit that is not a DNA sequence. Save files from a spreadsheet as
UTF-8 tab separated text; the byte-order mark Excel adds to UTF-8 files
is ignored.

Marker ids in the MLVAbank style `Name_<unit>bp_<size>bp_<n>u` carry
unit length, reference size and repeats in the name and need no values
in the columns.

## 4. Decide how alleles are numbered

Publications number alleles in one of three ways, and the scheme must say
which:

- **By the formula**: the allele is the repeat number, computed from the
  amplicon size, the reference size and the reference repeats. Leave
  `bins.tsv` out. This is the convention of MLVAbank and of most schemes
  designed on genome sequences.
- **By a size table**: each allele is a name assigned to a size or a size
  range, often with capillary sizing offsets built in, as in several
  Brucella and Salmonella schemes. Convert the paper's table with

      mlvatyper bins-from-table allele_sizes.csv --output schemes/organism_2015/bins.tsv

  The table needs locus, allele and size columns; column names are
  matched by common aliases or given with `--locus-column`,
  `--allele-column` and `--size-column`. Ranges may use any dash, as
  copied from a PDF (`100–110`), decimal sizes are rounded to whole base
  pairs, and rows without an allele are skipped with a warning; read the
  log for those. Alleles from the bin table take precedence over the
  formula and need not be numbers; sizes outside every bin are flagged
  OFF_TABLE and fall back to the formula.
- **By a table alone, for a locus with no single repeat unit**: a compound
  locus interleaves two repeat unit sizes, for example the 27 and 33 bp
  units of STTR3, so no one unit length can express its array and the
  formula cannot be used. The same holds when a locus has two states that
  do not lie on one repeat ladder, such as an insertion that is present in
  some strains and absent in others. Give the marker a `bins.tsv` and no
  `ref_repeats`. A size in the table gets that allele and is flagged
  BINNED; any other size gets no allele and OFF_TABLE, which is the
  honest result, since a formula that does not apply would give a wrong
  number. Cover as many sizes as the data support, state in the notes which
  genomes the table was built from, and say how a size that is not listed
  should be added.

When the publication gives no count at all, as in a gel based fingerprint,
the scheme may define one, but it must be a rule that needs no outside
value and it must be recorded as a convention in `notes` and in the comments
of `markers.tsv`. Alleles from such a locus compare only with each other and
not with any published scale.

Record the rounding convention of the paper in `scheme.yml`: `rounding:
nearest`, `up`, `down` or `intermediate` (the default, which reports
`n.5` for a value between two integers rather than choosing one).

## 5. Verify against the reference strain

Find the assembly accession of the strain the paper used as reference (a
GenBank `GCA_` or RefSeq `GCF_` accession with its version, from the NCBI
Assembly or Datasets pages) and write it into `scheme.yml`:

    assembly: GCF_000006945.2

Then run

    mlvatyper check-scheme --scheme schemes/organism_2015 --outdir check

The assembly is downloaded from NCBI into `~/.cache/mlvatyper/assemblies`
(or the directory in `MLVATYPER_CACHE`, or `--cache-dir`) the first time
and read from there afterwards. The sequence is never copied into the
scheme directory, so a scheme carries an accession and not a genome, and
the repository stays small. `--assembly GCF_...` checks against another
accession, `--reference ref.fasta` against a local file, `--no-download`
checks the files only, and `--fetch CP000000.1 --email you@example.org`
fetches a single nucleotide record for strains without an assembly entry.

`check/scheme_check.tsv` has one row per marker. Read the columns in this
order:

| column | expected | if not |
|---|---|---|
| `sts_ok` | 1 | the STS line is invalid; the note says why |
| `hash_word_11` | 1 | the primer has IUPAC codes in every 11-mer; it is still found by the rescue pass |
| `n_products` | 1 | 0: a primer is wrong, or the strain is not the one the primers were designed on; 2 or more: the locus is duplicated, and if that is real rather than a primer fault, list the marker in `known_multiple` |
| `size_ok` | 1 | the product differs from `ref_size`: a size convention mismatch (section 2) or a transcription error |
| `repeats_ok` | 1 | the tandem count differs from `ref_repeats`: the unit is wrong, in the wrong orientation, or the reference strain holds variant units |

The most common finding is a constant offset between the product and the
stated size across all markers, for example every product one base
longer: the paper measured the size by another convention. Correct the
sizes rather than the primers. A single marker off by a multiple of its
unit length means the downloaded assembly is not the strain the paper
used, or the strain carries a different allele; note it in `scheme.yml`.

### Checks that have caught errors

`check-scheme` confirms that the primers reproduce the stated sizes in one
assembly. These further checks have each found a mistake in a bundled scheme,
and are worth making for a new one.

- **Read the header of the assembly.** `check-scheme` cannot know that an
  accession is the strain the paper used. An earlier audit of two schemes
  named genomes it had not looked at, and two of its accessions were
  unrelated bacteria. Print the first line of each FASTA file and read it.
- **Reproduce the published reference profile**, not only the sizes: a run on
  the reference strain should return the paper's allele string, in the
  paper's locus order.
- **Verify the bibliographic details against the DOI registry.** A DOI
  written from memory pointed at a different paper by other authors, and a
  first author was wrong in a scheme's name. Look the DOI up and compare the
  authors, year, journal and pages with the paper in front of you.
- **Amplify more strains.** A locus keeps one flanking length, so its
  products in different strains should differ by whole multiples of the
  repeat unit. A locus that does not follow this has an indel in the flank,
  an insertion present in some lineages, a duplicated copy or a wrong unit;
  the test cannot tell these apart, so each exception needs an explanation
  in the notes and not just a mention. Use genomes from the whole species and
  record which subspecies they belong to.
- **Cross-check against a sibling scheme.** Where two schemes read the same
  array with different primers, run both on the same genomes: they should
  give the same allele wherever the array is the same, and the places they
  differ are informative. This is also how one scheme's counts can be taken
  over by another, by measuring the second scheme's product in a genome the
  first one counts.

When the paper gives no reference repeats, run the reference assembly
through mlvatyper and take the values from it:

    mlvatyper run --scheme schemes/organism_2015 --fasta ref.fasta --outdir ref_run \
        --prefix ref --reference-sample ref

`ref_run/marker_defs_derived.tsv` then holds `ref_size` and
`ref_repeats` for every marker that was found; copy them into
`markers.tsv` after reading the tandem counts in `ref_calls_long.tsv`,
where a count of 1 or 0 means the unit did not match and needs review
(a count of 0 is also flagged COUNT_DISAGREES).
When the paper gives no unit at all, `mlvatyper discover` on the same
run proposes one per marker with a category (see
[method.md](method.md)); accept category A and B units and check the
rest against the paper.

## 6. Record the provenance

`scheme.yml` must state where the values came from, because
`check-scheme` refuses a scheme without it:

    name: organism_2015
    organism: Organism species
    citation: Author A, Author B. Title. Journal 2015;12:345-356. doi:10.0000/example
    url: https://doi.org/10.0000/example (supplementary table S2)
    reference: strain XYZ, chromosome CP000000.1
    assembly: GCA_000000000.1
    retrieved: 2026-10-24
    known_multiple: Locus_7
    notes: sizes in the paper exclude the reverse primer; corrected by 21 bp. Locus_7 is duplicated in CP000000.1.

Write in `notes` every correction made in sections 2 to 5. The next
person to use the scheme, including you a year later, can then tell a
deliberate change from a mistake.

Add the scheme to [scheme_references.md](scheme_references.md), with its
publication, its DOI, the reference strain and a line on the basis of its
values; a test fails until every bundled scheme and every DOI it cites is
there.

`known_multiple` lists the markers that give more than one product in the
reference because the locus really is duplicated there, separated by
commas. Checks that require one product per marker skip them, so a
published locus can be kept with its ambiguity recorded rather than
dropped from the scheme. Calls for such a marker are flagged MULTIPLE.
Use it only when the duplication is established, not to silence a primer
that matches the wrong place.

A value is one line. `#` after a space starts a comment, so write a URL
fragment without a space before it (`page#table2`). `check-scheme` fails
when `citation`, `url` or `reference` is empty or a placeholder such as
`TBD`, and also when the value only explains its absence (`not recorded;
the DOI has to be added`): find the source before the scheme is used.

## 7. Type isolates and compare

    mlvatyper run --scheme schemes/organism_2015 --fasta 'isolates/*.fasta' \
        --outdir typed --prefix study --mlvabank

The report in `typed/study_report.html` summarises what was found and
flagged. Before the profiles are compared with published ones, type a
few isolates with known conventional MLVA results and check the
agreement per marker; a marker that disagrees by a constant has a wrong
reference value, and a marker that disagrees by varying amounts has a
wrong unit length or a size table with sizing offsets that the formula
does not reproduce.

## Checklist

- primers from the supplementary table, read back once
- sizes in the outer-edge convention, or corrected and noted
- unit sequences in forward primer orientation, variants separated by `;`
- reference size and repeats from the paper or from the reference assembly
- allele numbering decided: formula, `bins.tsv` or both, rounding stated, and any count that is a
  convention of the scheme recorded as one
- header of the downloaded assembly read, and the published reference profile reproduced
- DOI and bibliographic details checked against the registry
- `check-scheme` run against the reference assembly, every marker with one
  product of the stated size and repeat count, exceptions noted
- citation, url, reference, assembly accession and retrieval date in
  `scheme.yml`; no FASTA file in the scheme directory
- agreement with conventional MLVA checked on isolates with known profiles

Schemes that pass this list can be contributed to the repository under
`src/mlvatyper/schemes/`; the CI job then downloads the assembly of every
bundled scheme by its accession and requires one product of the reference
size for every marker at each change.
