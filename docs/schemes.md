# Building a scheme from a publication

A scheme is a directory with four files that mlvatyper reads through
`--scheme DIR`:

| file | holds | required |
|---|---|---|
| `primers.sts` | marker id, forward primer, reverse primer, expected size | yes |
| `markers.tsv` | repeat unit, unit length, reference size and repeats, export label | recommended |
| `bins.tsv` | amplicon size to allele table | only for schemes numbered by a size table |
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
  (section 5).
- `label`: the locus name used by the target database, for example the
  MLVAbank column name, when it differs from the id.

Marker ids in the MLVAbank style `Name_<unit>bp_<size>bp_<n>u` carry
unit length, reference size and repeats in the name and need no values
in the columns.

## 4. Decide how alleles are numbered

Publications number alleles in one of two ways, and the scheme must say
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
  `--allele-column` and `--size-column`. Alleles from the bin table take
  precedence over the formula; sizes outside every bin are flagged
  OFF_TABLE and fall back to the formula.

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
| `n_products` | 1 | 0: a primer is wrong, or the strain is not the one the primers were designed on; 2 or more: the locus is duplicated |
| `size_ok` | 1 | the product differs from `ref_size`: a size convention mismatch (section 2) or a transcription error |
| `repeats_ok` | 1 | the tandem count differs from `ref_repeats`: the unit is wrong, in the wrong orientation, or the reference strain holds variant units |

The most common finding is a constant offset between the product and the
stated size across all markers, for example every product one base
longer: the paper measured the size by another convention. Correct the
sizes rather than the primers. A single marker off by a multiple of its
unit length means the downloaded assembly is not the strain the paper
used, or the strain carries a different allele; note it in `scheme.yml`.

When the paper gives no reference repeats, run the reference assembly
through mlvatyper and take the values from it:

    mlvatyper run --scheme schemes/organism_2015 --fasta ref.fasta --outdir ref_run \
        --prefix ref --reference-sample ref

`ref_run/marker_defs_derived.tsv` then holds `ref_size` and
`ref_repeats` for every marker that was found; copy them into
`markers.tsv` after reading the tandem counts in `ref_calls_long.tsv`,
where a count of 1 or 0 means the unit did not match and needs review.
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
    notes: sizes in the paper exclude the reverse primer; corrected by 21 bp. Locus_7 is duplicated in CP000000.1.

Write in `notes` every correction made in sections 2 to 5. The next
person to use the scheme, including you a year later, can then tell a
deliberate change from a mistake.

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
- allele numbering decided: formula or `bins.tsv`, rounding stated
- `check-scheme --reference` run, every marker with one product of the
  stated size and repeat count, exceptions noted
- citation, url, reference, assembly accession and retrieval date in
  `scheme.yml`; no FASTA file in the scheme directory
- agreement with conventional MLVA checked on isolates with known profiles

Schemes that pass this list can be contributed to the repository under
`src/mlvatyper/schemes/`; the CI job then downloads the assembly of every
bundled scheme by its accession and requires one product of the reference
size for every marker at each change.
