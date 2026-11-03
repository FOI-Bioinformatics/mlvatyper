# Run report

Every run writes `<prefix>_report.html` next to the tables, unless
`--no-report` is given. The report is built only from files in the output
directory, so it can be rebuilt after editing marker definitions or
changing the rounding settings without a new search:

    mlvatyper summarize --outdir out_testM --prefix testM --rounding nearest ...
    mlvatyper report --outdir out_testM --prefix testM

Several prefixes can share an output directory: each report shows the
tables, parameters and marker definitions (`<prefix>_definitions.tsv`) of
its own prefix. `--title` replaces the default title.

Styles and figures are embedded and nothing is loaded over the network, so
the file opens on a machine without internet access and prints to PDF from
a browser. The pictures below are from a run of the bundled
`burkholderia_uren_2007` scheme on eight genomes; `examples/testM/testM_report.html`
is an older report with the first layout.

![Top of the report: scheme and publication, completeness, across-sample
summary and loci to review](images/report.png)

![Clustered heatmap of the allele profiles with the tree, and the pairwise
distance matrix](images/report_heatmap.png)

## Sections

1. **Header**: scheme, samples, markers, amplicons, reference sample,
   primer and definition files, search date and output directory, with
   links to the sections. Versions of mlvatyper and merPCR are in the line
   above.
2. **Summary**: per sample, the markers typed (an integer allele) and
   missing (no amplicon; a marker without definition is not missing), then
   the status code counts, followed by one line per flagged code naming
   the markers.
   For five or more samples a subsection **Across samples** follows: the mean
   and lowest share of loci with status OK, the samples below 75%, the loci
   to review, a per-locus table (typed, OK, missing, most frequent flag and
   its share, distinct alleles, Simpson index, lowest allele) and the missing
   products by reason. A locus is listed for review when one flag is on more
   than half of its typed calls or an allele is zero or below, which points at
   the scheme (reference values, repeat unit, patterns) more often than at the
   strains. The same tables are written as `<prefix>_sample_summary.tsv` and
   `<prefix>_locus_summary.tsv`.
3. **Allele profiles**: repeat numbers in STS order, shaded by status with
   a colour legend. The table is transposed when markers outnumber samples.
   The hyphen joined profile strings are in a collapsed block; a negative
   allele appears there in parentheses, `5-(-1)-3`.
4. **Figures**: amplicon length per marker against the expected size,
   residual per marker, and a status matrix of samples by markers (ordered
   like the heatmap below for five or more samples). For five or more samples
   there are six more: a **clustered heatmap of the allele profiles** (rows
   ordered by average linkage clustering, with the tree; each locus on its own
   colour scale; repeat numbers printed in the cells of small matrices), a
   **pairwise distance matrix** in the same order, an **allele spectrum** per
   locus, **call rate bars** per locus and per sample, a **minimum spanning tree** of
   the distinct profiles (identical profiles share a node, edge labels give the
   number of loci that differ, node colour is the share of OK calls, and the
   layout approximates path lengths along the tree) and a **classical
   multidimensional scaling** of the distances in two dimensions, whose axis
   labels give the share of the distance structure each axis keeps. The distance is the
   share of loci that differ among those typed in both samples; pairs with
   fewer than eight shared loci take the mean distance. Alleles are compared
   as categories, so the tree describes how many loci differ and not how far
   the repeat numbers are apart. It is a view of the data, not a phylogeny.
5. **Evidence per amplicon**: one collapsible table per sample with
   coordinates, length, expected size, residual, primer mismatches and
   indels, the search pass, the allele and the tandem array count. Rows
   carrying a flag are shaded.
6. **Markers without amplicon**: the reason established by the search, with
   its explanation.
7. **Repeat unit discovery**, when `discover` has run: the primary
   candidate per marker with its category and the cross-sample evidence.
8. **Methods**: a paragraph generated from the recorded parameters, ready
   to adapt for a manuscript, a table of software with citations, the
   command line, the input files with sizes and MD5 sums, and the full
   `run_info.json`.
9. **Appendix**: the marker definitions used for this prefix and the status
   code glossary.

## Layout conventions

- The header table names the scheme, its organism, its publication with the
  DOI and its reference strain when the run used a bundled scheme; the Methods
  paragraph cites the same publication.
- Tables are numbered in order, with the caption above the table; figures are
  numbered in order, with the caption below. Long per-sample tables are
  collapsed when there are more than 20 samples.
- Header rows and the first column stay in view when a wide table scrolls.
  Numbers are right aligned, and allele cells are centred.
- In the allele table, codes that stand in place of an allele are abbreviated
  (`ND` not detected, `CE` near contig end, `1P` one primer only, `NS` not
  searched, `NDef` no definition, `OFF` off table) so that columns keep one
  width; the key is above the table and the full code is the cell's tooltip.
- The status figure has a legend, grid lines and axis titles.
- Printing from a browser uses landscape A4, open sections and smaller
  tables; the navigation is left out.

## Reading the colours

The same colours are used in the allele table, the status figure and the
evidence table.

| colour | meaning |
|---|---|
| white | OK |
| yellow | several amplicons for one marker, or overlapping hits collapsed |
| orange | intermediate repeat number, or the tandem count disagrees with the length |
| blue | primer indel, or the amplicon contains N |
| green | allele taken from the bin table |
| grey | no amplicon, or no definition for the marker |
| red | outside every bin, negative, or censored |

A colour never replaces the status code: every flagged cell also carries
its code in the status table and in the evidence table, and the codes are
explained in the appendix of each report.

## What the report does not do

It shows what the run produced; it does not judge whether the calls are
correct. A profile should be compared with isolates typed by conventional
MLVA before it is used against published profiles, and the residual figure
is the first place to look when a marker's unit length or reference size is
in doubt.
