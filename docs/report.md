# Run report

Every run writes `<prefix>_report.html` next to the tables, unless
`--no-report` is given. The report is built only from files in the output
directory, so it can be rebuilt after editing marker definitions or
changing the rounding settings without a new search:

    mlvatyper report --outdir out_testM --prefix testM

Styles and figures are embedded and nothing is loaded over the network, so
the file opens on a machine without internet access and prints to PDF from
a browser. The example below is `examples/testM/testM_report.html`.

![Top of the report: metadata, section links, status summary and the allele
profile table](images/report.png)

## Sections

1. **Header**: scheme, samples, markers, amplicons, reference sample,
   primer and definition files, search date and output directory, with
   links to the sections. Versions of mlvatyper and merPCR are in the line
   above.
2. **Summary**: markers typed and missing per sample, then the status code
   counts, followed by one line per flagged code naming the markers.
3. **Allele profiles**: repeat numbers in STS order, shaded by status with
   a colour legend. The table is transposed when markers outnumber samples.
   The hyphen joined profile strings for MLVAbank are in a collapsed block.
4. **Figures**: amplicon length per marker against the expected size,
   residual per marker, and a status matrix of samples by markers.
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
9. **Appendix**: the marker definitions used and the status code glossary.

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
