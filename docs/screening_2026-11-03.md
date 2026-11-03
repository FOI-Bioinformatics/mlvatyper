# Screening evaluation: `burkholderia_uren_2007` on 40 genomes (2026-11-03)

This page records a first check of how a bundled scheme behaves on a set of
genomes it was not built on, and what that changed. It is meant as a method
for repeating the check with other schemes and with larger sets.

## Set and method

38 distinct genomes and 2 duplicates under other file names, all RefSeq
(GCF) except the K96243 control (GCA_000011545.1):

| Group | Genomes | Notes |
| --- | --- | --- |
| Control | 1 | K96243, the reference of the scheme |
| *B. pseudomallei*, complete | 20 | chosen to spread over countries |
| *B. pseudomallei*, draft | 12 | 4 to 445 contigs, N50 from 28 to 230 kb |
| *B. mallei* | 5 | 2 complete or chromosome level, 3 draft |
| Duplicates | 2 | same file as a genome above, to test reproducibility |

`mlvatyper --log-json log.jsonl run --scheme burkholderia_uren_2007 --fasta
'genomes/*.fasta' --jobs 4 --align none`. The headers of all FASTA files name
the intended organism. There is no external truth: the publication gives no
per-strain profiles. The checks are therefore internal: call rate, status
codes, agreement between the length based allele and the tandem count,
reproducibility, and a manual look at missing products. The analysis script
used reads `_status.tsv`, `_alleles.tsv`, `_calls_long.tsv`, `hits.tsv` and
`missing_markers.tsv`. The `--align mafft` comparison was not run because MAFFT
was not installed on the machine used.

## Findings from the first run, and what changed

| Finding | First run | After changes |
| --- | --- | --- |
| Alleles of 0 or below | 38 calls at 5 loci (3564k down to -7) | none |
| COUNT_DISAGREES | 39% of calls | 7% |
| Loci OK per genome (mean, of 32) | about 17 | 26 |
| Missing products | 61 | 47 |
| Run time, 40 genomes, 4 jobs | 19 s | 35 s |

1. **Alleles of 0 or below.** At five loci (2115k, 1500k, 3564k, 1788k, 1217k)
   the copy number measured in K96243 was lower than the arrays of other
   strains need. These five loci now use the copy number of the publication
   (15, 6, 30, 8, 15). With them, the shortest product seen at each locus
   gives an allele of 2 or more. The copy numbers of the other 27 loci were not
   contradicted by the screen and are unchanged. K96243 profiles made before
   this change differ at these five positions.
2. **COUNT_DISAGREES.** At 16 loci the arrays hold variant units that the
   repeat patterns do not follow, so the tandem count and the length based
   allele differ in most strains. The scheme lists them under
   `interrupted_arrays`: the tandem count is reported and the difference is not
   flagged. This removes a cross-check from half of the loci, which the scheme
   notes state. 2445k is now flagged in 52% of calls and 2356k and 140k in 45%
   and 38%; they are candidates for the same list.
3. **Aligned rescue.** `--rescue-indels` anchored on an exact match of one
   primer, so products with a mismatch in both primers were not found. It now
   accepts up to `-N` mismatches at the anchor. This recovers 14 products
   (3671k in 7 genomes, 2445k in 7 genomes, 2341k, 3145k) at sizes inside the
   expected window and changes none of the other 1,221 hits.
4. **NEAR_CONTIG_END.** `missing_markers.tsv` has a `detail` column that states
   where each primer lies relative to the contig ends.

## Result after the changes

- 40 of 40 samples processed, peak memory 236 MB. K96243 calls all 32 loci OK.
- Both duplicates give profiles and statuses identical to the original.
- An allele is reported at 31 of 32 loci on average (21 at the lowest). All 38
  distinct genomes have distinct profiles; 15 are typed at every locus.
- Per locus, 4 to 21 distinct alleles (mean Simpson index 0.83); no locus is
  fixed.
- INTERMEDIATE is raised for 81 calls (6%), mostly at 2518k, 1217k and 2971k.

## The 47 missing products

| Cause | Products | Remedy |
| --- | --- | --- |
| Primer near a contig end, usually the two primers on different contigs (an array that broke the assembly) | 17, all in draft assemblies; 11 in one genome (445 contigs) | none with this tool; type from reads or a better assembly |
| One primer found, the other too diverged or the product outside the size window | 20 | `-N`, `-M`, or an adjusted primer; 1217k, 2815k and 20k have an array of 1.5 to 1.7 kb against a 250 to 480 bp window |
| Neither primer found | 10 (3 in *B. mallei*) | sequence divergence or deletion; not checked against annotation |

In drafts, call rate falls weakly with contig count (rank correlation -0.29)
and rises with N50 (+0.30).

## Limits of this evaluation

- No external truth for alleles. A call can be internally consistent and
  still differ from the published allele name.
- 40 genomes, mostly from a few countries; rare lineages are not covered.
- `--align mafft` was not compared.
- The copy numbers taken from the publication at five loci rest on the
  assumption that the published value counts the same array; the screen shows
  that the measured value cannot be right there, not that the published one is.

## Repeating the check

Select genomes with `datasets summary genome taxon`, record a manifest (level,
contigs, N50, country), include duplicates and the reference, run as above and
compare per locus: call rate, share of calls with COUNT_DISAGREES or
INTERMEDIATE, alleles of 0 or below, and the cause of each missing product.
Treat a locus where more than half of the calls disagree with the tandem count,
or where an allele is 0 or below, as a scheme problem and not as strain
variation.
