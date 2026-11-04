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
`missing_markers.tsv`. The run was repeated with `--align mafft` (MAFFT 7.526):
the allele, status, residual and tandem count of all 1,235 calls with a product
are identical to the run without alignment, and the run time rises from 35 s to
63 s.

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
- The alignment step was compared on this set only; it does not feed the calls.
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

## All bundled schemes on 12 genomes each

The same method was applied to the 24 other published schemes, with up to 12
complete RefSeq genomes of the scheme's species or genus chosen at random from
NCBI (a fixed seed; the accession lists are not kept in the repository). This
is a first look: 12 genomes show gross problems, not rare ones, and the genomes
need not belong to the lineage the scheme was built on. Call rate is the
share of loci with status OK or BINNED; typed is the share of calls with an
integer allele.

| Scheme | Genomes | Loci | Mean OK | Lowest OK | Typed |
| --- | --- | --- | --- | --- | --- |
| `bacillus_mlva15_vanert_2007` | 12 | 15 | 85% | 60% | 89% |
| `bacillus_mlva25_lista_2006` | 12 | 25 | 85% | 20% | 94% |
| `bacillus_mlva31_mlvabank` | 12 | 31 | 90% | 58% | 92% |
| `bacillus_mlva8_keim_2000` | 12 | 8 | 86% | 62% | 91% |
| `blongum_mlva5_matamoros_2011` | 12 | 5 | 68% | 20% | 93% |
| `Bru_MLVA16` | 12 | 16 | 86% | 0% | 91% |
| `brucella_mlva15_lefleche_2006` | 12 | 15 | 89% | 0% | 91% |
| `coxiella_mlva17_arricau_bouvery_2006` | 12 | 17 | 93% | 71% | 98% |
| `francisella_johansson_2004` | 12 | 25 | 90% | 84% | 94% |
| `francisella_vogler_2008` | 12 | 11 | 98% | 91% | 98% |
| `Lpn_MLVA12` | 12 | 12 | 92% | 75% | 94% |
| `magalactiae_mlva4_mcauliffe_2008` | 11 | 4 | 57% | 25% | 70% |
| `mtb_miru15_supply_2006` | 12 | 15 | 91% | 80% | 99% |
| `mtb_miru24_supply_2006` | 12 | 24 | 93% | 88% | 98% |
| `nmeningitidis_mlva12_liao_2006` | 12 | 12 | 86% | 75% | 92% |
| `paeruginosa_mlva8_onteniente_2003` | 12 | 8 | 81% | 62% | 97% |
| `saureus_mlva14_pourcel_2009` | 12 | 14 | 92% | 64% | 97% |
| `saureus_mlva5_sabat_2003` | 12 | 5 | 80% | 60% | 82% |
| `saureus_mlva8_schouls_2009` | 12 | 8 | 90% | 75% | 100% |
| `salmonella_enteritidis_mlva9_malorny_2008` | 12 | 9 | 99% | 89% | 100% |
| `STm_ECDC` | 12 | 5 | 83% | 40% | 83% |
| `yersinia_klevytska_2001` | 12 | 42 | 97% | 93% | 98% |
| `yersinia_mlva25_lefleche_2001` | 12 | 25 | 96% | 88% | 99% |
| `yersinia_mlva26_li_2013` | 12 | 26 | 99% | 92% | 100% |

Loci that the summary lists for review (one flag on more than half of the typed
calls, or an allele below zero):

| Scheme | Locus | Reason |
| --- | --- | --- |
| `blongum_mlva5_matamoros_2011` | VNTR23 | NEGATIVE in 92% of calls; allele -2 |
| `Bru_MLVA16` | Bruce19 | INTERMEDIATE in 55% of calls |
| `mtb_miru15_supply_2006` | 580 | INTERMEDIATE in 100% of calls |
| `mtb_miru24_supply_2006` | 580 | INTERMEDIATE in 100% of calls |
| `paeruginosa_mlva8_onteniente_2003` | ms077 | INTERMEDIATE in 70% of calls |
| `saureus_mlva14_pourcel_2009` | Sa0906 | allele -1.5 |
| `saureus_mlva14_pourcel_2009` | Sa1425 | allele -2 |
| `yersinia_klevytska_2001` | M27 | allele -2 |

Of the loci listed here, three turned out not to be errors at all: Ft-M19 is
the published 30 bp *holarctica* deletion (all 12 short products are
*holarctica*, all 7 long ones are *novicida*, *tularensis* or *mediasiatica*),
Magal_VNTR17 carries the paper's own fractional copy number for an array with
72 to 83 per cent internal identity, and VNTR61_01 is the locus the Schouls
paper singles out for a heterogeneous upstream flank, which the authors binned
by hand. The schemes already recorded all three, so they now name them in
`explained_flags` and a screen lists them apart from the loci that need review.

Two remain, and neither can be fixed by changing a copy number. At MIRU 580 all
19 non-reference genomes sit on one grid and H37Rv alone is 24 bp off it, a
partial copy of the 77 bp unit. At Bruce19, over 31 *Brucella* genomes, only
three product sizes occur (163 bp in 16M alone, 178 bp and 184 bp); every size
is congruent modulo 3, and 16M differs from every other genome by an odd
multiple of 3, never by a multiple of the published 6 bp unit. In both the
reference sits off the grid the other strains share, which the `grid_offset`
mechanism cannot express, because a fractional `ref_repeats` assumes the
reference is on the grid and only moves the flag to the other strains. Both
need the published size-to-allele table as `bins.tsv`.

Two of these were resolved after the second pass: M27 and VNTR23 are now
numbered on the units their reference amplicon holds (19 and 7) rather than the
published counts, because the published counts imply a non-repeat flank longer
than the shortest product seen in other strains, which cannot be. Their alleles
are 13 and 4 above the published numbering, which the scheme notes record. The
rest of this section describes the state before that change.

These are candidates, not diagnoses. A negative allele means that the product
is shorter than the reference calibration allows (a reference count that is too
low, a flank deletion in these strains, or a wrong unit), and a flag on every
call means that the reference size or unit does not match the strains; each
needs checking against the publication and the genomes, as was done for the
U'Ren scheme. An allele of 0 is not listed, because it can be a state of the
scheme: at Ft-M20B it means no insertion, and MIRU 1955 has a single unit in
the reference, so 0 is a genome without it (1 of 12 here).

Two low values have other causes. The lowest OK for the two *Brucella* schemes
(0%) is *Brucella anthropi* (strain PBO), which the genus level taxon query
returned: it is not in the species the schemes were built for, and no product
was found at any locus. *M. agalactiae* has the lowest mean call rate (57%,
4 loci, 11 genomes), where one missing product weighs 25 percentage points.

## Checking the loci that were listed

Each locus listed above was examined by aligning the amplicons of the genomes
with different sizes against the reference amplicon and reading the scheme's
own notes. Three have a scheme cause and three are real variation (two loci of one genome count as one).

| Scheme and locus | Finding | Verdict |
| --- | --- | --- |
| `yersinia_klevytska_2001` M27 | The published copy number is 6 but the CO92 amplicon holds 19 units of 8 bp (the scheme notes record 6 against 19). The allele formula therefore places a product of 181 bp, which holds about 11 units, at allele -2, and allele 32 for a product of 453 bp. | Scheme cause: the published numbering cannot hold in other strains. |
| `blongum_mlva5_matamoros_2011` VNTR23 | The published copy number is 3 but NCC2705 holds about 7 copies of the 10 bp unit (CTGATATTCG). Eleven of 12 genomes give 324 to 334 bp, which is 40 to 50 bp shorter than NCC2705, and the size is below the lower end of the range of the paper's table (344 to 574 bp). The formula therefore gives -2 and -1. | Scheme cause, as for M27. The publication could be checked for the range. |
| `mtb_miru15_supply_2006` and `mtb_miru24_supply_2006` 580 | All 12 genomes give products of 98 + 77 n bp (252, 329, 406, 560, 637), whole units from a flank of 98 bp. H37Rv gives 353 bp, 24 bp more than three units, because its last unit is partial. The reference sits between two grid points, so every other genome is a half allele. | Scheme cause: the copy number is derived, not published (the scheme says so). The published allele table, as `bins.tsv`, is the fix. |
| `saureus_mlva14_pourcel_2009` Sa0906, Sa1425 | One genome (NAS_AN_157) gives a 599 bp product at Sa0906 and a 272 bp product at Sa1425, with two further loci flagged. The other 11 give sizes in steps of the unit. | Real variation or an atypical assembly: NEGATIVE is the right flag. |
| `paeruginosa_mlva8_onteniente_2003` ms077 | Seven genomes give 388 bp, 15 bp shorter than the 403 bp product and not a multiple of the 39 bp unit. The alignment shows a 15 bp deletion inside the array. | Real variation of an imperfect array: INTERMEDIATE is correct. |
| `Bru_MLVA16` Bruce19 | Products of 178 and 184 bp are 9 bp and 15 bp longer than the 169 bp product (unit 6 bp). The alignment shows a 15 bp insertion inside the array. | Real variation of an imperfect array. |

The three scheme causes share one pattern: the copy number published for the
reference does not match the array in the reference sequence (as at five loci of
the U'Ren scheme), so other strains get negative or half alleles. Keeping the
published number keeps the published numbering at the reference strain only;
using the physical count removes the negative alleles but gives numbers that do
not match the publication at the reference. This is a decision for the
maintainer and has not been changed here.

## Second pass: 20 genomes per scheme, with alignment and report

Each of the 24 other schemes was run on 20 RefSeq genomes of its species or
genus: 14 complete and 6 draft (scaffold or contig level), chosen at random with
a fixed seed, with `--jobs 4 --align mafft` and the HTML report, on one
laptop. Time is wall clock for the whole run including MAFFT and the report;
memory is the peak of the main process. Call rate is the share of loci with
status OK or BINNED, averaged over the complete and the draft genomes.

| Scheme | Genomes | Loci | Time (s) | Memory (MB) | OK, complete | OK, draft |
| --- | --- | --- | --- | --- | --- | --- |
| `bacillus_mlva15_vanert_2007` | 20 | 15 | 8 | 192 | 89% | 68% |
| `bacillus_mlva25_lista_2006` | 20 | 25 | 14 | 204 | 92% | 65% |
| `bacillus_mlva31_mlvabank` | 20 | 31 | 16 | 203 | 94% | 75% |
| `bacillus_mlva8_keim_2000` | 20 | 8 | 6 | 192 | 89% | 71% |
| `blongum_mlva5_matamoros_2011` | 20 | 5 | 4 | 153 | 67% | 83% |
| `Bru_MLVA16` | 20 | 16 | 8 | 170 | 60% | 89% |
| `brucella_mlva15_lefleche_2006` | 20 | 15 | 7 | 155 | 64% | 94% |
| `coxiella_mlva17_arricau_bouvery_2006` | 20 | 17 | 9 | 180 | 90% | 79% |
| `francisella_johansson_2004` | 20 | 25 | 10 | 195 | 75% | 70% |
| `francisella_vogler_2008` | 20 | 11 | 5 | 163 | 86% | 82% |
| `Lpn_MLVA12` | 20 | 12 | 8 | 168 | 88% | 71% |
| `magalactiae_mlva4_mcauliffe_2008` | 17 | 4 | 2 | 144 | 57% | 38% |
| `mtb_miru15_supply_2006` | 20 | 15 | 13 | 185 | 92% | 66% |
| `mtb_miru24_supply_2006` | 20 | 24 | 18 | 186 | 93% | 73% |
| `nmeningitidis_mlva12_liao_2006` | 20 | 12 | 6 | 161 | 79% | 74% |
| `paeruginosa_mlva8_onteniente_2003` | 20 | 8 | 10 | 249 | 88% | 52% |
| `saureus_mlva14_pourcel_2009` | 20 | 14 | 9 | 159 | 91% | 77% |
| `saureus_mlva5_sabat_2003` | 20 | 5 | 8 | 155 | 77% | 60% |
| `saureus_mlva8_schouls_2009` | 20 | 8 | 5 | 160 | 87% | 77% |
| `salmonella_enteritidis_mlva9_malorny_2008` | 20 | 9 | 6 | 180 | 99% | 89% |
| `STm_ECDC` | 20 | 5 | 5 | 191 | 81% | 67% |
| `yersinia_klevytska_2001` | 20 | 42 | 16 | 210 | 98% | 94% |
| `yersinia_mlva25_lefleche_2001` | 20 | 25 | 10 | 192 | 96% | 96% |
| `yersinia_mlva26_li_2013` | 20 | 26 | 10 | 193 | 100% | 93% |

- **Speed and memory.** 2 to 18 s and 140 to 250 MB for 20 genomes, so the
  time grows with the number of loci and not with genome size. A serial run
  (`--jobs 1`) of four schemes gave files identical to the parallel run,
  alignments included.
- **Draft assemblies lose products.** The call rate of the drafts is 5 to
  25 percentage points below that of the complete genomes for most schemes,
  mostly from NEAR_CONTIG_END (an array that broke the assembly). The two
  *Brucella* schemes are the exception: 5 of the 14 complete genomes belong to
  other species of the genus (*B. anthropi*, *B. intermedia*), which the
  primers do not amplify, so the comparison does not apply. For
  *B. longum* the complete genomes score lower than the drafts because of the
  negative alleles at VNTR23 (see above), not because of the assemblies.
- **A bug found and fixed.** Five of the 24 runs stopped in the report step
  with a numerical error (the schemes with 4 to 9 loci: `blongum_mlva5`,
  `magalactiae_mlva4`, `STm_ECDC`, `salmonella_enteritidis_mlva9`,
  `saureus_mlva5`). The distance between two samples needs 8 shared loci, which
  a panel of fewer than 8 loci cannot give, so every distance was undefined and
  the scaling failed. The threshold is now half of a small panel (at least 2),
  a panel with no shared loci gives equal distances, and a figure that cannot be
  drawn is left out with a warning and no longer stops the run. Tests added.
- **More loci to review with 20 genomes.** The additional ones are
  Ft-M19 (INTERMEDIATE in 63% of calls: 12 genomes give a product 30 bp shorter
  than the others, which is not a multiple of the 13 bp unit, so the locus is
  bimodal like Ft-M20B), Magal_VNTR17 and VNTR61_01 (half alleles because the
  sizes are not on the unit grid of the reference), and Bruce19, where all 15
  genomes are 15 or 21 bp longer than the reference 16M product (163 bp), so
  16M itself is off the grid of the other strains. These belong to the same
  class as the loci checked above: the allele formula is anchored on a reference
  whose product does not sit on the grid of most other strains.
