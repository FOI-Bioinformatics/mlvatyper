# Bundled schemes and their scientific references

Every bundled scheme reproduces a published MLVA or VNTR typing assay. This page names the
publication each one comes from, the reference strain its values are anchored on, and how its
reference values were obtained, because the last of these decides how far an allele can be
compared with the literature. `mlvatyper schemes` lists the schemes and `mlvatyper schemes NAME`
prints one scheme's loci, citation and notes; the `notes` field of each `scheme.yml` and the comments
at the end of its `markers.tsv` are the full account, and this page is the summary.

The DOIs were checked against the Crossref registry on 2026-11-03, and the bibliography at the end is
taken from those records. A test keeps this page in step with the bundled schemes: every scheme must
appear in it, and so must every DOI that a `scheme.yml` cites.

## How the reference values were obtained

An allele is `ref_repeats + (length - ref_size) / unit_length`, so it is only as comparable with a
publication as its reference size and count. The table gives one of these bases for each scheme:

- **Published**: the reference size, unit and repeat count are values in the publication, or in the
  database identifiers built from it, and a run on the reference strain reproduces them.
- **Measured**: the publication does not give the reference product or count, so they were measured in
  the reference assembly.
- **Derived**: counted from the reference assembly or from a size range, because the source table gives
  no counts. Confirm these against the published allele tables before comparing profiles.
- **Calibrated**: taken from another bundled scheme that counts the same array, through a genome both
  schemes can be run on.
- **Convention**: the publication gives no count, so the scheme defines one and says so. Alleles from
  such a locus compare only with each other.

## Schemes

| scheme | organism | loci | reference strain (assembly) | source | basis of the values |
|---|---|---|---|---|---|
| `bacillus_mlva15_vanert_2007` | *Bacillus anthracis* | 15 | *B. anthracis* Ames Ancestor (GCA_000008445.1) | [Van Ert et al. 2007](https://doi.org/10.1371/journal.pone.0000461) | Published values read from the MLVAbank locus identifiers, confirmed in the assembly. |
| `bacillus_mlva25_lista_2006` | *Bacillus anthracis* | 25 | *B. anthracis* Ames Ancestor (GCA_000008445.1) | [Lista et al. 2006](https://doi.org/10.1186/1471-2180-6-33) | Published size-to-repeat table (Table 3) as bins. The paper gives no primers, so those are MLVAbank's. Numbering differs from MLVAbank at five loci. |
| `bacillus_mlva31_mlvabank` | *Bacillus anthracis* | 31 | *B. anthracis* Ames Ancestor (GCA_000008445.1) | [Keim et al. 2000](https://doi.org/10.1128/JB.182.10.2928-2936.2000); [Lista et al. 2006](https://doi.org/10.1186/1471-2180-6-33); [Van Ert et al. 2007](https://doi.org/10.1371/journal.pone.0000461) | Values read from the MLVAbank identifiers; all 31 loci reproduced in the assembly. |
| `bacillus_mlva8_keim_2000` | *Bacillus anthracis* | 8 | *B. anthracis* Ames Ancestor (GCA_000008445.1) | [Keim et al. 2000](https://doi.org/10.1128/JB.182.10.2928-2936.2000) | Published values read from the MLVAbank locus identifiers, confirmed in the assembly. |
| `blongum_mlva5_matamoros_2011` | *Bifidobacterium longum subsp. longum* | 5 | *B. longum* NCC2705 (GCA_000007525.1) | [Matamoros et al. 2011](https://doi.org/10.1016/j.mimet.2011.10.005) | Published: units and motifs (Table 3) and the NCC2705 profile (Table 1), which a run reproduces at four loci. VNTR23 is numbered on the 7 units its amplicon holds, because the published 3 implies a flank longer than products seen in other strains; its alleles are 4 above the publication. |
| `Bru_MLVA16` | *Brucella* | 16 | *B. melitensis* 16M (GCF_000007125.1) | [Le Flèche et al. 2006](https://doi.org/10.1186/1471-2180-6-9); [Al Dahouk et al. 2007](https://doi.org/10.1016/j.mimet.2006.12.015) | Published: Table 1 sizes and 16M copy numbers, all 16 products match. Panels after Al Dahouk et al. No unit sequences, so no tandem cross-check. |
| `brucella_mlva15_lefleche_2006` | *Brucella melitensis* | 15 | *B. melitensis* 16M (GCA_000007125.1) | [Le Flèche et al. 2006](https://doi.org/10.1186/1471-2180-6-9) | Published, in the published locus order. Not interchangeable with MLVA-16, which adds Bruce19. |
| `burkholderia_uren_2007` | *Burkholderia pseudomallei* | 32 | *B. pseudomallei* K96243 (GCA_000011545.1) | [U'Ren et al. 2007](https://doi.org/10.1186/1471-2180-7-23) | Measured in K96243: the paper gives array sizes and ranges but not the K96243 product. Copy numbers are measured, except five (2115k, 1500k, 3564k, 1788k, 1217k) taken from the paper because the measured count gave alleles of 0 or below in other strains. |
| `coxiella_mlva17_arricau_bouvery_2006` | *Coxiella burnetii* | 17 | *C. burnetii* Nine Mile RSA 493 (GCA_000007765.2) | [Arricau-Bouvery et al. 2006](https://doi.org/10.1186/1471-2180-6-38) | Published sizes; 14 of 17 reproduced, and ms30, ms31 and ms36 carry measured sizes. |
| `francisella_johansson_2004` | *Francisella tularensis* | 25 | *F. tularensis* SCHU S4 (GCA_000008985.1) | [Johansson et al. 2004](https://doi.org/10.1128/JB.186.17.5808-5818.2004) | Published (Table 2): motifs, units and SCHU S4 counts for all 25 loci. Ft-M14 is duplicated and Ft-M8 is displaced in some lineages. |
| `francisella_vogler_2008` | *Francisella tularensis* | 11 | *F. tularensis* LVS (GCA_000009245.1) | [Vogler et al. 2009](https://doi.org/10.1111/j.1472-765X.2008.02484.x); [Johansson et al. 2004](https://doi.org/10.1128/JB.186.17.5808-5818.2004) | Ten loci take unit and count from Johansson by calibration through SCHU S4. Ft-M20B is a convention of this scheme, allele 0 without the insertion and (size - 306) / 15 with it. |
| `Lpn_MLVA12` | *Legionella pneumophila* | 12 | *L. pneumophila* Philadelphia-1 (GCF_000008485.1) | [Sobral et al. 2011](https://doi.org/10.1128/AEM.05556-11) | Published (Table 1); all 12 products match. Lpms31 gives fractional alleles in other genomes. |
| `mtb_miru15_supply_2006` | *Mycobacterium tuberculosis* | 15 | *M. tuberculosis* H37Rv (GCA_000195955.2) | [Supply et al. 2006](https://doi.org/10.1128/JCM.01392-06) | Derived from H37Rv, because the source table gives no copy numbers. Confirm against the published allele tables before comparing with MIRU codes. |
| `mtb_miru24_supply_2006` | *Mycobacterium tuberculosis* | 24 | *M. tuberculosis* H37Rv (GCA_000195955.2) | [Supply et al. 2006](https://doi.org/10.1128/JCM.01392-06) | Derived from H37Rv, because the source table gives no copy numbers. Confirm against the published allele tables before comparing with MIRU codes. |
| `magalactiae_mlva4_mcauliffe_2008` | *Mycoplasma agalactiae* | 4 | *M. agalactiae* PG2 (GCA_000063605.1) | [McAuliffe et al. 2008](https://doi.org/10.1186/1471-2180-8-193) | Published (Table 1): units and PG2 copy numbers, three of them fractional. Magal_VNTR19's reported length and span disagree, so its measured product is used. |
| `nmeningitidis_mlva12_liao_2006` | *Neisseria meningitidis* | 12 | *N. meningitidis* Z2491 (GCA_000009105.1) | [Liao et al. 2006](https://doi.org/10.1186/1471-2180-6-44) | Published (Tables 1 and 5), consistent across three genomes. NMTR12 is called from bins that cover the three sequenced genomes only. |
| `paeruginosa_mlva8_onteniente_2003` | *Pseudomonas aeruginosa* | 8 | *P. aeruginosa* PAO1 (GCA_000006765.1) | [Onteniente et al. 2003](https://doi.org/10.1128/JCM.41.11.4991-4997.2003) | Published (Table 1): units and PAO1 counts; the profile 11-12-5-8-7-12-14-45 is reproduced. |
| `salmonella_enteritidis_mlva9_malorny_2008` | *Salmonella enterica subsp. enterica serovar Enteritidis* | 9 | *S.* Enteritidis NCTC 13349 (P125109) (GCA_000009505.1) | [Malorny et al. 2008](https://doi.org/10.1186/1471-2180-8-84) | Published repeat numbers (Table 1); four sequence lengths measured. The paper's MLVA type is a string of arbitrary indices and is not reproduced. |
| `STm_ECDC` | *Salmonella enterica subsp. enterica serovar Typhimurium* | 5 | *S.* Typhimurium LT2 (GCF_000006945.2) | [Larsson et al. 2009](https://doi.org/10.2807/ese.14.15.19174-en) | Sizes are the LT2 products of the ECDC protocol; repeat numbers derived from its size ranges and confirmed by tandem counts. STTR3 from the published size-to-allele table (bins). |
| `saureus_mlva14_pourcel_2009` | *Staphylococcus aureus* | 14 | *S. aureus* Mu50 (GCA_000009665.1) | [Pourcel et al. 2009](https://doi.org/10.1128/JCM.00267-09); [Bogut et al. 2024](https://doi.org/10.2478/cipms-2024-0014) | Published (Table 1), with Mu50 sizes from Bogut et al.; checked on 31 clinical profiles. |
| `saureus_mlva5_sabat_2003` | *Staphylococcus aureus* | 5 | *S. aureus* NCTC 8325 (GCA_000013425.1) | [Sabat et al. 2003](https://doi.org/10.1128/JCM.41.4.1801-1804.2003) | Convention of this scheme: whole repeat units between the primer sites, since the paper is a gel fingerprint with no counts. |
| `saureus_mlva8_schouls_2009` | *Staphylococcus aureus* | 8 | *S. aureus* N315 (GCA_000009645.1) | [Schouls et al. 2009](https://doi.org/10.1371/journal.pone.0005082) | Published: the N315 profile of Figure 7, validated on seven genomes. VNTR61_01 has a heterogeneous flank. |
| `yersinia_klevytska_2001` | *Yersinia pestis* | 42 | *Y. pestis* CO92 (GCA_000009065.1) | [Klevytska et al. 2001](https://doi.org/10.1128/JCM.39.9.3179-3185.2001) | Published CO92 array column; three loci carry half-copy counts. M27 is numbered on the 19 units its amplicon holds, because the published 6 implies a flank longer than products seen in other strains; its alleles are 13 above the publication. A 42 locus survey, not a standard panel. |
| `yersinia_mlva25_lefleche_2001` | *Yersinia pestis* | 25 | *Y. pestis* CO92 (GCA_000009065.1) | [Le Flèche et al. 2001](https://doi.org/10.1186/1471-2180-1-2); [Pourcel et al. 2004](https://doi.org/10.1186/1471-2180-4-22) | Published unit, count and product for all 25 loci. |
| `yersinia_mlva26_li_2013` | *Yersinia pestis* | 26 | *Y. pestis* CO92 (GCA_000009065.1) | [Li et al. 2013](https://doi.org/10.1371/journal.pone.0066567) | Published (Table S2); five encoded sizes were measured and corrected, each corroborated by another bundled scheme. |
| `testM` | *Francisella tularensis* | 31 | a *holarctica* strain (none) | [Johansson et al. 2004](https://doi.org/10.1128/JB.186.17.5808-5818.2004) | No reference counts; anchored at run time with --reference-sample. Derived from the Johansson panel: 22 of its 31 markers share a Johansson primer pair and unit. |

## Related schemes

Several organisms have more than one bundled scheme, and the schemes within a group are not
interchangeable even where they read the same loci.

- ***Bacillus anthracis***. The four panels share Ames Ancestor and nest by locus: the MLVA-8 loci are in
  the 15, 25 and 31 locus panels, and MLVA-25 is the first 25 loci of the MLVAbank MLVA-31. The MLVA-8,
  MLVA-15 and MLVA-31 counts come from the MLVAbank identifiers. Lista et al. number vrrA, vrrB1, vrrB2,
  pXO2-at and bams34 differently, so a profile from `bacillus_mlva25_lista_2006` must not be mixed with
  the others.
- ***Brucella***. `Bru_MLVA16` and `brucella_mlva15_lefleche_2006` share their loci apart from Bruce19,
  which only the 16 locus panel includes, so their profiles are not interchangeable.
- ***Francisella tularensis***. `francisella_johansson_2004` and `francisella_vogler_2008` amplify the
  same VNTRs with different primers and names (Ft-M2 and Ft-M02), and Vogler splits Ft-M20 into Ft-M20A
  and Ft-M20B. Run on 22 genomes, the two give the same allele in all 145 comparisons at the nine loci
  other than Ft-M20. `testM` is built on the Johansson panel.
- ***Mycobacterium tuberculosis***. `mtb_miru15_supply_2006` is the discriminatory subset of the 24 locus
  set and gives the same values at the same loci.
- ***Salmonella***. `STm_ECDC` is for serovar Typhimurium and the Malorny scheme for Enteritidis. They
  share no locus.
- ***Staphylococcus aureus***. Three schemes read the *spa* repeat array with different primers, as
  VNTR24_01 (Schouls), Sa0122 (Pourcel) and `spa` (Sabat). The Pourcel array Sa1132 is VNTR63_01 in the
  Schouls scheme. Each scheme numbers an array on its own scale, so a count from one is not a count from
  another.
- ***Yersinia pestis***. Three overlapping panels. Six loci are shared by `yersinia_klevytska_2001` and
  `yersinia_mlva25_lefleche_2001` and numbered differently, and `yersinia_mlva26_li_2013` shares M loci
  with the Klevytska survey while differing from it at M15, M22, M23, M25 and M58.

## Where a value is not a published value

These are the places where a count rests on a convention, a measurement or a derivation and not on a
number printed in the paper. Each scheme's notes give the detail.

- `saureus_mlva5_sabat_2003`, all five loci. The paper is a gel fingerprint and publishes no counts, and
  its primers lie inside the repeat arrays, so the allele is the number of whole repeat units between the
  primer sites. It is not Sabat's pattern number.
- `francisella_vogler_2008`, Ft-M20B. No publication gives its count; the allele is 0 without the
  insertion and `(size - 306) / 15` with it, from a bins table.
- `nmeningitidis_mlva12_liao_2006`, NMTR12. The paper determined every allele by sequencing, so the bins
  cover the three sequenced genomes only, and any other size is OFF_TABLE.
- `mtb_miru15_supply_2006` and `mtb_miru24_supply_2006`, all loci, and the repeat numbers of
  `STm_ECDC`. Derived from the reference.
- `burkholderia_uren_2007`, all loci. Measured in K96243; copy numbers of five loci from the paper (see the scheme notes). Seventeen loci have arrays the repeat patterns cannot follow across strains and are not compared with the tandem count.
- `salmonella_enteritidis_mlva9_malorny_2008`. The repeat numbers are published, but the paper's own MLVA
  type is a string of arbitrary allele indices whose table is not printed, so it cannot be reproduced.

## Bibliography

1. Al Dahouk S, Le Flèche P, Nöckler K, Jacques I, Grayon M, Scholz HC, et al. Evaluation of Brucella MLVA typing for human brucellosis. Journal of Microbiological Methods 2007;69(1):137-145. doi:10.1016/j.mimet.2006.12.015
2. Arricau-Bouvery N, Hauck Y, Bejaoui A, Frangoulidis D, Bodier CC, Souriau A, et al. Molecular characterization of Coxiella burnetii isolates by infrequent restriction site-PCR and MLVA typing. BMC Microbiol 2006;6:38. doi:10.1186/1471-2180-6-38
3. Bogut A, Grzegorczyk A, Kuszneruk D. The Multiple-Locus Variable Number Tandem Repeat Analysis (MLVA) of Staphylococuss aureus clinical isolates recovered from the Provincial Specialist Hospital in Lublin, Poland. Current Issues in Pharmacy and Medical Sciences 2024;37(2):85-91. doi:10.2478/cipms-2024-0014
4. Johansson A, Farlow J, Larsson P, Dukerich M, Chambers E, Byström M, et al. Worldwide Genetic Relationships among Francisella tularensis Isolates Determined by Multiple-Locus Variable-Number Tandem Repeat Analysis. J Bacteriol 2004;186(17):5808-5818. doi:10.1128/JB.186.17.5808-5818.2004
5. Keim P, Price LB, Klevytska AM, Smith KL, Schupp JM, Okinaka R, et al. Multiple-Locus Variable-Number Tandem Repeat Analysis Reveals Genetic Relationships within Bacillus anthracis. J Bacteriol 2000;182(10):2928-2936. doi:10.1128/JB.182.10.2928-2936.2000
6. Klevytska AM, Price LB, Schupp JM, Worsham PL, Wong J, Keim P. Identification and Characterization of Variable-Number Tandem Repeats in the Yersinia pestis Genome. J Clin Microbiol 2001;39(9):3179-3185. doi:10.1128/JCM.39.9.3179-3185.2001
7. Larsson JT, Torpdahl M, Petersen RF, Sørensen G, Lindstedt BA, Nielsen EM. Development of a new nomenclature for Salmonella Typhimurium multilocus variable number of tandem repeats analysis (MLVA). Eurosurveillance 2009;14(15):19174. doi:10.2807/ese.14.15.19174-en
8. Le Flèche P, Hauck Y, Onteniente L, Prieur A, Denoeud F, Ramisse V, et al. A tandem repeats database for bacterial genomes: application to the genotyping of Yersinia pestis and Bacillus anthracis. BMC Microbiol 2001;1:2. doi:10.1186/1471-2180-1-2
9. Le Flèche P, Jacques I, Grayon M, Al Dahouk S, Bouchon P, Denoeud F, et al. Evaluation and selection of tandem repeat loci for a Brucella MLVA typing assay. BMC Microbiol 2006;6:9. doi:10.1186/1471-2180-6-9
10. Li Y, Cui Y, Cui B, Yan Y, Yang X, Wang H, et al. Features of Variable Number of Tandem Repeats in Yersinia pestis and the Development of a Hierarchical Genotyping Scheme. PLoS ONE 2013;8:e66567. doi:10.1371/journal.pone.0066567
11. Liao JC, Li CC, Chiou CS. Use of a multilocus variable-number tandem repeat analysis method for molecular subtyping and phylogenetic analysis of Neisseria meningitidis isolates. BMC Microbiol 2006;6:44. doi:10.1186/1471-2180-6-44
12. Lista F, Faggioni G, Valjevac S, Ciammaruconi A, Vaissaire J, le Doujet C, et al. Genotyping of Bacillus anthracis strains based on automated capillary 25-loci Multiple Locus Variable-Number Tandem Repeats Analysis. BMC Microbiol 2006;6:33. doi:10.1186/1471-2180-6-33
13. Malorny B, Junker E, Helmuth R. Multi-locus variable-number tandem repeat analysis for outbreak studies of Salmonella enterica serotype Enteritidis. BMC Microbiol 2008;8:84. doi:10.1186/1471-2180-8-84
14. Matamoros S, Savard P, Roy D. Genotyping of Bifidobacterium longum subsp. longum strains by multilocus variable number of tandem repeat analysis. Journal of Microbiological Methods 2011;87(3):378-380. doi:10.1016/j.mimet.2011.10.005
15. McAuliffe L, Churchward CP, Lawes JR, Loria G, Ayling RD, Nicholas RAJ. VNTR analysis reveals unexpected genetic diversity within Mycoplasma agalactiae, the main causative agent of contagious agalactia. BMC Microbiol 2008;8:193. doi:10.1186/1471-2180-8-193
16. Onteniente L, Brisse S, Tassios PT, Vergnaud G. Evaluation of the Polymorphisms Associated with Tandem Repeats for Pseudomonas aeruginosa Strain Typing. J Clin Microbiol 2003;41(11):4991-4997. doi:10.1128/JCM.41.11.4991-4997.2003
17. Pourcel C, André-Mazeaud F, Neubauer H, Ramisse F, Vergnaud G. Tandem repeats analysis for the high resolution phylogenetic analysis of Yersinia pestis. BMC Microbiol 2004;4:22. doi:10.1186/1471-2180-4-22
18. Pourcel C, Hormigos K, Onteniente L, Sakwinska O, Deurenberg RH, Vergnaud G. Improved Multiple-Locus Variable-Number Tandem-Repeat Assay for Staphylococcus aureus Genotyping, Providing a Highly Informative Technique Together with Strong Phylogenetic Value. J Clin Microbiol 2009;47(10):3121-3128. doi:10.1128/JCM.00267-09
19. Sabat A, Krzyszton-Russjan J, Strzalka W, Filipek R, Kosowska K, Hryniewicz W, et al. New Method for Typing Staphylococcus aureus Strains: Multiple-Locus Variable-Number Tandem Repeat Analysis of Polymorphism and Genetic Relationships of Clinical Isolates. J Clin Microbiol 2003;41(4):1801-1804. doi:10.1128/JCM.41.4.1801-1804.2003
20. Schouls LM, Spalburg EC, van Luit M, Huijsdens XW, Pluister GN, van Santen-Verheuvel MG, et al. Multiple-Locus Variable Number Tandem Repeat Analysis of Staphylococcus Aureus: Comparison with Pulsed-Field Gel Electrophoresis and spa-Typing. PLoS ONE 2009;4:e5082. doi:10.1371/journal.pone.0005082
21. Sobral D, Le Cann P, Gerard A, Jarraud S, Lebeau B, Loisy-Hamon F, et al. High-Throughput Typing Method To Identify a Non-Outbreak-Involved Legionella pneumophila Strain Colonizing the Entire Water Supply System in the Town of Rennes, France. Appl Environ Microbiol 2011;77(19):6899-6907. doi:10.1128/AEM.05556-11
22. Supply P, Allix C, Lesjean S, Cardoso-Oelemann M, Rüsch-Gerdes S, Willery E, et al. Proposal for Standardization of Optimized Mycobacterial Interspersed Repetitive Unit-Variable-Number Tandem Repeat Typing of Mycobacterium tuberculosis. J Clin Microbiol 2006;44(12):4498-4510. doi:10.1128/JCM.01392-06
23. U'Ren JM, Schupp JM, Pearson T, Hornstra H, Friedman CLC, Smith KL, et al. Tandem repeat regions within the Burkholderia pseudomallei genome and their application for high resolution genotyping. BMC Microbiol 2007;7:23. doi:10.1186/1471-2180-7-23
24. Van Ert MN, Easterday WR, Huynh LY, Okinaka RT, Hugh-Jones ME, Ravel J, et al. Global Genetic Population Structure of Bacillus anthracis. PLoS ONE 2007;2:e461. doi:10.1371/journal.pone.0000461
25. Vogler AJ, Birdsell D, Wagner DM, Keim P. An optimized, multiplexed multi-locus variable-number tandem repeat analysis system for genotyping Francisella tularensis. Letters in Applied Microbiology 2009;48(1):140-144. doi:10.1111/j.1472-765X.2008.02484.x
