# Validation, database PubMLST_2026-11-20

22 genomes: 15 complete genomes from NCBI and 7 draft assemblies from PubMLST.

Every genome is typed with the bundled database and compared with the PubMLST isolate record for the same strain. PubMLST leaves many of these fields empty, and an empty field is counted neither way.

| What | Agree | Compared |
|---|---|---|
| emm type | 21 | 21 |
| emm subtype | 18 | 19 |
| ST | 18 | 20 |
| emm cluster | 16 | 16 |
| emm pattern | 14 | 15 |
| FCT-region form | 15 | 15 |
| pilin type | – | 0 (PubMLST records none) |

PubMLST's public isolate records hold no value at all for pilin type, so that cannot be checked against them. It is checked instead against the worked examples in PubMLST's typing guide (MGAS6180 = pil028, Alab49, and isolate 6424 = pil088), which the tool reproduces.

**102 of 106 recorded values agree.** 4 differ for a reason traced to the source data (below) and **0 differ unexplained.**

### Differences traced to the source data

Each was investigated; none is a fault in the typing.

- **MGAS6180, emm subtype**: the tool says `emm28.0`, PubMLST records `emm28.4`. The assembly matches emm_28.0 exactly; emm_28.4, which the isolate record gives, differs from it by one base. The subtype of the deposited sequence is emm28.0.
- **MGAS8232, emm pattern**: the tool says `REA`, PubMLST records `A-C`. Measured, not predicted: the emm region of this genome holds mga-2 and SF1 only, with no SF3 or SF4 anywhere near it (checked down to 70% identity), which fits no pattern A-E, so the guide's rules give REA. A-C is the pattern predicted for emm18, which PubMLST also records as being seen with pattern D.
- **HSC5, ST**: the tool says `84`, PubMLST records `680`. PubMLST isolate 2548 lists gtr 2 (entered by hand in 2015); the deposited assembly NC_021807.1 matches gtr_22 exactly (100% identity and coverage). gtr_2 and gtr_22 differ by one base, and ST84 and ST680 differ only at gtr, so ST84 is the ST of this assembly.
- **i7422, ST**: the tool says `-`, PubMLST records `1778`. ST1778 was assigned after 31 December 2024, so it is not in the bundled free database; NOVEL is the right answer against it. Users with their own PubMLST key get this ST.

## How to read that number

The unit is a single field of a single genome, not a genome: 22 genomes with up to 7 comparable fields each is 154 possible values, and PubMLST records no value in 48 of them, which leaves 106. So it is not "102 genomes correct", and the denominator is decided by what PubMLST happens to have curated, not by what matters most.

**The fields are not independent of each other.** emm cluster and the predicted emm pattern are table lookups keyed on the emm type (`db.clusters[emm_type]` and `db.patterns[emm_type]`), so once the emm type is right they follow from the curated tables and cannot fail on their own. Those 31 values test that the tables were transcribed from the guide correctly, which is worth knowing, but they are not 31 further checks of the typing. The 75 emm type, emm subtype, ST and FCT-region values are the independent ones.

**emm pattern is not compared like with like.** PubMLST's public records hold no measured pattern for any of these genomes, so the tool's measured pattern is compared with the pattern PubMLST predicts for the emm type. They answer different questions, and MGAS8232 below is exactly that case.

**Both sides come from PubMLST**, though from different parts of it: the tool types against PubMLST's sequence definitions, and the expected values are PubMLST's curated isolate records. Agreement therefore shows the tool reproduces PubMLST's own assignments, including which hit is emm and which FCT form the loci make. It is not independent evidence that those assignments are biologically right.

**Not checked here at all:** the superantigen, resistance-gene, M1UK lineage and warning columns, because PubMLST's isolate records hold nothing to compare them with. M1UK is covered by the test suite instead, which writes the published SNPs into a genome and expects M1UK back.

**All nine FCT forms are exercised, but not all against PubMLST.** FCT-1 to FCT-6 are the ones that agree above. FCT-9 and FCT-10 are called on isolates PubMLST records no FCT form for, so they count as coverage rather than agreement, and FCT-9 has three further genomes of its own below. FCT-3/4 is not a form a genome carries: it is the answer when pilA34 and pilB34 are found but neither regulator is, which the test suite covers by masking rofA out of an FCT-4 genome.

**The sample leans on reference strains:** 15 finished genomes and only 7 draft assemblies, so fragmented real-world assemblies are under-represented. The test suite covers a contig break through emm separately.

## Compared with PubMLST

`-` in an expected column means PubMLST has no value for it.

| Genome | Source | emm type | emm subtype | ST | emm cluster | emm pattern | FCT form | pilin type |
|---|---|---|---|---|---|---|---|---|
| `MGAS5005` | NCBI NC_007297.2 | 1 ✓ | emm1.0 ✓ | 28 ✓ | A-C3 ✓ | A-C ✓ | FCT-2 ✓ | – |
| `MGAS6180` | NCBI NC_007296.2 | 28 ✓ | emm28.0 (PubMLST emm28.4, explained) | 52 ✓ | E4 ✓ | E ✓ | FCT-4 ✓ | – |
| `Alab49` | NCBI CP003068.1 | 53 ✓ | emm53.0 ✓ | 11 ✓ | D4 ✓ | D ✓ | FCT-3 ✓ | – |
| `SF370` | NCBI NC_002737.2 | 1 ✓ | emm1.6 ✓ | 28 ✓ | A-C3 ✓ | A-C ✓ | FCT-2 ✓ | – |
| `MGAS10394` | NCBI NC_006086.1 | 6 ✓ | emm6.4 ✓ | 382 ✓ | clade_Y_M6 ✓ | A-C ✓ | FCT-1 ✓ | – |
| `Manfredo` | NCBI NC_009332.1 | 5 ✓ | emm5.0 ✓ | 99 ✓ | clade_Y_M5 ✓ | A-C ✓ | FCT-3 ✓ | – |
| `SSI-1` | NCBI NC_004606.1 | 3 ✓ | emm3.1 ✓ | 15 ✓ | A-C5 ✓ | A-C ✓ | FCT-3 ✓ | – |
| `MGAS8232` | NCBI NC_003485.1 | 18 ✓ | emm18.19 ✓ | 42 ✓ | clade_Y_M18 ✓ | REA (PubMLST A-C, explained) | FCT-3 ✓ | – |
| `MGAS315` | NCBI NC_004070.1 | 3 ✓ | emm3.1 ✓ | 15 ✓ | A-C5 ✓ | A-C ✓ | FCT-3 ✓ | – |
| `MGAS9429` | NCBI NC_008021.1 | 12 ✓ | emm12.0 ✓ | 36 ✓ | A-C4 ✓ | A-C ✓ | FCT-4 ✓ | – |
| `MGAS15252` | NCBI NC_017040.1 | 59 ✓ | emm59.0 ✓ | 172 ✓ | E6 ✓ | D ✓ | FCT-4 ✓ | – |
| `MGAS1882` | NCBI NC_017053.1 | 59 ✓ | emm59.0 ✓ | 172 ✓ | E6 ✓ | D ✓ | FCT-4 ✓ | – |
| `HSC5` | NCBI NC_021807.1 | 14 ✓ | emm14.3 ✓ | 84 (PubMLST 680, explained) | clade_Y_M14 ✓ | A-C ✓ | FCT-4 ✓ | – |
| `MGAS10750` | NCBI NC_008024.1 | 4 ✓ | emm4.0 ✓ | – | E1 ✓ | E ✓ | FCT-5 ✓ | – |
| `MGAS10270` | NCBI NC_008022.1 | 2 ✓ | emm2.0 ✓ | 55 ✓ | E4 ✓ | E ✓ | FCT-6 ✓ | – |
| `i4078` | PubMLST isolate 4078 | 22 ✓ | – | 1235 ✓ | – | – | – | – |
| `i4115` | PubMLST isolate 4115 | 44 ✓ | emm44.2 ✓ | 1279 ✓ | – | – | – | – |
| `i6424` | PubMLST isolate 6424 | – | – | 117 ✓ | – | – | – | – |
| `i6537` | PubMLST isolate 6537 | 11 ✓ | emm11.0 ✓ | – | – | – | – | – |
| `i6541` | PubMLST isolate 6541 | 53 ✓ | emm53.0 ✓ | 1477 ✓ | – | – | – | – |
| `i6549` | PubMLST isolate 6549 | 18 ✓ | emm18.21 ✓ | 1482 ✓ | – | – | – | – |
| `i7422` | PubMLST isolate 7422 | 75 ✓ | – | - (PubMLST 1778, explained) | E6 ✓ | – | – | – |

## Full output of the tool

Every column of `results_PubMLST_2026-11-20.tsv`, as the tool wrote it.

| Sample | emm | emm status | emmL | emmCG | emm cluster | emm pattern | ST | ALLELE1 | ALLELE2 | ALLELE3 | ALLELE4 | ALLELE5 | ALLELE6 | ALLELE7 | MLST status | M1 lineage | Superantigens | Resistance genes | FCT type | Warnings | Database |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| MGAS5005 | emm1.0 | PERFECT | - | - | A-C3 | A | 28 | gki(4) | gtr(3) | murI(4) | mutS(4) | recP(4) | xpt(2) | yqiL(4) | PERFECT | M1global (0 M1UK SNPs) | smeZ;speA;speG;speJ | - | FCT-2 (pil001) | - | PubMLST_2026-11-20 |
| MGAS6180 | emm28.0 | PERFECT | ~emmL203.4 (97.2%) | - | E4 | E | 52 | gki(11) | gtr(6) | murI(14) | mutS(5) | recP(9) | xpt(17) | yqiL(19) | PERFECT | - | smeZ;speC;speG;speJ;speK | - | FCT-4 (pil028) | - | PubMLST_2026-11-20 |
| Alab49 | emm53.0 | PERFECT | emmL134.2;~emmL156.1 (98.9%) | - | D4 | D | 11 | gki(2) | gtr(6) | murI(1) | mutS(2) | recP(2) | xpt(2) | yqiL(2) | PERFECT | - | smeZ;speC;speG;speL;speM | - | FCT-3 (pil033) | - | PubMLST_2026-11-20 |
| SF370 | emm1.6 | PERFECT | - | - | A-C3 | A | 28 | gki(4) | gtr(3) | murI(4) | mutS(4) | recP(4) | xpt(2) | yqiL(4) | PERFECT | M1global (0 M1UK SNPs) | smeZ;speC;speG;speH;speI;speJ | - | FCT-2 (pil001) | M1UK: 1 of 27 SNP sites not found | PubMLST_2026-11-20 |
| MGAS10394 | emm6.4 | PERFECT | - | - | clade_Y_M6 | A | 382 | gki(5) | gtr(52) | murI(5) | mutS(5) | recP(5) | xpt(4) | yqiL(3) | PERFECT | - | smeZ;speA;speC;speG;speK | msrD(2);mefA(1) | FCT-1 (pil006) | superantigen not called, nearest known allele under 90% identity (a new variant, or a disrupted gene): speI 87.4% id / 59% of the gene, speH 84.3% id / 98% of the gene | PubMLST_2026-11-20 |
| Manfredo | emm5.0 | PERFECT | - | - | clade_Y_M5 | A | 99 | gki(33) | gtr(30) | murI(7) | mutS(5) | recP(5) | xpt(26) | yqiL(3) | PERFECT | - | smeZ;speC;speG;speH;speI | - | FCT-3 (pil005) | - | PubMLST_2026-11-20 |
| SSI-1 | emm3.1 | PERFECT | - | - | A-C5 | A | 15 | gki(2) | gtr(6) | murI(8) | mutS(5) | recP(2) | xpt(3) | yqiL(2) | PERFECT | - | smeZ;speA;speG;speK;ssa | - | FCT-3 (pil003) | - | PubMLST_2026-11-20 |
| MGAS8232 | emm18.19 | PERFECT | - | - | clade_Y_M18 | REA | 42 | gki(6) | gtr(5) | murI(7) | mutS(1) | recP(5) | xpt(6) | yqiL(3) | PERFECT | - | smeZ;speA;speC;speG;speL;speM | - | FCT-3 (pil018) | emm pattern markers mga-2;SF1 fit no pattern A-E: REA; predicted for this emm type: A-C/D | PubMLST_2026-11-20 |
| MGAS315 | emm3.1 | PERFECT | - | - | A-C5 | A | 15 | gki(2) | gtr(6) | murI(8) | mutS(5) | recP(2) | xpt(3) | yqiL(2) | PERFECT | - | smeZ;speA;speG;speK;ssa | - | FCT-3 (pil003) | - | PubMLST_2026-11-20 |
| MGAS9429 | emm12.0 | PERFECT | - | - | A-C4 | A | 36 | gki(5) | gtr(2) | murI(2) | mutS(6) | recP(6) | xpt(2) | yqiL(2) | PERFECT | - | smeZ;speC;speG;speH;speI | - | FCT-4 (pil012) | - | PubMLST_2026-11-20 |
| MGAS15252 | emm59.0 | PERFECT | emmL174.1 | - | E6 | D | 172 | gki(56) | gtr(24) | murI(39) | mutS(7) | recP(30) | xpt(2) | yqiL(33) | PERFECT | - | smeZ;speG;speJ | - | FCT-4 | pilA nt80_13 + pilB nt80_28: no pilin type listed | PubMLST_2026-11-20 |
| MGAS1882 | emm59.0 | PERFECT | emmL174.1 | - | E6 | D | 172 | gki(56) | gtr(24) | murI(39) | mutS(7) | recP(30) | xpt(2) | yqiL(33) | PERFECT | - | smeZ;speG;speJ;speK | - | FCT-4 | pilA nt80_13 + pilB nt80_28: no pilin type listed | PubMLST_2026-11-20 |
| HSC5 | emm14.3 | PERFECT | emm51.0 | - | clade_Y_M14 | B | 84 | gki(20) | gtr(22) | murI(2) | mutS(3) | recP(27) | xpt(18) | yqiL(12) | PERFECT | - | smeZ;speG;speL;speM | - | FCT-4 (pil014) | - | PubMLST_2026-11-20 |
| MGAS10750 | emm4.0 | PERFECT | emmL236.3;emmL156.0 | - | E1 | E | 1072 | gki(5) | gtr(11) | murI(8) | mutS(5) | recP(152) | xpt(2) | yqiL(1) | PERFECT | - | smeZ;speC;ssa | ermA(1) | FCT-5 (pil004) | - | PubMLST_2026-11-20 |
| MGAS10270 | emm2.0 | PERFECT | emmL149.2 | - | E4 | E | 55 | gki(11) | gtr(9) | murI(1) | mutS(9) | recP(2) | xpt(3) | yqiL(4) | PERFECT | - | speC;speG;speK | - | FCT-6 (pil002) | - | PubMLST_2026-11-20 |
| i4078 | emm22.0 | PERFECT | ~emmL149.3 (99.4%);emmL156.0 | - | E4 | E | 1235 | gki(188) | gtr(8) | murI(1) | mutS(1) | recP(1) | xpt(3) | yqiL(4) | PERFECT | - | smeZ;speG;speL;speM;ssa | tetM(~35) | FCT-4 (pil066) | - | PubMLST_2026-11-20 |
| i4115 | emm44.2 | PERFECT | ~emmL203.4 (97.2%) | - | E3 | E | 1279 | gki(4) | gtr(2) | murI(3) | mutS(6) | recP(14) | xpt(3) | yqiL(171) | PERFECT | - | speG;speH;speJ | mefA(8);msrD(7);catA(6);tetM(26) | FCT-4 | pilA34(~7) not exact: no nt80 cluster, pilin type not called | PubMLST_2026-11-20 |
| i6424 | emm81.0 | PERFECT | ~emmL156.1 (98.3%);emmL164.3 | - | E6 | D | 117 | gki(38) | gtr(2) | murI(2) | mutS(23) | recP(39) | xpt(2) | yqiL(25) | PERFECT | - | speG;speH;speJ | - | FCT-10 (pil088) | - | PubMLST_2026-11-20 |
| i6537 | emm11.0 | PERFECT | emmL202.1 | - | E6 | E | - | gki(3) | gtr(4) | murI(6) | mutS(7) | recP(1) | xpt(5) | yqiL(-) | MISSING | - | smeZ;speC;speG | - | FCT-4 (pil011) | MLST MISSING: yqiL(-) not found | PubMLST_2026-11-20 |
| i6541 | emm53.0 | PERFECT | emmL134.2;~emmL156.1 (98.9%) | - | D4 | D | 1477 | gki(2) | gtr(6) | murI(1) | mutS(2) | recP(186) | xpt(2) | yqiL(2) | PERFECT | - | smeZ;speC;speG;speL;speM | - | FCT-3 (pil033) | - | PubMLST_2026-11-20 |
| i6549 | emm18.21 | PERFECT | ~emmL156.1 (98.9%);~emmL156.1 (98.3%);~emmL167.0 (98.9%) | - | clade_Y_M18 | D | 1482 | gki(1) | gtr(2) | murI(2) | mutS(130) | recP(35) | xpt(3) | yqiL(109) | PERFECT | - | smeZ;speG;speJ | - | FCT-3 (pil053) | emm pattern markers not exact (new variant?): mga-2 | PubMLST_2026-11-20 |
| i7422 | emm75.0 | PERFECT | emmL170.0 | - | E6 | E | - | gki(3) | gtr(2) | murI(1) | mutS(37) | recP(12) | xpt(8) | yqiL(7) | NOVEL | - | smeZ;speC;speG;speJ | - | FCT-9 (pil025) | MLST NOVEL: every allele is an exact match but the combination 3/2/1/37/12/8/7 is not a known ST | PubMLST_2026-11-20 |
## Form coverage (not agreement)

PubMLST records an FCT-region form for only 134 of its isolates, **none of which has contigs**, and it holds a single FCT-9 isolate with no assembly. So for some forms there is no curated value anywhere to compare against. The genomes below, found by typing 500 complete NCBI genomes, instead show that the rule engine reaches the form and that the loci it saw match that form's definition. **These are not counted in the agreement figure above**, because there is nothing to agree with.

| Genome | emm | Form called | Pilin type | Why it is evidence |
|---|---|---|---|---|
| `NZ_CP066541.1` | emm92.0 | FCT-9 (pil019) ✓ | pil019 | rofA and pilB9 exact, no pilA hit at all even at 70% identity |
| `NZ_CP035455.1` | emm197.0 | FCT-9 (pil065) ✓ | pil065 | rofA and pilB9 exact, no pilA hit at all even at 70% identity |
| `NZ_CP033335.1` | emm75.0 | FCT-9 (pil025) ✓ | pil025 | rofA and pilB9 exact; the only pilA-like hits cover 8% of the gene, so pilA is absent |

### What 500 complete genomes look like

Every complete *S. pyogenes* genome in NCBI that the search returned (500 of 776) was typed in one run, as a check that the tool does not fall over and to see the spread of forms. It found **45 distinct pilin types**. Full table: [sweep_500_genomes.tsv](sweep_500_genomes.tsv).

| Form | Genomes |
|---|---|
| FCT-4 | 177 |
| FCT-2 | 134 |
| FCT-9 | 86 |
| FCT-3 | 55 |
| FCT-5 | 39 |
| FCT-1 | 7 |
| FCT-6 | 2 |

FCT-9 is 17% of these genomes but only 1 of PubMLST's 134 curated FCT isolates. That curated set is old and uneven (it still uses the retired FCT-7 and FCT-8 labels and has no FCT-3/4 or FCT-10), so the difference is in what has been curated, not in the calls: every FCT-9 genome checked has rofA and pilB9 as exact matches and no pilA at all. FCT-10 did not appear in this set, and FCT-3/4 cannot: it is what the tool answers when the regulator is missing, which the test suite covers by masking rofA.


## Genomes

Listed in [expected.tsv](expected.tsv) with the PubMLST isolate each expected value comes from; the files and their checksums are in [genomes/MANIFEST.tsv](genomes/MANIFEST.tsv). Rebuild the set with `python scripts/fetch_validation_genomes.py` and re-run this record with `python scripts/validate.py`.

