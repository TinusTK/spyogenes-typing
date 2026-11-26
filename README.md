# spyogenes-typing

I built this tool for public health and epidemiology use, to provide rapid,
reproducible typing of *Streptococcus pyogenes* (group A *Streptococcus*,
GAS) from assembled genomes, against PubMLST's database and following
PubMLST's own typing rules:

- **emm** type and subtype (e.g. `emm1.0`), plus emm-like genes (**emmL**)
  and *S. dysgalactiae* emm genes (**emmCG**)
- **emm cluster** (Sanderson-Smith 2014) and **emm pattern** (A-E, REA)
- **7-gene MLST** sequence type (ST) and the allele of every locus
- **M1UK lineage** for emm1: M1UK, M1intermediate, M1DK or M1global
- **superantigen** genes and **resistance** genes
- **FCT-region** form and **pilin type** (the T-antigen pilus region)
- per-row **warnings** and the **database version** on every row
- one tabular report per run, ready for Excel and downstream reporting

Database: [PubMLST *S. pyogenes*](#database), emm data curated by the CDC.
The bundled copy holds **only data submitted up to 31 December 2024**, which
is the part PubMLST gives without an account, refreshed through a reviewed
pull request. For current data, add your own PubMLST key with `--setup`.

> **Prefer no command line?**
> [Genomics Workflow Microbials (GWM)](https://github.com/TinusTK/GenomicsWorkflowMicrobials)
> is a free desktop app for Windows, Mac and Linux with a point-and-click
> interface. *S. pyogenes* typing is built in: choose a folder of reads or
> assemblies, press Run and read the results as tables. Nothing to install
> by hand, and no server needed: it all runs locally on your computer.

## How settled each part is

PubMLST marks several of its *S. pyogenes* schemes as still developing, and
this tool reports them anyway, so it is worth knowing which results rest on
a settled scheme and which on a provisional one. These flags are read from
PubMLST itself, not assigned here, and the database updater reports it when
they change:

| What it reports | PubMLST scheme | PubMLST flags |
|---|---|---|
| ST, alleles, MLST status | MLST (1) | **established** (no flags) |
| emm, emmL, emmCG, emm pattern markers | Emm typing (4) | *in development* |
| Resistance genes | Resistance (5) | *in development* |
| FCT type (regulator plus pilA / pilB) | FCT-pilA (2), FCT-pilB (3), FCT-other (6) | *experimental*, *in development* |
| pilL alleles (in the `--hits` file only; they do not change the FCT call) | FCT-pilL (7) | *experimental*, *in development* |

Four more columns come from outside PubMLST's schemes: **emm cluster** and
the predicted **emm pattern** are curated by hand from PubMLST's typing
guide, **superantigens** from NCBI RefSeq genomes, and the **M1UK** SNPs
from Lynskey et al. 2019.

PubMLST also warns that its partial typing loci "may not be 100%
inclusive", so a locus that does not match is a possible new variant rather
than an error. The tool says so in the Warnings column instead of guessing.

## Output

One row per genome, one table per run, so several genomes share one header:

```bash
spyogenes-typing MGAS5005.fasta MGAS6180.fasta Alab49.fasta > typing.tsv
spyogenes-typing *.fasta.gz --outfile typing.tsv
```

| Sample | emm | emm status | emmL | emmCG | emm cluster | emm pattern | ST | MLST status | M1 lineage | FCT type | Warnings |
|---|---|---|---|---|---|---|---|---|---|---|---|
| MGAS5005 | emm1.0 | PERFECT | - | - | A-C3 | A | 28 | PERFECT | M1global (0 M1UK SNPs) | FCT-2 (pil001) | - |
| MGAS6180 | emm28.0 | PERFECT | ~emmL203.4 (97.2%) | - | E4 | E | 52 | PERFECT | - | FCT-4 (pil028) | - |
| Alab49 | emm53.0 | PERFECT | emmL134.2;~emmL156.1 (98.9%) | - | D4 | D | 11 | PERFECT | - | FCT-3 (pil033) | - |

Those are the headline columns; the full report has 22. `emmL` and `emmCG`
are always there, `-` when the genome has none: an emm-like gene is listed
as emmL and never reported as the emm type, which is what the emmPrime rule
below settles.

**Warnings is the column to read first.** Every other column can put
something in it, so one glance tells you whether any part of the row needs
checking; `-`, as in all three rows above, means nothing did. The text names
the column it is about, and several are joined with `; `. Real examples:

| Genome | Warnings |
|---|---|
| MGAS8232 | `emm pattern markers mga-2;SF1 fit no pattern A-E: REA; predicted for this emm type: A-C/D` |
| MGAS10394 | `superantigen not called, nearest known allele under 90% identity (a new variant, or a disrupted gene): speI 87.4% id / 59% of the gene, speH 84.3% id / 98% of the gene` |
| a genome with a locus missing | `MLST MISSING: gki(-) not found` |
| a genome with most MLST genes missing | `only 3 of 7 S. pyogenes MLST loci found at >= 95% identity` |
| a new emm subtype (emm column `~emm63.2`) | `emm NOVEL: 1 SNP from emm63.2 (99.4%); cluster and pattern from emm63` |
| a genome cut through emm | `emm PARTIAL: emm cut by a contig end (54% of the allele), subtype unknown` |

Here is one whole row, laid out downwards (MGAS6180):

| Column | Value |
|---|---|
| Sample | `MGAS6180` |
| emm | `emm28.0` |
| emm status | `PERFECT` |
| emmL | `~emmL203.4 (97.2%)` |
| emmCG | `-` |
| emm cluster | `E4` |
| emm pattern | `E` |
| ST | `52` |
| ALLELE1..7 | `gki(11)` `gtr(6)` `murI(14)` `mutS(5)` `recP(9)` `xpt(17)` `yqiL(19)` |
| MLST status | `PERFECT` |
| M1 lineage | `-` |
| Superantigens | `smeZ;speC;speG;speJ;speK` |
| Resistance genes | `-` |
| FCT type | `FCT-4 (pil028)` |
| Warnings | `-` |
| Database | `PubMLST_2026-11-20` |

Empty values are `-`. MGAS6180 and Alab49 are the worked examples in
PubMLST's own typing guide, and the tool reproduces both. What every column
means is in [The report](#the-report) below.

## Install

Needs Python 3.9+ and BLAST+ (`blastn`, `makeblastdb`).

```bash
mamba env create -f environment.yml
conda activate spyogenes-typing
pip install .
```

## Input

Assembled genomes (contigs or complete), FASTA, GenBank or EMBL, plain or
`.gz` / `.bz2` / `.zip`. Several files give several rows under one header.
`--fofn FILE` reads the file names from a file (one per line, an optional
TAB and a sample name).

## The report

`Sample` is the file name without its extension (or `--label` / the
name in `--fofn`). Empty values are `-`.

| Column | What it is | Meaning |
|---|---|---|
| emm | The emm allele: type and subtype together, e.g. `emm1.0` (type 1, subtype 1.0). `~emm63.2` when there is no exact match: the closest allele, with differences (emm status `NOVEL`, details in Warnings). Taken from the 180 bp 5' region PubMLST holds, curated by the CDC. | The main label for GAS in surveillance and outbreak work: isolates with the same emm type may be linked, and emm types are tracked over time and place. |
| emm status | `PERFECT` exact allele; `NOVEL` no exact allele, the nearest one (≥95% identity) is shown with a `~`, e.g. `~emm63.2` (a new subtype, see Warnings); `PARTIAL` emm cut by a contig end, only the type is shown (e.g. `emm1`); `MULTICOPY` the same allele at more than one emm locus; `MIXED` different emm alleles at more than one emm locus (mixed culture or contamination?); `MISSING` no emm found. | How sure the emm call is. `PERFECT` can be reported as is; `NOVEL` and `PARTIAL` need a look; `MIXED` suggests the sample is not a single strain. |
| emmL | emm-like genes (e.g. *mrp*, *enn*) found away from the emm locus, as `emmL134.2`; when not exact, the closest allele with a `~` and its identity, e.g. `~emmL203.4 (97.2%)`. | Part of the strain's gene layout, not its type. Shown so an emm-like gene is never taken for the emm type. |
| emmCG | *S. dysgalactiae* (stG / stC) emm genes found away from the emm locus, as `STG653.0`; `~STG653.0 (98.1%)` when not exact. | Often a sign of genes taken up from *S. dysgalactiae*; worth noting, and together with the warning that most MLST genes were not found it may mean the genome is not *S. pyogenes*. |
| emm cluster | The emm cluster of the emm type (Sanderson-Smith et al. 2014), as PubMLST writes it: E1-E6 (clade X), D1-D5 and A-C1 to A-C5 (clade Y), singletons `clade_X_M164`, `clade_Y_M6`, outliers `outlier_M55`. Types in one cluster share M-protein binding properties: used for vaccine design (cross-protection) and to group types in epidemiology. | Groups emm types whose M proteins behave alike, so it says which other types an isolate may resemble for immunity and vaccine coverage. |
| emm pattern | How the emm and emm-like genes are laid out (Bessen; PubMLST). `A`, `B`, `C` (throat specialists), `D` (skin specialists), `E` (generalists), `REA` (rearranged). Measured from marker hits next to emm (below). When the markers cannot all be seen (a short contig), the pattern known for the emm type is shown as `A-C (predicted)`. | Hints at where the strain tends to cause infection: throat (A-C), skin (D) or both (E). |
| ST | MLST sequence type. | The strain's genetic background, from the 7 MLST genes. Isolates with the same emm type but different STs are probably not closely related. |
| ALLELE1..7 | MLST alleles, e.g. `gki(4)`. `~4` nearest allele with differences, `4?` cut by a contig end, `-` missing, `4,5` two alleles. `--combine-alleles` puts them in one ALLELES column (`gki(4);gtr(3);...`); `--no-alleles` leaves them out. | Which MLST gene caused a missing or new ST, when there is one. |
| MLST status | `PERFECT` all exact with an ST; `NOVEL` all exact, new combination; `PARTIAL` some alleles not exact; `MULTICOPY`; `MIXED`; `MISSING` a locus not found; `BAD` many alleles not exact; `NONE` nothing found. Same words as mlst-typing. | How sure the ST is. `PERFECT` can be reported as is; `NOVEL` may be a new ST to submit to PubMLST; the others need a look. |
| M1 lineage | emm1 only. Counts the 27 M1UK SNPs (Lynskey et al. 2019) and the M1DK SNPs: `M1UK (27/27 SNPs)`, `M1intermediate (13/27 M1UK SNPs)`, `M1DK`, `M1global (0 M1UK SNPs)`. | For emm1, whether the isolate belongs to M1UK, the lineage that makes more SpeA toxin and has spread widely since it emerged in England in the 2010s. |
| Superantigens | Superantigen genes found (≥90% identity): speA, speC, speG, speH, speI, speJ, speK, speL, speM, smeZ, ssa. `(partial)` when cut by a contig end. | The toxin genes the strain carries. Superantigens are linked to scarlet fever and streptococcal toxic shock; speA in particular to invasive disease. |
| Resistance genes | PubMLST resistance alleles, e.g. `ermB(1)`; `~1` with differences, `1?` partial. | Predicted antibiotic resistance, e.g. erm genes for macrolides and clindamycin, tet genes for tetracycline. A prediction from the genome, not a lab susceptibility result. |
| FCT type | FCT-region form and pilin type, e.g. `FCT-4 (pil028)`. The form comes from the regulator (rofA or nra) and the pilA / pilB loci found; the pilin type from the pilA and pilB nt80 clusters of exact alleles. | The pilus the strain makes to stick to host tissue. The pilin is the T antigen of the older serological T-typing, so the pilin type can be set roughly against T-typing data. |
| Warnings | Anything to check, from any column above, separated by `; `. Covers emm (NOVEL, PARTIAL, MIXED, an emmL at the emm locus, no emmPrime), emm cluster and pattern (none listed for the type, markers that fit no pattern, measured against predicted), MLST (which locus is missing, mixed, or not exact), M1 lineage (SNP sites not found), superantigens and resistance genes (present but under the 90% calling limit), FCT (pilin not exact, no pilin type listed, hits that fit no form) and the MLST genes: when 5 or fewer of the 7 *S. pyogenes* MLST genes are found at 95% identity or more, it says `only n of 7 S. pyogenes MLST loci found at >= 95% identity` (another species such as *S. dysgalactiae*, or a poor assembly; this counts MLST genes only and is not a species identification: a close relative can still pass). `-` means nothing needs checking. | The column to read first: `-` means the row can be used as it is; anything else names the column to check and why. |
| Database | The database version, e.g. `PubMLST_2026-11-20`. | Which data the calls came from, so results can be compared and repeated later. |

### How the calls are made

All rules come from PubMLST's guide *Rules-based assignment of emm and
pilin genotypes for group A streptococci* (Feb 2024).

- **Which hit is emm?** A genome often has several emm-like genes. The
  short `emmPrime` marker sits about 50 nt upstream of emm; the hit just
  downstream of it is emm, the rest are emm-like. If an emm-like (emmL) or
  emmCG allele sits at the emm locus, that is the emm type (with a
  warning). Without emmPrime only an allele of the `emm` locus is called
  emm, never an emm-like one.
- **emm pattern** from the marker hits within 10 kb of emm:
  A = mga-1 + SF1; B = mga-1 + 2 × SF1; C = mga-1 + SF1 + SF3;
  D = mga-2 + SF1 + SF3 + SF4; E = mga-2 + SF2 + SF3 + SF4. Other
  combinations are REA, with a warning (a new marker variant is possible:
  PubMLST's marker lists are not complete).
- **FCT-region form**: FCT-1 rofA + pilA1/pilB1; FCT-2 rofA + pilA2/pilB2;
  FCT-3 nra + pilA34/pilB34; FCT-4 rofA + pilA34/pilB34; FCT-3/4
  pilA34/pilB34 alone; FCT-5 rofA + pilA5/pilB5; FCT-6 rofA + pilA6/pilB6;
  FCT-9 rofA + pilB9; FCT-10 rofA + pilA10/pilB10. PubMLST files the
  guide's pilA10 alleles under pilA5, so the report names them as
  PubMLST does.
- A pilin allele with no nt80 cluster in the guide gives the form without a
  pilin type, with a warning.

## Other output

- `--json FILE`: the same results as JSON.
- `--hits FILE`: every place something hit, with the closest allele,
  identity, coverage, contig and position.

## Database

The bundled `database/` holds data from
[PubMLST *S. pyogenes*](https://pubmlst.org/organisms/streptococcus-pyogenes)
submitted up to 31 December 2024 (the part PubMLST gives without an account),
the curated tables in `database/tables/` and the references in
`database/references/`.

For newer data use your own PubMLST account and key:

```bash
spyogenes-typing --setup          # asks for the key, downloads your own database
spyogenes-typing --update         # refresh it
spyogenes-typing --update-check   # is PubMLST newer? (JSON, downloads nothing)
spyogenes-typing --setup --remove
```

Your key is kept in `~/.config/spyogenes-typing/credentials.json`
(readable only by you) and your database in
`~/.local/share/spyogenes-typing/database`. The database used is, in
order: `--db`, `$SPYOGENES_TYPING_DB`, your own, the bundled one.

Rebuild the bundled database: `python scripts/update_database.py --db database`
(`--references` also rebuilds the NCBI parts).

## Validation

Every result is checked against PubMLST's own isolate records, and the whole
record is in [validation/](validation) so you can see it rather than take my
word for it: [VALIDATION.md](validation/VALIDATION.md) holds the comparison
and the tool's full output for every genome, [expected.tsv](validation/expected.tsv)
holds what PubMLST records for each strain and which isolate it came from,
and [genomes/](validation/genomes) holds the genomes themselves with their
checksums.

**22 genomes** against database `PubMLST_2026-11-20`: 15 complete genomes
from NCBI (including PubMLST's own worked examples MGAS6180 and Alab49) and
7 draft assemblies from PubMLST. On every field PubMLST records for these
strains, the tool reproduces PubMLST's own assignment, and each of the four
differences is traced to the source data rather than the typing:

| What | Agree | Compared | What that field tests |
|---|---|---|---|
| emm type | 21 | 21 | which hit is emm (the emmPrime rule) and the allele |
| emm subtype | 18 | 19 | the same, to the subtype |
| ST | 18 | 20 | the MLST calls and the profile lookup |
| FCT-region form | **15** | **15** | the regulator plus pilA / pilB rule |
| emm cluster | 16 | 16 | the curated table, keyed on the emm type |
| emm pattern | 14 | 15 | measured here against PubMLST's *predicted* |

**Read that as 102 of 106 field values, not 102 genomes**, and note that the
last two rows are table lookups keyed on the emm type, so once the emm type
is right they follow from the curated tables and cannot fail on their own.
The 75 emm type, subtype, ST and FCT values are the independent ones.
PubMLST records no pilin type, superantigen, resistance or M1UK value for
these isolates, so those columns are not checked here at all: the pilin type
is checked against the worked examples in the typing guide (MGAS6180 =
pil028, isolate 6424 = pil088) and M1UK by the test suite, which writes the
published SNPs into a genome and expects M1UK back.

**All nine FCT forms are exercised, but not all of them against PubMLST.**
FCT-1 to FCT-6 are the 15 that agree above. FCT-9 and FCT-10 are called here
too, on isolates PubMLST records no FCT value for, so they count as coverage
rather than agreement. FCT-3/4 is not a form a genome carries at all: it is
the answer when pilA34 and pilB34 are found but neither regulator is, so the
test suite covers it by masking rofA out of an FCT-4 genome.

**Form coverage is recorded separately**, because agreement and coverage are
different things. PubMLST has curated an FCT-region form for only 134
isolates, none with contigs, so for some forms there is nothing to compare
against. All 500 complete *S. pyogenes* genomes in NCBI were typed in one run
instead ([sweep_500_genomes.tsv](validation/sweep_500_genomes.tsv)): that
found **45 distinct pilin types** and seven of the nine forms, including 86
FCT-9 genomes, three of which are kept here as coverage cases. No FCT-10 came
up in those 500, and FCT-3/4 could not. None of that is counted in the 106.

[VALIDATION.md](validation/VALIDATION.md) sets out all of this, including
what each difference turned out to be and what the figure does and does not
show. Re-run it with `python scripts/validate.py`, which exits non-zero if
any difference is not one of the explained ones.

## Licence

MIT (see `LICENSE`). The bundled PubMLST data and the curated tables keep
the terms of their own sources, listed below.

## Sources and credits

- [PubMLST *S. pyogenes*](https://pubmlst.org/organisms/streptococcus-pyogenes)
  typing databases ([Jolley et al. 2018](https://doi.org/10.12688/wellcomeopenres.14826.1),
  Wellcome Open Res 3:124), emm data curated by the
  [CDC *Streptococcus* laboratory](https://www.cdc.gov/strep-lab/php/group-a-strep/index.html).
  Please cite PubMLST.
- The typing rules: PubMLST's
  [Rules-based assignment of emm and pilin genotypes for group A streptococci](https://dev.pubmlst.org/static/organisms/streptococcus-pyogenes/Guide_emm_pilin_genotypes_Feb2024.pdf)
  (Feb 2024), which the curated tables in `database/tables/` are typed from.
- emm clusters: [Sanderson-Smith et al. 2014](https://doi.org/10.1093/infdis/jiu260),
  J Infect Dis 210:1325-38.
- emm patterns and pilin types: [Bessen et al. 2024](https://doi.org/10.1128/mbio.00693-24),
  mBio 15:e00693-24.
- M1UK SNPs: [Lynskey et al. 2019](https://doi.org/10.1016/S1473-3099%2819%2930446-3),
  Lancet Infect Dis 19:1209-18; positions as in
  [github.com/boasvdp/assembly_snptyper](https://github.com/boasvdp/assembly_snptyper)
  (MIT licence).
- Superantigen gene sequences, from the annotated CDS of these NCBI RefSeq
  genomes: [MGAS5005](https://www.ncbi.nlm.nih.gov/nuccore/NC_007297.2),
  [SF370](https://www.ncbi.nlm.nih.gov/nuccore/NC_002737.2),
  [MGAS315](https://www.ncbi.nlm.nih.gov/nuccore/NC_004070.1),
  [MGAS8232](https://www.ncbi.nlm.nih.gov/nuccore/NC_003485.1),
  [MGAS10394](https://www.ncbi.nlm.nih.gov/nuccore/NC_006086.1),
  [MGAS10270](https://www.ncbi.nlm.nih.gov/nuccore/NC_008022.1),
  [MGAS10750](https://www.ncbi.nlm.nih.gov/nuccore/NC_008024.1),
  [MGAS15252](https://www.ncbi.nlm.nih.gov/nuccore/NC_017040.1) and
  [MGAS6180](https://www.ncbi.nlm.nih.gov/nuccore/NC_007296.2), plus
  [L29565.1](https://www.ncbi.nlm.nih.gov/nuccore/L29565.1) for *ssa*.
