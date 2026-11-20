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
- a species check, per-row **warnings** and the **database version** on
  every row
- one tabular report per run, ready for Excel and downstream reporting

Database: [PubMLST *S. pyogenes*](https://pubmlst.org/organisms/streptococcus-pyogenes),
emm data curated by the CDC. Refreshed through a reviewed pull request.

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
| ST, alleles, MLST status, species check | MLST (1) | **established** (no flags) |
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

| Sample | emm | emm status | emmL | emmCG | emm cluster | emm pattern | ST | MLST status | M1 lineage | FCT type |
|---|---|---|---|---|---|---|---|---|---|---|
| MGAS5005 | emm1.0 | PERFECT | - | - | A-C3 | A | 28 | PERFECT | M1global (0 M1UK SNPs) | FCT-2 (pil001) |
| MGAS6180 | emm28.0 | PERFECT | emmL203.4 (97.2%) | - | E4 | E | 52 | PERFECT | - | FCT-4 (pil028) |
| Alab49 | emm53.0 | PERFECT | emmL134.2;emmL156.1 (98.9%) | - | D4 | D | 11 | PERFECT | - | FCT-3 (pil033) |

Those are the headline columns; the full report has 23. `emmL` and `emmCG`
are always there, `-` when the genome has none: an emm-like gene is listed
as emmL and never reported as the emm type, which is what the emmPrime rule
below settles. Here is one whole row, laid out downwards (MGAS6180):

| Column | Value |
|---|---|
| Sample | `MGAS6180` |
| Species check | `S. pyogenes` |
| emm | `emm28.0` |
| emm status | `PERFECT` |
| emmL | `emmL203.4 (97.2%)` |
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

| Column | What it is |
|---|---|
| Species check | `S. pyogenes` when at least 6 of the 7 MLST loci are found at 95% identity or more. Otherwise `not S. pyogenes? (n/7 MLST loci)`: the genome may be another species (e.g. *S. dysgalactiae*) or a poor assembly. |
| emm | The emm allele: type and subtype together, e.g. `emm1.0` (type 1, subtype 1.0). Taken from the 180 bp 5' region PubMLST holds, curated by the CDC. |
| emm status | `PERFECT` exact allele; `NOVEL` no exact allele, the nearest one (≥95% identity) is shown (a new subtype, see Warnings); `PARTIAL` emm cut by a contig end, only the type is shown (e.g. `emm1`); `MULTICOPY` the same allele at more than one emm locus; `MIXED` different emm alleles at more than one emm locus (mixed culture or contamination?); `MISSING` no emm found. |
| emmL | emm-like genes (e.g. *mrp*, *enn*) found away from the emm locus, as `emmL134.2`; the identity is added when not exact. |
| emmCG | *S. dysgalactiae* (stG / stC) emm genes found away from the emm locus. |
| emm cluster | The emm cluster of the emm type (Sanderson-Smith et al. 2014), as PubMLST writes it: E1-E6 (clade X), D1-D5 and A-C1 to A-C5 (clade Y), singletons `clade_X_M164`, `clade_Y_M6`, outliers `outlier_M55`. Types in one cluster share M-protein binding properties: used for vaccine design (cross-protection) and to group types in epidemiology. |
| emm pattern | How the emm and emm-like genes are laid out (Bessen; PubMLST). `A`, `B`, `C` (throat specialists), `D` (skin specialists), `E` (generalists), `REA` (rearranged). Measured from marker hits next to emm (below). When the markers cannot all be seen (a short contig), the pattern known for the emm type is shown as `A-C (predicted)`. |
| ST | MLST sequence type. |
| ALLELE1..7 | MLST alleles, e.g. `gki(4)`. `~4` nearest allele with differences, `4?` cut by a contig end, `-` missing, `4,5` two alleles. `--combine-alleles` puts them in one ALLELES column (`gki(4);gtr(3);...`); `--no-alleles` leaves them out. |
| MLST status | `PERFECT` all exact with an ST; `NOVEL` all exact, new combination; `PARTIAL` some alleles not exact; `MULTICOPY`; `MIXED`; `MISSING` a locus not found; `BAD` many alleles not exact; `NONE` nothing found. Same words as mlst-typing. |
| M1 lineage | emm1 only. Counts the 27 M1UK SNPs (Lynskey et al. 2019) and the M1DK SNPs: `M1UK (27/27 SNPs)`, `M1intermediate (13/27 M1UK SNPs)`, `M1DK`, `M1global (0 M1UK SNPs)`. |
| Superantigens | Superantigen genes found (≥90% identity): speA, speC, speG, speH, speI, speJ, speK, speL, speM, smeZ, ssa. `(partial)` when cut by a contig end. |
| Resistance genes | PubMLST resistance alleles, e.g. `ermB(1)`; `~1` with differences, `1?` partial. |
| FCT type | FCT-region form and pilin type, e.g. `FCT-4 (pil028)`. The form comes from the regulator (rofA or nra) and the pilA / pilB loci found; the pilin type from the pilA and pilB nt80 clusters of exact alleles. |
| Warnings | Anything to check, separated by `; `. |
| Database | The database version, e.g. `PubMLST_2026-11-20`. |

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

The bundled `database/` holds PubMLST data submitted up to 31 December
2024 (the part PubMLST gives without an account), the curated tables in
`database/tables/` and the references in `database/references/`.

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
7 draft assemblies from PubMLST. **102 of 106 recorded values agree, with
nothing unexplained:**

| What | Agree | Compared |
|---|---|---|
| emm type | 21 | 21 |
| emm subtype | 18 | 19 |
| ST | 18 | 20 |
| emm cluster | 16 | 16 |
| emm pattern | 14 | 15 |
| FCT-region form | **15** | **15** |

The four differences were each traced to the source data, not to the typing:
two are PubMLST isolate records that disagree with the deposited genome by a
single base (HSC5's ST comes down to one base in *gtr*; MGAS6180's subtype to
one base in *emm*), one is a measured pattern where PubMLST only records a
predicted one (MGAS8232 genuinely has no SF3 or SF4 beside its emm), and one
is an ST assigned after the free data cut-off. Each is set out in full in the
record.

PubMLST's public isolate records hold no pilin type at all, so that column
cannot be checked against them; it is checked against the worked examples in
the typing guide instead (MGAS6180 = pil028, isolate 6424 = pil088).

Re-run it yourself with `python scripts/validate.py`, which exits non-zero if
any difference is not one of the explained ones.

## Licence

MIT (see `LICENSE`). The bundled PubMLST data and the curated tables keep
the terms of their own sources, listed below.

## Sources and credits

- PubMLST *S. pyogenes* typing databases (Jolley et al. 2018, Wellcome Open
  Res 3:124), emm data curated by the CDC. Please cite PubMLST.
- emm clusters: Sanderson-Smith et al. 2014, J Infect Dis 210:1325-38.
- emm patterns and pilin types: Bessen et al. 2024, mBio.
- M1UK SNPs: Lynskey et al. 2019, Lancet Infect Dis 19:1209-18; positions
  as in github.com/boasvdp/assembly_snptyper (MIT licence).
- Superantigen gene sequences: NCBI RefSeq genomes MGAS5005, SF370,
  MGAS315, MGAS8232, MGAS10394, MGAS10270 and L29565.1 (ssa).
