# spyogenes-typing

Typing of *Streptococcus pyogenes* (group A *Streptococcus*, GAS) from
assembled genomes, using PubMLST's data and PubMLST's own typing rules:

- **emm** type and subtype, emm-like genes (emmL) and *S. dysgalactiae*
  emm genes (emmCG)
- **emm cluster** (Sanderson-Smith 2014) and **emm pattern** (A-E)
- **MLST** (7 loci, ST)
- **M1 lineage** for emm1: M1UK, M1intermediate, M1DK or M1global
- **superantigen** genes and **resistance** genes
- **FCT-region** form and pilin type (T-antigen pilus region)

One row per genome, in one table. Runs on the command line, in Galaxy and
in GWM.

```bash
spyogenes-typing genome.fasta > typing.tsv
spyogenes-typing *.fasta.gz --outfile typing.tsv
```

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
