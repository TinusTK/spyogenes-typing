#!/usr/bin/env python3
"""
spyogenes-typing - Streptococcus pyogenes (group A Streptococcus) typing from
PubMLST: emm type and subtype, emm-like genes, emm cluster and pattern, MLST,
M1UK lineage, superantigen genes, resistance genes and FCT-region form.

    spyogenes-typing genome.fasta [more.fasta ...] > typing.tsv

One row per genome; several genomes share one header. See the README for
what every column means and how to read it.

The emm, emm pattern, emm cluster, pilin and FCT rules are PubMLST's own,
from its guide 'Rules-based assignment of emm and pilin genotypes for group A
streptococci' (Feb 2024). Needs only Python 3.8+ and BLAST+ (blastn).
"""

import argparse
import bz2
import csv
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections import defaultdict

VERSION = "1.0.0.dev0"
EXE = "spyogenes-typing"
URL = "https://github.com/TinusTK/spyogenes-typing"
ENV_DB = "SPYOGENES_TYPING_DB"
SEP = "__"
BLAST_LONG = "spyogenes"
BLAST_SHORT = "spyogenes_short"
BLAST_FMT = "6 stitle slen length nident qseqid qstart qend sstart send qlen sstrand"

MLST_SCHEME = "mlst"
EMM_SCHEME = "emm"
SAG_SCHEME = "superantigens"
RES_SCHEME = "resistance"
EMM_FAMILY = ("emm", "emmL", "emmCG")
PATTERN_LOCI = ("mga-1", "mga-2", "SF1", "SF2", "SF3", "SF4")
PRIME_GAP = 150        # emmPrime sits ~50 nt upstream of the emm hit
EMM_WINDOW = 10000     # mga / SF markers counted within this distance of emm
MARKER_ID = 90         # short markers: >90% identity counts (guide p4: new variant)
GENE_ID = 90           # superantigen, resistance and pilin genes: minimum identity
NEAR_ID = 75           # ... and the floor for saying "something is there" instead
GENE_COV = 0.6         # ... and minimum share of the gene covered
EDGE = 3               # a hit this close to a contig end is cut by the edge
SPECIES_LOCI = 6       # MLST loci at >= --minid needed to pass the species check

# emm pattern from the marker hits next to emm (guide p4): marker -> count;
# every other marker must be absent
PATTERNS = (
    ("A", {"mga-1": 1, "SF1": 1}),
    ("B", {"mga-1": 1, "SF1": 2}),
    ("C", {"mga-1": 1, "SF1": 1, "SF3": 1}),
    ("D", {"mga-2": 1, "SF1": 1, "SF3": 1, "SF4": 1}),
    ("E", {"mga-2": 1, "SF2": 1, "SF3": 1, "SF4": 1}),
)

PERFECT, NOVEL, PARTIAL = "PERFECT", "NOVEL", "PARTIAL"
MULTICOPY, MIXED, MISSING = "MULTICOPY", "MIXED", "MISSING"


class TypingError(Exception):
    pass


class Log:
    quiet = False
    debug = False

    @classmethod
    def msg(cls, text):
        if not cls.quiet:
            print(text, file=sys.stderr)

    @classmethod
    def dbg(cls, text):
        if cls.debug:
            print(f"DEBUG: {text}", file=sys.stderr)


def err(text):
    raise TypingError(text)


# --------------------------------------------------------------------------
# Database
# --------------------------------------------------------------------------

def config_dir():
    """Where the setup answer and the user's PubMLST key are kept."""
    base = os.environ.get("XDG_CONFIG_HOME") or os.path.join(os.path.expanduser("~"), ".config")
    return os.path.join(base, EXE)


def user_database():
    """The user's own database, built by `--setup` with their PubMLST key."""
    base = (os.environ.get("XDG_DATA_HOME")
            or os.path.join(os.path.expanduser("~"), ".local", "share"))
    return os.path.join(base, EXE, "database")


def find_database(arg_db, own=True):
    """--db, then $SPYOGENES_TYPING_DB, then the user's own database (built
    with --setup), then the bundled one: database/ next to or above this
    script, or spyogenes_typing_database/ beside the installed module."""
    here = os.path.dirname(os.path.realpath(__file__))
    mine = user_database()
    if not (own and os.path.exists(os.path.join(mine, "blast", BLAST_LONG + ".nsq"))):
        mine = None
    for cand in (arg_db, os.environ.get(ENV_DB), mine,
                 os.path.join(here, "database"),
                 os.path.join(here, "..", "database"),
                 os.path.join(here, "spyogenes_typing_database")):
        if cand and os.path.isdir(cand):
            return os.path.realpath(cand)
    return os.path.realpath(arg_db or os.environ.get(ENV_DB) or os.path.join(here, "database"))


def db_version(db):
    try:
        with open(os.path.join(db, "VERSION.txt")) as fh:
            return fh.readline().strip() or "unknown"
    except OSError:
        return "unknown"


def read_table(path):
    """Rows of a TSV with '#' comment lines and a header."""
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


class Database:
    def __init__(self, path):
        self.path = path
        need = [os.path.join(path, "blast", BLAST_LONG + ".nsq"),
                os.path.join(path, "blast", BLAST_SHORT + ".nsq"),
                os.path.join(path, "pubmlst", MLST_SCHEME, MLST_SCHEME + ".txt")]
        need += [os.path.join(path, "tables", t) for t in
                 ("emm_clusters.tsv", "emm_patterns.tsv", "pilin_nt80.tsv",
                  "pilin_types.tsv", "fct_regions.tsv")]
        for f in need:
            if not os.path.exists(f):
                err(f"Database incomplete or not found: {path} (missing {os.path.relpath(f, path)})")
        self.version = db_version(path)
        self.blast_long = os.path.join(path, "blast", BLAST_LONG)
        self.blast_short = os.path.join(path, "blast", BLAST_SHORT)
        # MLST profiles
        with open(os.path.join(path, "pubmlst", MLST_SCHEME, MLST_SCHEME + ".txt")) as fh:
            rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
        self.mlst_loci = [c for c in rows[0][1:] if c not in ("clonal_complex", "species", "CC")]
        n = len(self.mlst_loci)
        self.profiles = {"/".join(r[1:n + 1]): r[0] for r in rows[1:]}
        # curated tables (PubMLST guide, Feb 2024)
        t = os.path.join(path, "tables")
        self.clusters = {r["emm_type"].upper(): r["emm_cluster"]
                         for r in read_table(f"{t}/emm_clusters.tsv")}
        self.patterns = defaultdict(list)
        for r in read_table(f"{t}/emm_patterns.tsv"):
            if r["pattern"] not in self.patterns[r["emm_type"].upper()]:
                self.patterns[r["emm_type"].upper()].append(r["pattern"])
        self.nt80 = {(r["locus"], r["allele"]): r["nt80_cluster"]
                     for r in read_table(f"{t}/pilin_nt80.tsv")}
        self.pilin_types = {(r["pilA_nt80"], r["pilB_nt80"]): (r["pilin_type"], r["fct_form"])
                            for r in read_table(f"{t}/pilin_types.tsv")}
        self.fct_rules = [(r["fct_form"], r["regulator"], frozenset(r["pilin_loci"].split("+")))
                          for r in read_table(f"{t}/fct_regions.tsv")]
        m1 = os.path.join(path, "references", "m1_snps.tsv")
        self.m1_snps = read_table(m1) if os.path.exists(m1) else []
        # flanks as regex pieces; a neighbouring SNP site may hold either base
        sites = {int(x["position"]): x for x in self.m1_snps}
        for x in self.m1_snps:
            p, left, right = int(x["position"]), x["left"], x["right"]
            x["left_re"] = [f"[{sites[q]['ref']}{sites[q]['alt']}]" if q in sites else c
                            for q, c in zip(range(p - len(left), p), left)]
            x["right_re"] = [f"[{sites[q]['ref']}{sites[q]['alt']}]" if q in sites else c
                             for q, c in zip(range(p + 1, p + 1 + len(right)), right)]


# --------------------------------------------------------------------------
# Input: FASTA / GenBank / EMBL, plain or .gz / .bz2 / .zip (as mlst-typing)
# --------------------------------------------------------------------------

def read_bytes(path):
    with open(path, "rb") as fh:
        data = fh.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    elif data[:3] == b"BZh":
        data = bz2.decompress(data)
    elif data[:4] == b"PK\x03\x04":
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            files = [n for n in z.namelist() if not n.endswith("/")]
            if not files:
                err(f"Zip file has no files in it: {path}")
            data = z.read(files[0])
    return data.decode("utf-8", errors="replace")


def parse_sequences(text):
    """Return a list of (id, sequence) from FASTA, GenBank or EMBL text."""
    lines = text.splitlines()
    first = next((l for l in lines if l.strip()), "")
    seqs = []
    if first.startswith(">"):
        name, parts = None, []
        for line in lines:
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(parts)))
                name, parts = (line[1:].split() or ["seq"])[0], []
            elif name is not None:
                parts.append("".join(line.split()))
        if name is not None:
            seqs.append((name, "".join(parts)))
    elif first.startswith("LOCUS") or first.startswith("ID "):
        genbank = first.startswith("LOCUS")
        name, parts, inseq = None, [], False
        for line in lines:
            if genbank and line.startswith("LOCUS"):
                name = (line.split()[1:2] or ["seq"])[0]
            elif not genbank and line.startswith("ID "):
                name = line[5:].split(";")[0].strip() or "seq"
            elif (genbank and line.startswith("ORIGIN")) or (not genbank and line.startswith("SQ ")):
                inseq, parts = True, []
            elif line.startswith("//"):
                if inseq:
                    seqs.append((name or f"seq{len(seqs) + 1}", "".join(parts)))
                inseq = False
            elif inseq:
                parts.append(re.sub(r"[^A-Za-z*-]", "", line))
        if inseq:
            seqs.append((name or f"seq{len(seqs) + 1}", "".join(parts)))
    else:
        err("Input is not FASTA, GenBank or EMBL")
    seqs = [(n, s.upper()) for n, s in seqs if s]
    if not seqs:
        err("Sequence contains no data")
    return seqs


def write_fasta(seqs, path):
    """Contigs are renamed c1, c2, ... so odd FASTA ids cannot upset BLAST."""
    with open(path, "w") as fh:
        for i, (_, seq) in enumerate(seqs, 1):
            fh.write(f">c{i}\n")
            for j in range(0, len(seq), 80):
                fh.write(seq[j:j + 80] + "\n")


def sample_name(path):
    base = os.path.basename(path)
    base = re.sub(r"\.(gz|bz2|zip)$", "", base, flags=re.I)
    return re.sub(r"\.(fasta|fas|fa|fna|ffn|fsa|contigs|gbk|gb|gbff|embl)$", "", base, flags=re.I)


def revcom(dna):
    return dna[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


# --------------------------------------------------------------------------
# BLAST hits and regions
# --------------------------------------------------------------------------

class Hit:
    """One BLAST HSP of a database allele on a contig (the contig is the query)."""
    __slots__ = ("scheme", "locus", "allele", "slen", "alen", "nident", "contig",
                 "start", "end", "slo", "shi", "qlen", "minus")

    def __init__(self, f):
        _, self.scheme, self.locus, self.allele = f[0].split(SEP)
        self.slen, self.alen, self.nident = int(f[1]), int(f[2]), int(f[3])
        self.contig = f[4]
        q1, q2, s1, s2 = int(f[5]), int(f[6]), int(f[7]), int(f[8])
        self.start, self.end = min(q1, q2), max(q1, q2)
        self.slo, self.shi = min(s1, s2), max(s1, s2)
        self.qlen = int(f[9])
        self.minus = f[10] == "minus"

    @property
    def pid(self):                # % identity over the aligned part
        return 100.0 * self.nident / self.alen

    @property
    def cov(self):                # share of the allele in the alignment
        return (self.shi - self.slo + 1) / self.slen

    @property
    def full(self):
        return self.slo == 1 and self.shi == self.slen

    @property
    def exact(self):
        return self.full and self.alen == self.slen and self.nident == self.slen

    @property
    def edge(self):
        """The alignment stops at a contig end (within EDGE bp) on a side
        where part of the allele is missing: the gene is cut by the edge."""
        miss5, miss3 = self.slo - 1, self.slen - self.shi
        left, right = (miss3, miss5) if self.minus else (miss5, miss3)
        return bool((left and self.start - 1 <= EDGE) or (right and self.qlen - self.end <= EDGE))

    @property
    def quality(self):
        """Sort key: exact first, then full length, then fewest differences."""
        diffs = (self.alen - self.nident) + (self.slen - (self.shi - self.slo + 1))
        return (self.exact, self.full, -diffs, self.nident)

    @property
    def name(self):
        return f"{self.locus}_{self.allele}"

    @property
    def diff_text(self):
        """Differences from the allele in words: '1 SNP', '2 SNPs, 3 bp indel'."""
        gaps = (self.alen - (self.shi - self.slo + 1)) + (self.alen - (self.end - self.start + 1))
        snps = self.alen - self.nident - gaps
        parts = [f"{snps} SNP{'s' if snps != 1 else ''}"] if snps else []
        if gaps:
            parts.append(f"{gaps} bp indel")
        return ", ".join(parts) or "0 SNPs"


def run_blast(fasta, db, short, threads, outfile):
    if short:   # 19-54 bp signature markers
        opts = ["-task", "blastn-short", "-word_size", "7", "-evalue", "1",
                "-perc_identity", "85"]
    else:
        opts = ["-task", "megablast", "-word_size", "20", "-evalue", "1e-10",
                "-perc_identity", "80"]
    cmd = ["blastn", "-query", fasta, "-db", db, "-out", outfile, "-dust", "no",
           "-num_threads", str(threads), "-max_target_seqs", "100000",
           "-outfmt", BLAST_FMT, *opts]
    Log.dbg("Running: " + " ".join(cmd))
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        err(f"blastn failed: {' '.join(cmd)}\n{r.stderr.strip()}")
    with open(outfile) as fh:
        return [Hit(l.rstrip("\n").split("\t")) for l in fh if l.strip()]


class Region:
    """Overlapping hits at one place on a contig; best = the closest allele."""

    def __init__(self, hit):
        self.contig, self.start, self.end = hit.contig, hit.start, hit.end
        self.hits = [hit]

    @property
    def best(self):
        return self.hits[0]

    def overlaps(self, h):
        if h.contig != self.contig:
            return False
        ov = min(self.end, h.end) - max(self.start, h.start) + 1
        return ov >= 0.5 * min(self.end - self.start + 1, h.end - h.start + 1)


def regions(hits):
    """Group hits by place, best hit first; regions in contig order."""
    out = []
    for h in sorted(hits, key=lambda h: h.quality, reverse=True):
        for r in out:
            if r.overlaps(h):
                r.hits.append(h)
                break
        else:
            out.append(Region(h))
    return sorted(out, key=lambda r: (int(r.contig[1:]), r.start))


def call(hit, minid):
    """Status of the closest allele at a region: PERFECT (exact), NOVEL (full
    length, >= minid), PARTIAL (cut by a contig edge, aligned part >= minid),
    or None (too far from any allele to name)."""
    if hit.exact:
        return PERFECT
    if hit.pid >= minid and hit.cov >= 0.95:
        return NOVEL
    if hit.pid >= minid and hit.edge and hit.cov >= 0.5:
        return PARTIAL
    return None


def gene_regions(hits, minid=GENE_ID, mincov=GENE_COV):
    """Regions of whole genes (superantigens, resistance, pilins) whose
    closest allele passes the identity and coverage limits."""
    return [r for r in regions(hits)
            if r.best.pid >= minid and (r.best.cov >= mincov or r.best.edge)]


def near_misses(hits, minid=GENE_ID, mincov=GENE_COV, floor=NEAR_ID):
    """Places where a gene is clearly present but too far from any known
    allele to be called: identity between floor and minid, or a short piece
    of one. Reported so the user knows the tool left something out, instead
    of a silent absence. Only hits BELOW the calling limit appear here; a
    hit that was called is in its own column already."""
    out = []
    for r in regions(hits):
        b = r.best
        if b.pid >= minid and (b.cov >= mincov or b.edge):
            continue                      # called, not a near miss
        if b.pid < floor or b.cov < 0.5:
            continue                      # too little to say anything about
        out.append(b)
    return out


def near_text(kind, hits):
    """One warning line for the near misses of one gene set."""
    if not hits:
        return []
    items = ", ".join(f"{h.locus} {h.pid:.1f}% id / {100 * h.cov:.0f}% of the gene"
                      for h in sorted(hits, key=lambda h: -h.pid)[:6])
    more = f" and {len(hits) - 6} more" if len(hits) > 6 else ""
    return [f"{kind} not called, nearest known allele under {GENE_ID:g}% identity "
            f"(a new variant, or a disrupted gene): {items}{more}"]


def gene_label(h):
    """'ermB(1)' exact, 'ermB(~1)' full length with differences, 'ermB(1?)' partial."""
    if h.exact:
        return f"{h.locus}({h.allele})"
    if h.cov >= 0.95:
        return f"{h.locus}(~{h.allele})"
    return f"{h.locus}({h.allele}?)"


# --------------------------------------------------------------------------
# emm type, emm-like genes, emm cluster and emm pattern
# --------------------------------------------------------------------------

def emm_type_of(h):
    """'1' from emm_1.0, '134' from emmL_134.2, 'STG653' from emmCG_STG653.0."""
    return h.allele.split(".")[0].upper()


def emm_label(h):
    """How the allele at the emm locus is written: emm1.0, emm134.2, STG653.0."""
    return h.allele if h.locus == "emmCG" else f"emm{h.allele}"


def tied_types(region):
    """emm types of every hit as good as the best one (a piece of emm cut by a
    contig end often fits many subtypes equally well)."""
    return list(dict.fromkeys(emm_type_of(h) for h in region.hits
                              if h.quality == region.best.quality))


def other_label(h, minid):
    """emm-like and emmCG hits away from the emm locus: emmL134.2,
    'emmL134.2 (98.3%)' when not exact, emm hits as 'emm59.0'."""
    name = h.allele if h.locus == "emmCG" else f"{h.locus}{h.allele}"
    if h.exact:
        return name
    return f"{name} ({h.pid:.1f}%{', partial' if h.cov < 0.95 else ''})"


def near_prime(region, primes):
    """True when an emmPrime hit sits just upstream of the region, same strand."""
    b = region.best
    for p in primes:
        if p.contig != b.contig or p.minus != b.minus:
            continue
        gap = (b.start - p.end) if not b.minus else (p.start - b.end)
        if 0 <= gap <= PRIME_GAP:
            return True
    return False


class EmmResult:
    def __init__(self):
        self.emm = "-"
        self.status = MISSING
        self.types = []          # emm type(s) for cluster / pattern
        self.emmL = []
        self.emmCG = []
        self.cluster = "-"
        self.pattern = "-"
        self.region = None       # the emm locus region (for the marker window)
        self.warnings = []


def type_emm(long_hits, short_hits, db, o):
    r = EmmResult()
    fam = [h for h in long_hits if h.scheme == EMM_SCHEME and h.locus in EMM_FAMILY]
    regs = [g for g in regions(fam) if g.best.pid >= 80 and g.best.cov >= 0.5]
    primes = [h for h in short_hits if h.locus == "emmPrime" and h.exact]
    at_prime = [g for g in regs if near_prime(g, primes)]
    others = [g for g in regs if g not in at_prime]
    if at_prime:
        emm_regs = at_prime
    else:
        # no emmPrime next to any hit: only an emm-locus allele may be emm
        emm_regs = [g for g in regs if g.best.locus == "emm" and call(g.best, o.minid)]
        others = [g for g in regs if g not in emm_regs]
        if emm_regs:
            r.warnings.append("emmPrime not found next to emm: emm taken from the emm locus hit")
        elif regs:
            r.warnings.append("no emm at an emmPrime site; emm-like hits only")
    emm_types = {t for g in emm_regs for t in tied_types(g)}
    for g in others:
        c = call(g.best, o.minid)
        if c is None:
            continue
        if c == PARTIAL and emm_types & set(tied_types(g)):
            continue            # the other piece of an emm cut by a contig end
        (r.emmCG if g.best.locus == "emmCG" else r.emmL).append(other_label(g.best, o.minid))

    if not emm_regs:
        r.status = MISSING
    else:
        calls = [(g, call(g.best, o.minid)) for g in emm_regs]
        named = [(g, c) for g, c in calls if c]
        if not named:
            g = emm_regs[0]
            r.status = NOVEL
            r.warnings.append(f"possible new emm type: nearest {emm_label(g.best)} at "
                              f"{g.best.pid:.1f}% identity (below {o.minid:g}%)")
            r.region = g
        else:
            labels = []
            for g, c in named:
                lab = emm_label(g.best)
                if lab not in labels:
                    labels.append(lab)
            r.region = named[0][0]
            r.emm = ",".join(labels)
            r.types = list(dict.fromkeys(emm_type_of(g.best) for g, _ in named))
            if len(named) > 1:
                r.status = MULTICOPY if len(labels) == 1 else MIXED
                if r.status == MIXED:
                    r.warnings.append(f"{len(named)} emm loci with different alleles: {r.emm}")
            else:
                g, c = named[0]
                r.status = c
                if c == NOVEL:
                    r.warnings.append(f"emm NOVEL: {g.best.diff_text} from {emm_label(g.best)} "
                                      f"({g.best.pid:.1f}%); cluster and pattern "
                                      f"from emm{emm_type_of(g.best)}")
                elif c == PARTIAL:
                    types = tied_types(g)
                    if g.best.locus != "emmCG":
                        r.emm = ",".join(f"emm{t}" for t in types)
                    r.types = types
                    r.warnings.append(f"emm PARTIAL: emm cut by a contig end ({100 * g.best.cov:.0f}% "
                                      f"of the allele), subtype unknown")
            for g, _ in named:
                if g.best.locus == "emmL":
                    r.warnings.append(f"emm type from emmL{g.best.allele} at the emm locus")
    # cluster: Sanderson-Smith 2014, as PubMLST lists it
    if r.types:
        r.cluster, cluster_warn = emm_cluster(r.types, db)
        r.warnings += cluster_warn
    r.pattern = emm_pattern(r, short_hits, db, o)
    return r


def emm_cluster(types, db):
    """(cluster column, warnings) for the emm type(s) called. The guide's
    cluster list does not cover every type, so a type it misses is said in
    Warnings rather than left as a bare '-' that reads like a failure."""
    cluster = ",".join(db.clusters.get(t, "-") for t in types)
    missing = [t for t in types if t not in db.clusters]
    warnings = []
    if missing:
        warnings.append("no emm cluster listed in the typing guide for "
                        + ", ".join(f"emm{t}" for t in missing)
                        + " (Sanderson-Smith 2014 does not cover every type)")
    if len(missing) == len(types):
        cluster = "-"
    return cluster, warnings


def marker_counts(short_hits, region):
    """Number of places each mga / SF marker hits within EMM_WINDOW of emm,
    whether any hit is not exact, and whether the window runs off the contig."""
    b = region.best
    lo, hi = b.start - EMM_WINDOW, b.end + EMM_WINDOW
    near = [h for h in short_hits if h.locus in PATTERN_LOCI and h.contig == b.contig
            and lo <= h.start and h.end <= hi and h.full and h.pid >= MARKER_ID]
    counts, inexact = {}, []
    for locus in PATTERN_LOCI:
        regs = regions([h for h in near if h.locus == locus])
        if regs:
            counts[locus] = len(regs)
            inexact += [g.best.locus for g in regs if not g.best.exact]
    cut = lo < 1 or hi > b.qlen
    return counts, inexact, cut


def emm_pattern(r, short_hits, db, o):
    """Measured pattern (A-E or REA) from the markers next to emm (guide p4);
    the type's predicted pattern (guide p5) when it cannot be measured."""
    predicted = []
    for t in r.types:
        predicted += [p for p in db.patterns.get(t, []) if p not in predicted]
    pred = "/".join(predicted)
    if r.region is None:
        return "-"
    counts, inexact, cut = marker_counts(short_hits, r.region)
    found = ";".join(f"{k}x{v}" if v > 1 else k for k, v in counts.items()) or "none"
    measured = next((p for p, need in PATTERNS if counts == need), None)
    if inexact:
        r.warnings.append(f"emm pattern markers not exact (new variant?): {','.join(inexact)}")
    if measured is None and r.status in (PERFECT, NOVEL, PARTIAL, MULTICOPY):
        mga = "mga-1" in counts or "mga-2" in counts
        if mga and not cut:
            measured = "REA"
            extra = f"; predicted for this emm type: {pred}" if pred else ""
            r.warnings.append(f"emm pattern markers {found} fit no pattern A-E: REA{extra}")
    if measured:
        group = "A-C" if measured in "ABC" else measured
        if predicted and group not in predicted and measured != "REA":
            r.warnings.append(f"measured emm pattern {measured} differs from predicted {pred}")
        return measured
    if pred:
        r.warnings.append(f"emm pattern not measurable (markers: {found}): predicted from emm type")
        return f"{pred} (predicted)"
    return "-"


# --------------------------------------------------------------------------
# MLST (allele calls and status words as mlst-typing)
# --------------------------------------------------------------------------

def collapse_copies(c):
    parts = c.split(",")
    return parts[0] if len(parts) > 1 and len(set(parts)) == 1 else c


def mlst_warnings(st, code, loci, status):
    """Say in Warnings why the ST is not a plain match, so the column that
    lists everything worth checking covers MLST too (the status column alone
    does not say which locus is at fault)."""
    bad = lambda f: ", ".join(f"{l}({c})" for l, c in zip(loci, code) if f(c))
    if status == PERFECT:
        return []
    if status == MULTICOPY:
        where = bad(lambda c: collapse_copies(c) != c)
        tail = "; the ST still stands" if st != "-" else ""
        return [f"MLST MULTICOPY: identical copies of {where}{tail}"]
    if status == NOVEL:
        return [f"MLST NOVEL: every allele is an exact match but the combination "
                f"{'/'.join(code)} is not a known ST"]
    if status == MIXED:
        return [f"MLST MIXED: more than one allele at {bad(lambda c: ',' in collapse_copies(c))}"
                " (contamination, or two strains in one assembly)"]
    if status == MISSING:
        return [f"MLST MISSING: {bad(lambda c: c == '-')} not found"]
    if status == "NONE":
        return ["MLST NONE: no MLST locus found"]
    near = bad(lambda c: c.startswith("~"))
    part = bad(lambda c: c.endswith("?"))
    bits = [x for x in (f"{near} not an exact match" if near else "",
                        f"{part} cut by a contig end" if part else "") if x]
    return [f"MLST {status}: " + "; ".join(bits)] if bits else [f"MLST {status}"]


def type_mlst(long_hits, db, o):
    """(ST, [allele call per locus], status, loci found at >= minid, warnings)."""
    code, good = [], 0
    for locus in db.mlst_loci:
        regs = regions([h for h in long_hits if h.scheme == MLST_SCHEME and h.locus == locus])
        exact = sorted({g.best.allele for g in regs if g.best.exact}, key=float)
        if exact:
            # identical copies stay "5,5" so MULTICOPY can be seen
            n = sum(1 for g in regs if g.best.exact)
            c = ",".join(exact) if len(exact) > 1 else ",".join([exact[0]] * n)
        else:
            c = "-"
            for g in regs:
                st = call(g.best, o.minid)
                if st == NOVEL:
                    c = f"~{g.best.allele}"
                    break
                if st == PARTIAL and c == "-":
                    c = f"{g.best.allele}?"
        if exact or any(g.best.pid >= o.minid and g.best.cov >= 0.95 for g in regs):
            good += 1
        code.append(c)
    st = db.profiles.get("/".join(collapse_copies(c) for c in code), "-")
    status = mlst_status(st, code)
    return st, code, status, good, mlst_warnings(st, code, db.mlst_loci, status)


def mlst_status(st, code):
    multi = any(collapse_copies(c) != c for c in code)
    if st != "-":
        return MULTICOPY if multi else PERFECT
    if all(c == "-" for c in code):
        return "NONE"
    if any("," in collapse_copies(c) for c in code):
        return MIXED
    if multi:
        return MULTICOPY
    if not any(re.search(r"[^\d,]", c) for c in code):
        return NOVEL            # every allele exact, new combination
    if any(c == "-" for c in code):
        return MISSING
    score = 100 - 30 * sum(c.startswith("~") for c in code) / len(code) \
        - 50 * sum(c.endswith("?") for c in code) / len(code)
    return PARTIAL if score >= 70 else "BAD"


# --------------------------------------------------------------------------
# M1 lineage (emm1 only): M1UK / M1DK SNPs on MGAS5005 (NC_007297.2)
# --------------------------------------------------------------------------

def base_at(genome, rc, left, right):
    """The base between two flanks (lists of regex pieces), found once in the
    genome (either strand); flanks are shortened when a difference near the
    site stops the longer match."""
    for k in (30, 20, 14):
        pat = re.compile("".join(left[-k:]) + "([ACGT])" + "".join(right[:k]))
        found = [m.group(1) for m in pat.finditer(genome)]
        found += [m.group(1) for m in pat.finditer(rc)]
        if len(found) == 1:
            return found[0]
        if len(found) > 1:
            return None
    return None


def type_m1(seqs, emm, db):
    if "1" not in emm.types or not db.m1_snps:
        return "-", []
    genome = "N".join(s for _, s in seqs)
    rc = revcom(genome)
    tally = defaultdict(lambda: [0, 0, 0])      # lineage -> [alt, found, total]
    for s in db.m1_snps:
        t = tally[s["lineage"]]
        t[2] += 1
        b = base_at(genome, rc, s["left_re"], s["right_re"])
        if b is None:
            continue
        t[1] += 1
        if b == s["alt"]:
            t[0] += 1
    warnings = []
    for lin, (alt, found, total) in tally.items():
        if found < total:
            warnings.append(f"{lin}: {total - found} of {total} SNP sites not found")
    uk_alt, uk_found, uk_total = tally.get("M1UK", [0, 0, 0])
    dk_alt, dk_found, dk_total = tally.get("M1DK", [0, 0, 0])
    if uk_found < 0.9 * uk_total:
        return "-", warnings + ["M1 lineage not called: too few M1UK SNP sites found"]
    if uk_alt == uk_found:
        return f"M1UK ({uk_alt}/{uk_total} SNPs)", warnings
    if uk_alt:
        return f"M1intermediate ({uk_alt}/{uk_total} M1UK SNPs)", warnings
    if dk_total and dk_found >= 0.9 * dk_total and dk_alt == dk_found:
        return f"M1DK ({dk_alt}/{dk_total} SNPs)", warnings
    return "M1global (0 M1UK SNPs)", warnings


# --------------------------------------------------------------------------
# Superantigens, resistance genes, FCT-region form and pilin type
# --------------------------------------------------------------------------

def type_superantigens(long_hits):
    sag = [h for h in long_hits if h.scheme == SAG_SCHEME]
    genes = {}
    for g in gene_regions(sag):
        b = g.best
        name = b.locus if b.cov >= 0.95 else f"{b.locus}(partial)"
        if b.locus not in genes or "(" in genes[b.locus]:
            genes[b.locus] = name
    called = set(genes)
    near = [h for h in near_misses(sag) if h.locus not in called]
    return [genes[k] for k in sorted(genes, key=str.lower)], near_text("superantigen", near)


def type_resistance(long_hits):
    res = [h for h in long_hits if h.scheme == RES_SCHEME]
    calls = [gene_label(g.best) for g in gene_regions(res)]
    called = {g.best.locus for g in gene_regions(res)}
    near = [h for h in near_misses(res) if h.locus not in called]
    return calls, near_text("resistance gene", near)


def type_fct(long_hits, short_hits, db):
    """FCT-region form from rofA / nra and the pilA / pilB loci (guide p15),
    pilin type from the pilA and pilB nt80 clusters (guide p8-14)."""
    warnings = []
    regs = [h.locus for h in short_hits
            if h.locus in ("rofA", "nra") and h.full and h.pid >= MARKER_ID]
    regulators = sorted(set(regs))
    pila = [g.best for g in gene_regions([h for h in long_hits if h.scheme == "fct_pila"])]
    pilb = [g.best for g in gene_regions([h for h in long_hits if h.scheme == "fct_pilb"])]
    if not (pila or pilb or regulators):
        return "-", warnings
    if len(pila) > 1 or len(pilb) > 1:
        warnings.append("more than one pilA or pilB: "
                        + ",".join(h.name for h in pila + pilb))
    nt80a = {h.name: db.nt80.get((h.locus, h.allele)) if h.exact else None for h in pila}
    nt80b = {h.name: db.nt80.get((h.locus, h.allele)) if h.exact else None for h in pilb}
    loci = set()
    for h in pila:
        # PubMLST files the guide's pilA10 alleles under pilA5 (nt80_2)
        is10 = h.locus == "pilA5" and (nt80a[h.name] == "nt80_2"
                                       or any(b.locus == "pilB10" for b in pilb))
        loci.add("pilA10" if is10 else h.locus)
    loci |= {h.locus for h in pilb}
    tries = [r for r in regulators] or [""]
    forms = []
    for reg in tries:
        forms += [f for f, rr, need in db.fct_rules if rr == reg and need == loci]
    forms = list(dict.fromkeys(forms))
    if len(regulators) > 1:
        warnings.append("both rofA and nra found")
    for h in pila + pilb:
        if not h.exact:
            warnings.append(f"{gene_label(h)} not exact: no nt80 cluster, pilin type not called")
    if not forms:
        found = ",".join(regulators + sorted(loci)) or "none"
        warnings.append(f"FCT hits {found} fit no FCT-region form")
        return "-", warnings
    form = ",".join(forms)      # "FCT-3/4" is itself one form
    # pilin type: the pair of nt80 clusters
    a = [v for v in nt80a.values()]
    b = [v for v in nt80b.values()]
    if len(a) <= 1 and len(b) == 1 and all(a + b):
        key = (a[0] if a else "", b[0])
        if key in db.pilin_types:
            ptype, pform = db.pilin_types[key]
            if pform and pform not in forms:
                warnings.append(f"{ptype} is listed as {pform}, hits give {form}")
            return f"{form} ({ptype})", warnings
        warnings.append(f"pilA {key[0] or 'none'} + pilB {key[1]}: no pilin type listed")
    else:
        missing = [n for n, v in {**nt80a, **nt80b}.items() if v is None]
        exact_missing = [n for n in missing if any(h.name == n and h.exact for h in pila + pilb)]
        if exact_missing:
            warnings.append(f"no nt80 cluster listed for {','.join(exact_missing)}: "
                            "pilin type not called")
    return form, warnings


# --------------------------------------------------------------------------
# One genome
# --------------------------------------------------------------------------

def join(items, sep=";"):
    return sep.join(items) if items else "-"


def type_genome(path, label, db, o, tmp):
    seqs = parse_sequences(read_bytes(path))
    names = [n for n, _ in seqs]
    fna = os.path.join(tmp, "genome.fna")
    write_fasta(seqs, fna)
    long_hits = run_blast(fna, db.blast_long, False, o.threads, os.path.join(tmp, "long.bls"))
    short_hits = run_blast(fna, db.blast_short, True, o.threads, os.path.join(tmp, "short.bls"))
    Log.dbg(f"{len(long_hits)} long and {len(short_hits)} short hits")

    emm = type_emm(long_hits, short_hits, db, o)
    st, code, mlst_st, good, mlst_warn = type_mlst(long_hits, db, o)
    m1, m1_warn = type_m1(seqs, emm, db)
    sags, sag_warn = type_superantigens(long_hits)
    res, res_warn = type_resistance(long_hits)
    fct, fct_warn = type_fct(long_hits, short_hits, db)
    n = len(db.mlst_loci)
    species = ("S. pyogenes" if good >= SPECIES_LOCI
               else f"not S. pyogenes? ({good}/{n} MLST loci)")
    warnings = (emm.warnings + mlst_warn + m1_warn + fct_warn
                + sag_warn + res_warn)
    if good < SPECIES_LOCI:
        warnings.insert(0, f"only {good} of {n} S. pyogenes MLST loci found at "
                           f">= {o.minid:g}% identity: check the species")
    row = {
        "Sample": label,
        "Species check": species,
        "emm": emm.emm,
        "emm status": emm.status,
        "emmL": join(emm.emmL),
        "emmCG": join(emm.emmCG),
        "emm cluster": emm.cluster,
        "emm pattern": emm.pattern,
        "ST": st,
        "alleles": [f"{l}({c})" for l, c in zip(db.mlst_loci, code)],
        "MLST status": mlst_st,
        "M1 lineage": m1,
        "Superantigens": join(sags),
        "Resistance genes": join(res),
        "FCT type": fct,
        "Warnings": join(warnings, "; "),
        "Database": db.version,
    }
    hits = detail_rows(label, names, long_hits, short_hits, o) if o.hits else []
    return row, hits


def detail_rows(label, names, long_hits, short_hits, o):
    """--hits: the closest allele at every place something hit."""
    out = []
    groups = defaultdict(list)
    for h in long_hits:
        if h.pid >= 80:
            key = h.scheme if h.scheme != EMM_SCHEME else "emm"
            groups[("long", key)].append(h)
    for h in short_hits:
        if h.full and h.pid >= MARKER_ID:
            groups[("short", h.locus)].append(h)
    for (_, _), hs in sorted(groups.items()):
        for g in regions(hs):
            b = g.best
            out.append([label, b.scheme, b.locus, b.allele, f"{b.pid:.2f}",
                        f"{100 * b.cov:.1f}", names[int(b.contig[1:]) - 1], b.start, b.end,
                        "-" if b.minus else "+", call(b, o.minid) or "-"])
    return out


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

HEAD = ["Sample", "Species check", "emm", "emm status", "emmL", "emmCG",
        "emm cluster", "emm pattern", "ST"]
TAIL = ["MLST status", "M1 lineage", "Superantigens", "Resistance genes",
        "FCT type", "Warnings", "Database"]
HITS_HEAD = ["Sample", "Scheme", "Locus", "Allele", "Identity", "Coverage",
             "Contig", "Start", "End", "Strand", "Call"]


def header(db, o):
    if o.no_alleles:
        mid = []
    elif o.combine_alleles:
        mid = ["ALLELES"]
    else:
        mid = [f"ALLELE{i + 1}" for i in range(len(db.mlst_loci))]
    return HEAD + mid + TAIL


def cells(row, o):
    if o.no_alleles:
        mid = []
    elif o.combine_alleles:
        mid = [";".join(row["alleles"])]
    else:
        mid = row["alleles"]
    return [row[c] for c in HEAD] + mid + [row[c] for c in TAIL]


def tsv(fh, values):
    fh.write("\t".join(str(v).replace("\t", " ") for v in values) + "\n")


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------

def get_args(argv):
    p = argparse.ArgumentParser(
        prog=EXE,
        description="Streptococcus pyogenes typing from PubMLST: emm, emm cluster and "
                    "pattern, MLST, M1UK, superantigens, resistance and FCT",
        epilog=f"Homepage: {URL}")
    p.add_argument("files", nargs="*", metavar="genome",
                   help="assembled genome: FASTA/GenBank/EMBL, may be .gz/.bz2/.zip")
    g = p.add_argument_group("general")
    g.add_argument("--version", action="store_true", help="print version and exit")
    g.add_argument("--check", action="store_true", help="check dependencies and exit")
    g.add_argument("-q", "--quiet", action="store_true", help="no stderr messages")
    g.add_argument("--debug", action="store_true", help="verbose debug output")
    g.add_argument("--threads", type=int, default=1, help="BLAST threads [1]")
    g.add_argument("--fofn", default="", help="file of input filenames (optional TAB label)")
    g.add_argument("--label", default="", help="sample name for a single input file")
    g = p.add_argument_group("output")
    g.add_argument("--outfile", default="", help="write the report here [stdout]")
    g.add_argument("--combine-alleles", action="store_true",
                   help="MLST alleles in one ALLELES column (gki(4);gtr(3);...)")
    g.add_argument("--no-alleles", action="store_true",
                   help="no MLST allele columns: only ST and MLST status")
    g.add_argument("--json", default="", help="also write the results as JSON")
    g.add_argument("--hits", default="", metavar="FILE",
                   help="also write every hit (closest allele per place) to FILE")
    g = p.add_argument_group("scoring")
    g.add_argument("--minid", type=float, default=95,
                   help="%% identity of a full-length allele to name it (emm, MLST) [95]")
    g = p.add_argument_group("database")
    g.add_argument("--db", default="",
                   help=f"database folder [${ENV_DB}, your own (--setup) or the bundled one]")
    g = p.add_argument_group("up-to-date PubMLST data with your own key (command line only)")
    g.add_argument("--setup", action="store_true",
                   help="add or replace your PubMLST key and download your own database")
    g.add_argument("--update", action="store_true",
                   help="refresh your own database with the saved key")
    g.add_argument("--remove", action="store_true",
                   help="with --setup: forget the key and your own database")
    g.add_argument("--key-file", default="", metavar="FILE",
                   help="with --setup: read the key from FILE ('-' = stdin), ask nothing")
    g.add_argument("--oauth-start", action="store_true",
                   help="with --setup --key-file: OAuth step 1 (prints approve_url as JSON)")
    g.add_argument("--oauth-finish", default="", metavar="CODE",
                   help="with --setup: OAuth step 2 with PubMLST's verification code")
    g.add_argument("--update-check", action="store_true",
                   help="ask PubMLST whether anything is newer than your database (JSON)")
    g.add_argument("--setup-status", action="store_true",
                   help="print where setup stands as JSON and exit")
    o = p.parse_args(argv)
    if o.combine_alleles and o.no_alleles:
        p.error("--combine-alleles and --no-alleles cannot be used together")
    return o


def load_setup():
    """The own-key setup lives in spyogenes_typing_setup.py; the Galaxy tool
    ships without it and always passes --db."""
    try:
        import spyogenes_typing_setup
    except ImportError:
        return None
    return spyogenes_typing_setup


def run_setup(o):
    su = load_setup()
    if su is None:
        print("ERROR: setup is not part of this install; "
              "it comes with the command-line version (pip install)", file=sys.stderr)
        return 1
    if (o.remove or o.key_file or o.oauth_start or o.oauth_finish) and not o.setup:
        print("ERROR: --remove, --key-file, --oauth-start and --oauth-finish "
              "go with --setup", file=sys.stderr)
        return 1
    if o.oauth_start and not o.key_file:
        print("ERROR: --oauth-start reads the client key and client secret "
              "from --key-file FILE", file=sys.stderr)
        return 1
    paths = {"config": config_dir(), "db": user_database()}
    try:
        if o.setup_status:
            db = find_database(o.db)
            print(json.dumps({**su.status(paths), "database": db,
                              "database_version": db_version(db),
                              "version": VERSION}, indent=1))
            return 0
        if o.remove:
            return su.remove(paths)
        if o.oauth_start:
            return su.oauth_start(paths, o.key_file)
        bundled = find_database(o.db, own=False)
        if o.update_check:
            return su.check_newer(paths, bundled)
        if o.oauth_finish:
            return su.oauth_finish(paths, bundled, o.oauth_finish)
        if o.key_file:
            return su.setup_with_key(paths, bundled, o.key_file)
        if o.update:
            return su.update(paths, bundled)
        return su.setup(paths, bundled)
    except KeyboardInterrupt:
        print("\nStopped. Nothing was changed.", file=sys.stderr)
        return 1
    except (su.SetupError, RuntimeError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


def first_run_due(o):
    """True only for a first run typed by hand: a terminal on both ends, the
    default database, and no setup answer saved yet."""
    if o.db or os.environ.get(ENV_DB) or o.version or o.check:
        return False
    try:
        if not (sys.stdin.isatty() and sys.stderr.isatty()):
            return False
    except (AttributeError, ValueError):
        return False
    return not os.path.exists(os.path.join(config_dir(), "setup.json"))


def main(argv=None):
    o = get_args(argv)
    Log.quiet, Log.debug = o.quiet, o.debug
    if (o.setup or o.update or o.remove or o.key_file or o.setup_status
            or o.oauth_start or o.oauth_finish or o.update_check):
        return run_setup(o)
    if first_run_due(o):
        su = load_setup()
        if su is not None:
            su.first_run({"config": config_dir(), "db": user_database()},
                         find_database("", own=False))
    dbpath = find_database(o.db)
    if o.version:
        print(f"{EXE} {VERSION} (database {db_version(dbpath)})")
        return 0
    try:
        return run(o, dbpath)
    except TypingError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


def run(o, dbpath):
    Log.msg(f"This is {EXE} {VERSION} running on {sys.platform} with Python "
            f"{sys.version.split()[0]}")
    exe = shutil.which("blastn")
    if not exe:
        err("Could not find 'blastn' - please install BLAST+")
    Log.msg(f"Found 'blastn' => {exe}")
    db = Database(dbpath)
    Log.msg(f"Database: {dbpath} ({db.version})")
    if o.check:
        Log.msg("OK.")
        return 0

    files, labels = list(o.files), {}
    if o.fofn:
        files = []
        with open(o.fofn) as fh:
            for line in fh:
                if not line.strip():
                    continue
                path, _, lab = line.rstrip("\n").partition("\t")
                files.append(path)
                if lab:
                    labels[path] = lab
    if not files:
        err("Please give one or more assembled genomes (FASTA/GenBank/EMBL, may be .gz)")
    if o.label and len(files) > 1:
        err("--label only works with a single input file")
    for path in files:
        if os.path.isdir(path):
            err(f"'{path}' is a directory, not a file")
        if not os.access(path, os.R_OK):
            err(f"Unable to read from '{path}'")

    out = open(o.outfile, "w") if o.outfile else sys.stdout
    hits_fh = open(o.hits, "w") if o.hits else None
    rows = []
    try:
        tsv(out, header(db, o))
        if hits_fh:
            tsv(hits_fh, HITS_HEAD)
        with tempfile.TemporaryDirectory(prefix="spyogenes-typing_") as tmp:
            for i, path in enumerate(files, 1):
                # names from --fofn (Galaxy dataset names) lose their extension too
                label = o.label or sample_name(labels.get(path) or path)
                Log.msg(f"[{i}/{len(files)}] {label}")
                row, hits = type_genome(path, label, db, o, tmp)
                tsv(out, cells(row, o))
                out.flush()
                for h in hits:
                    tsv(hits_fh, h)
                rows.append(row)
    finally:
        if o.outfile:
            out.close()
        if hits_fh:
            hits_fh.close()
    if o.json:
        with open(o.json, "w") as fh:
            json.dump([json_row(r, db) for r in rows], fh, indent=2)
            fh.write("\n")
    Log.msg("Please cite PubMLST (Jolley et al. 2018, https://pubmed.ncbi.nlm.nih.gov/30345391)")
    Log.msg("Done.")
    return 0


def json_row(r, db):
    out = {k: v for k, v in r.items() if k != "alleles"}
    out["alleles"] = dict(zip(db.mlst_loci, (a[a.index("(") + 1:-1] for a in r["alleles"])))
    out["Warnings"] = r["Warnings"].split("; ") if r["Warnings"] != "-" else []
    return out


if __name__ == "__main__":
    sys.exit(main())
