#!/usr/bin/env python3
"""
spyogenes_typing_update.py - build or refresh database/ for spyogenes-typing.
Run it as scripts/update_database.py (same options):

    python scripts/update_database.py --db database --summary CHANGES_DB.md
    python scripts/update_database.py --db database --references   # also NCBI parts

What goes in the database folder
--------------------------------
  schemes.tsv        PubMLST schemes to download (kept, edited by hand)
  pubmlst/<scheme>/  <locus>.tfa allele FASTA files (+ <scheme>.txt, the STs,
                     for schemes with profiles: MLST)
  references/        superantigens.fa and m1_snps.tsv (built from NCBI with
                     --references, or when missing; carried over otherwise)
  tables/            emm_clusters.tsv, emm_patterns.tsv, fct_regions.tsv
                     (curated from the PubMLST typing guide; never downloaded,
                     always carried over)
  blast/             spyogenes (loci of 60 bp and longer, plus the
                     superantigen genes) and spyogenes_short (the short
                     signature loci) BLAST indexes
  manifest.tsv       one row per scheme (counts, last_updated, restricted)
  VERSION.txt        one line, e.g. PubMLST_2026-11-20

Data access
-----------
Without authentication PubMLST only returns data submitted up to 31 Dec 2024
(the standard, freely available data), and that is all this script fetches by
itself. Newer data needs the user's own PubMLST account and key:
`spyogenes-typing --setup` (spyogenes_typing_setup.py) sets SIGNER below and
runs this updater on a private copy of the database, outside the repository.

Exit codes: 0 = updated (changes written), 3 = no changes, other = error.
"""

import argparse
import csv
import datetime as dt
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request

NO_CHANGES = 3
USER_AGENT = "spyogenes-typing-updater"
SEP = "__"                      # BLAST titles: spy__<scheme>__<locus>__<allele>
SHORT = 60                      # loci whose longest allele is shorter go to the short index
BLAST_LONG = "spyogenes"
BLAST_SHORT = "spyogenes_short"
SAG_SCHEME = "superantigens"    # scheme name used for the superantigen genes in BLAST titles
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


# --------------------------------------------------------------------------
# HTTP
# --------------------------------------------------------------------------

# Requests with the user's own PubMLST key, only set by spyogenes_typing_setup:
# an object with wants(url) -> bool, apply(url) -> (url, extra headers) and
# renew() -> get a new OAuth session token. None = no login (free data).
SIGNER = None


class NotFound(RuntimeError):
    pass


def http_get(url, retries=4, pause=0.2):
    """Text of url. HTTP 404 raises NotFound at once (no retries)."""
    signed = SIGNER is not None and SIGNER.wants(url)
    renewed = False
    attempt = 0
    while True:
        try:
            target, extra = SIGNER.apply(url) if signed else (url, {})
            req = urllib.request.Request(target, headers={"User-Agent": USER_AGENT, **extra})
            with urllib.request.urlopen(req, timeout=300) as resp:
                data = resp.read()
            time.sleep(pause)  # be gentle with the public APIs
            return data.decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise NotFound(f"GET {url}: not found") from e
            # signed request refused: the session token (12 h) may have run
            # out, so renew it once; after that the account has no access
            if signed and e.code in (401, 403):
                if e.code == 401 and not renewed:
                    renewed = True
                    SIGNER.renew()
                    continue
                raise RuntimeError(f"GET {url} refused (HTTP {e.code})") from e
            attempt += 1
            if attempt >= retries:
                raise RuntimeError(f"GET {url} failed: {e}") from e
            time.sleep(5 * attempt)
        except Exception as e:  # noqa: BLE001
            attempt += 1
            if attempt >= retries:
                raise RuntimeError(f"GET {url} failed: {e}") from e
            time.sleep(5 * attempt)


def http_json(url):
    return json.loads(http_get(url))


# --------------------------------------------------------------------------
# PubMLST schemes
# --------------------------------------------------------------------------

def read_scheme_list(path):
    """schemes.tsv -> list of dicts (scheme, source, api, description)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t") + ["", "", ""]
            rows.append(dict(scheme=f[0], source=f[1], api=f[2], description=f[3]))
    names = [r["scheme"] for r in rows]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        raise RuntimeError(f"duplicate scheme names in {path}: {sorted(dup)}")
    if SAG_SCHEME in names:
        raise RuntimeError(f"{path}: the scheme name {SAG_SCHEME!r} is reserved")
    return rows


def normalise_fasta(text, locus):
    """Rewrite headers to >{locus}_{allele_id}; return (fasta, ids)."""
    out, ids = [], []
    prefix = locus + "_"
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            raw = line[1:].split()[0]
            allele = raw[len(prefix):] if raw.startswith(prefix) else raw.rsplit("_", 1)[-1]
            ids.append(allele)
            out.append(f">{locus}_{allele}")
        else:
            out.append(line.upper())
    return ("\n".join(out) + "\n") if out else "", ids


def limited(info):
    message = str(info.get("message", ""))
    return "restricted" in message or "2024-12-31" in message


def fetch_scheme(row, outdir):
    """Download one scheme into outdir/<scheme>/. Returns manifest dict.

    Loci without alleles (PubMLST answers 404) are skipped. Schemes with a
    primary key (MLST) also get their profiles as <scheme>.txt (ST + loci).
    """
    name, url = row["scheme"], row["api"].rstrip("/")
    info = http_json(url)
    loci = [u.rstrip("/").split("/")[-1] for u in info["loci"]]
    base = url.split("/schemes/")[0]
    sdir = os.path.join(outdir, name)
    os.makedirs(sdir)

    n_alleles, kept = 0, []
    for locus in loci:
        try:
            fasta, ids = normalise_fasta(http_get(f"{base}/loci/{locus}/alleles_fasta"), locus)
        except NotFound:
            print(f"  {locus}: no alleles, skipped", file=sys.stderr)
            continue
        if not ids:
            print(f"  {locus}: no alleles, skipped", file=sys.stderr)
            continue
        with open(os.path.join(sdir, f"{locus}.tfa"), "w") as fh:
            fh.write(fasta)
        n_alleles += len(ids)
        kept.append(locus)

    n_prof = 0
    if info.get("primary_key_field"):
        # profiles_csv is tab-delimited despite the name; keep ST + loci only
        prof = http_get(info.get("profiles_csv") or f"{url}/profiles_csv")
        rows = [r for r in csv.reader(io.StringIO(prof), delimiter="\t") if r]
        header = rows[0] if rows else ["ST"] + loci
        if header[0] != "ST" or any(l not in header for l in loci):
            raise RuntimeError(f"{name}: unexpected profiles header {header[:12]}")
        order = [l for l in header if l in loci]
        idx = [header.index(l) for l in order]
        with open(os.path.join(sdir, f"{name}.txt"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["ST"] + order)
            for r in rows[1:]:
                if not r[0]:
                    continue
                w.writerow([r[0]] + [r[i] if i < len(r) else "" for i in idx])
                n_prof += 1

    return dict(scheme=name, source=row["source"], loci=len(kept), alleles=n_alleles,
                profiles=n_prof, last_updated=info.get("last_updated") or "",
                restricted=limited(info), description=row["description"], api=url)


# --------------------------------------------------------------------------
# NCBI references (--references): superantigen genes and M1 lineage SNPs
# --------------------------------------------------------------------------

# Superantigen genes taken from the annotated CDS of complete genomes. Each
# distinct sequence of a gene is kept once, named after the first genome it
# was seen in.
SAG_GENOMES = ["NC_007297.2",   # MGAS5005 (M1)
               "NC_002737.2",   # SF370 (M1)
               "NC_004070.1",   # MGAS315 (M3)
               "NC_003485.1",   # MGAS8232 (M18)
               "NC_006086.1",   # MGAS10394 (M6)
               "NC_008022.1"]   # MGAS10270 (M2)
SAG_GENES = ["speA", "speC", "speG", "speH", "speI", "speJ", "speK", "speL", "speM", "smeZ"]
SSA_RECORD = "L29565.1"         # ssa (streptococcal superantigen), Mollick et al. 1993

# M1 lineage SNPs on MGAS5005 (NC_007297.2), position REF ALT. M1UK: 27 SNPs
# (Lynskey et al. 2019, Lancet Infect Dis); M1DK: Johannesen et al. 2023.
# Positions as listed in assembly_snptyper (github.com/boasvdp/assembly_snptyper,
# data/M1UK.vcf and data/M1DK.vcf, MIT licence).
M1_REFERENCE = "NC_007297.2"
M1_SNPS = {
    "M1UK": """115646 C T;116162 A C;116163 C A;250832 T C;513254 G T;528360 A T;
563631 G A;613633 T C;626494 G A;661707 G A;730823 C T;784467 T C;819098 G A;
923079 G A;942633 G T;983438 G C;1082253 C T;1238124 G A;1238673 G A;
1251193 G A;1373176 C A;1407497 C T;1446116 C T;1535209 A G;1702540 C T;
1734749 G A;1828734 G A""",
    "M1DK": """50024 A G;600192 T C;655648 C T;751431 G A;798142 C T;874083 G A;
890776 G A;908687 C T;948344 A C;967436 A G;1045871 G A;1058475 G A;
1062670 G A;1272188 T C;1484362 C A""",
}
FLANK = 50


def efetch_genbank(acc):
    return http_get(f"{EUTILS}?db=nuccore&id={acc}&rettype=gbwithparts&retmode=text", pause=0.5)


def genbank_sequence(text):
    seq = text.split("\nORIGIN", 1)[1].split("\n//")[0]
    return re.sub(r"[^acgtnACGTN]", "", seq).upper()


def revcom(dna):
    return dna[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def genbank_cds(text):
    """(gene, sequence) of every complete CDS with a /gene name and a simple
    location (start..end, maybe complement). Joined or partial CDS are skipped."""
    seq = genbank_sequence(text)
    feats = text.split("\nFEATURES", 1)[1].split("\nORIGIN", 1)[0]
    out = []
    for block in re.split(r"\n     (?=\S)", feats):
        if not block.startswith("CDS "):
            continue
        loc = block.split("\n", 1)[0][3:].strip()
        gene = re.search(r'/gene="([^"]+)"', block)
        if not gene or "/pseudo" in block:
            continue
        m = re.fullmatch(r"(complement\()?(\d+)\.\.(\d+)\)?", loc)
        if not m:
            continue
        s = seq[int(m.group(2)) - 1:int(m.group(3))]
        out.append((gene.group(1), revcom(s) if m.group(1) else s))
    return out


def build_references(refdir):
    """Write references/superantigens.fa and references/m1_snps.tsv from NCBI."""
    os.makedirs(refdir)
    seen, records = {}, []          # gene -> set of sequences
    m1_seq = None
    # NCBI spells some names in lower case (spej, speh, SSA): match without case
    canonical = {g.lower(): g for g in SAG_GENES}
    for acc in SAG_GENOMES:
        print(f"  NCBI {acc}", file=sys.stderr)
        text = efetch_genbank(acc)
        if acc == M1_REFERENCE:
            m1_seq = genbank_sequence(text)
        for gene, s in genbank_cds(text):
            gene = canonical.get(gene.lower())
            if gene and s not in seen.setdefault(gene, set()):
                seen[gene].add(s)
                records.append((gene, acc, s))
    print(f"  NCBI {SSA_RECORD}", file=sys.stderr)
    ssa = [s for g, s in genbank_cds(efetch_genbank(SSA_RECORD)) if g.lower() == "ssa"]
    if not ssa:
        raise RuntimeError(f"{SSA_RECORD}: no CDS named ssa")
    records.append(("ssa", SSA_RECORD, ssa[0]))
    missing = [g for g in SAG_GENES if g not in seen]
    if missing:
        raise RuntimeError(f"superantigen genes not found in the reference genomes: {missing}")
    with open(os.path.join(refdir, "superantigens.fa"), "w") as fh:
        for gene, acc, s in sorted(records, key=lambda r: (r[0], SAG_GENOMES.index(r[1])
                                                           if r[1] in SAG_GENOMES else 99)):
            fh.write(f">{gene}_{acc}\n")
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + "\n")

    # M1 lineage SNPs: a flank of FLANK bases each side, so the typing step can
    # find each site in an assembly by BLAST and read the base in the middle.
    with open(os.path.join(refdir, "m1_snps.tsv"), "w") as fh:
        fh.write(f"# M1 lineage SNPs on {M1_REFERENCE} (MGAS5005); flanks are "
                 f"{FLANK} bp each side of the site\n")
        fh.write("lineage\tposition\tref\talt\tleft\tright\n")
        for lineage, block in M1_SNPS.items():
            for item in block.replace("\n", "").split(";"):
                pos, ref, alt = item.split()
                p = int(pos)
                if m1_seq[p - 1] != ref:
                    raise RuntimeError(f"{lineage} {pos}: reference base is "
                                       f"{m1_seq[p - 1]}, expected {ref}")
                fh.write(f"{lineage}\t{p}\t{ref}\t{alt}\t"
                         f"{m1_seq[p - 1 - FLANK:p - 1]}\t{m1_seq[p:p + FLANK]}\n")
    return {g: len(v) for g, v in seen.items()} | {"ssa": 1}


# --------------------------------------------------------------------------
# BLAST databases
# --------------------------------------------------------------------------

def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def build_blastdbs(pubmlst_dir, ref_dir, blast_dir):
    """Two BLAST indexes: loci whose longest allele is SHORT bp or more (plus
    the superantigen genes) and the short signature loci on their own, which
    need other blastn settings. Titles: spy__<scheme>__<locus>__<allele>."""
    os.makedirs(blast_dir, exist_ok=True)
    fa = {BLAST_LONG: os.path.join(blast_dir, BLAST_LONG + ".fa"),
          BLAST_SHORT: os.path.join(blast_dir, BLAST_SHORT + ".fa")}
    out = {k: open(v, "w") for k, v in fa.items()}
    try:
        for scheme in sorted(os.listdir(pubmlst_dir)):
            sdir = os.path.join(pubmlst_dir, scheme)
            if not os.path.isdir(sdir):
                continue
            for tfa in sorted(f for f in os.listdir(sdir) if f.endswith(".tfa")):
                locus = tfa[:-4]
                recs = list(read_fasta(os.path.join(sdir, tfa)))
                index = BLAST_SHORT if max(len(s) for _, s in recs) < SHORT else BLAST_LONG
                for name, s in recs:
                    allele = name[len(locus) + 1:]
                    out[index].write(f">spy{SEP}{scheme}{SEP}{locus}{SEP}{allele}\n{s}\n")
        sag = os.path.join(ref_dir, "superantigens.fa")
        for name, s in read_fasta(sag):
            gene, acc = name.split("_", 1)
            out[BLAST_LONG].write(f">spy{SEP}{SAG_SCHEME}{SEP}{gene}{SEP}{acc}\n{s}\n")
    finally:
        for fh in out.values():
            fh.close()
    for name, path in fa.items():
        subprocess.run(["makeblastdb", "-in", path, "-dbtype", "nucl",
                        "-title", name, "-out", os.path.join(blast_dir, name)],
                       check=True, stdout=subprocess.DEVNULL)
        os.remove(path)  # the index holds the sequences; keep the repo smaller


# --------------------------------------------------------------------------
# Comparison
# --------------------------------------------------------------------------

def read_manifest(path):
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return {r["scheme"]: r for r in csv.DictReader(fh, delimiter="\t")}


def folder_digest(d):
    """Content of every file under folder d, for change detection."""
    out = {}
    for root, _, files in os.walk(d):
        for f in sorted(files):
            p = os.path.join(root, f)
            with open(p) as fh:
                out[os.path.relpath(p, d)] = fh.read()
    return out


def compare(old_pub, new_pub):
    old = set(os.listdir(old_pub)) if os.path.isdir(old_pub) else set()
    new = set(os.listdir(new_pub))
    changed = sorted(s for s in old & new
                     if folder_digest(os.path.join(old_pub, s))
                     != folder_digest(os.path.join(new_pub, s)))
    return sorted(new - old), sorted(old - new), changed


def write_summary(path, version, man, added, removed, changed, failed, refs):
    lines = [f"# Database update {version}", ""]
    if any(m["restricted"] for m in man):
        lines += ["PubMLST data is the **standard, freely available** set: alleles "
                  "and profiles submitted up to **31 Dec 2024** (no authentication).", ""]
    lines += [f"- Schemes: {len(man)} "
              f"({sum(m['loci'] for m in man)} loci, "
              f"{sum(m['alleles'] for m in man)} alleles, "
              f"{sum(m['profiles'] for m in man)} STs)",
              f"- Added: {', '.join(added) or 'none'}",
              f"- Removed: {', '.join(removed) or 'none'}",
              f"- Changed: {len(changed)}" + (f" ({', '.join(changed)})" if changed else ""),
              f"- References (superantigens, M1 SNPs): {refs}"]
    if failed:
        lines += ["", "**Failed (kept the previous copy):**", ""]
        lines += [f"- {n}: {e}" for n, e in failed]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def old_row(old_man, name):
    """Manifest row of the previous database, typed like a fresh one."""
    m = dict(old_man[name])
    for k in ("loci", "alleles", "profiles"):
        m[k] = int(m[k])
    m["restricted"] = m.get("restricted") == "True"
    return m


def read_version(db):
    try:
        with open(os.path.join(db, "VERSION.txt")) as fh:
            return fh.readline().strip()
    except OSError:
        return ""


# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(prog="update_database.py",
                                 description=__doc__.split("\n\n")[0])
    ap.add_argument("--db", default="database", help="database folder [database]")
    ap.add_argument("--schemes", help="scheme list [<db>/schemes.tsv]")
    ap.add_argument("--only", help="comma-separated scheme names (for testing)")
    ap.add_argument("--references", action="store_true",
                    help="rebuild references/ from NCBI (done anyway when it is missing)")
    ap.add_argument("--version-label", help="VERSION.txt line [PubMLST_<YYYY-MM-DD>]")
    ap.add_argument("--summary", help="write a Markdown change summary here")
    args = ap.parse_args(argv)

    db = os.path.abspath(args.db)
    rows = read_scheme_list(args.schemes or os.path.join(db, "schemes.tsv"))
    want = set(args.only.split(",")) if args.only else None
    if want and want - {r["scheme"] for r in rows}:
        sys.exit(f"--only: unknown schemes {sorted(want - {r['scheme'] for r in rows})}")
    fetch = [r for r in rows if want is None or r["scheme"] in want]
    if shutil.which("makeblastdb") is None:
        sys.exit("makeblastdb not found (install BLAST+)")

    old_man = read_manifest(os.path.join(db, "manifest.tsv"))
    old_pub, old_ref = os.path.join(db, "pubmlst"), os.path.join(db, "references")
    tmp = tempfile.mkdtemp(prefix="spydb_", dir=db)
    new_pub, new_ref = os.path.join(tmp, "pubmlst"), os.path.join(tmp, "references")
    os.makedirs(new_pub)
    man, failed = [], []
    try:
        n = len(fetch) + 1
        for i, row in enumerate(fetch, 1):
            name, old = row["scheme"], os.path.join(old_pub, row["scheme"])
            print(f"[{i}/{n}] {name}", file=sys.stderr)
            try:
                man.append(fetch_scheme(row, new_pub))
            except Exception as e:  # noqa: BLE001
                failed.append((name, str(e)[:200]))
                print(f"ERROR: {name}: {e}", file=sys.stderr)
                shutil.rmtree(os.path.join(new_pub, name), ignore_errors=True)
                if os.path.isdir(old):
                    shutil.copytree(old, os.path.join(new_pub, name))
                    if name in old_man:
                        man.append(old_row(old_man, name))
        # --only: carry over every listed scheme that was not fetched
        listed = {r["scheme"] for r in rows}
        if want is not None and os.path.isdir(old_pub):
            done = set(os.listdir(new_pub))
            for s in sorted(set(os.listdir(old_pub)) & listed - done):
                shutil.copytree(os.path.join(old_pub, s), os.path.join(new_pub, s))
                if s in old_man:
                    man.append(old_row(old_man, s))
        if not os.listdir(new_pub):
            raise RuntimeError("no scheme could be downloaded")

        print(f"[{n}/{n}] references", file=sys.stderr)
        refs = "kept"
        if args.references or not os.path.isdir(old_ref):
            counts = build_references(new_ref)
            refs = "rebuilt (" + ", ".join(f"{g} {c}" for g, c in sorted(counts.items())) + ")"
        else:
            shutil.copytree(old_ref, new_ref)

        added, removed, changed = compare(old_pub, new_pub)
        ref_changed = (not os.path.isdir(old_ref)
                       or folder_digest(old_ref) != folder_digest(new_ref))
        if not (added or removed or changed or ref_changed) \
                and os.path.isdir(os.path.join(db, "blast")):
            print("No changes.", file=sys.stderr)
            return NO_CHANGES

        build_blastdbs(new_pub, new_ref, os.path.join(tmp, "blast"))
        pub_changed = bool(added or removed or changed) or not read_version(db)
        version = (args.version_label or f"PubMLST_{dt.date.today().isoformat()}") \
            if pub_changed else read_version(db)

        man.sort(key=lambda m: m["scheme"])
        with open(os.path.join(tmp, "manifest.tsv"), "w", newline="") as fh:
            cols = ["scheme", "source", "loci", "alleles", "profiles",
                    "last_updated", "restricted", "description", "api"]
            w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(man)
        with open(os.path.join(tmp, "VERSION.txt"), "w") as fh:
            fh.write(version + "\n")

        # swap in (schemes.tsv and tables/ are never touched)
        for part in ("pubmlst", "references", "blast", "manifest.tsv", "VERSION.txt"):
            dst = os.path.join(db, part)
            if os.path.isdir(dst):
                shutil.rmtree(dst)
            elif os.path.exists(dst):
                os.remove(dst)
            shutil.move(os.path.join(tmp, part), dst)
        if args.summary:
            write_summary(args.summary, version, man, added, removed, changed, failed, refs)
        print(f"Database {version}: {len(man)} schemes, "
              f"+{len(added)} -{len(removed)} ~{len(changed)}, references {refs.split(' ')[0]}, "
              f"{len(failed)} failed", file=sys.stderr)
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
