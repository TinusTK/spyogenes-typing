#!/usr/bin/env python3
"""Fetch the validation genomes and write validation/genomes/MANIFEST.tsv
and validation/expected.tsv.

    python scripts/fetch_validation_genomes.py [--out validation]

Two sources, both public:

  * **Complete genomes from NCBI** (nuccore, FASTA). PubMLST's own isolate
    records hold the expected emm type and FCT-region form for these named
    strains, but almost none of them carry a sequence accession, so each
    strain is matched to its PubMLST isolate record by name. The isolate id
    that was matched is written into expected.tsv, so every expected value
    can be traced back.
  * **Draft assemblies from PubMLST** (`/isolates/<id>/contigs_fasta`), for
    the isolates that have contigs. Their own record gives the expected
    values directly.

Nothing here is needed to run or test the tool; it only rebuilds the
validation set. The genomes are committed, so this script rarely runs.
"""

import argparse
import gzip
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ISOLATES = "https://rest.pubmlst.org/db/pubmlst_spyogenes_isolates"
UA = "spyogenes-typing-validation"

# Complete NCBI genomes, with the PubMLST isolate name they are recorded under.
# "pubmlst" is the isolate name to look up in the PubMLST isolate database.
NCBI = [
    ("MGAS5005", "NC_007297.2", "MGAS5005"),
    ("MGAS6180", "NC_007296.2", "MGAS6180"),
    ("Alab49", "CP003068.1", "ALAB49"),
    ("SF370", "NC_002737.2", "SF370"),
    ("MGAS10394", "NC_006086.1", "MGAS10394"),
    ("Manfredo", "NC_009332.1", "Manfredo"),
    ("SSI-1", "NC_004606.1", "SSI-1"),
    ("MGAS8232", "NC_003485.1", "MGAS8232"),
    ("MGAS315", "NC_004070.1", "MGAS315"),
    ("MGAS9429", "NC_008021.1", "MGAS9429"),
    ("MGAS15252", "NC_017040.1", "MGAS15252"),
    ("MGAS1882", "NC_017053.1", "MGAS1882"),
    ("HSC5", "NC_021807.1", "HSC5"),
    ("MGAS10750", "NC_008024.1", "MGAS10750"),
    ("MGAS10270", "NC_008022.1", "MGAS10270"),
]

# PubMLST isolates that have contigs (found by trying the isolates that carry
# an emm type or an FCT-region form; most older records have no sequence data).
PUBMLST_IDS = [4078, 4115, 6424, 6537, 6541, 6549, 7422]

FIELDS = ("isolate", "emm_type", "emm_subtype", "emm_cluster",
          "emm_pattern_measured", "emm_pattern_predicted", "FCT_region",
          "pilin_type")


def get(url, tries=4, data=None):
    """Bytes of url. With data= the request is a JSON POST (PubMLST's search
    endpoint only answers POST)."""
    body = json.dumps(data).encode() if data is not None else None
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, data=body,
                headers={"User-Agent": UA,
                         **({"Content-Type": "application/json"} if body else {})})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(3 + 3 * i)
        except Exception:
            time.sleep(3 + 3 * i)
    return None


def write_gz(path, data):
    with gzip.open(path, "wb", compresslevel=9) as fh:
        fh.write(data)


def fasta_stats(data):
    text = data.decode("utf-8", "replace")
    seqs = [l for l in text.splitlines() if l.startswith(">")]
    length = sum(len(l.strip()) for l in text.splitlines() if not l.startswith(">"))
    first = seqs[0][1:].strip() if seqs else ""
    return len(seqs), length, first


def record(isolate_id):
    raw = get(f"{ISOLATES}/isolates/{isolate_id}")
    if raw is None:
        return None
    rec = json.loads(raw)
    prov = rec.get("provenance", {})
    out = {k: prov.get(k) for k in FIELDS}
    out["ST"] = next((s.get("fields", {}).get("ST")
                      for s in rec.get("schemes", [])
                      if s.get("description") == "MLST"), None)
    out["id"] = isolate_id
    out["has_contigs"] = bool(rec.get("sequence_bin"))
    return out


def find_by_name(name):
    """The PubMLST isolate record whose isolate name matches (case-insensitive)."""
    raw = get(f"{ISOLATES}/isolates/search?page_size=100",
              data={"field.isolate": name})
    if raw is None:
        return None
    ids = [int(u.rsplit("/", 1)[1]) for u in json.loads(raw).get("isolates", [])]
    # Several records can share a strain name; prefer one that carries the
    # fields the validation compares, so a bare duplicate does not win.
    best = None
    for i in ids:
        rec = record(i)
        if not rec or (rec.get("isolate") or "").strip().lower() != name.lower():
            continue
        filled = sum(1 for k in ("emm_type", "FCT_region", "pilin_type", "ST")
                     if rec.get(k))
        if best is None or filled > best[0]:
            best = (filled, rec)
    return best[1] if best else None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="validation")
    ap.add_argument("--skip-download", action="store_true",
                    help="keep the genomes already in place, refresh the tables only")
    a = ap.parse_args(argv)
    gdir = os.path.join(a.out, "genomes")
    os.makedirs(gdir, exist_ok=True)

    manifest, expected = [], []

    for sample, acc, pubmlst_name in NCBI:
        path = os.path.join(gdir, f"{sample}.fasta.gz")
        if not (a.skip_download and os.path.exists(path)):
            print(f"NCBI {acc} ({sample})", file=sys.stderr)
            data = get(f"{EUTILS}?db=nuccore&id={acc}&rettype=fasta&retmode=text")
            if not data or not data.startswith(b">"):
                sys.exit(f"could not download {acc}")
            write_gz(path, data)
            time.sleep(0.4)
        with gzip.open(path, "rb") as fh:
            data = fh.read()
        n, length, desc = fasta_stats(data)
        manifest.append(dict(file=f"{sample}.fasta.gz", source=f"NCBI {acc}",
                             sequences=n, length_bp=length,
                             sha256=hashlib.sha256(data).hexdigest(),
                             description=desc))
        rec = find_by_name(pubmlst_name)
        if rec is None:
            print(f"  WARNING: no PubMLST isolate named {pubmlst_name}", file=sys.stderr)
            rec = {}
        expected.append(dict(
            sample=sample, genome=f"{sample}.fasta.gz", accession=acc,
            pubmlst_isolate_id=rec.get("id", ""),
            pubmlst_isolate_name=rec.get("isolate", ""),
            emm_type=rec.get("emm_type") or "", emm_subtype=rec.get("emm_subtype") or "",
            emm_cluster=rec.get("emm_cluster") or "", ST=rec.get("ST") or "",
            FCT_region=rec.get("FCT_region") or "",
            pilin_type=rec.get("pilin_type") or "",
            emm_pattern_predicted=rec.get("emm_pattern_predicted") or "",
            source="complete genome from NCBI; expected values from the PubMLST "
                   "isolate record of the same strain"))

    for i in PUBMLST_IDS:
        rec = record(i)
        if not rec or not rec["has_contigs"]:
            print(f"  WARNING: isolate {i} has no contigs, skipped", file=sys.stderr)
            continue
        sample = f"i{i}"
        path = os.path.join(gdir, f"{sample}.fasta.gz")
        if not (a.skip_download and os.path.exists(path)):
            print(f"PubMLST isolate {i} ({rec.get('isolate')})", file=sys.stderr)
            data = get(f"{ISOLATES}/isolates/{i}/contigs_fasta")
            if not data or not data.startswith(b">"):
                sys.exit(f"could not download contigs of isolate {i}")
            write_gz(path, data)
            time.sleep(0.4)
        with gzip.open(path, "rb") as fh:
            data = fh.read()
        n, length, desc = fasta_stats(data)
        manifest.append(dict(file=f"{sample}.fasta.gz",
                             source=f"PubMLST isolate {i}",
                             sequences=n, length_bp=length,
                             sha256=hashlib.sha256(data).hexdigest(),
                             description=f"{rec.get('isolate')} (draft assembly, "
                                         f"{n} contigs)"))
        expected.append(dict(
            sample=sample, genome=f"{sample}.fasta.gz", accession="",
            pubmlst_isolate_id=i, pubmlst_isolate_name=rec.get("isolate") or "",
            emm_type=rec.get("emm_type") or "", emm_subtype=rec.get("emm_subtype") or "",
            emm_cluster=rec.get("emm_cluster") or "", ST=rec.get("ST") or "",
            FCT_region=rec.get("FCT_region") or "",
            pilin_type=rec.get("pilin_type") or "",
            emm_pattern_predicted=rec.get("emm_pattern_predicted") or "",
            source="draft assembly and expected values from the same PubMLST "
                   "isolate record"))

    def table(path, rows, cols):
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    table(os.path.join(gdir, "MANIFEST.tsv"), manifest,
          ["file", "source", "sequences", "length_bp", "sha256", "description"])
    table(os.path.join(a.out, "expected.tsv"), expected,
          ["sample", "genome", "accession", "pubmlst_isolate_id",
           "pubmlst_isolate_name", "emm_type", "emm_subtype", "emm_cluster",
           "ST", "emm_pattern_predicted", "FCT_region", "pilin_type", "source"])
    print(f"{len(manifest)} genomes, {len(expected)} expected rows", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
