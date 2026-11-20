"""Tests for spyogenes-typing.

Unit tests need nothing. Run tests need blastn and the bundled database/.
tests/data/*_loci.fasta.gz hold the typed loci (and the M1 SNP sites) cut
out of three complete genomes from NCBI: MGAS5005 (NC_007297.2),
MGAS6180 (NC_007296.2) and Alab49 (CP003068.1). MGAS6180 and Alab49 are
the worked examples in PubMLST's typing guide (Feb 2024).
"""

import csv
import gzip
import io
import json
import os
import shutil
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")
DB = os.path.join(ROOT, "database")
sys.path.insert(0, os.path.join(ROOT, "tool"))

import spyogenes_typing as st  # noqa: E402

TOOL = [sys.executable, os.path.join(ROOT, "tool", "spyogenes_typing.py"), "--db", DB]
HAVE_BLAST = shutil.which("blastn") is not None
HAVE_DB = os.path.isdir(os.path.join(DB, "blast"))
needs_blast = pytest.mark.skipif(not (HAVE_BLAST and HAVE_DB),
                                 reason="needs blastn and database/")


def genome(name):
    return os.path.join(DATA, f"{name}_loci.fasta.gz")


def run(*args, ok=True):
    r = subprocess.run([*TOOL, "-q", *args], capture_output=True, text=True)
    if ok:
        assert r.returncode == 0, r.stderr
    return r


def rows(text):
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def hit(locus="emm", allele="1.0", slen=180, alen=180, nident=180, start=101,
        end=280, s1=1, s2=180, qlen=10000, strand="plus", scheme="emm", contig="c1"):
    return st.Hit([f"spy__{scheme}__{locus}__{allele}", slen, alen, nident, contig,
                   start, end, s1, s2, qlen, strand])


# --------------------------------------------------------------------------
# Unit tests
# --------------------------------------------------------------------------

def test_parse_fasta_genbank_and_names():
    assert st.parse_sequences(">a x\nacgt\nAC\n>b\nGG\n") == [("a", "ACGTAC"), ("b", "GG")]
    gb = "LOCUS       X1   6 bp\nORIGIN\n        1 acgtac\n//\n"
    assert st.parse_sequences(gb) == [("X1", "ACGTAC")]
    with pytest.raises(st.TypingError):
        st.parse_sequences("hello\n")
    assert st.sample_name("/x/ERR123.contigs.fa.gz") == "ERR123.contigs"
    assert st.sample_name("s1.fasta") == "s1"


def test_compressed_input(tmp_path):
    p = tmp_path / "g.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">c\nACGT\n")
    assert st.parse_sequences(st.read_bytes(str(p))) == [("c", "ACGT")]


def test_hit_exact_full_and_edge():
    h = hit()
    assert h.exact and h.full and not h.edge
    # first 80 bp of the allele missing, alignment starts at the contig start
    h = hit(alen=100, nident=100, start=1, end=100, s1=81, s2=180)
    assert not h.full and h.edge
    # same on the minus strand: the missing 5' end lies past the contig end
    h = hit(alen=100, nident=100, start=9901, end=10000, s1=180, s2=81, strand="minus")
    assert h.edge
    # minus strand, 3' end missing: it would lie left of the hit, inside the contig
    h = hit(alen=100, nident=100, start=9901, end=10000, s1=100, s2=1, strand="minus")
    assert not h.edge
    # alignment stops inside the contig: a different gene, not an edge
    h = hit(alen=100, nident=100, start=500, end=599, s1=81, s2=180)
    assert not h.edge


def test_call_words():
    assert st.call(hit(), 95) == st.PERFECT
    assert st.call(hit(nident=176), 95) == st.NOVEL
    assert st.call(hit(nident=160), 95) is None
    assert st.call(hit(alen=100, nident=100, start=1, end=100, s1=81, s2=180), 95) == st.PARTIAL


def test_emm_labels():
    assert st.emm_label(hit(allele="1.0")) == "emm1.0"
    assert st.emm_label(hit(locus="emmL", allele="134.2")) == "emm134.2"
    assert st.emm_label(hit(locus="emmCG", allele="STG653.0")) == "STG653.0"
    assert st.emm_type_of(hit(locus="emmCG", allele="STG653.0")) == "STG653"


def test_near_prime_strand_and_distance():
    emm = st.Region(hit(start=1000, end=1179, strand="minus"))
    prime = hit(locus="emmPrime", allele="2", slen=19, alen=19, nident=19,
                start=1227, end=1245, s1=19, s2=1, strand="minus")
    assert st.near_prime(emm, [prime])
    far = hit(locus="emmPrime", allele="2", slen=19, alen=19, nident=19,
              start=1600, end=1618, s1=19, s2=1, strand="minus")
    assert not st.near_prime(emm, [far])
    wrong_side = hit(locus="emmPrime", allele="2", slen=19, alen=19, nident=19,
                     start=900, end=918, s1=19, s2=1, strand="minus")
    assert not st.near_prime(emm, [wrong_side])


def test_pattern_rules_from_the_guide():
    rules = dict(st.PATTERNS)
    assert rules["A"] == {"mga-1": 1, "SF1": 1}
    assert rules["B"] == {"mga-1": 1, "SF1": 2}
    assert rules["C"] == {"mga-1": 1, "SF1": 1, "SF3": 1}
    assert rules["D"] == {"mga-2": 1, "SF1": 1, "SF3": 1, "SF4": 1}
    assert rules["E"] == {"mga-2": 1, "SF2": 1, "SF3": 1, "SF4": 1}


def test_mlst_status_words():
    assert st.mlst_status("28", ["4", "3", "4", "4", "4", "2", "4"]) == "PERFECT"
    assert st.mlst_status("28", ["4,4", "3", "4", "4", "4", "2", "4"]) == "MULTICOPY"
    assert st.mlst_status("-", ["4", "3", "4", "4", "4", "2", "999"]) == "NOVEL"
    assert st.mlst_status("-", ["4,5", "3", "4", "4", "4", "2", "4"]) == "MIXED"
    assert st.mlst_status("-", ["-", "3", "4", "4", "4", "2", "4"]) == "MISSING"
    assert st.mlst_status("-", ["~4", "3", "4", "4", "4", "2", "4"]) == "PARTIAL"
    assert st.mlst_status("-", ["-"] * 7) == "NONE"


def test_base_at_with_neighbouring_snp():
    left = list("ACGTTGCAAGTCCATGACTGATCGGATCCA")
    right = ["[CT]"] + list("GGATTACAGCTTGACCATGGTACCGATTA")
    g = "".join(left) + "G" + "T" + "".join(right[1:])
    assert st.base_at(g, st.revcom(g), left, right) == "G"
    assert st.base_at(st.revcom(g), g, left, right) == "G"


@needs_blast
def test_database_tables_load():
    db = st.Database(DB)
    assert db.mlst_loci == ["gki", "gtr", "murI", "mutS", "recP", "xpt", "yqiL"]
    assert db.clusters["1"] == "A-C3"
    assert db.patterns["18"] == ["A-C", "D"]
    assert db.pilin_types[("nt80_13", "nt80_34")][0] == "pil028"
    assert len(db.m1_snps) == 42


# --------------------------------------------------------------------------
# Run tests
# --------------------------------------------------------------------------

EXPECTED = {
    "MGAS5005": {"emm": "emm1.0", "emm status": "PERFECT", "emm cluster": "A-C3",
                 "emm pattern": "A", "ST": "28", "MLST status": "PERFECT",
                 "M1 lineage": "M1global (0 M1UK SNPs)",
                 "Superantigens": "smeZ;speA;speG;speJ", "FCT type": "FCT-2 (pil001)"},
    "MGAS6180": {"emm": "emm28.0", "emm status": "PERFECT", "emm cluster": "E4",
                 "emm pattern": "E", "ST": "52", "M1 lineage": "-",
                 "Superantigens": "smeZ;speC;speG;speJ;speK", "FCT type": "FCT-4 (pil028)"},
    "Alab49": {"emm": "emm53.0", "emm status": "PERFECT", "emmL": "emmL134.2;emmL156.1 (98.9%)",
               "emm cluster": "D4", "emm pattern": "D", "ST": "11",
               "Superantigens": "smeZ;speC;speG;speL;speM", "FCT type": "FCT-3 (pil033)"},
}


@needs_blast
@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_reference_genomes(name):
    r = rows(run(genome(name)).stdout)
    assert len(r) == 1
    row = r[0]
    assert row["Sample"] == f"{name}_loci"
    assert row["Species check"] == "S. pyogenes"
    for k, v in EXPECTED[name].items():
        assert row[k] == v, k
    assert row["Database"] == st.db_version(DB)


@needs_blast
def test_allele_layouts_and_files(tmp_path):
    g = genome("MGAS5005")
    default = rows(run(g).stdout)[0]
    assert [default[f"ALLELE{i}"] for i in range(1, 8)] == \
        ["gki(4)", "gtr(3)", "murI(4)", "mutS(4)", "recP(4)", "xpt(2)", "yqiL(4)"]
    comb = rows(run("--combine-alleles", g).stdout)[0]
    assert comb["ALLELES"] == "gki(4);gtr(3);murI(4);mutS(4);recP(4);xpt(2);yqiL(4)"
    none = run("--no-alleles", "--label", "x", g).stdout
    assert "ALLELE" not in none.splitlines()[0] and rows(none)[0]["Sample"] == "x"
    bad = run("--no-alleles", "--combine-alleles", g, ok=False)
    assert bad.returncode != 0
    js, hits = tmp_path / "o.json", tmp_path / "h.tsv"
    run("--json", str(js), "--hits", str(hits), g, genome("Alab49"))
    data = json.load(open(js))
    assert [d["ST"] for d in data] == ["28", "11"]
    assert data[0]["alleles"]["gki"] == "4"
    h = rows(open(hits).read())
    assert any(x["Locus"] == "emmPrime" and x["Call"] == "PERFECT" for x in h)


@needs_blast
def test_m1uk_and_split_emm(tmp_path):
    """M1UK SNPs written into MGAS5005 give M1UK; emm cut in two by a
    contig break gives PARTIAL with the type only, and no stray emm-like hit."""
    seqs = st.parse_sequences(st.read_bytes(genome("MGAS5005")))
    db = st.Database(DB)
    out = []
    for name, s in seqs:
        a = int(name.split("_")[-2])
        s = list(s)
        for x in db.m1_snps:
            p = int(x["position"]) - a
            if x["lineage"] == "M1UK" and 0 <= p < len(s):
                assert s[p] == x["ref"]
                s[p] = x["alt"]
        out.append((name, "".join(s)))
    uk = tmp_path / "uk.fa"
    uk.write_text("".join(f">{n}\n{s}\n" for n, s in out))
    assert rows(run(str(uk)).stdout)[0]["M1 lineage"] == "M1UK (27/27 SNPs)"
    # split the contig holding emm1 (NC_007297.2 1683318-1683497, minus strand)
    split = []
    for name, s in seqs:
        a = int(name.split("_")[-2])
        cut = 1683400 - a
        if 0 < cut < len(s):
            split += [(name + "_a", s[:cut]), (name + "_b", s[cut:])]
        else:
            split.append((name, s))
    sp = tmp_path / "split.fa"
    sp.write_text("".join(f">{n}\n{s}\n" for n, s in split))
    row = rows(run(str(sp)).stdout)[0]
    assert (row["emm"], row["emm status"], row["emmL"]) == ("emm1", "PARTIAL", "-")
    assert row["emm cluster"] == "A-C3"


@needs_blast
def test_not_spyogenes(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(">x\n" + "ACGT" * 5000 + "\n")
    row = rows(run(str(p)).stdout)[0]
    assert row["Species check"].startswith("not S. pyogenes")
    assert row["emm status"] == "MISSING" and row["MLST status"] == "NONE"
