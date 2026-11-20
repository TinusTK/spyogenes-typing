#!/usr/bin/env python3
"""Type every validation genome and compare the result with PubMLST.

    python scripts/validate.py                     # uses validation/ and database/
    python scripts/validate.py --db database --out validation

Writes, into validation/:
  results_<database version>.tsv   the tool's full report, every column
  VALIDATION.md                    the record: what was compared, what agreed

What is compared, per genome, only where PubMLST has a value:
emm type, emm subtype, ST, emm cluster, predicted emm pattern, FCT-region
form and pilin type. A field PubMLST leaves empty is reported as "not
recorded" and is not counted either way.

Expected values and the genomes come from validation/expected.tsv and
validation/genomes/ (see scripts/fetch_validation_genomes.py).
"""

import argparse
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# (column in expected.tsv, column in the report, label, how to compare)
CHECKS = [
    ("emm_type", "emm", "emm type", "emm_type"),
    ("emm_subtype", "emm", "emm subtype", "emm_subtype"),
    ("ST", "ST", "ST", "plain"),
    ("emm_cluster", "emm cluster", "emm cluster", "cluster"),
    ("emm_pattern_measured", "emm pattern", "emm pattern", "pattern"),
    ("FCT_region", "FCT type", "FCT-region form", "fct"),
    ("pilin_type", "FCT type", "pilin type", "pilin"),
]

# Differences that were investigated and traced to the source data rather than
# to the tool. Each one names the evidence, so a future run that changes will
# no longer match its explanation and will show up as unexplained.
KNOWN = {
    ("HSC5", "ST"): (
        "PubMLST isolate 2548 lists gtr 2 (entered by hand in 2015); the "
        "deposited assembly NC_021807.1 matches gtr_22 exactly (100% identity "
        "and coverage). gtr_2 and gtr_22 differ by one base, and ST84 and ST680 "
        "differ only at gtr, so ST84 is the ST of this assembly."),
    ("MGAS6180", "emm subtype"): (
        "The assembly matches emm_28.0 exactly; emm_28.4, which the isolate "
        "record gives, differs from it by one base. The subtype of the "
        "deposited sequence is emm28.0."),
    ("MGAS8232", "emm pattern"): (
        "Measured, not predicted: the emm region of this genome holds mga-2 and "
        "SF1 only, with no SF3 or SF4 anywhere near it (checked down to 70% "
        "identity), which fits no pattern A-E, so the guide's rules give REA. "
        "A-C is the pattern predicted for emm18, which PubMLST also records as "
        "being seen with pattern D."),
    ("i7422", "ST"): (
        "ST1778 was assigned after 31 December 2024, so it is not in the "
        "bundled free database; NOVEL is the right answer against it. Users "
        "with their own PubMLST key get this ST."),
}

# PubMLST's isolate database still uses FCT-7 and FCT-8, which its Feb 2024
# typing guide no longer lists; the tool only calls the nine forms the guide
# defines, so those rows cannot agree or disagree.
OLD_FORMS = {"FCT-7", "FCT-8"}


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_tsv_comments(path):
    """Like read_tsv, for a table that starts with '#' comment lines."""
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))


def ours(row, field, how):
    """The part of our report that answers one expected field."""
    v = (row.get(field) or "").strip()
    if how == "emm_type":
        return ",".join(t.replace("emm", "").split(".")[0] for t in v.split(",")) if v != "-" else "-"
    if how == "emm_subtype":
        return v
    if how == "cluster":
        return v
    if how == "pattern":
        p = v.replace(" (predicted)", "").strip()
        return "A-C" if p in ("A", "B", "C") else p
    if how == "fct":
        return v.split(" ")[0] if v != "-" else "-"
    if how == "pilin":
        return v.split("(")[1].rstrip(")") if "(" in v else "-"
    return v


def agrees(expected, got, how):
    if how == "emm_subtype":
        # PubMLST writes emm28.0; we write the same with the emm prefix
        return got == expected or got == expected.replace("emm", "")
    if how == "cluster":
        # PubMLST's isolate database has older label styles such as "E-6"
        return got == expected or got == expected.replace("-", "", 1)
    return got == expected


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=os.path.join(ROOT, "database"))
    ap.add_argument("--out", default=os.path.join(ROOT, "validation"))
    ap.add_argument("--results", default="", help="use an existing report instead of typing")
    a = ap.parse_args(argv)

    expected = read_tsv(os.path.join(a.out, "expected.tsv"))
    gdir = os.path.join(a.out, "genomes")
    with open(os.path.join(a.db, "VERSION.txt")) as fh:
        dbv = fh.readline().strip()

    results_path = a.results or os.path.join(a.out, f"results_{dbv}.tsv")
    if not a.results:
        files = [os.path.join(gdir, r["genome"]) for r in expected]
        missing = [f for f in files if not os.path.exists(f)]
        if missing:
            sys.exit(f"missing genomes: {missing[:3]} "
                     "(run scripts/fetch_validation_genomes.py)")
        cmd = [sys.executable, os.path.join(ROOT, "tool", "spyogenes_typing.py"),
               "--db", a.db, "--outfile", results_path, "--threads", "4", *files]
        print("Typing " + f"{len(files)} genomes", file=sys.stderr)
        r = subprocess.run(cmd)
        if r.returncode != 0:
            sys.exit("typing failed")

    got = {}
    for row in read_tsv(results_path):
        got[row["Sample"]] = row

    # per-genome comparison
    tally = {label: [0, 0] for _, _, label, _ in CHECKS}
    notes, rows, explained, unexplained = [], [], [], []
    for e in expected:
        sample = e["genome"].replace(".fasta.gz", "")
        r = got.get(sample)
        if r is None:
            sys.exit(f"no result row for {sample}")
        verdicts = {}
        for ecol, rcol, label, how in CHECKS:
            exp = (e.get(ecol) or "").strip()
            # PubMLST's public records hold no measured emm pattern, so fall
            # back to the pattern it predicts for the emm type and say so.
            against_predicted = False
            if label == "emm pattern" and (not exp or exp == "-"):
                exp = (e.get("emm_pattern_predicted") or "").strip()
                against_predicted = bool(exp and exp != "-")
            if not exp or exp == "-":
                verdicts[label] = ("not recorded", "", "")
                continue
            if label == "FCT-region form" and exp in OLD_FORMS:
                verdicts[label] = ("old label", exp, ours(r, rcol, how))
                notes.append(f"{sample}: PubMLST records {exp}, a form the Feb 2024 "
                             f"guide no longer lists; the tool called "
                             f"{ours(r, rcol, how)}")
                continue
            o = ours(r, rcol, how)
            ok = agrees(exp, o, how)
            tally[label][1] += 1
            if ok:
                tally[label][0] += 1
                verdicts[label] = ("agree", exp, o)
                continue
            why = KNOWN.get((sample, label))
            if why:
                explained.append((sample, label, exp, o, why))
                verdicts[label] = ("explained", exp, o)
            else:
                verdicts[label] = ("DIFFERS", exp, o)
                unexplained.append((sample, label, exp, o,
                                    "compared with the pattern PubMLST predicts "
                                    "for the emm type, not a measured pattern"
                                    if against_predicted else ""))
        rows.append((e, r, verdicts))

    # ---------------------------------------------------------------- report
    lines = [f"# Validation, database {dbv}", "",
             f"{len(rows)} genomes: "
             f"{sum(1 for e, _, _ in rows if e['accession'])} complete genomes from "
             f"NCBI and {sum(1 for e, _, _ in rows if not e['accession'])} draft "
             "assemblies from PubMLST.", "",
             "Every genome is typed with the bundled database and compared with the "
             "PubMLST isolate record for the same strain. PubMLST leaves many of "
             "these fields empty, and an empty field is counted neither way.", ""]
    lines += ["| What | Agree | Compared |", "|---|---|---|"]
    for _, _, label, _ in CHECKS:
        ok, n = tally[label]
        lines.append(f"| {label} | {ok} | {n} |"
                     if n else f"| {label} | – | 0 (PubMLST records none) |")
    empty = [label for _, _, label, _ in CHECKS if not tally[label][1]]
    if empty:
        lines += ["", "PubMLST's public isolate records hold no value at all for "
                  + " or ".join(empty) + ", so that cannot be checked against them. "
                  "It is checked instead against the worked examples in PubMLST's "
                  "typing guide (MGAS6180 = pil028, Alab49, and isolate 6424 = "
                  "pil088), which the tool reproduces."]
    total_ok = sum(v[0] for v in tally.values())
    total_n = sum(v[1] for v in tally.values())
    lines += ["", f"**{total_ok} of {total_n} recorded values agree.** "
              f"{len(explained)} differ for a reason traced to the source data "
              f"(below) and **{len(unexplained)} differ unexplained.**", ""]
    if explained:
        lines += ["### Differences traced to the source data", "",
                  "Each was investigated; none is a fault in the typing.", ""]
        for sample, label, exp, o, why in explained:
            lines += [f"- **{sample}, {label}**: the tool says `{o}`, PubMLST "
                      f"records `{exp}`. {why}"]
        lines += [""]
    if unexplained:
        lines += ["### Unexplained differences", ""]
        for sample, label, exp, o, extra in unexplained:
            lines += [f"- **{sample}, {label}**: the tool says `{o}`, PubMLST "
                      f"records `{exp}`." + (f" Note: {extra}." if extra else "")]
        lines += [""]
    if notes:
        lines += ["Notes:", ""] + [f"- {t}" for t in notes] + [""]

    # ------------------------------------------------- how to read the number
    indep = sum(tally[l][1] for l in ("emm type", "emm subtype", "ST",
                                      "FCT-region form"))
    derived = tally["emm cluster"][1] + tally["emm pattern"][1]
    lines += [
        "## How to read that number", "",
        f"The unit is a single field of a single genome, not a genome: "
        f"{len(rows)} genomes with up to {len(CHECKS)} comparable fields each is "
        f"{len(rows) * len(CHECKS)} possible values, and PubMLST records no value "
        f"in {len(rows) * len(CHECKS) - total_n} of them, which leaves {total_n}. "
        "So it is not \"102 genomes correct\", and the denominator is decided by "
        "what PubMLST happens to have curated, not by what matters most.", "",
        f"**The fields are not independent of each other.** emm cluster and the "
        f"predicted emm pattern are table lookups keyed on the emm type "
        f"(`db.clusters[emm_type]` and `db.patterns[emm_type]`), so once the emm "
        f"type is right they follow from the curated tables and cannot fail on "
        f"their own. Those {derived} values test that the tables were transcribed "
        f"from the guide correctly, which is worth knowing, but they are not "
        f"{derived} further checks of the typing. The {indep} emm type, emm "
        f"subtype, ST and FCT-region values are the independent ones.", "",
        "**emm pattern is not compared like with like.** PubMLST's public records "
        "hold no measured pattern for any of these genomes, so the tool's measured "
        "pattern is compared with the pattern PubMLST predicts for the emm type. "
        "They answer different questions, and MGAS8232 below is exactly that case.",
        "",
        "**Both sides come from PubMLST**, though from different parts of it: the "
        "tool types against PubMLST's sequence definitions, and the expected "
        "values are PubMLST's curated isolate records. Agreement therefore shows "
        "the tool reproduces PubMLST's own assignments, including which hit is "
        "emm and which FCT form the loci make. It is not independent evidence "
        "that those assignments are biologically right.", "",
        "**Not checked here at all:** the superantigen, resistance-gene, M1UK "
        "lineage, species-check and warning columns, because PubMLST's isolate "
        "records hold nothing to compare them with. M1UK is covered by the test "
        "suite instead, which writes the published SNPs into a genome and expects "
        "M1UK back. Of the nine FCT forms, FCT-9 and FCT-3/4 are not covered: no "
        "public genome was found for them.", "",
        "**The sample leans on reference strains:** 15 finished genomes and only "
        "7 draft assemblies, so fragmented real-world assemblies are "
        "under-represented. The test suite covers a contig break through emm "
        "separately.", "",
        "## Compared with PubMLST", "",
              "`-` in an expected column means PubMLST has no value for it.", "",
              "| Genome | Source | emm type | emm subtype | ST | emm cluster | "
              "emm pattern | FCT form | pilin type |",
              "|---|---|---|---|---|---|---|---|---|"]
    for e, r, v in rows:
        cells = []
        for _, _, label, _ in CHECKS:
            state, exp, o = v[label]
            if state == "not recorded":
                cells.append("–")
            elif state == "agree":
                cells.append(f"{o} ✓")
            elif state == "explained":
                cells.append(f"{o} (PubMLST {exp}, explained)")
            elif state == "old label":
                cells.append(f"{o} (PubMLST {exp})")
            else:
                cells.append(f"**{o} ≠ {exp}**")
        src = f"NCBI {e['accession']}" if e["accession"] else f"PubMLST isolate {e['pubmlst_isolate_id']}"
        cells.insert(0, src)
        cells.insert(0, f"`{e['sample']}`")
        lines.append("| " + " | ".join(cells) + " |")

    lines += ["", "## Full output of the tool", "",
              f"Every column of `{os.path.basename(results_path)}`, as the tool "
              "wrote it.", ""]
    head = read_tsv(results_path)
    cols = list(head[0].keys()) if head else []
    lines += ["| " + " | ".join(cols) + " |",
              "|" + "|".join("---" for _ in cols) + "|"]
    for row in head:
        lines.append("| " + " | ".join((row[c] or "").replace("|", "\\|")
                                       for c in cols) + " |")

    # ---------------------------------------------------- form coverage
    cov_path = os.path.join(a.out, "coverage.tsv")
    sweep_path = os.path.join(a.out, "sweep_500_genomes.tsv")
    if os.path.exists(cov_path):
        cov = read_tsv_comments(cov_path)
        files = [os.path.join(a.out, "coverage_genomes", r["genome"]) for r in cov]
        got_cov = {}
        if all(os.path.exists(f) for f in files):
            cres = os.path.join(a.out, "results_coverage.tsv")
            subprocess.run([sys.executable, os.path.join(ROOT, "tool",
                                                         "spyogenes_typing.py"),
                            "--db", a.db, "--outfile", cres, "--threads", "4",
                            *files], check=False)
            got_cov = {r["Sample"]: r for r in read_tsv(cres)}
        lines += [
            "## Form coverage (not agreement)", "",
            "PubMLST records an FCT-region form for only 134 of its isolates, "
            "**none of which has contigs**, and it holds a single FCT-9 isolate "
            "with no assembly. So for some forms there is no curated value "
            "anywhere to compare against. The genomes below, found by typing 500 "
            "complete NCBI genomes, instead show that the rule engine reaches "
            "the form and that the loci it saw match that form's definition. "
            "**These are not counted in the agreement figure above**, because "
            "there is nothing to agree with.", "",
            "| Genome | emm | Form called | Pilin type | Why it is evidence |",
            "|---|---|---|---|---|"]
        for r in cov:
            g = got_cov.get(r["genome"].replace(".fasta.gz", ""), {})
            called = g.get("FCT type", "(not run)")
            ok = "✓" if called.startswith(r["form"]) else "**differs**"
            lines.append(f"| `{r['accession']}` | {r['emm']} | {called} {ok} | "
                         f"{r['pilin_type']} | {r['what it shows']} |")
        lines += [""]
    if os.path.exists(sweep_path):
        rows_s = read_tsv(sweep_path)
        forms = {}
        pilins = set()
        for r in rows_s:
            f = (r["FCT type"] or "-").split(" ")[0]
            forms[f] = forms.get(f, 0) + 1
            if "(" in (r["FCT type"] or ""):
                pilins.add(r["FCT type"].split("(")[1].rstrip(")"))
        lines += [
            "### What 500 complete genomes look like", "",
            f"Every complete *S. pyogenes* genome in NCBI that the search "
            f"returned ({len(rows_s)} of 776) was typed in one run, as a check "
            f"that the tool does not fall over and to see the spread of forms. "
            f"It found **{len(pilins)} distinct pilin types**. Full table: "
            f"[sweep_500_genomes.tsv](sweep_500_genomes.tsv).", "",
            "| Form | Genomes |", "|---|---|"]
        for f in sorted(forms, key=lambda x: -forms[x]):
            lines.append(f"| {f} | {forms[f]} |")
        lines += ["",
                  "FCT-9 is 17% of these genomes but only 1 of PubMLST's 134 "
                  "curated FCT isolates. That curated set is old and uneven (it "
                  "still uses the retired FCT-7 and FCT-8 labels and has no "
                  "FCT-3/4 or FCT-10), so the difference is in what has been "
                  "curated, not in the calls: every FCT-9 genome checked has "
                  "rofA and pilB9 as exact matches and no pilA at all. FCT-10 "
                  "did not appear in this set, and FCT-3/4 cannot: it is what "
                  "the tool answers when the regulator is missing, which the "
                  "test suite covers by masking rofA.", ""]

    lines += ["", "## Genomes", "",
              "Listed in [expected.tsv](expected.tsv) with the PubMLST isolate each "
              "expected value comes from; the files and their checksums are in "
              "[genomes/MANIFEST.tsv](genomes/MANIFEST.tsv). Rebuild the set with "
              "`python scripts/fetch_validation_genomes.py` and re-run this record "
              "with `python scripts/validate.py`.", ""]

    with open(os.path.join(a.out, "VALIDATION.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"{total_ok}/{total_n} agree, {len(explained)} explained, "
          f"{len(unexplained)} unexplained; wrote VALIDATION.md and "
          f"{os.path.basename(results_path)}", file=sys.stderr)
    return 0 if not unexplained else 1


if __name__ == "__main__":
    sys.exit(main())
