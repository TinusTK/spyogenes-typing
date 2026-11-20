# Validation

See [VALIDATION.md](VALIDATION.md) for the validation record: what was
compared, what agreed, and the differences that were traced to the source
data.

- [expected.tsv](expected.tsv) — one row per genome: what PubMLST records
  for that strain, and the PubMLST isolate the values come from.
- [genomes/](genomes) — every genome used, with
  [MANIFEST.tsv](genomes/MANIFEST.tsv) (checksums, length, description).
- `results_<database version>.tsv` — the tool's full report for all of them.

Rebuild the set and the record:

```bash
python scripts/fetch_validation_genomes.py   # genomes + expected.tsv
python scripts/validate.py                   # types them, writes VALIDATION.md
```

`validate.py` exits non-zero if any difference is not one of the explained
ones recorded in its `KNOWN` table, so an unexpected change is noticed.
