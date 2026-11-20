#!/usr/bin/env bash
# Assemble a Galaxy-installable tool directory with the database bundled:
#   bash scripts/bundle_tool.sh [outdir]   -> <outdir>/spyogenes_typing/
# The Galaxy files (tool.xml, test-data/, .shed.yml) live in galaxywrapper/,
# the tool script in tool/. Lint and test this bundle, never galaxywrapper/
# (the script and database/ only sit next to tool.xml here, and
# <required_files> must match it).
# Only the database parts the typing reads are bundled: the BLAST indexes,
# the MLST profiles, the curated tables and the M1 SNPs (the allele FASTA
# files are inside the BLAST indexes).
# Left out on purpose: tool/spyogenes_typing_setup.py and
# tool/spyogenes_typing_update.py (the user's own PubMLST key and database).
# They are command-line only; the Galaxy tool always uses the bundled database.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$ROOT/build}"
DEST="$OUT/spyogenes_typing"
DB="$ROOT/database"
[ -f "$DB/blast/spyogenes.nsq" ] || {
  echo "database/ is empty - run scripts/update_database.py first" >&2; exit 1; }
rm -rf "$DEST"
mkdir -p "$DEST/database/pubmlst/mlst" "$DEST/database/references"
cp "$ROOT"/galaxywrapper/tool.xml "$ROOT"/galaxywrapper/.shed.yml "$ROOT"/tool/spyogenes_typing.py "$DEST"/
if [ -f "$ROOT/LICENSE" ]; then cp "$ROOT/LICENSE" "$DEST"/; fi
cp -r "$ROOT"/galaxywrapper/test-data "$DEST"/
cp -r "$DB"/blast "$DB"/tables "$DEST"/database/
cp "$DB"/pubmlst/mlst/mlst.txt "$DEST"/database/pubmlst/mlst/
cp "$DB"/references/m1_snps.tsv "$DEST"/database/references/
cp "$DB"/VERSION.txt "$DB"/manifest.tsv "$DEST"/database/
echo "Bundled tool at $DEST (database $(head -1 "$DB/VERSION.txt"), $(du -sh "$DEST" | cut -f1))"
