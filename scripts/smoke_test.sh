#!/usr/bin/env bash
# Offline end-to-end check on the tiny bundled structure: run -> validate -> reproduce.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="$(mktemp -d)"
FIX=tests/fixtures/tiny_test.pdb
pdbmi run --structure "$FIX" --chain A --mutation A50V --data-dir "$OUT" --output "$OUT/report.json"
pdbmi validate "$OUT/report.json"
pdbmi reproduce "$OUT/report.json" --structure "$FIX"
pdbmi dag > /dev/null
echo "SMOKE TEST PASSED"
