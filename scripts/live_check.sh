#!/usr/bin/env bash
# Needs network. Real structure (hen egg-white lysozyme 1AKI): offline measurements, then a 3-position scan.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-data/results/live}"
mkdir -p "$OUT"
pdbmi run --pdb 1AKI --chain A --mutation S50V --output "$OUT/1aki-S50V.json"
pdbmi validate "$OUT/1aki-S50V.json"
pdbmi reproduce "$OUT/1aki-S50V.json"
pdbmi scan --pdb 1AKI --chain A --positions 50-52 --out "$OUT/scan"
echo "LIVE CHECK PASSED (remote predictors are separate: see README 'Remote services')"
