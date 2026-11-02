#!/bin/bash
# Package C engine probe runner (WSL). Usage: bash -l run_probes.sh <inputs dir> <output dir>
# Copies the inputs into <output dir>/work, runs the three macros with discus_suite, keeps the logs.
set -u
SRC="$1"
OUT="$2"
mkdir -p "$OUT/work"
cp "$SRC"/*.stru "$SRC"/*.mac "$OUT/work/"
cd "$OUT/work" || exit 1
sha256sum *.stru *.mac > "$OUT/input_sha256.txt"
{ command -v discus_suite; sha256sum "$(command -v discus_suite)"; uname -a; date -u +%Y-%m-%dT%H:%M:%SZ; } > "$OUT/engine.txt"
: > "$OUT/exit_codes.txt"
for m in occ_bcc disp_sc mixed_bcc; do
  timeout 300 discus_suite -noautorun -macro "$m.mac" > "$OUT/$m.log" 2>&1
  echo "$m exit=$?" >> "$OUT/exit_codes.txt"
done
cat "$OUT/exit_codes.txt"
