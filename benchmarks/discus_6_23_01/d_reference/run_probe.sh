#!/bin/bash
# Package D reference-scattering probe runner (WSL). Usage: bash -l run_probe.sh <inputs dir> <output dir>
# Copies the inputs into <output dir>/work (LF line endings), runs every *.mac with discus_suite, keeps logs.
set -u
SRC="$1"
OUT="$2"
mkdir -p "$OUT/work"
cp "$SRC"/*.stru "$SRC"/*.mac "$OUT/work/"
cd "$OUT/work" || exit 1
sed -i 's/\r$//' *.stru *.mac
sha256sum *.stru *.mac > "$OUT/input_sha256.txt"
{ command -v discus_suite; sha256sum "$(command -v discus_suite)"; uname -a; date -u +%Y-%m-%dT%H:%M:%SZ; echo "OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}"; } > "$OUT/engine.txt"
: > "$OUT/exit_codes.txt"
for m in *.mac; do
  r="${m%.mac}"
  timeout 600 discus_suite -noautorun -macro "$m" > "$OUT/$r.log" 2>&1
  echo "$r exit=$?" >> "$OUT/exit_codes.txt"
done
sha256sum *.h5 > "$OUT/output_sha256.txt" 2>/dev/null
cat "$OUT/exit_codes.txt"
