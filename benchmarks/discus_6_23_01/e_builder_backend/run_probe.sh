#!/bin/bash
# Builder-backend probe runner (WSL). Usage: bash -l run_probe.sh <inputs dir> <output dir>
# Copies every run directory into <output dir>/work (LF line endings) and runs its run.mac with
# discus_suite in the same environment as PDFRunner (no DISPLAY, PGPLOT=none, stdin "exit").
set -u
SRC="$1"
OUT="$2"
mkdir -p "$OUT/work"
OUT=$(realpath -m "$OUT")   # absolute: the script cds into $OUT/work
cp -r "$SRC"/* "$OUT/work/"
cd "$OUT/work" || exit 1
find . -type f \( -name '*.mac' -o -name '*.stru' -o -name '*.cll' \) -exec sed -i 's/\r$//' {} +
find . -type f \( -name '*.mac' -o -name '*.stru' -o -name '*.cll' \) | sort | xargs sha256sum > "$OUT/input_sha256.txt"
{ command -v discus_suite; sha256sum "$(command -v discus_suite)"; uname -a; date -u +%Y-%m-%dT%H:%M:%SZ; echo "OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}"; } > "$OUT/engine.txt"
: > "$OUT/exit_codes.txt"
for d in */; do
  r="${d%/}"
  ( cd "$r" && echo exit | env -u DISPLAY PGPLOT=none timeout 900 discus_suite -noautorun -macro run.mac > "../../$r.log" 2>&1 )
  echo "$r exit=$?" >> "$OUT/exit_codes.txt"
done
find . -type f \( -name '*.h5' -o -name 'result.pdf' \) | sort | xargs sha256sum > "$OUT/output_sha256.txt" 2>/dev/null
cat "$OUT/exit_codes.txt"
