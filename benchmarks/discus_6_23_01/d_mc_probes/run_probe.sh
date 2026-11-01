#!/bin/bash
# Package D MMC probe runner (WSL). Usage: bash -l run_probe.sh <inputs dir> <output dir>
# Runs every line "name omp" of runs.txt (omp "-" = OMP_NUM_THREADS unset). Logs, exit codes, hashes.
set -u
SRC="$1"
OUT="$2"
mkdir -p "$OUT/work"
cp "$SRC"/*.stru "$SRC"/*.mac "$SRC"/runs.txt "$OUT/work/"
cd "$OUT/work" || exit 1
sed -i 's/\r$//' *.stru *.mac runs.txt
sha256sum *.stru *.mac runs.txt > "$OUT/input_sha256.txt"
{ command -v discus_suite; sha256sum "$(command -v discus_suite)"; uname -a; nproc; date -u +%Y-%m-%dT%H:%M:%SZ; } > "$OUT/engine.txt"
: > "$OUT/exit_codes.txt"
while read -r r omp; do
  [ -z "$r" ] && continue
  start=$(date +%s.%N)
  if [ "$omp" = "-" ]; then
    env -u OMP_NUM_THREADS timeout 900 discus_suite -noautorun -macro "$r.mac" > "$OUT/$r.log" 2>&1
  else
    OMP_NUM_THREADS="$omp" timeout 900 discus_suite -noautorun -macro "$r.mac" > "$OUT/$r.log" 2>&1
  fi
  code=$?
  end=$(date +%s.%N)
  echo "$r exit=$code omp=$omp wall=$(echo "$end - $start" | bc)" >> "$OUT/exit_codes.txt"
done < runs.txt
sha256sum *_out.stru > "$OUT/output_sha256.txt" 2>/dev/null
cat "$OUT/exit_codes.txt"
