#!/usr/bin/env bash
# S1 step 2: run selected DISCUS 3D probe cases (C0..C5) in WSL2.
#
# Usage (from spike_discus_3d):
#   bash run_cases.sh <fresh-output-dir> [C0 C1 ... C5]
# Without case names the default set is C0 C1 C2 C3 C4.
#
# Safety model:
#   * Refuses to write into an occupied output directory (use --overwrite to
#     force; this is the ONLY way to clobber previous evidence).
#   * Validates case names BEFORE writing any files.
#   * Generates and processes ONLY the selected case macros.
#   * Treats a missing engine, nonzero DISCUS exit, missing/empty/malformed
#     HDF5 output, or a failed inspector as a NONZERO script exit.
#   * Expected HDF5 names are derived from each generated macro's `outf`
#     declarations (C0/C1/C2/C5 -> 3 files each; C3 -> 2; C4 -> 1).
#   * Saves per-case logs, timing, inspector JSON and checksums in the output
#     directory. The optional archive is uniquely named and never overwrites
#     files in the source tree. Large HDF5 files are retained in the output
#     directory and are NOT committed to Git.
set -u
set -o pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"

# ---- argument parsing ------------------------------------------------------
OVERWRITE=0
ARGS=()
for a in "$@"; do
  case "$a" in
    --overwrite) OVERWRITE=1 ;;
    *) ARGS+=("$a") ;;
  esac
done

if [ "${#ARGS[@]}" -lt 1 ]; then
  echo "Usage: bash run_cases.sh <fresh-output-dir> [C0 C1 ... C5]" >&2
  echo "  (add --overwrite to clobber an occupied output directory)" >&2
  exit 2
fi

RUN="${ARGS[0]}"
if [ "${#ARGS[@]}" -gt 1 ]; then
  CASES=("${ARGS[@]:1}")
else
  CASES=(C0 C1 C2 C3 C4)
fi

# Validate case names BEFORE writing anything.
for c in "${CASES[@]}"; do
  case "$c" in
    C0|C1|C2|C3|C4|C5) ;;
    *) echo "Unknown case '$c' (choose C0 C1 C2 C3 C4 C5)" >&2; exit 2 ;;
  esac
done

# Refuse an occupied/colliding run directory (never silently delete evidence).
if [ -e "$RUN" ] && [ ! -d "$RUN" ]; then
  echo "Refusing: '$RUN' exists and is not a directory." >&2
  exit 3
fi
if [ -d "$RUN" ] && [ -n "$(ls -A "$RUN" 2>/dev/null)" ]; then
  if [ "$OVERWRITE" -ne 1 ]; then
    echo "Refusing: output directory '$RUN' is not empty." >&2
    echo "Delete it or pass --overwrite to clobber previous evidence." >&2
    exit 3
  fi
  # Explicit overwrite: remove prior contents so stale files cannot be mistaken
  # for fresh results.
  rm -rf "$RUN"
fi
mkdir -p "$RUN"

PY="${PYTHON:-python3}"
ENGINE="${DISCUS_ENGINE:-discus_suite}"

# ---- stage shared inputs ---------------------------------------------------
cp "$HERE/../examples/gaas.cll" "$RUN/gaas.cll"
cp "$HERE/macros/mono2.cll" "$RUN/mono2.cll"

# ---- generate ONLY the selected case macros --------------------------------
# Use literal \n in the disorder text; convert to real newlines after generation.
gen() { # $1=case  $2=seed  $3=disorder-lines (with literal \n)
  sed -e "s/@CASE@/$1/g" -e "s/@SEED@/$2/g" -e "s|@DISORDER@|$3|" \
      "$HERE/macros/case.mac.in" > "$RUN/$1.mac"
  sed -i 's/\\n/\n/g' "$RUN/$1.mac"
}
for c in "${CASES[@]}"; do
  case "$c" in
    C0) gen C0 0 "# perfect crystal: no disorder" ;;
    C1) gen C1 1 "seed 1\n  replace GA,AS,all,0.15" ;;
    C2) gen C2 2 "seed 2\n  replace GA,AS,all,0.15" ;;
    C5) gen C5 1 "seed 1\n  replace GA,AS,all,0.15"
        sed -i 's/^  xray$/  xray table:inter/' "$RUN/C5.mac" ;;
    C3|C4) cp "$HERE/macros/$c.mac" "$RUN/$c.mac" ;;
  esac
done

# ---- record staged hashes --------------------------------------------------
cd "$RUN" || exit 1
sha256sum gaas.cll mono2.cll "${CASES[@]}.mac" | tee hashes.txt

# ---- run each selected case; fail fast with a clear summary ---------------
FAILED=()
for c in "${CASES[@]}"; do
  echo "=== $c $(date -u +%FT%TZ)"
  if ! command -v "$ENGINE" >/dev/null 2>&1; then
    echo "ERROR: engine '$ENGINE' not found on PATH; aborting." >&2
    FAILED+=("$c:engine-missing")
    break
  fi
  /usr/bin/time -v -o "$c.time.txt" env -u DISPLAY PGPLOT=none "$ENGINE" -noautorun -macro "$c.mac" < /dev/null > "$c.log" 2>&1
  rc=$?
  echo "exit $rc ; $(grep -E 'Elapsed|Maximum resident' "$c.time.txt" 2>/dev/null | tr '\n' ' ')"
  if [ "$rc" -ne 0 ]; then
    echo "ERROR: $c engine exited $rc (see $c.log)" >&2
    grep -n -i -E "error|unknown|warn|\*\*\*\*" "$c.log" | head -20 >&2
    FAILED+=("$c:engine-exit-$rc")
    break
  fi
  # Expected HDF5 names come from THIS macro's outf declarations.
  expected="$(grep -oE 'outf[[:space:]]+[^[:space:]]+' "$c.mac" | awk '{print $2}')"
  if [ -z "$expected" ]; then
    echo "ERROR: $c macro declares no 'outf' output; aborting." >&2
    FAILED+=("$c:no-outf")
    break
  fi
  for f in $expected; do
    if [ ! -e "$f" ] || [ ! -s "$f" ]; then
      echo "ERROR: $c expected output '$f' missing or empty; aborting." >&2
      FAILED+=("$c:missing-$f")
      break
    fi
    if ! "$PY" "$HERE/inspect_output.py" "$f" --json "${f%.h5}.inspect.json" >/dev/null 2>&1; then
      echo "ERROR: $c inspector failed for '$f'; aborting." >&2
      FAILED+=("$c:inspect-fail-$f")
      break
    fi
    if ! grep -q '"format": "hdf5"' "${f%.h5}.inspect.json"; then
      echo "ERROR: $c output '$f' is not valid HDF5; aborting." >&2
      FAILED+=("$c:not-hdf5-$f")
      break
    fi
  done
  [ "${#FAILED[@]}" -gt 0 ] && break
done

# ---- summary ---------------------------------------------------------------
if [ "${#FAILED[@]}" -gt 0 ]; then
  echo "=== FAILED: $(printf '%s; ' "${FAILED[@]}")" >&2
  exit 1
fi
echo "=== ALL SELECTED CASES SUCCEEDED: ${CASES[*]}"
ls -l *.h5 2>/dev/null | tee h5_list.txt
sha256sum *.h5 2>/dev/null | tee -a hashes.txt

# Optional, uniquely named archive (never overwrites the source tree).
ARCHIVE="$RUN/probe_3d_$(date -u +%Y%m%dT%H%M%SZ).tgz"
tar czf "$ARCHIVE" --exclude='*.h5' -C "$(dirname "$RUN")" "$(basename "$RUN")" 2>/dev/null \
  && echo "Archive (no big .h5): $ARCHIVE" \
  || echo "Archive creation skipped."
