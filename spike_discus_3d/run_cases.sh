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
#   * --overwrite is guarded: it refuses to delete protected paths (home,
#     repo root, script dir, /, ., ..) and symlink aliases that resolve
#     into protected locations.
#   * Validates case names BEFORE writing any files.
#   * Generates and processes ONLY the selected case macros.
#   * Treats a missing engine, nonzero DISCUS exit, missing/empty/malformed
#     HDF5 output, or a failed inspector as a NONZERO script exit.
#   * Expected HDF5 names are derived from each generated macro's `outf`
#     declarations (C0/C1/C2/C5 -> 3 files each; C3 -> 2; C4 -> 1).
#   * Saves per-case logs, timing, inspector JSON and checksums in the output
#     directory. The optional archive is uniquely named, created OUTSIDE the
#     archived tree, and never overwrites files in the source tree. Large HDF5
#     files are retained in the output directory and are NOT committed to Git.
set -u
set -o pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$HERE/.." && pwd)"

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

# ---- resolve and validate the output directory -----------------------------
# Convert to an absolute path early.
if [ ! -e "$RUN" ]; then
  # Path does not exist yet: resolve the parent and append the basename.
  parent="$(dirname "$RUN")"
  if [ -z "$parent" ] || [ "$parent" = "." ]; then
    RUN_ABS="$PWD/$(basename "$RUN")"
  else
    parent_abs="$(cd "$parent" 2>/dev/null && pwd)" || {
      echo "Refusing: cannot resolve parent directory of '$RUN'." >&2
      exit 3
    }
    RUN_ABS="$parent_abs/$(basename "$RUN")"
  fi
else
  # Path exists: resolve it.
  RUN_ABS="$(cd "$RUN" 2>/dev/null && pwd)" || {
    echo "Refusing: cannot resolve path '$RUN'." >&2
    exit 3
  }
fi

# Guard destructive --overwrite: refuse protected paths.
if [ "$OVERWRITE" -eq 1 ]; then
  # Normalize: strip trailing slashes.
  norm="$RUN_ABS"
  while [ "${norm%/}" != "$norm" ]; do norm="${norm%/}"; done
  # Refuse the root.
  if [ "$norm" = "/" ]; then
    echo "Refusing --overwrite: cannot target /." >&2
    exit 3
  fi
  # Refuse the home directory.
  if [ "$norm" = "$HOME" ]; then
    echo "Refusing --overwrite: cannot target the home directory." >&2
    exit 3
  fi
  # Refuse the repository root.
  if [ "$norm" = "$REPO_ROOT" ]; then
    echo "Refusing --overwrite: cannot target the repository root." >&2
    exit 3
  fi
  # Refuse the script/source directory.
  if [ "$norm" = "$HERE" ]; then
    echo "Refusing --overwrite: cannot target the script/source directory." >&2
    exit 3
  fi
  # Refuse if the path is an ancestor of the repo root.
  case "$REPO_ROOT/" in
    "$norm"/*) echo "Refusing --overwrite: path is an ancestor of the repository." >&2; exit 3 ;;
  esac
  # Refuse if the path is INSIDE the repo root (any subdirectory, including
  # .git, src, tests, etc.).
  case "$norm" in
    "$REPO_ROOT"/*) echo "Refusing --overwrite: path is inside the repository." >&2; exit 3 ;;
  esac
  # Refuse if the resolved path is a symlink pointing into a protected location.
  if [ -L "$RUN_ABS" ]; then
    target="$(readlink -f "$RUN_ABS" 2>/dev/null)" || target=""
    case "$target" in
      "$HOME"*) echo "Refusing --overwrite: symlink resolves into protected home directory." >&2; exit 3 ;;
      "$REPO_ROOT"*) echo "Refusing --overwrite: symlink resolves into the repository." >&2; exit 3 ;;
      "$HERE"*) echo "Refusing --overwrite: symlink resolves into the script directory." >&2; exit 3 ;;
    esac
  fi
fi

# Refuse an occupied/colliding run directory (never silently delete evidence).
if [ -e "$RUN_ABS" ] && [ ! -d "$RUN_ABS" ]; then
  echo "Refusing: '$RUN_ABS' exists and is not a directory." >&2
  exit 3
fi
if [ -d "$RUN_ABS" ] && [ -n "$(ls -A "$RUN_ABS" 2>/dev/null)" ]; then
  if [ "$OVERWRITE" -ne 1 ]; then
    echo "Refusing: output directory '$RUN_ABS' is not empty." >&2
    echo "Delete it or pass --overwrite to clobber previous evidence." >&2
    exit 3
  fi
  # Explicit overwrite: remove prior contents so stale files cannot be mistaken
  # for fresh results.
  echo "WARNING: --overwrite: deleting contents of '$RUN_ABS'"
  rm -rf "$RUN_ABS"
fi
mkdir -p "$RUN_ABS"

PY="${PYTHON:-python3}"
ENGINE="${DISCUS_ENGINE:-discus_suite}"

# ---- stage shared inputs ---------------------------------------------------
cp "$REPO_ROOT/examples/gaas.cll" "$RUN_ABS/gaas.cll"
cp "$HERE/macros/mono2.cll" "$RUN_ABS/mono2.cll"

# ---- generate ONLY the selected case macros --------------------------------
# Use literal \n in the disorder text; convert to real newlines after generation.
gen() { # $1=case  $2=seed  $3=disorder-lines (with literal \n)
  sed -e "s/@CASE@/$1/g" -e "s/@SEED@/$2/g" -e "s|@DISORDER@|$3|" \
      "$HERE/macros/case.mac.in" > "$RUN_ABS/$1.mac"
  sed -i 's/\\n/\n/g' "$RUN_ABS/$1.mac"
}
for c in "${CASES[@]}"; do
  case "$c" in
    C0) gen C0 0 "# perfect crystal: no disorder" ;;
    C1) gen C1 1 "seed 1\n  replace GA,AS,all,0.15" ;;
    C2) gen C2 2 "seed 2\n  replace GA,AS,all,0.15" ;;
    C5) gen C5 1 "seed 1\n  replace GA,AS,all,0.15"
        sed -i 's/^  xray$/  xray table:inter/' "$RUN_ABS/C5.mac" ;;
    C3|C4) cp "$HERE/macros/$c.mac" "$RUN_ABS/$c.mac" ;;
  esac
done

# ---- record staged hashes --------------------------------------------------
cd "$RUN_ABS" || exit 1
# Build a proper array of macro filenames from the validated case names.
MACROS=()
for c in "${CASES[@]}"; do MACROS+=("$c.mac"); done
sha256sum gaas.cll mono2.cll "${MACROS[@]}" > hashes.txt || {
  echo "ERROR: failed to compute SHA-256 manifest; aborting." >&2
  exit 1
}

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

# Optional, uniquely named archive created OUTSIDE the archived tree.
ARCHIVE_BASE="$(date -u +%Y%m%dT%H%M%SZ)"
PARENT_DIR="$(dirname "$RUN_ABS")"
TMP_ARCHIVE="$(mktemp "$PARENT_DIR/.probe_3d_${ARCHIVE_BASE}_XXXXXX.tgz")" || {
  echo "WARNING: could not create temp archive file." >&2
  exit 1
}
if tar czf "$TMP_ARCHIVE" --exclude='*.h5' -C "$PARENT_DIR" "$(basename "$RUN_ABS")" 2>/dev/null; then
  ARCHIVE="$RUN_ABS/probe_3d_${ARCHIVE_BASE}.tgz"
  # Do not silently overwrite an existing archive with the same timestamp.
  if [ -e "$ARCHIVE" ]; then
    mv "$TMP_ARCHIVE" "${ARCHIVE}.tmp_renamed"
    echo "WARNING: archive '$ARCHIVE' already exists; saved as '${ARCHIVE}.tmp_renamed'"
  else
    mv "$TMP_ARCHIVE" "$ARCHIVE"
    echo "Archive (no big .h5): $ARCHIVE"
  fi
else
  rm -f "$TMP_ARCHIVE"
  echo "WARNING: archive creation failed." >&2
fi
