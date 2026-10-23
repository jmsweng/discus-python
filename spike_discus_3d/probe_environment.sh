#!/usr/bin/env bash
# S1 step 0 — OBSERVATION ONLY. Run inside Ubuntu/WSL2 (same OS as DISCUS).
# Collects facts about the installed DISCUS_SUITE and its documentation.
# It does not install anything, use sudo, or run a diffraction calculation.
# Usage: bash probe_environment.sh [outdir]   (default: ./probe_out_<utc>)
set -u
OUT="${1:-probe_out_$(date -u +%Y%m%dT%H%M%SZ)}"
mkdir -p "$OUT"; cd "$OUT" || exit 1
log() { echo "== $*" | tee -a summary.txt; }

log "date/host"; { date -u; uname -a; cat /etc/os-release | head -4; grep -i microsoft /proc/version || echo "not WSL?"; nproc; free -m | head -2; } 2>&1 | tee -a summary.txt
log "which discus_suite"; BIN=$(command -v discus_suite || true); echo "${BIN:-NOT FOUND}" | tee -a summary.txt
if [ -n "$BIN" ]; then
  ls -l "$BIN" | tee -a summary.txt; sha256sum "$(readlink -f "$BIN")" | tee -a summary.txt; file "$(readlink -f "$BIN")" | tee -a summary.txt
  log "dynamic libs (hdf5 / x11 / etc.)"; ldd "$BIN" 2>&1 | tee ldd.txt | grep -i -E "hdf|netcdf|x11|pgplot|png" | tee -a summary.txt || true
fi
log "python/venv"; for p in "$HOME/projects/discuspdf_starter/.venv/bin/python" python3; do
  [ -x "$(command -v $p 2>/dev/null || echo $p)" ] && { $p -c "import sys;print(sys.executable, sys.version.split()[0])"; $p -c "import numpy,h5py;print('numpy',numpy.__version__,'h5py',h5py.__version__)" 2>&1 | tail -1; } | tee -a summary.txt
done

log "banner (headless, empty macro)"
printf 'exit\n' > banner.mac
env -u DISPLAY PGPLOT=none timeout 60 discus_suite -noautorun -macro banner.mac < /dev/null > banner.log 2>&1; echo "exit code $?" | tee -a summary.txt
head -40 banner.log | tee -a summary.txt

log "documentation files"
find /usr/local/share /usr/local -maxdepth 3 \( -iname "*man*.pdf" -o -iname "*discus*" -o -iname "*suite*" -o -iname "*kuplot*" \) 2>/dev/null | head -50 | tee docfiles.txt | tee -a summary.txt
log "keyword hits in manuals (uses pdftotext if installed; otherwise lists PDFs only)"
if command -v pdftotext >/dev/null; then
  while read -r f; do
    if [[ "$f" == *.pdf ]]; then
      b=$(basename "$f" .pdf); pdftotext -layout "$f" "$b.txt" 2>/dev/null
      echo "-- $f -> $b.txt ($(wc -l < "$b.txt") lines)" | tee -a summary.txt
      grep -n -i -E "hdf|[.]h5|lots|bragg|3d|three.dimensional|mrc|ccp4|diffuse|delta.?pdf|reciprocal" "$b.txt" | head -150 > "kw_$b.txt"
    fi
  done < docfiles.txt
else echo "pdftotext not installed; copy the PDFs listed above back instead (do NOT sudo-install without owner approval)" | tee -a summary.txt; fi

log "in-suite help for the fourier menu (discovery; errors are data, keep them)"
cat > help_fourier.mac <<'MAC'
discus
fourier
help
exit
exit
exit
MAC
env -u DISPLAY PGPLOT=none timeout 120 discus_suite -noautorun -macro help_fourier.mac < /dev/null > help_fourier.log 2>&1; echo "exit code $?" | tee -a summary.txt
wc -l help_fourier.log | tee -a summary.txt; head -80 help_fourier.log
log "done -> $(pwd)"; echo "Return the whole directory (tar czf) plus any manual PDFs not converted."
