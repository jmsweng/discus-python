#!/usr/bin/env bash
# S1 step 2: run C0 (perfect), C1 and C2 (random Ga->As antisites, seeds 1 and 2) in WSL.
# Usage (from spike_discus_3d):  bash run_cases.sh [rundir]     default ~/probe_3d
# Observation-only apart from writing into <rundir>. No sudo, no installs.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; RUN="${1:-$HOME/probe_3d}"; mkdir -p "$RUN"
PY="${PYTHON:-python}"; shift 2>/dev/null; CASES="${*:-C0 C1 C2 C3 C4}"
cp "$HERE/../examples/gaas.cll" "$HERE/macros/mono2.cll" "$RUN/"
for c in $CASES; do rm -f "$RUN"/${c}_*.h5 "$RUN"/${c}_*.json "$RUN"/${c}.log; done
gen() { # case seed disorder-lines
  sed -e "s/@CASE@/$1/g" -e "s/@SEED@/$2/g" -e "s|@DISORDER@|$3|" "$HERE/macros/case.mac.in" > "$RUN/$1.mac"; }
for c in $CASES; do case $c in
  C0) gen C0 0 "# perfect crystal: no disorder";;
  C1) gen C1 1 "seed 1\n  replace GA,AS,all,0.15";;
  C2) gen C2 2 "seed 2\n  replace GA,AS,all,0.15";;
  C5) gen C5 1 "seed 1\n  replace GA,AS,all,0.15"; sed -i 's/^  xray$/  xray table:inter/' "$RUN/C5.mac";;
  C3|C4) cp "$HERE/macros/$c.mac" "$RUN/$c.mac";;
esac; done
sed -i 's/\\n/\n/g' "$RUN"/C*.mac
cd "$RUN" || exit 1
sha256sum gaas.cll mono2.cll C*.mac | tee hashes.txt
for c in $CASES; do
  echo "=== $c $(date -u +%FT%TZ)"
  /usr/bin/time -v -o "$c.time.txt" env -u DISPLAY PGPLOT=none discus_suite -noautorun -macro "$c.mac" < /dev/null > "$c.log" 2>&1
  echo "exit $? ; $(grep -E 'Elapsed|Maximum resident' $c.time.txt | tr '\n' ' ')"
  grep -n -i -E "error|unknown|warn|\*\*\*\*" "$c.log" | head -20
done
ls -l *.h5 2>/dev/null | tee h5_list.txt
for f in $(for c in $CASES; do ls ${c}_*.h5 2>/dev/null; done); do [ -e "$f" ] && "$PY" "$HERE/inspect_output.py" "$f" --json "${f%.h5}.inspect.json" > /dev/null && echo "inspected $f"; done
sha256sum *.h5 2>/dev/null | tee -a hashes.txt
tar czf "$HERE/probe_3d.tgz" --exclude="*.h5" -C "$(dirname "$RUN")" "$(basename "$RUN")" 2>/dev/null
echo "Small results tarball (no big .h5 data): $HERE/probe_3d.tgz  -- send me this."
