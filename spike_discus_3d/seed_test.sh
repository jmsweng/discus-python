#!/usr/bin/env bash
# Why does `seed 1` not reproduce the same crystal? Builds ONLY the disordered crystal (no Fourier), saves it,
# and compares md5 of the saved structure between repeated runs. Observation only; writes to ~/seed_test.
set -u; RUN="${1:-$HOME/seed_test}"; mkdir -p "$RUN"; cd "$RUN" || exit 1
cp "$(dirname "$(readlink -f "$0")")/../examples/gaas.cll" .
mk() { # name  seed-position(before|after)
  { echo "discus"; [ "$2" = before ] && echo "seed 1"; echo "read"; echo "  cell gaas.cll,8,8,8"
    [ "$2" = after ] && echo "seed 1"; echo "replace GA,AS,all,0.15"; echo "save"; echo "  format keyword"
    echo "  outfile $1.stru"; echo "  run"; echo "exit"; echo "exit"; echo "exit"; } > "$1.mac"; }
mk before before; mk after after
echo "variant                      run1_md5                          run2_md5                          same?"
for v in before after; do for th in default 1; do
  for r in 1 2; do
    if [ $th = 1 ]; then export OMP_NUM_THREADS=1; else unset OMP_NUM_THREADS; fi
    rm -f $v.stru; env -u DISPLAY PGPLOT=none discus_suite -noautorun -macro $v.mac < /dev/null > $v.$th.$r.log 2>&1
    m[$r]=$(md5sum < $v.stru 2>/dev/null | cut -c1-32); done
  printf "seed %-6s threads=%-7s %s %s %s\n" $v $th "${m[1]}" "${m[2]}" "$([ "${m[1]}" = "${m[2]}" ] && echo SAME || echo DIFFERENT)"
done; done | tee seed_result.txt
tar czf "$(dirname "$(readlink -f "$0")")/seed_test.tgz" -C "$RUN" seed_result.txt *.mac *.log
