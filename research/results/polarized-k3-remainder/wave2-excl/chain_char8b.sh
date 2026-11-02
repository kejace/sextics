#!/bin/bash
# Remaining degree-8 character sectors (full target chars 1,2,3,5; K3 target chars 0,1,2,3,5); residual dim/degree only.
cd /tmp/hodge/results/polarized-k3-remainder/wave2-excl
export PYTHONPATH=/home/kejace/devel/sextics-work
run_sector() { t=$1; b=$2
  echo "=== target=$t char=$b $(date)"
  nice -n 15 python3 relax.py 8 65521 --target=$t --char=$b --squares --out=relax-$t-d8-char$b-p65521 2>&1 | grep -v "engine v2" | grep "squares-summary\|\[result\]" | cut -c1-500
  timeout 1200 nice -n 15 python3 survivors.py relax-$t-d8-char$b-p65521 2>&1 | grep -v "engine v2" | grep -E "params|c[0-9]+ =|residual|DIM=" | head -60
}
for b in 1 2 3 5; do run_sector full $b; done
for b in 0 1 2 3 5; do run_sector k3 $b; done
echo "=== CHAR8B DONE $(date)"
