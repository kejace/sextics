#!/bin/bash
# Degree-8 relaxation restricted to mu6-homogeneous octics, character by character, full target then K3 target.
cd /tmp/hodge/results/polarized-k3-remainder/wave2-excl
export PYTHONPATH=/home/kejace/devel/sextics-work
for t in full k3; do for b in 0 1 2 3 5; do
  echo "=== target=$t char=$b $(date)"
  nice -n 15 python3 relax.py 8 65521 --target=$t --char=$b --squares --out=relax-$t-d8-char$b-p65521 2>&1 | grep -v "engine v2" | grep "squares-summary\|\[result\]" | cut -c1-400
  nice -n 15 python3 survivors.py relax-$t-d8-char$b-p65521 --primes 2>&1 | grep -v "engine v2" | tail -30
done; done
echo "=== CHAR8 DONE $(date)"
