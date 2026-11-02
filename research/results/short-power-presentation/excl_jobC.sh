#!/bin/bash
# wave-2 excl lens: degree-7 relaxations (even, odd), prime index 0; msolve first, Singular second
cd /tmp/hodge/results/short-power-presentation
export W2JOB=excl-spp-relax-d7
for spec in "7 e" "7 o"; do
  set -- $spec
  echo "=== d=$1 pars=$2 $(date)"
  python3 excl_relax.py $1 $2 --threads 32 2>&1 | grep -v "hodge engine"
done
echo "=== done $(date)"
