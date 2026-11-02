#!/bin/bash
# wave-2 excl lens: low-degree relaxations (d=5,6; even/odd/mixed), prime index 0
cd /tmp/hodge/results/short-power-presentation
export W2JOB=excl-spp-relax-low
for spec in "5 o" "5 m" "6 e" "6 o" "6 m"; do
  set -- $spec
  echo "=== d=$1 pars=$2 $(date)"
  python3 excl_relax.py $1 $2 --threads 8 2>&1 | grep -v "hodge engine"
done
echo "=== second prime (pidx 1) for d=6"
for spec in "6 e" "6 o" "6 m"; do
  set -- $spec
  echo "=== d=$1 pars=$2 pidx=1 $(date)"
  python3 excl_relax.py $1 $2 --threads 8 --pidx 1 2>&1 | grep -v "hodge engine"
done
echo "=== done $(date)"
