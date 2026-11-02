#!/bin/bash
# wave-2 excl lens: low-degree relaxations, prime index 0 (and index 1 for d=6)
cd /tmp/hodge/results/short-power-presentation
export W2JOB=excl-spp-relax-low
for spec in "5 o" "5 m" "6 m" "6 e --pidx 1" "6 o --pidx 1" "6 m --pidx 1"; do
  echo "=== d/pars: $spec  $(date)"
  python3 excl_relax.py $spec --threads 16 2>&1 | grep -v "hodge engine"
done
echo "=== done $(date)"
