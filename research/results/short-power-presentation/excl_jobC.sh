#!/bin/bash
# wave-2 excl lens (instance B): degree-7 pure-parity relaxations, prime 1073741789, msolve (32 threads) then Singular
cd /tmp/hodge/results/short-power-presentation
export W2JOB=exclB-spp-relax-d7
for spec in "7 e" "7 o"; do
  echo "=== d/pars: $spec  $(date)"
  python3 excl_relax.py $spec --threads 32 2>&1 | grep -v "hodge engine"
done
echo "=== done $(date)"
