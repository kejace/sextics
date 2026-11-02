#!/bin/bash
# wave-2 excl lens (instance B): mixed-parity relaxations d=6 (two primes) then d=7
cd /tmp/hodge/results/short-power-presentation
export W2JOB=exclB-spp-relax-mixed
for spec in "6 m --pidx 0" "6 m --pidx 1" "7 m --pidx 0"; do
  echo "=== d/pars: $spec  $(date)"
  python3 excl_relax.py $spec --threads 32 --no-singular 2>&1 | grep -v "hodge engine"
done
echo "=== done $(date)"
