#!/bin/bash
# E1 annihilator sweep, second prime, both parities (re-verification d<=8 even, d<=6 odd; extension to d=9)
export OMP_NUM_THREADS=1 FLINT_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /tmp/hodge/results/quadratic-cubic-pullback/wave2-excl
P=1073741827
for d in 1 2 3 4 5 6 7 8 9; do
  for par in 1 0; do
    [ -f e1k_d${d}_par${par}_p${P}.json ] && continue
    nice -n 15 python3 e1k.py $d $par $P --basis 2>&1 | grep -v "engine v2"
  done
done
echo CHAIN_A_DONE
