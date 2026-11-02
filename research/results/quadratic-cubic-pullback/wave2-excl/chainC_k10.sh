#!/bin/bash
# d = 10 E1 kernels, both parities, second prime (basis + (D,N) pairs), then first prime for cross-check
export OMP_NUM_THREADS=1 FLINT_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /tmp/hodge/results/quadratic-cubic-pullback/wave2-excl
for d in 10 11 12; do for par in 0 1; do nice -n 15 python3 dual.py $d $par 1073741827 2>&1 | grep -v "engine v2"; done; done
for P in 1073741827 2147483647; do
  for par in 0 1; do
    [ -f e1k_d10_par${par}_p${P}.json ] && continue
    nice -n 15 python3 e1k.py 10 $par $P --basis 2>&1 | grep -v "engine v2"
  done
done
echo CHAIN_C_DONE
