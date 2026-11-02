#!/bin/bash
# Sector-by-sector Koszul relaxation at degree d (semi-invariant denominators), gauge-fixed m-primary test.
# usage: bash sectors.sh d p threads   (requires kernel_exact_d{d}.json from kernel_multi.py)
d=$1; p=$2; T=$3
cd /tmp/hodge/results/cyclic-triple-residual/wave2-excl
for c0 in 0 2 4; do for c1 in 0 1 2 3 4 5; do
  SKIP_FORMS=1 python3 relax0.py $d $p $T sector $c0,$c1 > sector_d${d}_${c0}${c1}.out 2> sector_d${d}_${c0}${c1}.err
  tail -1 sector_d${d}_${c0}${c1}.out
done; done
