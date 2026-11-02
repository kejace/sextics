#!/bin/bash
# Hessian relaxation per character sector (semi-invariant denominators), gauge-fixed m-primary test, all 18 sectors.
# usage: bash sectors_hess.sh d p threads
d=$1; p=$2; T=$3
cd /tmp/hodge/results/cyclic-triple-residual/wave2-excl
for c0 in 0 2 4; do for c1 in 0 1 2 3 4 5; do
  SKIP_FORMS=1 SECTOR=$c0,$c1 python3 relax_run.py $d $p $T > hsector_d${d}_${c0}${c1}.out 2> hsector_d${d}_${c0}${c1}.err
  grep "gauge-fixed" hsector_d${d}_${c0}${c1}.out | sed "s/^/sector ($c0,$c1): /"
done; done
