#!/bin/bash
# Jacobian-type degrees: full E2-reduced varpro over J_{d,chi}, all four sectors, d=5..7
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /tmp/hodge/results/reflection-coupled-residual
for d in 5 6 7; do for chi in "0 0" "1 1" "0 1" "1 0"; do set -- $chi
  python3 w2_search.py --d $d --chi $1 $2 --family jac --mode full --restarts 300 --procs 4 --seed 11 --tag jac-d$d-chi$1$2-full
done; done
echo ALLDONE
