#!/bin/bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /tmp/hodge/results/reflection-coupled-residual
for d in 5 6 7; do for chi in "0 0" "1 1" "0 1" "1 0"; do set -- $chi
  python3 w2_search.py --d $d --chi $1 $2 --family jac --mode ob --restarts 3000 --procs 4 --seed 21 --tag jac-d$d-chi$1$2-ob
done; done
echo ALLDONE
