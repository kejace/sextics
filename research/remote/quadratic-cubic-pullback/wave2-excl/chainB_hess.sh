#!/bin/bash
# Hessian relaxation of polar families: d=5 (mixed), d=6 (even, odd, mixed), d=7 (.ms generation; Macaulay skipped automatically)
export OMP_NUM_THREADS=1 FLINT_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /tmp/hodge/results/quadratic-cubic-pullback/wave2-excl
P=1073741827
nice -n 15 python3 hess.py 5 mixed $P 2>&1 | grep -v "engine v2"
for mode in even odd mixed; do nice -n 15 python3 hess.py 6 $mode $P 2>&1 | grep -v "engine v2"; done
for mode in even odd mixed; do nice -n 15 python3 hess.py 7 $mode $P 2>&1 | grep -v "engine v2"; done
echo CHAIN_B_DONE
