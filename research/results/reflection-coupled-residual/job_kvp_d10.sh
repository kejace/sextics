#!/bin/bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 rcr2.py e1 --d 10 --sec all > e1_d10.out 2>&1
run8 () { for k in 0 1 2 3; do python3 kvp.py --d $1 --sec $2 --family $3 --restarts $4 --seed $(($5 + k)) --opt lm --tag d10 > kvp_d$1_$2_$3_s$(($5 + k)).out 2>&1 & done; for k in 4 5 6 7; do python3 kvp.py --d $1 --sec $2 --family $3 --restarts $4 --seed $(($5 + k)) --opt lbfgs --maxiter 800 --tag d10 > kvp_d$1_$2_$3_s$(($5 + k)).out 2>&1 & done; wait; }
for s in 00 01 10 11; do run8 10 $s JK 15 800; done
echo KVP_D10_DONE
