#!/bin/bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 rcr2.py e1 --d 9 --sec all > e1_d9.out 2>&1
run8 () { for k in 0 1 2 3; do python3 kvp.py --d $1 --sec $2 --family $3 --restarts $4 --seed $(($5 + k)) --opt lm --tag d9 > kvp_d$1_$2_$3_s$(($5 + k)).out 2>&1 & done; for k in 4 5 6 7; do python3 kvp.py --d $1 --sec $2 --family $3 --restarts $4 --seed $(($5 + k)) --opt lbfgs --maxiter 800 --tag d9 > kvp_d$1_$2_$3_s$(($5 + k)).out 2>&1 & done; wait; }
for s in 00 01 10 11; do python3 rcr2.py famrelax --d 9 --sec $s --family JK --restarts 1500 --seed 203 > famrelax_d9_$s.out 2>&1; done
for s in 00 01 10 11; do run8 9 $s JK 40 600; done
for s in e o; do run8 9 $s JK 20 700; done
echo KVP_D9_DONE
