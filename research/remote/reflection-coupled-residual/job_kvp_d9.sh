#!/bin/bash
# d = 9: E1 spaces per sector (cached: 00 and 10 carry 21 non-Jacobian classes each; 01, 11 are J-only), then Koszul varpro
# over kernel + J_9 families; 8 parallel seeds per family (4 LM, 4 L-BFGS with 900 iterations)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 rcr2.py e1 --d 9 --sec all > e1_d9.out 2>&1
run8 () { for k in 0 1 2 3; do python3 kvp.py --d $1 --sec $2 --family $3 --restarts $4 --seed $(($5 + k)) --opt lm --maxiter 150 --tag d9 > kvp_d$1_$2_$3_s$(($5 + k)).out 2>&1 & done; for k in 4 5 6 7; do python3 kvp.py --d $1 --sec $2 --family $3 --restarts $4 --seed $(($5 + k)) --opt lbfgs --maxiter 900 --tag d9 > kvp_d$1_$2_$3_s$(($5 + k)).out 2>&1 & done; wait; }
run8 9 00 JK 25 600
run8 9 10 JK 25 620
run8 9 01 J 15 640
run8 9 11 J 15 660
for s in e o; do run8 9 $s JK 10 680; done
echo KVP_D9_DONE
