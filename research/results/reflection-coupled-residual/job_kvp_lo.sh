#!/bin/bash
# Koszul-reduced full-E2 variable projection over the complete E1-admissible families at d = 7, 8; 8 parallel seeds per family
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
run8 () {  # run8 <d> <sec> <family> <restarts-per-seed> <seedbase>
  for k in 0 1 2 3 4 5 6 7; do python3 kvp.py --d $1 --sec $2 --family $3 --restarts $4 --seed $(($5 + k)) --tag lo > kvp_d$1_$2_$3_s$(($5 + k)).out 2>&1 & done; wait
}
for s in 00 01 10 11; do run8 7 $s J 50 100; done
for s in e o; do run8 7 $s J 30 200; done
run8 8 00 JK 80 300
for s in 01 10 11; do run8 8 $s J 40 400; done
for s in e o; do run8 8 $s JK 30 500; done
echo KVP_LO_DONE
