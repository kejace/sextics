#!/bin/bash
# chain3: kernel-family joint search at the higher degrees d = 9..12 (sector (0,0) first, then (0,1)).  8 CPUs.
cd /tmp/hodge/results/one-reflection-skew
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
S="python3 w2_search.py --procs 8"
$S --job ker-d9-r0s0 --d 9 --rho 0 --sigma 0 --family kernel --restarts 600 --method lm --seed0 9000 --time-limit 10800 >> chain3.out 2>> chain3.err
$S --job ker-d10-r0s0 --d 10 --rho 0 --sigma 0 --family kernel --restarts 300 --method lm --seed0 10000 --time-limit 10800 >> chain3.out 2>> chain3.err
$S --job ker-d9-r0s1 --d 9 --rho 0 --sigma 1 --family kernel --restarts 300 --method lm --seed0 9100 --time-limit 7200 >> chain3.out 2>> chain3.err
$S --job ker-d11-r0s0 --d 11 --rho 0 --sigma 0 --family kernel --restarts 150 --method lm --seed0 11000 --time-limit 10800 >> chain3.out 2>> chain3.err
$S --job ker-d12-r0s0 --d 12 --rho 0 --sigma 0 --family kernel --restarts 80 --method lm --seed0 12000 --time-limit 10800 >> chain3.out 2>> chain3.err
$S --job ker-d10-r0s1 --d 10 --rho 0 --sigma 1 --family kernel --restarts 200 --method lm --seed0 10100 --time-limit 7200 >> chain3.out 2>> chain3.err
$S --job ker-d9-r1s0 --d 9 --rho 1 --sigma 0 --family kernel --restarts 300 --method lm --seed0 9200 --time-limit 7200 >> chain3.out 2>> chain3.err
$S --job ker-d9-r1s1 --d 9 --rho 1 --sigma 1 --family kernel --restarts 300 --method lm --seed0 9300 --time-limit 7200 >> chain3.out 2>> chain3.err
echo CHAIN3 DONE >> chain3.out
