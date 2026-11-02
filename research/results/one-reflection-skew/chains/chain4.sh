#!/bin/bash
# chain4: half-target A+ (sigma-transport suffices), controls on the solved neighbour one-reflection, and the
# Jacobian-type families at d = 7, 8 (D in J_d only).  8 CPUs.
cd /tmp/hodge/results/one-reflection-skew
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
S="python3 w2_search.py --procs 8"
$S --job aplus-d7-r0s0 --d 7 --rho 0 --sigma 0 --family kernel --target Aplus --restarts 800 --method lm --seed0 17000 --time-limit 5400 >> chain4.out 2>> chain4.err
$S --job aplus-d7-r0s1 --d 7 --rho 0 --sigma 1 --family kernel --target Aplus --restarts 400 --method lm --seed0 17100 --time-limit 3600 >> chain4.out 2>> chain4.err
$S --job aplus-d8-r0s0 --d 8 --rho 0 --sigma 0 --family kernel --target Aplus --restarts 500 --method lm --seed0 18000 --time-limit 7200 >> chain4.out 2>> chain4.err
$S --inst one-reflection --job ctrl-onerefl-d7 --d 7 --rho 0 --sigma 0 --family kernel --restarts 800 --method lm --seed0 27000 --time-limit 5400 >> chain4.out 2>> chain4.err
$S --inst one-reflection --job ctrl-onerefl-d8 --d 8 --rho 0 --sigma 0 --family kernel --restarts 500 --method lm --seed0 28000 --time-limit 7200 >> chain4.out 2>> chain4.err
for d in 7 8; do for r in 0 1; do for s in 0 1; do
  $S --job jac-d${d}-r${r}s${s} --d $d --rho $r --sigma $s --family jacobian --restarts 300 --method lm --seed0 $((100*d+10*r+s)) --time-limit 2700 >> chain4.out 2>> chain4.err
done; done; done
$S --job aplus-d9-r0s0 --d 9 --rho 0 --sigma 0 --family kernel --target Aplus --restarts 300 --method lm --seed0 19000 --time-limit 7200 >> chain4.out 2>> chain4.err
$S --job aplus-d6-r0s0-full --d 6 --rho 0 --sigma 0 --family kernel --target Aplus --restarts 300 --method both --seed0 16000 --time-limit 1800 >> chain4.out 2>> chain4.err
echo CHAIN4 DONE >> chain4.out
