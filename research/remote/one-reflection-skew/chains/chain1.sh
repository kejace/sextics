#!/bin/bash
# chain1: dims table, Jacobian-type degrees (5-var d=5,6, all sectors), optimizer calibration on an accepted
# witness, then the kernel-family joint search at d=7 and d=8 (all sectors).  8 CPUs.
cd /tmp/hodge/results/one-reflection-skew
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
S="python3 w2_search.py --procs 8"
python3 w2_dims.py --dmin 4 --dmax 12 --dmax6 9 > w2-dims.out 2> w2-dims.err
for d in 5 6; do for r in 0 1; do for s in 0 1; do
  $S --job jac-d${d}-r${r}s${s} --d $d --rho $r --sigma $s --family jacobian --restarts 160 --method both --seed0 $((100*d+10*r+s)) --time-limit 1800 >> chain1.out 2>> chain1.err
done; done; done
$S --inst family-hesse-web --job calib-hesse-d6 --d 6 --rho 0 --sigma 0 --family kernel --restarts 400 --method both --seed0 5000 --time-limit 3600 >> chain1.out 2>> chain1.err
$S --inst xtype-strict-rational-113 --job calib-xtype-d6 --d 6 --rho 0 --sigma 1 --family kernel --restarts 200 --method lm --seed0 5500 --time-limit 1800 >> chain1.out 2>> chain1.err
$S --job ker-d7-r0s0 --d 7 --rho 0 --sigma 0 --family kernel --restarts 2000 --method lm --seed0 7000 --time-limit 7200 >> chain1.out 2>> chain1.err
$S --job ker-d7-r0s1 --d 7 --rho 0 --sigma 1 --family kernel --restarts 600 --method lm --seed0 7100 --time-limit 3600 >> chain1.out 2>> chain1.err
$S --job ker-d7-r1s0 --d 7 --rho 1 --sigma 0 --family kernel --restarts 600 --method lm --seed0 7200 --time-limit 3600 >> chain1.out 2>> chain1.err
$S --job ker-d7-r1s1 --d 7 --rho 1 --sigma 1 --family kernel --restarts 600 --method lm --seed0 7300 --time-limit 3600 >> chain1.out 2>> chain1.err
$S --job ker-d8-r0s0 --d 8 --rho 0 --sigma 0 --family kernel --restarts 1200 --method lm --seed0 8000 --time-limit 10800 >> chain1.out 2>> chain1.err
$S --job ker-d8-r0s1 --d 8 --rho 0 --sigma 1 --family kernel --restarts 400 --method lm --seed0 8100 --time-limit 5400 >> chain1.out 2>> chain1.err
$S --job ker-d8-r1s0 --d 8 --rho 1 --sigma 0 --family kernel --restarts 400 --method lm --seed0 8200 --time-limit 5400 >> chain1.out 2>> chain1.err
$S --job ker-d8-r1s1 --d 8 --rho 1 --sigma 1 --family kernel --restarts 400 --method lm --seed0 8300 --time-limit 5400 >> chain1.out 2>> chain1.err
$S --job ker-d7-r0s0-lbfgs --d 7 --rho 0 --sigma 0 --family kernel --restarts 800 --method lbfgs --seed0 7900 --time-limit 7200 >> chain1.out 2>> chain1.err
echo CHAIN1 DONE >> chain1.out
