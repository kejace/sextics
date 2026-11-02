#!/bin/bash
# chain5: optimizer calibration at larger K (hesse-web d=8: the needle set contains D*·S_2, a 6-dim family),
# the solved-neighbour control (one-reflection d=7,8), the half-target A+ at d=7,8, and deflated-LM diversity runs.
cd /tmp/hodge/results/one-reflection-skew
LOCK=chains/chain5.lock; if [ -e "$LOCK" ]; then echo "chain5: lock present ($(cat $LOCK)); exiting"; exit 0; fi; echo "$$ $(date)" > "$LOCK"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
S="python3 w2_search.py --procs 4"
$S --inst family-hesse-web --job calib-hesse-d8 --d 8 --rho 0 --sigma 0 --family kernel --restarts 160 --method lm --seed0 5800 --time-limit 3600 >> chain5.out 2>> chain5.err
$S --inst one-reflection --job ctrl-onerefl-d7 --d 7 --rho 0 --sigma 0 --family kernel --restarts 600 --method lm --seed0 27000 --time-limit 5400 >> chain5.out 2>> chain5.err
$S --job aplus-d7-r0s0 --d 7 --rho 0 --sigma 0 --family kernel --target Aplus --restarts 600 --method lm --seed0 17000 --time-limit 5400 >> chain5.out 2>> chain5.err
$S --job ker-d7-r0s0-defl --d 7 --rho 0 --sigma 0 --family kernel --restarts 400 --method defl --defl-runs 10 --seed0 7700 --time-limit 5400 >> chain5.out 2>> chain5.err
$S --inst one-reflection --job ctrl-onerefl-d8 --d 8 --rho 0 --sigma 0 --family kernel --restarts 400 --method lm --seed0 28000 --time-limit 7200 >> chain5.out 2>> chain5.err
$S --job aplus-d8-r0s0 --d 8 --rho 0 --sigma 0 --family kernel --target Aplus --restarts 400 --method lm --seed0 18000 --time-limit 7200 >> chain5.out 2>> chain5.err
$S --job ker-d8-r0s0-defl --d 8 --rho 0 --sigma 0 --family kernel --restarts 300 --method defl --defl-runs 10 --seed0 8700 --time-limit 7200 >> chain5.out 2>> chain5.err
echo CHAIN5 DONE >> chain5.out
