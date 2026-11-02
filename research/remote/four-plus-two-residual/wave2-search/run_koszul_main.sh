#!/bin/bash
# Main kernel+J joint search, Koszul (E2-only) formulation, LM optimizer; resumable. Order: d=6 then 7 then 8.
cd /tmp/hodge/results/four-plus-two-residual
S="python3 wave2-search/fpt_search.py --procs 8 --threads 1 --formulation koszul --method lm --maxiter 200"
$S --job kK6 --family kernelJ --degrees 6 --restarts 50 --batches 40 --seed0 11000
$S --job kK7 --family kernelJ --degrees 7 --restarts 50 --batches 40 --seed0 12000
$S --job kK8 --family kernelJ --degrees 8 --restarts 25 --batches 40 --seed0 13000
echo "[koszul-main] done"
