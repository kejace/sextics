#!/bin/bash
# Structured sub-families with the Koszul formulation (LM). Resumable.
cd /tmp/hodge/results/four-plus-two-residual
S="python3 wave2-search/fpt_search.py --procs 8 --threads 1 --formulation koszul --method lm --maxiter 200"
for V in X0,X1,Y0,Y1 X0,X1,Y0,X2 X0,X1,Y0,Y2 X0,X1,Y1,X2 X0,X1,Y1,Y2 X0,Y0,Y1,X2 X0,Y0,Y1,Y2 X1,Y0,Y1,X2 X1,Y0,Y1,Y2 X0,X1,X2,Y2 Y0,Y1,X2,Y2 X0,Y0,X2,Y2 X0,Y1,X2,Y2 X1,Y0,X2,Y2 X1,Y1,X2,Y2; do
  $S --job kvars4 --family vars:$V --degrees 6,7,8,9,10 --restarts 25 --batches 20 --seed0 14000
done
for V in X0,X1,Y0 X0,X1,Y1 X0,Y0,Y1 X1,Y0,Y1 X0,X1,X2 X0,X1,Y2 Y0,Y1,X2 Y0,Y1,Y2 X0,Y0,X2 X0,Y0,Y2 X1,Y1,X2 X1,Y1,Y2 X0,Y1,X2 X0,Y1,Y2 X1,Y0,X2 X1,Y0,Y2 X0,X2,Y2 X1,X2,Y2 Y0,X2,Y2 Y1,X2,Y2; do
  $S --job kvars3 --family vars:$V --degrees 6,7,8,9,10,11,12 --restarts 25 --batches 10 --seed0 15000
done
$S --job kJw --family J --degrees 6,7,8,9 --restarts 25 --batches 20 --seed0 16000
for P in X2 Y2 "X2*Y2" "X2**2" "Y2**2" "X2**2+Y2**2" "X2**2-Y2**2"; do
  $S --job kprod --family "prod:$P" --degrees 7,8,9 --restarts 25 --batches 12 --seed0 17000
done
python3 wave2-search/fpt_search.py --procs 8 --threads 1 --job hess --family hess --degrees 6,7,8,9,10 --restarts 50 --batches 20 --seed0 5000
echo "[koszul-struct] done"
