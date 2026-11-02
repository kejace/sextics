#!/bin/bash
# Structured sub-families (small K): pencils/grids are exhaustive, others get many restarts. Resumable.
cd /tmp/hodge/results/four-plus-two-residual
S="python3 wave2-search/fpt_search.py --procs 8 --threads 1"
# cones with vertex a plane (3 core variables) and the core line (4 core variables)
for V in X0,X1,Y0 X0,X1,Y1 X0,Y0,Y1 X1,Y0,Y1 X0,X1,X2 X0,X1,Y2 Y0,Y1,X2 Y0,Y1,Y2 X0,Y0,X2 X0,Y0,Y2 X1,Y1,X2 X1,Y1,Y2 X0,Y1,X2 X0,Y1,Y2 X1,Y0,X2 X1,Y0,Y2 X0,X2,Y2 X1,X2,Y2 Y0,X2,Y2 Y1,X2,Y2; do
  $S --job vars3 --family vars:$V --degrees 6,7,8,9,10 --restarts 20 --batches 10 --seed0 3000
done
for V in X0,X1,Y0,Y1 X0,X1,Y0,X2 X0,X1,Y0,Y2 X0,X1,Y1,X2 X0,X1,Y1,Y2 X0,Y0,Y1,X2 X0,Y0,Y1,Y2 X1,Y0,Y1,X2 X1,Y0,Y1,Y2 X0,X1,X2,Y2 Y0,Y1,X2,Y2 X0,Y0,X2,Y2 X0,Y1,X2,Y2 X1,Y0,X2,Y2 X1,Y1,X2,Y2; do
  $S --job vars4 --family vars:$V --degrees 6,7,8 --restarts 20 --batches 8 --seed0 4000
done
# Jacobian-type w-space and its Hessian relaxation
$S --job hess --family hess --degrees 6,7,8,9,10 --restarts 50 --batches 20 --seed0 5000
$S --job Jw --family J --degrees 6,7,8 --restarts 20 --batches 10 --seed0 6000
# products with fixed binary factors
for P in X2 Y2 "X2*Y2" "X2**2" "Y2**2" "X2**2+Y2**2" "X2**2-Y2**2"; do
  $S --job prod --family "prod:$P" --degrees 7,8 --restarts 20 --batches 8 --seed0 7000
done
echo "[structured] all done"
