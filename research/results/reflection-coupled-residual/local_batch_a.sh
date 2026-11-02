#!/bin/bash
# local light batch (one process at a time): E1 at d=6, sparse-kvp controls, d=8 timing
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
cd /tmp/hodge/results/reflection-coupled-residual
timeout 120 python3 rcr2.py e1 --d 6 --sec all 2>&1 | grep "e1space" | cut -c1-200
echo "== control Fermat d6 sparse"
HODGE_RCR_INSTANCE='custom:fermatX03Y03' HODGE_RCR_F='X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6' HODGE_RCR_A='X0**3*Y0**3' timeout 150 python3 kvp.py --d 6 --sec 00 --family J --restarts 6 --seed 2 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|control" | cut -c1-250
echo "== control fermat-pair d4 sparse"
HODGE_RCR_INSTANCE='fermat-pair-rational-d4-strict' timeout 100 python3 kvp.py --d 4 --sec 00 --family K --restarts 12 --seed 3 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|control" | cut -c1-250
echo "== our d8 sec00 JK sparse (timing, 3 restarts)"
timeout 170 python3 kvp.py --d 8 --sec 00 --family JK --restarts 3 --seed 11 --tag local 2>&1 | grep -v "engine v2" | grep "kvp\]\|confirm\]\|candidate" | cut -c1-250
echo LOCAL_BATCH_A_DONE
