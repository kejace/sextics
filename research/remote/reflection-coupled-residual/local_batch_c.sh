#!/bin/bash
# local batch C: d=9 E1 spaces (sector by sector), then degree-8 sectors with L-BFGS variable projection (graceful timeouts)
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
cd /tmp/hodge/results/reflection-coupled-residual
n=0; while ps -eo comm,args | awk '$1=="python3" && /kvp.py --d/' | grep -q . && [ $n -lt 30 ]; do sleep 2; n=$((n+1)); done
for s in 00 01 10 11; do timeout 300 python3 rcr2.py e1 --d 9 --sec $s 2>&1 | grep "e1space\|coker" | cut -c1-260; done
for s in 01 10 11; do timeout 300 python3 kvp.py --d 8 --sec $s --family J --restarts 3 --seed 19 --opt lbfgs --maxiter 500 --tag localC 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|stopped" | cut -c1-250; done
timeout 300 python3 kvp.py --d 8 --sec 00 --family JK --restarts 2 --seed 20 --opt lbfgs --maxiter 500 --tag localC 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|stopped" | cut -c1-250
echo LOCAL_BATCH_C_DONE
