#!/bin/bash
# local batch F: long L-BFGS run on the degree-9 sector-10 family (21 kernel classes + J_9)
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
cd /tmp/hodge/results/reflection-coupled-residual
n=0; while ps -eo args | grep -E '^python3 (kvp\.py|rcr2\.py e1) ' | grep -q . && [ $n -lt 240 ]; do sleep 5; n=$((n+1)); done
timeout 1000 python3 kvp.py --d 9 --sec 10 --family JK --restarts 1 --seed 51 --opt lbfgs --maxiter 900 --tag localF 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|stopped" | cut -c1-250
echo LOCAL_BATCH_F_DONE
