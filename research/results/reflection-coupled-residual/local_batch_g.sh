#!/bin/bash
# local batch G: 1500-iteration L-BFGS on the degree-9 sector-00 family from seed 41 (keeps the best vector for warm starts)
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
cd /tmp/hodge/results/reflection-coupled-residual
n=0; while ps -eo args | grep -E '^python3 (kvp\.py|rcr2\.py e1) ' | grep -q . && [ $n -lt 360 ]; do sleep 5; n=$((n+1)); done
timeout 3000 python3 kvp.py --d 9 --sec 00 --family JK --restarts 1 --seed 41 --opt lbfgs --maxiter 1500 --tag localG 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|stopped" | cut -c1-250
echo LOCAL_BATCH_G_DONE
