#!/bin/bash
# local batch E: does the degree-9 sector-00 landscape keep descending? one long L-BFGS run (900 iterations) from the seed that reached 0.018
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
cd /tmp/hodge/results/reflection-coupled-residual
n=0; while ps -eo args | grep -E '^python3 (kvp\.py|rcr2\.py e1) ' | grep -q . && [ $n -lt 120 ]; do sleep 5; n=$((n+1)); done
timeout 1000 python3 kvp.py --d 9 --sec 00 --family JK --restarts 1 --seed 41 --opt lbfgs --maxiter 900 --tag localE 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|stopped" | cut -c1-250
echo LOCAL_BATCH_E_DONE
