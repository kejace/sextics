#!/bin/bash
# local batch D: Koszul varpro over the degree-9 sector-00 family (21 kernel classes + J_9), one restart per process
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
cd /tmp/hodge/results/reflection-coupled-residual
n=0; while ps -eo args | grep -E '^python3 (kvp\.py|rcr2\.py e1) ' | grep -q . && [ $n -lt 60 ]; do sleep 5; n=$((n+1)); done
for seed in 41 42; do timeout 420 python3 kvp.py --d 9 --sec 00 --family JK --restarts 1 --seed $seed --opt lbfgs --maxiter 300 --tag localD 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|stopped" | cut -c1-250; done
timeout 420 python3 kvp.py --d 9 --sec 00 --family JK --restarts 1 --seed 43 --opt lm --maxiter 60 --tag localD 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|stopped" | cut -c1-250
for seed in 44 45; do timeout 420 python3 kvp.py --d 9 --sec 10 --family JK --restarts 1 --seed $seed --opt lbfgs --maxiter 300 --tag localD 2>&1 | grep -v "engine v2" | grep "kvp\] d=\|summary\|confirm\]\|candidate\|stopped" | cut -c1-250; done
echo LOCAL_BATCH_D_DONE
