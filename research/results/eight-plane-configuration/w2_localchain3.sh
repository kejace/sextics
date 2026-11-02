#!/bin/bash
cd /tmp/hodge/results/eight-plane-configuration
while pgrep -f "[w]2_localchain1.sh" > /dev/null; do sleep 15; done
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -u w2_jvarpro.py --family jacobian --d 7 --par 0 --restarts 12 --seed 31 --nproc 1 --stage1-only 1 --stage1-nfev 200 --tag local-stage1b-d7-p0 > local-stage1b-d7-p0.out 2>&1
python3 -u w2_continue.py local-stage1b-d7-p0 0 7 600 > local-continue-d7-p0.out 2>&1
