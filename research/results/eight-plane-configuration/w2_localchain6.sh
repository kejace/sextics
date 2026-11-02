#!/bin/bash
cd /tmp/hodge/results/eight-plane-configuration
while [ ! -f w2_Cproj_d8.npy ]; do sleep 30; done
sleep 60
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
nice -n 5 python3 -u w2_jvarpro.py --family jacobian --d 8 --par 0 --restarts 20 --seed 41 --nproc 1 --stage1-only 1 --stage1-nfev 200 --tag local-stage1-d8-p0 > local-stage1-d8-p0.out 2>&1
nice -n 5 python3 -u w2_jvarpro.py --family jacobian --d 8 --par 1 --restarts 20 --seed 41 --nproc 1 --stage1-only 1 --stage1-nfev 200 --tag local-stage1-d8-p1 > local-stage1-d8-p1.out 2>&1
