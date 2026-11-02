#!/bin/bash
cd /tmp/hodge/results/eight-plane-configuration
while pgrep -f w2_localchain1.sh > /dev/null; do sleep 15; done
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -u w2_jvarpro.py --family jacobian --d 6 --par 0 --restarts 60 --seed 21 --nproc 1 --max-nfev 120 --tag local-search-d6-p0 > local-search-d6-p0.out 2>&1
python3 -u w2_jvarpro.py --family jacobian --d 6 --par 1 --restarts 60 --seed 21 --nproc 1 --max-nfev 120 --tag local-search-d6-p1 > local-search-d6-p1.out 2>&1
python3 -u w2_jvarpro.py --family jacobian --d 7 --par 0 --restarts 16 --seed 21 --nproc 1 --max-nfev 150 --stage1 1 --tag local-search-d7-p0 > local-search-d7-p0.out 2>&1
python3 -u w2_jvarpro.py --family jacobian --d 7 --par 1 --restarts 16 --seed 21 --nproc 1 --max-nfev 150 --stage1 1 --tag local-search-d7-p1 > local-search-d7-p1.out 2>&1
