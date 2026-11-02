#!/bin/bash
cd /tmp/hodge/results/eight-plane-configuration
while [ ! -f w2_Cproj_d8.npy ]; do sleep 30; done
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
nice -n 10 timeout 5400 python3 -u w2_hessmixed.py 7 > hessmixed7.out 2>&1
timeout 1200 python3 w2_sing_hilb.py hessmixed-quadrics-d7.ms > hessmixed7-hilb.out 2>&1
