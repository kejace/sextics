#!/bin/bash
cd /tmp/hodge/results/eight-plane-configuration
while pgrep -f "[w]2_timage.py 7" > /dev/null; do sleep 20; done
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
nice -n 10 timeout 7200 python3 -u w2_timage.py 8 1 > timage8.out 2>&1
