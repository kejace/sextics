#!/bin/bash
cd /tmp/hodge/results/klein-coupled-character
# wait until runA reaches d=10, then kill it and start the reduced d=10 run
until grep -q "d=10 parD=0 k=0" runA.out 2>/dev/null || ! pgrep -f runA_kernel.py > /dev/null; do sleep 10; done
pkill -f runA_kernel.py; sleep 2
python3 runA2_d10.py > runA2.out 2>&1
echo A2_DONE >> runA2.out
