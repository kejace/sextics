#!/bin/bash
cd /tmp/hodge/results/klein-coupled-character
while pgrep -f runB_jacobian.py > /dev/null; do sleep 5; done
python3 runC_lift_served.py > runC.out 2>&1
python3 runD_scans.py > runD.out 2>&1
echo CHAIN_DONE >> runD.out
