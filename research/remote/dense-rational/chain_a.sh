#!/bin/bash
# wait for direct 1..7 to finish, then run direct 8, then direct 9
cd /tmp/hodge/results/dense-rational
while pgrep -f "e1kernel.py direct 1 2" > /dev/null; do sleep 5; done
python3 e1kernel.py direct 8 > e1_direct_8.log 2>&1
python3 e1kernel.py direct 9 > e1_direct_9.log 2>&1
