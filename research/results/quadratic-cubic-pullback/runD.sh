#!/bin/bash
ulimit -v 14000000
cd /tmp/hodge/results/quadratic-cubic-pullback
while pgrep -f "run_list2.py ours listC" > /dev/null; do sleep 30; done
exec python3 run_list2.py ours listD_ours.txt
