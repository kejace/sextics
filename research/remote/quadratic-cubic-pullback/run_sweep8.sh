#!/bin/bash
ulimit -v 14000000
cd /tmp/hodge/results/quadratic-cubic-pullback
# wait for delta_served to finish (keep <= 2 solver processes)
while pgrep -f delta_served.py > /dev/null; do sleep 20; done
exec python3 sweep.py 8 8 ours
