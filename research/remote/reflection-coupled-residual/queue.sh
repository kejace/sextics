#!/bin/bash
# wait for my current solver processes, then run jobs with concurrency 2
cd /tmp/hodge/results/reflection-coupled-residual
while pgrep -f "batch.py cands_d4" > /dev/null || pgrep -f "fermat_control.py" > /dev/null; do sleep 20; done
cat "$1" | xargs -P 2 -I{} bash -c "{}"
echo "QUEUE $1 DONE" >> queue.log
