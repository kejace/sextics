#!/bin/bash
cd /tmp/hodge/results/polarized-k3-remainder
export PYTHONPATH=/home/kejace/devel/sextics-work
python3 -m hodge.sweep --instance polarized-k3-remainder --candidates wave2-excl/cands_d8d10_xy.jsonl \
  --out wave2-excl/sweep-d8d10.jsonl --workers 8 --dmax 10 --timeout 14400 --witness-dir wave2-excl/witnesses 2>&1 | grep -v "engine v2"
echo SWEEP-DONE; date
