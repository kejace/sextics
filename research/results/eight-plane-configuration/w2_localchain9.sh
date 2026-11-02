#!/bin/bash
cd /tmp/hodge/results/eight-plane-configuration
until grep -q '"event": "end"' local-stage1-d8-p1.jsonl 2>/dev/null; do sleep 30; done
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
nice -n 5 python3 -u w2_continue.py local-stage1-d8-p1 1 8 800 > local-continue-d8-p1.out 2>&1
