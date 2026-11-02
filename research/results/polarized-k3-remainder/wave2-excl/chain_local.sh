#!/bin/bash
# wait for the running d=8 exact test, then run the K3-target d=7 relaxation locally (one process at a time, nice 19)
while kill -0 254171 2>/dev/null; do sleep 20; done
cd /tmp/hodge/results/polarized-k3-remainder/wave2-excl
echo "=== K3 d=7 start $(date)"
nice -n 19 python3 relax.py 7 65521 --target=k3 --exact-ranks --squares --out=relax-k3-d7-p65521-local 2>&1 | grep -v "engine v2"
echo "=== K3 d=7 done $(date)"
