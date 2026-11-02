#!/bin/sh
while kill -0 3898750 2>/dev/null; do sleep 30; done
cd /tmp/hodge/results/eight-plane-configuration && python3 run_kernel10.py > run_kernel10.out 2>&1
