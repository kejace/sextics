#!/bin/sh
cd /tmp/hodge/results/diagonal-permutation
while pgrep -f "D8_images" >/dev/null; do sleep 10; done
python3 run_D8.py > run_D8.out 2>&1
python3 e1kernel.py 9 9 > e1kernel_d9.out 2>&1
