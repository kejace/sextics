#!/bin/bash
cd /tmp/hodge/results/reflection-coupled-residual
while pgrep -f "python3 fermat_control.py" > /dev/null; do sleep 10; done
timeout 5400 python3 e2lin.py kernel_d8_X03Y03.json 0 > e2lin_d8.out 2>&1
echo QUEUE4 DONE >> queue.log
