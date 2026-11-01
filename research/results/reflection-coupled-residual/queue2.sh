#!/bin/bash
cd /tmp/hodge/results/reflection-coupled-residual
while pgrep -f "k8.py 8" > /dev/null || pgrep -f "calib.py" > /dev/null; do sleep 10; done
( timeout 3600 python3 sweep_fast.py 'X0**3*Y0**3' 9 11 A > sweep_fast_A.out 2>&1 ) &
( timeout 1800 python3 fermat_control.py > fermat_control.out 2>&1 ) &
wait
echo QUEUE2 DONE >> queue.log
