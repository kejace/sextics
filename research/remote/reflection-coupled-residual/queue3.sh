#!/bin/bash
cd /tmp/hodge/results/reflection-coupled-residual
while pgrep -f "python3 k8.py 8" > /dev/null; do sleep 10; done
K8=$(python3 -c "import json;print(json.load(open('kernel_d8_X03Y03.json'))[0])")
echo "K8 = $K8" > k8_expr.txt
( timeout 7200 python3 testD.py E1-kernel-d8 "$K8" > testD_K8.out 2>&1 ) &
while pgrep -f "python3 calib.py" > /dev/null; do sleep 10; done
( timeout 1800 python3 fermat_control.py > fermat_control.out 2>&1 ) &
wait
( timeout 3600 python3 sweep_fast.py 'X0**3*Y0**3' 9 11 A > sweep_fast_A.out 2>&1 ) &
wait
echo QUEUE3 DONE >> queue.log
