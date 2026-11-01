#!/bin/bash
cd /tmp/hodge/results/reflection-coupled-residual
while pgrep -f "python3 e2lin.py" > /dev/null || pgrep -f "python3 fermat_control.py" > /dev/null; do sleep 15; done
{
python3 e1rank.py 8 'X0**3*Y1**3' instance A-X0Y1
python3 e1rank.py 8 'X0*X1*X2*Y0**2*Y1' instance A-mixed
python3 e1rank.py 8 'X0**3*Y0**3' "$(cat F_coupling1.txt)" F-coupling1
python3 e1rank.py 8 'X0**3*Y0**3' "$(cat F_randsym.txt)" F-randsym
python3 e1rank.py 8 'X0**3*Y2**3' instance A-companion
} > e1rank_d8.out 2>&1
echo QUEUE5 DONE >> queue.log
