#!/bin/bash
# follow-up chain (lens excl, lane cycle-q2): resume sectors d=6, then sectors d=7, then cones d=7..9
cd /tmp/hodge/results/cycle-q2
export W2JOB=excl-cycle-q2-chain2
python3 wave2-excl/sectors.py --deg 6 --workers 6 --random 2 >> wave2-excl/sectors6.out 2>&1
python3 wave2-excl/sectors.py --deg 7 --workers 8 --random 2 > wave2-excl/sectors7.out 2>&1
for d in 7 8 9; do python3 wave2-excl/cones.py --deg $d --workers 8 --random 1 > wave2-excl/cones$d.out 2>&1; done
echo CHAIN2 DONE >> wave2-excl/sectors7.out
