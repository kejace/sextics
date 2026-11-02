#!/bin/bash
# (1) unions of 3 and 4 sectors; (2) fast remaining re-verifications (deg <= 6); (3) msolve leading ideal of the 62-variable
# relaxation ideal; (4) remaining re-verifications (deg 10/12).  Resumable pieces: reverify skips finished tags.
cd /tmp/hodge/results/cycle-q2
export W2JOB=excl-cycle-q2-gb
python3 wave2-excl/subsets6.py --size 3 --batch 1000 > wave2-excl/subsets6_3.out 2>&1
python3 wave2-excl/subsets6.py --size 4 --batch 3000 > wave2-excl/subsets6_4.out 2>&1
python3 wave2-excl/reverify.py --workers 8 --max-d 6 >> wave2-excl/reverify.out 2>&1
python3 wave2-excl/gb_check.py wave2-excl/hess6_p1000000007.ms --threads 24 --gmode 1 > wave2-excl/hess6gb.out 2>&1
python3 wave2-excl/reverify.py --workers 8 >> wave2-excl/reverify.out 2>&1
echo GB CHAIN DONE >> wave2-excl/hess6gb.out
