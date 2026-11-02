#!/bin/bash
# dims-only survey of E1-admissible sectors and coordinate-cone sub-families, degrees 8..12 (no element tests, no relaxation for d>=10)
cd /tmp/hodge/results/cycle-q2
export W2JOB=excl-cycle-q2-survey W2_NO_ELEMENTS=1
# second-prime confirmation of the degree-6 relaxation structure
python3 wave2-excl/hess6.py --deg 6 --prime 998244353 > wave2-excl/hess6_p2.out 2>&1
for s in 1 2; do python3 wave2-excl/subsets6.py --size $s --prime 998244353 --batch 1000 > wave2-excl/subsets6_${s}_p2.out 2>&1; done
for d in 8 9; do python3 wave2-excl/sectors.py --deg $d --workers 6 --no-farkas --random 0 --cones > wave2-excl/survey$d.out 2>&1; done
export W2_NO_RELAX=1
for d in 10 11 12; do python3 wave2-excl/sectors.py --deg $d --workers 6 --no-farkas --random 0 --cones > wave2-excl/survey$d.out 2>&1; done
echo SURVEY DONE >> wave2-excl/survey12.out
