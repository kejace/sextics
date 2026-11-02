#!/bin/bash
# degree-11 Klein-invariant cone families (5-dim each): rational basis from 3 primes, fixed-D screens (1 basis element + 1 random)
cd /tmp/hodge/results/reflection-coupled-residual/wave2-excl
while kill -0 250447 2>/dev/null; do sleep 30; done     # wait for the current (0,0) msolve only
nice -n 10 timeout 7000 python3 conescreen.py 11 Y0,Y1 00 3 1 1 > conescreen_d11_noY0Y1.out 2>&1
nice -n 10 timeout 7000 python3 conescreen.py 11 Y1,Y2 00 3 1 1 > conescreen_d11_noY1Y2.out 2>&1
echo CONES DONE >> chain3_cones.log
