#!/bin/bash
cd /tmp/hodge/results/dense-rational
P=1073741827
while pgrep -f "polar_d7_p1073741827.ms" > /dev/null; do sleep 15; done
echo "exit (detached run finished)" >> polar_d7_msolve.log
timeout 21600 msolve -v 2 -t 8 -f polar_d7_pd1_slice_p$P.ms -o polar_d7_pd1_slice.out > polar_d7_pd1_slice_msolve.log 2>&1; echo "exit $?" >> polar_d7_pd1_slice_msolve.log
timeout 21600 msolve -v 2 -t 8 -f polar_d7_pd0_homog_p$P.ms -o polar_d7_pd0_homog.out > polar_d7_pd0_homog_msolve.log 2>&1; echo "exit $?" >> polar_d7_pd0_homog_msolve.log
timeout 21600 msolve -v 2 -t 8 -f polar_d7_pd1_homog_p$P.ms -o polar_d7_pd1_homog.out > polar_d7_pd1_homog_msolve.log 2>&1; echo "exit $?" >> polar_d7_pd1_homog_msolve.log
echo done > chain_g.done
timeout 10800 python3 e1kernel2.py 10 0 > e1_direct2_10_par0.log 2>&1; echo done > chain_g2.done
for i in 0 1 2; do timeout 7200 python3 e2norm.py e1kernel_d10_par0_reps.json $i >> e2norm_d10.log 2>&1; done; echo done > chain_g3.done
timeout 3600 python3 kstruct10.py e1kernel_d10_par0_basis.json > kstruct10.log 2>&1
timeout 10800 python3 e1kernel2.py 10 1 > e1_direct2_10_par1.log 2>&1
