#!/bin/bash
cd /tmp/hodge/results/dense-rational
while pgrep -f "polar_d7_p1073741827.ms" > /dev/null; do sleep 10; done
P=1073741827
timeout 3600 msolve -v 2 -t 8 -f polar_d7_pd1_slice_p$P.ms -o polar_d7_pd1_slice.out > polar_d7_pd1_slice_msolve.log 2>&1; echo "exit $?" >> polar_d7_pd1_slice_msolve.log
timeout 5400 msolve -v 2 -t 8 -f polar_d7_pd0_homog_p$P.ms -o polar_d7_pd0_homog.out > polar_d7_pd0_homog_msolve.log 2>&1; echo "exit $?" >> polar_d7_pd0_homog_msolve.log
timeout 5400 msolve -v 2 -t 8 -f polar_d7_pd1_homog_p$P.ms -o polar_d7_pd1_homog.out > polar_d7_pd1_homog_msolve.log 2>&1; echo "exit $?" >> polar_d7_pd1_homog_msolve.log
echo done > chain_g.done
bash /tmp/hodge/results/dense-rational/chain_g2_tail.sh
bash /tmp/hodge/results/dense-rational/chain_g3_tail.sh
