#!/bin/bash
# characteristic-0 msolve runs for the degenerate-polar locus on the fourfold X: 'on' (triple points of X) and 'off' (exchange directions)
cd /tmp/hodge/results/klein-coupled-character/w2c
for v in on off; do for k in 0 1 2 3 4 5; do
  ( timeout 7000 msolve -v 0 -t 1 -P 1 -f sigmaQ_X_${v}_c$k.ms -o sigmaQP_X_${v}_c$k.out; echo "XQ $v $k exit $?" >> run_XQ_remote.log ) &
done; done
wait
echo ALLDONE >> run_XQ_remote.log
