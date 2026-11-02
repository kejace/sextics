#!/bin/bash
cd /tmp/hodge/results/cyclic-triple-residual
echo "[$(date +%T)] d6 even k=5 restricted (vectorized)"
nice -n 19 python3 -u w2_macblocks.py 6 0 5 1073741827 restrict 2>&1 | grep --line-buffered -v "engine v2" | head -2
echo "[$(date +%T)] d6 even point identification"
nice -n 19 python3 -u w2_point_d6.py 40 2>&1 | grep --line-buffered -v "engine v2"
for par in 0 1; do
  echo "[$(date +%T)] d7 par $par k=3 restricted"
  nice -n 19 python3 -u w2_macblocks.py 7 $par 3 1073741827 restrict 2>&1 | grep --line-buffered -v "engine v2" | head -2
done
echo "[$(date +%T)] LOCAL_CHAIN_DONE"
