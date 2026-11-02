#!/bin/bash
# d=7 rigorous relaxation exclusion (restricted to the complement of the F-gauge locus), degree-3 Macaulay blocks, both parities
cd /tmp/hodge/results/cyclic-triple-residual
for par in 0 1; do
  nice -n 19 python3 -u w2_macblocks.py 7 $par 3 1073741827 restrict 2>&1 | grep -v "engine v2" | head -2
done
echo LOCAL_D7_DONE
