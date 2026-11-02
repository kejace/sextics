#!/bin/bash
# Full parametric sector systems at degree d for all 18 sectors (modular evidence).  usage: bash bilin_sectors.sh d p threads timeout_s
d=$1; p=$2; T=$3; TMO=${4:-3600}
cd /tmp/hodge/results/cyclic-triple-residual/wave2-excl
for c0 in 0 2 4; do for c1 in 0 1 2 3 4 5; do
  timeout $((TMO+600)) python3 bilin.py $d $p $T $c0,$c1 $TMO > bilin_d${d}_${c0}${c1}.out 2> bilin_d${d}_${c0}${c1}.err
  tail -1 bilin_d${d}_${c0}${c1}.out
done; done
