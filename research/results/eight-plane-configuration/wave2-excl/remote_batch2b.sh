#!/bin/bash
# wave-2 excl batch 2b (remote): d = 8 Hessian relaxation on the larger structured subfamilies (msolve, 8 threads each)
cd /tmp/hodge/results/eight-plane-configuration
TAG=excl-epc-batch2b
for par in 0 1; do
  for fam in arr grad blockX blockY coefX; do
    python3 wave2-excl/w2x_hess8_local.py $par $fam 8 7200 >> wave2-excl/$TAG-hsub8.out 2>&1
  done
done
