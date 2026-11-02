#!/bin/bash
# wave-2 excl batch 2b (remote): d = 8 Hessian relaxation on the larger structured subfamilies (msolve, 8 threads each; dreg 4 feasible ones)
cd /tmp/hodge/results/eight-plane-configuration
TAG=excl-epc-batch2b
for par in 0 1; do
  for fam in coefX arr grad; do
    python3 wave2-excl/w2x_hess8_local.py $par $fam 8 14400 >> wave2-excl/$TAG-hsub8.out 2>&1
  done
done
python3 wave2-excl/w2x_served.py 6 2 $TAG > wave2-excl/$TAG-served6.out 2>&1
python3 wave2-excl/w2x_served.py 7 1 $TAG > wave2-excl/$TAG-served7.out 2>&1
for par in 0 1; do for s in 1 2 3 4; do python3 wave2-excl/w2x_hess8_local.py $par hamrand:$s 4 600 >> wave2-excl/$TAG-hsub8.out 2>&1; done; done
