#!/bin/bash
# wave-2 excl batch 2 (remote): sampling d=7,8,9 + exact Y-partials + d=8 subfamily Hessian relaxations with msolve
cd /tmp/hodge/results/eight-plane-configuration
TAG=excl-epc-batch2
python3 wave2-excl/w2x_koszul.py exact $TAG partials > wave2-excl/$TAG-partials-exact.out 2>&1
python3 wave2-excl/w2x_jsample.py 7 8 $TAG > wave2-excl/$TAG-jsample7.out 2>&1
python3 wave2-excl/w2x_jsample.py 8 8 $TAG > wave2-excl/$TAG-jsample8.out 2>&1
for par in 0 1; do
  for fam in arr grad blockX blockY coefX coefB0 ham:K0 ham:H0 ham:diag0 ham:Q0Q1 ham:Q0^2; do
    python3 wave2-excl/w2x_hess_sub.py 8 $par $fam 2147483587 $TAG >> wave2-excl/$TAG-hsub8.out 2>&1
  done
done
for f in wave2-excl/hsub_d8_par*_p2147483587.ms; do
  b=$(basename $f); bash wave2-excl/run_msolve.sh $b 8 $TAG >> wave2-excl/$TAG-hsub8-msolve.out 2>&1
done
python3 wave2-excl/w2x_jsample.py 9 8 $TAG > wave2-excl/$TAG-jsample9.out 2>&1
