#!/bin/bash
cd /tmp/hodge/results/dense-rational
export PYTHONPATH=/home/kejace/devel/sextics-work
P2=2147483587
for d in 11 12; do for par in 0 1; do python3 wave2_excl_e1kernel.py $d $par $P2 --basis >> wave2-excl/e1kernel_sweep_d1112_p$P2.log 2>&1; done; done
echo done > wave2-excl/chain_k1112.done
