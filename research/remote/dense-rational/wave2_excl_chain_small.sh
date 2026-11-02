#!/bin/bash
# degree-6 Hessian relaxation in all three sectors (tiny msolve) + E1-kernel sweep d=1..10 both parities with prime 2147483587
cd /tmp/hodge/results/dense-rational
export PYTHONPATH=/home/kejace/devel/sextics-work
P=1073741789
for sec in even odd mixed; do bash wave2_excl_run_quadsys.sh 6 $sec $P 8 excl-d6-$sec; done
P2=2147483587
for d in 1 2 3 4 5 6 7 8 9; do for par in 0 1; do python3 wave2_excl_e1kernel.py $d $par $P2 >> wave2-excl/e1kernel_sweep_p$P2.log 2>&1; done; done
for par in 0 1; do python3 wave2_excl_e1kernel.py 10 $par $P2 --basis >> wave2-excl/e1kernel_sweep_p$P2.log 2>&1; done
echo done > wave2-excl/chain_small.done
