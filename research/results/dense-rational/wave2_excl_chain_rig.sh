#!/bin/bash
# Rigorous chain: E1-kernel sweep d=7..10 (both parities, p=2147483587, bases at d=10), then msolve on structured families
# (degree-7 mixed sector: pairs, gradient, xonly, yonly; degree-8 even/odd: blockdiag/cross, weuler, gradient, xonly, yonly).
cd /tmp/hodge/results/dense-rational
export PYTHONPATH=/home/kejace/devel/sextics-work OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 FLINT_NUM_THREADS=4
P2=2147483587
for d in 7 8 9; do for par in 0 1; do [ -f wave2-excl/e1kernel_d${d}_par${par}_p$P2.json ] || python3 wave2_excl_e1kernel.py $d $par $P2 >> wave2-excl/e1kernel_sweep_p$P2.log 2>&1; done; done
for par in 0 1; do [ -f wave2-excl/e1kernel_d10_par${par}_p${P2}_basis.json ] || python3 wave2_excl_e1kernel.py 10 $par $P2 --basis >> wave2-excl/e1kernel_sweep_p$P2.log 2>&1; done
echo kernels-done > wave2-excl/chain_rig.stage
P=1073741789
for fam in gradient xonly yonly pair-0-3 pair-0-1 pair-3-4 pair-0-4 pair-1-3 pair-1-2 pair-4-5 pair-0-2 pair-1-4 pair-1-5 pair-2-3 pair-2-4 pair-2-5 pair-3-5 pair-0-5; do
  bash wave2_excl_run_fam.sh 7 mixed $fam $P 8 rig-d7-mixed-$fam >> wave2-excl/fam_d7_mixed.log 2>&1
done
echo d7fam-done > wave2-excl/chain_rig.stage
for fam in blockdiag weuler gradient xonly yonly; do bash wave2_excl_run_fam.sh 8 even $fam $P2 8 rig-d8-even-$fam >> wave2-excl/fam_d8.log 2>&1; done
for fam in cross weuler gradient xonly yonly; do bash wave2_excl_run_fam.sh 8 odd $fam $P2 8 rig-d8-odd-$fam >> wave2-excl/fam_d8.log 2>&1; done
echo done > wave2-excl/chain_rig.done
