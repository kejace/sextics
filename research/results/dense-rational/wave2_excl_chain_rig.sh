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
# rigorous served-target test by E1 alone (cheap) for ~60 named denominators
python3 wave2_excl_e1dual.py $P2 8 > wave2-excl/e1dual_p$P2.log 2>&1
fams="gradient xonly yonly"
for i in 0 1 2 3 4 5; do for j in 0 1 2 3 4 5; do [ $i -lt $j ] && fams="$fams pair-$i-$j"; done; done
for i in 0 1 2 3 4 5; do for j in 0 1 2 3 4 5; do for k in 0 1 2 3 4 5; do [ $i -lt $j ] && [ $j -lt $k ] && fams="$fams comp-$i-$j-$k"; done; done; done
for fam in $fams; do
  [ -f wave2-excl/fam_d7_mixed_${fam}_p${P}.g1.out ] && continue
  bash wave2_excl_run_fam.sh 7 mixed $fam $P 8 rig-d7-mixed-$fam >> wave2-excl/fam_d7_mixed.log 2>&1
done
echo d7fam-done > wave2-excl/chain_rig.stage
for fam in blockdiag weuler gradient xonly yonly; do bash wave2_excl_run_fam.sh 8 even $fam $P2 8 rig-d8-even-$fam >> wave2-excl/fam_d8.log 2>&1; done
for fam in cross weuler gradient xonly yonly; do bash wave2_excl_run_fam.sh 8 odd $fam $P2 8 rig-d8-odd-$fam >> wave2-excl/fam_d8.log 2>&1; done
echo d8fam-done > wave2-excl/chain_rig.stage
# modular / numerical screens (not proofs)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLINT_NUM_THREADS=1
python3 wave2_excl_timage.py $P2 8 > wave2-excl/timage_p$P2.log 2>&1
python3 wave2_excl_deform.py $P2 8 > wave2-excl/deform_p$P2.log 2>&1
for par in 0 1; do
  b=wave2-excl/e1kernel_d10_par${par}_p${P2}_basis.json; [ -f $b ] && python3 wave2_excl_kstruct.py $b > wave2-excl/kstruct_d10_par$par.log 2>&1
  r=wave2-excl/e1kernel_d10_par${par}_p${P2}_reps.json; [ -f $r ] && python3 wave2_excl_e2reps.py $r 8 > wave2-excl/e2reps_d10_par$par.log 2>&1
done
echo screens-done > wave2-excl/chain_rig.stage
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
python3 wave2_excl_lm_d8.py even 2 0 > wave2-excl/lm_d8_prep.log 2>&1
for s in 1 2 3; do python3 wave2_excl_lm7_wrap.py 300 $s mixed > wave2-excl/lm7_mixed_s$s.log 2>&1 & done
for s in 1 2; do python3 wave2_excl_lm_d8.py even 100 $s > wave2-excl/lm_d8_even_s$s.log 2>&1 & done
for s in 1 2; do python3 wave2_excl_lm_d8.py odd 100 $s > wave2-excl/lm_d8_odd_s$s.log 2>&1 & done
wait
echo done > wave2-excl/chain_rig.done
