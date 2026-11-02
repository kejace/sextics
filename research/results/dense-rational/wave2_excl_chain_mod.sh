#!/bin/bash
# modular screens after the kernel sweep: served-target spaces, Fermat first-order deformation, kernel structure, E2 on degree-10 reps
cd /tmp/hodge/results/dense-rational
export PYTHONPATH=/home/kejace/devel/sextics-work OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
P2=2147483587
python3 wave2_excl_timage.py $P2 12 > wave2-excl/timage_p$P2.log 2>&1
python3 wave2_excl_deform.py $P2 12 > wave2-excl/deform_p$P2.log 2>&1
for par in 0 1; do
  b=wave2-excl/e1kernel_d10_par${par}_p${P2}_basis.json
  [ -f $b ] && python3 wave2_excl_kstruct.py $b > wave2-excl/kstruct_d10_par$par.log 2>&1
  r=wave2-excl/e1kernel_d10_par${par}_p${P2}_reps.json
  [ -f $r ] && python3 wave2_excl_e2reps.py $r 10 > wave2-excl/e2reps_d10_par$par.log 2>&1
done
echo done > wave2-excl/chain_mod.done
