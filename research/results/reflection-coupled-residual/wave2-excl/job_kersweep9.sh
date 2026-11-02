#!/bin/bash
# fixed-D rigorous tests of the degree-9 non-Jacobian E1-kernel representatives (all characters with nonzero kernel) + random combinations
cd /tmp/hodge/results/reflection-coupled-residual/wave2-excl
export W2LIB_SUFFIX=.r2
python3 e1char.py 9 9 2147483587 2147483579 > e1char_d9_rerun.out 2>&1
for c in 00 10 01 11; do
  if [ -f e1reps_d9_chi${c}_Q.json ]; then python3 kersweep.py 9 $c 6 >> kersweep_d9.out 2>&1; fi
done
echo KERSWEEP9 DONE >> job_kersweep9.log
