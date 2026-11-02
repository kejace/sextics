#!/bin/bash
# numerical screens: degree-7 relaxation (mixed sector = the only one not decided by msolve; even/odd as cross-checks), degree-8 pure sectors
cd /tmp/hodge/results/dense-rational
export PYTHONPATH=/home/kejace/devel/sextics-work OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
python3 wave2_excl_lm_d8.py even 2 0 > wave2-excl/lm_d8_prep.log 2>&1      # builds the cached complement of J_16^-
for s in 1 2 3 4 5 6; do python3 wave2_excl_lm7_wrap.py 400 $s mixed > wave2-excl/lm7_mixed_s$s.log 2>&1 & done
python3 wave2_excl_lm7_wrap.py 200 11 even > wave2-excl/lm7_even_s11.log 2>&1 &
python3 wave2_excl_lm7_wrap.py 200 12 odd > wave2-excl/lm7_odd_s12.log 2>&1 &
for s in 1 2 3 4; do python3 wave2_excl_lm_d8.py even 150 $s > wave2-excl/lm_d8_even_s$s.log 2>&1 & done
for s in 1 2 3 4; do python3 wave2_excl_lm_d8.py odd 150 $s > wave2-excl/lm_d8_odd_s$s.log 2>&1 & done
wait
echo done > wave2-excl/chain_lm.done
