#!/bin/bash
cd /tmp/hodge/results/klein-coupled-character/wave2-excl
for spec in "d7_even_k1 7 0 1" "d7_even_k2 7 0 2" "d7_even_k4 7 0 4" "d7_odd_k0 7 1 0"; do
  set -- $spec
  nice -n 10 timeout 14400 python3 -u w2_radical_exact.py $1 $2 $3 $4 3 > local_radical_exact_$1.out 2>&1
done
echo ALLX_DONE >> local_radical_exact_d7_odd_k0.out
