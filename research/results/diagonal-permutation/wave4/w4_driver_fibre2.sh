#!/bin/bash
# second-prime sweep with p = 1000000009 (= 1 mod 6); 2^31-1 gave empty msolve output
cd /tmp/hodge/results/diagonal-permutation/wave4
ulimit -v 6000000
for mp in "3/7 1000000009" "7/3 1000000009" "-5/3 1000000009"; do
  set -- $mp; mu=$1; pr=$2; tag=$(echo $mu | tr '/' 'o' | tr '-' 'm')
  python3 w4_hess1.py fibre $mu $pr 1 > hess1_fibre_${tag}_p${pr}.log 2>&1
  python3 w4_smooth_h0.py fibre $mu $pr 1 > smooth_fibre_${tag}_p${pr}.log 2>&1
  python3 w4_aut_fibre.py $mu $pr 1 > autfib_${tag}_p${pr}.log 2>&1
  echo "done $mu $pr"
done
