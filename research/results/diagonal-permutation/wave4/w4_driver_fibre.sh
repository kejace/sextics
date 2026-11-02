#!/bin/bash
# robustness sweep: rank-1 locus + automorphism group + smoothness/h0 of the fibre for several (mu, prime)
cd /tmp/hodge/results/diagonal-permutation/wave4
ulimit -v 6000000
for mp in "3/7 2147483647" "2/9 1073741827" "-5/3 1073741827" "1 1073741827" "0 1073741827" "7/3 2147483647"; do
  set -- $mp; mu=$1; pr=$2; tag=$(echo $mu | tr '/' 'o' | tr '-' 'm')
  python3 w4_hess1.py fibre $mu $pr 1 > hess1_fibre_${tag}_p${pr}.log 2>&1
  python3 w4_smooth_h0.py fibre $mu $pr 1 > smooth_fibre_${tag}_p${pr}.log 2>&1
  python3 w4_aut_fibre.py $mu $pr 1 > autfib_${tag}_p${pr}.log 2>&1
  echo "done $mu $pr"
done
