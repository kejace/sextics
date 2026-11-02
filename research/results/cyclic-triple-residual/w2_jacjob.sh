#!/bin/bash
# rigorous relaxation exclusion (block Macaulay ranks, then msolve leading ideal as a cross-check) for d = 5,6,7, both
# parities; then the numerical Jacobian-type searches
cd /tmp/hodge/results/cyclic-triple-residual
for d in 5 6 7; do for par in 0 1; do
  python3 w2_macblocks.py $d $par 3 > w2_macblocks_d${d}_p${par}.out 2>&1
done; done
for d in 6 7; do for par in 0 1; do
  python3 w2_macaulay.py $d $par > w2_mac_d${d}_p${par}.out 2>&1
  timeout 7200 msolve -t 8 -g 1 -f w2-quadrics-d${d}-par${par}.ms -o w2-quadrics-d${d}-par${par}.lm > w2_msolve_d${d}_p${par}.out 2>&1
  python3 w2_hilb.py w2-quadrics-d${d}-par${par}.lm $d $par >> w2_mac_d${d}_p${par}.out 2>&1
done; done
python3 w2_jac.py 8 1.0 > w2_jac.out 2> w2_jac.err
