#!/bin/bash
# rigorous relaxation exclusion (block Macaulay ranks on the complement of the F-gauge locus) for d = 6, 7, both parities,
# msolve cross-check at d = 6, then the numerical Jacobian-type searches (Fermat calibration first)
cd /tmp/hodge/results/cyclic-triple-residual
for par in 0 1; do
  python3 w2_macblocks.py 7 $par 3 1073741827 restrict > w2_macblocks_d7_p${par}_k3.out 2>&1
  python3 w2_macblocks.py 6 $par 5 1073741827 restrict > w2_macblocks_d6_p${par}_k5.out 2>&1
done
for par in 0 1; do
  python3 w2_macaulay.py 6 $par > w2_mac_d6_p${par}.out 2>&1
  timeout 3600 msolve -t 8 -g 1 -f w2-quadrics-d6-par${par}.ms -o w2-quadrics-d6-par${par}.lm > w2_msolve_d6_p${par}.out 2>&1
  python3 w2_hilb.py w2-quadrics-d6-par${par}.lm 6 $par >> w2_mac_d6_p${par}.out 2>&1
done
python3 w2_jac.py 8 1.0 > w2_jac.out 2> w2_jac.err
# second prime for the d=7 block ranks (independent confirmation)
for par in 0 1; do python3 w2_macblocks.py 7 $par 3 2147483647 restrict > w2_macblocks_d7_p${par}_k3_p2.out 2>&1; done
