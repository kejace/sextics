#!/bin/bash
# Positive control of the relaxation pipeline on the Fermat sextic (same A = X0^3 Y0^3): D0 = X0^6+X1^6+X2^6 is a known
# certificate denominator there (DiamondJM control) and lies in J_6 (g = (X0,X1,X2,0,0,0)/6), so the degree-6 quadratic
# system on W must have a NONTRIVIAL solution cone (dim >= 1).  A zero-dimensional result would reveal a bug.
cd /tmp/hodge/results/reflection-coupled-residual/wave2-excl
export W2LIB_F='X0**6 + X1**6 + X2**6 + Y0**6 + Y1**6 + Y2**6' W2LIB_TAG=fermat
P=2147483587
python3 hess.py export 6 $P 00 01 10 11 > control-fermat/export_d6.out 2>&1
python3 gbcheck.py control-fermat/hess_d6_W00011011_p${P}.ms 2 1200 8 > control-fermat/gbcheck_d6.out 2>&1
echo CONTROL DONE >> control-fermat/control.log
