#!/bin/bash
# Control on the same F: the served companion target A' = X0^3 Y2^3 (PR #150, existential). E1 kernels by character d = 1..9.
cd /tmp/hodge/results/reflection-coupled-residual/wave2-excl
export W2LIB_A='X0**3*Y2**3' W2LIB_TAG=companion W2LIB_SUFFIX=.companion
python3 e1char.py 1 9 2147483587 > control-companion/e1char.out 2>&1
echo COMPANION DONE >> control-companion/control.log
