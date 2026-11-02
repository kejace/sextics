#!/bin/bash
# Exp B driver: reflection searches (X itself, Fermat control, six coordinate pencils). p = 2147483647 (msolve needs p < 2^31).
cd /tmp/hodge/results/reflection-coupled-residual/construct
export MSOLVE_THREADS=16
P=2147483647
FEXPR='X0**6 + 5*X0**3*X1*Y0*Y1/7 + 5*X0*X1*Y0**3*Y1/7 + 5*X0*X1*Y0**2*Y2**2/7 + 10*X0*X1*Y1**2*Y2**2/7 + 15*X0*X2**3*Y0*Y1/7 + X1**6 + 10*X1**3*X2*Y0*Y1/7 + 5*X1*X2*Y0*Y1**3/7 + X2**6 + Y0**6 + Y1**6 + Y2**6'
run() { name=$1; shift; echo "=== $name $(date)"; python3 reflsearch.py "$name" $P 7 "$@" > "$name.out" 2>&1; tail -2 "$name.out"; }
run ctrl_fermat5 'X0**6+X1**6+X2**6+Y0**6+Y1**6'
run fib_Y0Y1 "$FEXPR" 'Y0=l*Y1'
run fib_Y0Y2 "$FEXPR" 'Y0=l*Y2'
run fib_Y1Y2 "$FEXPR" 'Y1=l*Y2'
run fib_X0X1 "$FEXPR" 'X0=l*X1'
run fib_X0X2 "$FEXPR" 'X0=l*X2'
run fib_X1X2 "$FEXPR" 'X1=l*X2'
run full_X "$FEXPR"
echo "=== ALL DONE $(date)"
