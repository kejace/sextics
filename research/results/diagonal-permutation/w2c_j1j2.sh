#!/bin/bash
export NTHREADS=4
for c in 0 1 2 3 4 5; do python3 w2c_reflection.py $c offX > reflection_offX_v$c.log 2>&1; done
for c in 0 1 2 3 4 5; do python3 w2c_reflection.py $c onX > reflection_onX_v$c.log 2>&1; done
for c in X0 X1 X2 Y0 Y1 Y2; do timeout 2400 python3 w2c_hessrank.py $c > hessrank_$c.log 2>&1; done
