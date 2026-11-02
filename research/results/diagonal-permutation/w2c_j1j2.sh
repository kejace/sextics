#!/bin/bash
# run J1 (6 charts) and J2 (6 charts) sequentially with 16 threads each; outputs named after the job
export NTHREADS=16
for c in X0 X1 X2 Y0 Y1 Y2; do python3 w2c_hessrank.py $c > hessrank_$c.log 2>&1; done
for c in 0 1 2 3 4 5; do python3 w2c_reflection.py $c > reflection_v$c.log 2>&1; done
