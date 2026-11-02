#!/bin/bash
export NTHREADS=4
python3 w2c_fibre_tests.py J2f > fibre_J2f.log 2>&1
for c in X0 X1 X2 Y0 Y1 Y2; do timeout 2400 python3 w2c_hessrank.py $c 1 > hessrank1_$c.log 2>&1; done
python3 w2c_fibre_tests.py J1f > fibre_J1f.log 2>&1
