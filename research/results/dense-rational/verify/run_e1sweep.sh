#!/bin/bash
cd /tmp/hodge/results/dense-rational/verify
P=1000003
for d in 10 11 12 13 14 15 16 17 18; do for par in 0 1; do python3 e1sweep_v.py $P dual $d $par; done; done
for d in 1 2 3 4 5 6 7; do for par in 0 1; do python3 e1sweep_v.py $P direct $d $par; done; done
for d in 8 9; do for par in 0 1; do python3 e1sweep_v.py $P direct $d $par; done; done
echo "SWEEP DONE"
