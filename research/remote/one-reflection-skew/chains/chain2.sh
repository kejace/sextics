#!/bin/bash
# chain2: structured candidates (fixed-D served analysis) and coordinate-cone families via varpro, d = 6..12.  8 CPUs.
cd /tmp/hodge/results/one-reflection-skew
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
S="python3 w2_search.py --procs 8"
# cone families: kernel ∩ C[vars] (vertex a plane / a line), sector (0,0), LM
for d in 6 7 8 9 10 11 12; do
  for cone in X0,X1,X2 X0,X2,Y0 X0,X1,Y0 X0,X2,Y1 X1,X2,Y0 X0,Y0,Y1 X2,Y0,Y1 X1,Y0,Y1 X0,X1,X2,Y0 X0,X1,X2,Y1 X0,X2,Y0,Y1 X0,X1,Y0,Y1 X1,X2,Y0,Y1; do
    tag=$(echo $cone | tr -d ',')
    $S --job cone-${tag}-d${d} --d $d --rho 0 --sigma 0 --family cone:$cone --restarts 120 --method lm --seed0 $((1000*d)) --time-limit 1500 >> chain2.out 2>> chain2.err
  done
done
OPENBLAS_NUM_THREADS=4 python3 w2_struct.py --job w2-struct --dmin 6 --dmax 10 >> chain2.out 2>> chain2.err
OPENBLAS_NUM_THREADS=4 python3 w2_struct.py --job w2-struct-hi --dmin 11 --dmax 12 --max-per-degree 60 >> chain2.out 2>> chain2.err
echo CHAIN2 DONE >> chain2.out
