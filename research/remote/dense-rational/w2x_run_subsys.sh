#!/bin/bash
# Macaulay-degree-3 emptiness tests for the Y-degree strata of the degree-7 Hessian relaxation (VM excl agent)
cd /tmp/hodge/results/dense-rational
export PYTHONPATH=/home/kejace/devel/sextics-work
sec=$1
if [ "$sec" = even ]; then ms=wave2-excl/quadsys_d7_even_p1073741789.ms; sets="X0,X2 Y1,X2 X0,Y1 X0 X2 Y1"; else ms=wave2-excl/quadsys_d7_odd_p1073741789.ms; sets="Y0,Y2 X1,Y2 Y0,X1 Y0 Y2 X1"; fi
for z in $sets; do
  echo "=== zero $z"; python3 w2x_subsys.py $ms --zero $z --macaulay 3 2>&1 | grep -v NOTE
done
echo done
