#!/bin/bash
# degree-6 Hessian relaxation, three sectors, msolve -g 1 leading ideal + ltcheck (VM agent)
cd /tmp/hodge/results/dense-rational
export PYTHONPATH=/home/kejace/devel/sextics-work
P=1073741789
for sec in even odd mixed; do
  ms=wave2-excl/quadsys_d6_${sec}_p$P.ms
  t0=$(date +%s)
  timeout 7200 msolve -v 1 -t 2 -g 1 -f $ms -o ${ms%.ms}.g1.out > ${ms%.ms}.msolve.log 2>&1; rc=$?
  t1=$(date +%s); echo "msolve exit $rc after $((t1-t0)) s" >> ${ms%.ms}.msolve.log
  python3 wave2_excl_ltcheck.py $ms ${ms%.ms}.g1.out vm-d6-$sec
done
echo done > wave2-excl/d6.done
