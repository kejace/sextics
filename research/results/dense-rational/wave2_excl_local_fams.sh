#!/bin/bash
# local, light: degree-7 mixed-sector families decided at Macaulay degree 3 (<= 63 unknowns), one msolve process at a time
cd /tmp/hodge/results/dense-rational
export PYTHONPATH=/home/kejace/devel/sextics-work
P=1073741789
fams="gradient xonly yonly"
for i in 0 1 2 3 4 5; do for j in 0 1 2 3 4 5; do [ $i -lt $j ] && fams="$fams pair-$i-$j"; done; done
for i in 0 1 2 3 4 5; do for j in 0 1 2 3 4 5; do for k in 0 1 2 3 4 5; do [ $i -lt $j ] && [ $j -lt $k ] && fams="$fams comp-$i-$j-$k"; done; done; done
for fam in $fams; do
  ms=wave2-excl/fam_d7_mixed_${fam}_p$P.ms
  [ -f ${ms%.ms}_ltcheck.json ] && continue
  timeout 120 python3 wave2_excl_quadsys_fam.py 7 mixed $fam $P > /dev/null 2>&1 || { echo "[$fam] generator failed"; continue; }
  t0=$(date +%s); timeout 300 msolve -v 1 -t 4 -g 1 -f $ms -o ${ms%.ms}.g1.out > ${ms%.ms}.msolve.log 2>&1; rc=$?; t1=$(date +%s)
  echo "msolve exit $rc after $((t1-t0)) s" >> ${ms%.ms}.msolve.log
  res=$(timeout 60 python3 wave2_excl_ltcheck.py $ms ${ms%.ms}.g1.out local-d7-mixed-$fam 2>/dev/null | grep -o '"projectively_empty_mod_p": [a-z]*')
  echo "[$fam] rc=$rc $((t1-t0))s $res"
  python3 - <<PY
import sys; sys.path.insert(0,'/tmp/hodge/results/dense-rational')
from wave2_excl_common import log2
log2(exp='msolve-run', job='local-d7-mixed-$fam', d=7, sector='mixed', family='$fam', prime=$P, threads=4, msolve_exit=$rc, seconds=$((t1-t0)), file='$ms', note='run locally (small system)')
PY
done 2>&1 | grep -v NOTE
echo ALL-DONE
