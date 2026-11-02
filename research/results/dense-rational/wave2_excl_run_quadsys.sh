#!/bin/bash
# Usage: bash wave2_excl_run_quadsys.sh d sector prime threads jobname   (runs in the lane dir on the remote)
d=$1; sector=$2; p=$3; T=$4; job=$5
cd /tmp/hodge/results/dense-rational
export PYTHONPATH=/home/kejace/devel/sextics-work
ms=wave2-excl/quadsys_d${d}_${sector}_p${p}.ms
python3 wave2_excl_quadsys.py $d $sector $p || exit 1
t0=$(date +%s)
timeout 72000 msolve -v 2 -t $T -g 1 -f $ms -o ${ms%.ms}.g1.out > ${ms%.ms}.msolve.log 2>&1
rc=$?
t1=$(date +%s)
echo "msolve exit $rc after $((t1-t0)) s" >> ${ms%.ms}.msolve.log
python3 wave2_excl_ltcheck.py $ms ${ms%.ms}.g1.out $job
python3 - <<PY
import json,sys
sys.path.insert(0,'/tmp/hodge/results/dense-rational')
from wave2_excl_common import log2
log2(exp='msolve-run', job='$job', d=$d, sector='$sector', prime=$p, threads=$T, msolve_exit=$rc, seconds=$((t1-t0)), file='$ms')
PY
