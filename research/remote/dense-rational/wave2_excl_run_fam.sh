#!/bin/bash
# Usage: bash wave2_excl_run_fam.sh d sector family prime threads jobname
d=$1; sector=$2; fam=$3; p=$4; T=$5; job=$6
cd /tmp/hodge/results/dense-rational
export PYTHONPATH=/home/kejace/devel/sextics-work
ms=wave2-excl/fam_d${d}_${sector}_${fam}_p${p}.ms
python3 wave2_excl_quadsys_fam.py $d $sector $fam $p || exit 1
t0=$(date +%s)
timeout 36000 msolve -v 1 -t $T -g 1 -f $ms -o ${ms%.ms}.g1.out > ${ms%.ms}.msolve.log 2>&1
rc=$?; t1=$(date +%s)
echo "msolve exit $rc after $((t1-t0)) s" >> ${ms%.ms}.msolve.log
python3 wave2_excl_ltcheck.py $ms ${ms%.ms}.g1.out $job
python3 - <<PY
import sys; sys.path.insert(0,'/tmp/hodge/results/dense-rational')
from wave2_excl_common import log2
log2(exp='msolve-run', job='$job', d=$d, sector='$sector', family='$fam', prime=$p, threads=$T, msolve_exit=$rc, seconds=$((t1-t0)), file='$ms')
PY
