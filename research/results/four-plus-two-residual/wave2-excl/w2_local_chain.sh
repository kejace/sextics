#!/bin/bash
# Sequential, single-process (niced) local fallback while the remote box is saturated.
# Each job logs to wave2-excl/remote/<job>.jsonl (merged by w2_pull.py --nopull).
cd /tmp/hodge/results/four-plus-two-residual || exit 1
export PYTHONPATH=/home/kejace/devel/sextics-work
run() { echo "[chain] $(date +%T) start $*"; nice -n 19 timeout "$1" python3 "wave2-excl/$2" "${@:3}" > "wave2-excl/remote/$3.out" 2>&1; echo "[chain] $(date +%T) exit $? $2"; }
run 7200  w2_symsweep.py excl-symsweep-local-k6-8 1 6 8 2500
run 2400  w2_sum.py      excl-sum-local 1
run 3600  w2_struct.py   excl-struct-local 1
run 10800 w2_symsweep.py excl-symsweep-local-k9-12 1 9 12 2500
echo "[chain] all done $(date +%T)"
