#!/bin/bash
# Pipeline D: <s0,r1>-invariant families (sector 0000) at d = 8, 10, 12 and reflection-adapted cone families (u/v).
cd /tmp/hodge/results/quadric-contact-residual
export PYTHONPATH=/home/kejace/devel/sextics-work
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-search.pipeD.log.jsonl
run() { echo "=== $(date '+%F %T') $*"; python3 "$@" 2>&1 | grep -v 'hodge engine v2'; echo "=== exit $? $(date '+%F %T')"; }
run w2_varpro.py --d 8  --chi 0000 --sym --restarts 600 --per-dim 30 --min-restarts 50 --procs 8 --seed 21000 --tag sym --time-limit 1800
run w2_varpro.py --d 10 --chi 0000 --sym --restarts 500 --per-dim 8  --min-restarts 50 --procs 8 --seed 22000 --tag sym --time-limit 3000
run w2_uvcones.py --dmin 8 --dmax 12 --restarts 40 --procs 8 --time-limit 1200
run w2_varpro.py --d 12 --chi 0000 --sym --restarts 400 --per-dim 4  --min-restarts 50 --procs 8 --seed 23000 --tag sym --time-limit 3600
echo "PIPELINE D DONE $(date)"
