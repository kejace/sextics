#!/bin/bash
# Pipeline C: structured geometry-suggested candidates (exact sector solves) + cone families (varpro).
cd /tmp/hodge/results/quadric-contact-residual
export PYTHONPATH=/home/kejace/devel/sextics-work
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-search.pipeC.log.jsonl
run() { echo "=== $(date '+%F %T') $*"; python3 "$@" 2>&1 | grep -v 'hodge engine v2'; echo "=== exit $? $(date '+%F %T')"; }
run w2_struct.py 8
run w2_struct.py 10
run w2_struct.py 12
run w2_cones.py --dmin 6 --dmax 12 --restarts 60 --procs 8
echo "PIPELINE C DONE $(date)"
