#!/bin/bash
# Pipeline B: character-MIXED families (never searched before): Hessian relaxation + full-E2 Koszul varpro.
cd /tmp/hodge/results/quadric-contact-residual
export PYTHONPATH=/home/kejace/devel/sextics-work
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-search.pipeB.log.jsonl
run() { echo "=== $(date '+%F %T') $*"; python3 "$@" 2>&1 | grep -v 'hodge engine v2'; echo "=== exit $? $(date '+%F %T')"; }
run w2_hess.py   --d 7 --parity 0 --restarts 1500 --procs 8 --seed 11 --tag hd7e
run w2_hess.py   --d 7 --parity 1 --restarts 1500 --procs 8 --seed 12 --tag hd7o
run w2_varpro.py --d 7 --parity 0 --restarts 300 --per-dim 0 --procs 8 --seed 6000 --tag w2m --time-limit 3600
run w2_varpro.py --d 7 --parity 1 --restarts 300 --per-dim 0 --procs 8 --seed 6100 --tag w2m --time-limit 3600
run w2_hess.py   --d 7 --mixed    --restarts 1000 --procs 8 --seed 13 --tag hd7m
run w2_varpro.py --d 7 --mixed    --restarts 200 --per-dim 0 --procs 8 --seed 6150 --tag w2m --time-limit 3600
run w2_hess.py   --d 8 --parity 0 --restarts 1500 --procs 8 --seed 14 --tag hd8e
run w2_varpro.py --d 8 --parity 0 --restarts 200 --per-dim 0 --procs 8 --seed 6200 --tag w2m --time-limit 5400
run w2_hess.py   --d 8 --parity 1 --restarts 1000 --procs 8 --seed 15 --tag hd8o
run w2_varpro.py --d 8 --parity 1 --restarts 200 --per-dim 0 --procs 8 --seed 6300 --tag w2m --time-limit 5400
run w2_hess.py   --d 9 --parity 0 --restarts 800 --procs 8 --seed 16 --tag hd9e
run w2_varpro.py --d 9 --parity 0 --restarts 120 --per-dim 0 --procs 8 --seed 6400 --tag w2m --time-limit 7200
echo "PIPELINE B DONE $(date)"
