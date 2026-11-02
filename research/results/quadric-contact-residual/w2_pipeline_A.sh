#!/bin/bash
# Pipeline A: full-E2 Koszul varpro over every character sector, degrees 9..12, then second-seed passes.
cd /tmp/hodge/results/quadric-contact-residual
export PYTHONPATH=/home/kejace/devel/sextics-work
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-search.pipeA.log.jsonl
run() { echo "=== $(date '+%F %T') $*"; python3 "$@" 2>&1 | grep -v 'hodge engine v2'; echo "=== exit $? $(date '+%F %T')"; }
run w2_varpro.py --d 8  --chi all --restarts 300 --per-dim 10 --min-restarts 24 --procs 8 --seed 12000 --tag w2b --time-limit 1800
run w2_varpro.py --d 9  --chi all --restarts 120 --per-dim 6 --min-restarts 16 --procs 8 --seed 2000 --tag w2a --time-limit 1500
run w2_varpro.py --d 10 --chi all --restarts 160 --per-dim 4 --min-restarts 16 --procs 8 --seed 3000 --tag w2a --time-limit 2400
run w2_varpro.py --d 11 --chi all --restarts 160 --per-dim 3 --min-restarts 16 --procs 8 --seed 4000 --tag w2a --time-limit 2400
run w2_varpro.py --d 12 --chi all --restarts 160 --per-dim 2 --min-restarts 16 --procs 8 --seed 5000 --tag w2a --time-limit 3000
# second passes with fresh seeds
run w2_varpro.py --d 9  --chi all --restarts 200 --per-dim 8 --min-restarts 24 --procs 8 --seed 13000 --tag w2b --time-limit 1800
run w2_varpro.py --d 10 --chi all --restarts 200 --per-dim 5 --min-restarts 24 --procs 8 --seed 14000 --tag w2b --time-limit 2400
run w2_varpro.py --d 12 --chi 0000,0001,1100,1101 --restarts 200 --per-dim 3 --min-restarts 24 --procs 8 --seed 15000 --tag w2b --time-limit 3600
run w2_varpro.py --d 11 --chi all --restarts 200 --per-dim 4 --min-restarts 24 --procs 8 --seed 16000 --tag w2b --time-limit 2400
run w2_varpro.py --d 8  --chi all --restarts 400 --per-dim 15 --min-restarts 30 --procs 8 --seed 17000 --tag w2c --time-limit 2400
echo "PIPELINE A DONE $(date)"
