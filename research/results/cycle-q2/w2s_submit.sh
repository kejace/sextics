#!/bin/bash
# Local stand-in for `hodge.rjob submit` (we are on the VM; ssh-to-self is unavailable).
# usage: w2s_submit.sh <name> <cpus> <command...>   -> runs in /tmp/hodge/results/cycle-q2, log /tmp/hodge/jobs/<name>.log
name=$1; cpus=$2; shift 2
cmd="$*"
mkdir -p /tmp/hodge/jobs /tmp/hodge/results/cycle-q2
cd /tmp/hodge/results/cycle-q2 || exit 1
export PYTHONPATH=/home/kejace/devel/sextics-work OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FLINT_NUM_THREADS=1
nohup setsid bash -c "$cmd" > /tmp/hodge/jobs/$name.log 2>&1 < /dev/null &
echo $! > /tmp/hodge/jobs/$name.pid
echo "{\"name\": \"$name\", \"lane\": \"cycle-q2\", \"pid\": $!, \"cpus\": $cpus, \"cmd\": \"$cmd\", \"ts\": \"$(date -Is)\"}" >> /tmp/hodge/jobs-registry.jsonl
echo "submitted $name pid $! log /tmp/hodge/jobs/$name.log"
