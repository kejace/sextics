#!/bin/bash
# local job launcher mirroring hodge.rjob conventions: jobsub.sh <name> <cpus> <command...>
name=$1; cpus=$2; shift 2
mkdir -p /tmp/hodge/jobs /tmp/hodge/results/polarized-k3-remainder
cd /tmp/hodge/results/polarized-k3-remainder
export PYTHONPATH=/home/kejace/devel/sextics-work OMP_NUM_THREADS=$cpus OPENBLAS_NUM_THREADS=$cpus FLINT_NUM_THREADS=$cpus
nohup setsid bash -c "$*" > /tmp/hodge/jobs/$name.log 2>&1 < /dev/null &
echo $! > /tmp/hodge/jobs/$name.pid
echo "submitted $name pid $! log /tmp/hodge/jobs/$name.log"
echo "{\"name\": \"$name\", \"lane\": \"polarized-k3-remainder\", \"pid\": \"$!\", \"cmd\": \"$*\", \"ts\": \"$(date +%FT%T)\"}" >> /tmp/hodge/jobs-registry.jsonl
