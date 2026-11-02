#!/bin/bash
# usage: submit.sh <name> <cpus> <command...>   (local replacement for hodge.rjob submit on the VM)
name=$1; cpus=$2; shift 2
mkdir -p /tmp/hodge/jobs
cd /tmp/hodge/results/symmetric-sextic-residual
export PYTHONPATH=/home/kejace/devel/sextics-work OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLINT_NUM_THREADS=1 W2_NPROC=$cpus
nohup setsid bash -c "$*" > /tmp/hodge/jobs/$name.log 2>&1 < /dev/null &
echo $! > /tmp/hodge/jobs/$name.pid
echo "submitted $name pid $! cpus $cpus"
