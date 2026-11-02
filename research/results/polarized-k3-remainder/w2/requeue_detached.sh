#!/bin/bash
# Detached re-queue of the wave-2 remote jobs for lane polarized-k3-remainder (they obey the runner's capacity guard).
# Each waits up to 12 h for capacity, then runs; outputs land in the remote lane dir (w2-*.jsonl, *.out); collect with
#   cd /tmp/hodge/results/polarized-k3-remainder && python3 w2_merge.py && python3 w2_report.py
cd /home/kejace/devel/sextics-work || exit 1
L=polarized-k3-remainder; W=43200; O=/tmp/hodge/results/$L/w2
q() { name=$1; cpus=$2; shift 2; setsid nohup python3 -m hodge.rjob submit --wait $W --lane $L --name $name --cpus $cpus -- "$@" > $O/submit-$name.out 2>&1 < /dev/null & }
q w2-k3-hess-tensor-d8 2 'python3 w2_hess_tensor.py 8 1500 101 --tag r1 --maxiter 150 && python3 w2_hess_tensor.py 8 1500 102 --deflate-P8 --tag r2 --maxiter 200 && python3 w2_hess_tensor.py 8 400 103 --mixed --tag r3 --maxiter 150 && python3 w2_hess_tensor.py 8 400 104 --mixed --deflate-P8 --tag r4 --maxiter 200'
q w2-k3-sweep-a 4 'python3 w2_sweep.py w2-candidates.json a --dmax 12 --e1only-above 10'
q w2-k3-full-d8p0-s1 8 'python3 w2_full_e2.py 8 0 600 1 --maxiter 60 --tag s1 --procs 8'
echo "detached submissions launched at $(date)"
