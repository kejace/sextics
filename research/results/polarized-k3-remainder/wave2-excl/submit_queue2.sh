#!/bin/bash
# Second submission queue: starts after queue 1 has submitted everything; keeps <= 3 of my jobs running.
cd /home/kejace/devel/sextics-work
L=polarized-k3-remainder
until grep -q "all submitted" /tmp/hodge/results/$L/wave2-excl/submit_queue.log 2>/dev/null; do sleep 120; done
myrunning() { python3 -m hodge.rjob status 2>/dev/null | grep -c "RUNNING excl-polarized-k3-remainder"; }
submit() { name=$1; cpus=$2; shift 2
  while [ "$(myrunning)" -ge 3 ]; do sleep 120; done
  python3 -m hodge.rjob submit --lane $L --name $name --cpus $cpus --wait 14400 -- "$@" 2>&1 | grep -v "engine v2"; echo "[queue2] $name submitted $(date)"; }
submit excl-polarized-k3-remainder-k3d7exact 1 bash wave2-excl/job_relax_k3_d7_exact.sh
submit excl-polarized-k3-remainder-relaxd8 1 bash wave2-excl/job_relax_d8.sh
submit excl-polarized-k3-remainder-sweep 8 bash wave2-excl/job_sweep_d8d10.sh
submit excl-polarized-k3-remainder-k3kernel 2 "cd wave2-excl && python3 k3kernel_uv.py remote 4"
echo "[queue2] all submitted $(date)"
