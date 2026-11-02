#!/bin/bash
# Sequentially submit remote jobs, waiting for capacity (rjob refuses when load>150 or <60GB free).
cd /home/kejace/devel/sextics-work
L=polarized-k3-remainder
submit() { name=$1; cpus=$2; shift 2; python3 -m hodge.rjob submit --lane $L --name $name --cpus $cpus --wait 14400 -- "$@" 2>&1 | grep -v "engine v2"; echo "[queue] $name rc=$? $(date)"; }
submit excl-polarized-k3-remainder-relaxfull 1 bash wave2-excl/job_relax_full.sh
sleep 30
submit excl-polarized-k3-remainder-relaxk3 1 bash wave2-excl/job_relax_k3.sh
sleep 30
submit excl-polarized-k3-remainder-e1table 2 "cd wave2-excl && python3 e1table.py 2147483629 9 e1table-p2147483629.json"
echo "[queue] all submitted $(date)"
