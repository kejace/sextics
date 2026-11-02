#!/bin/bash
# submit chain5 (calibration at K~90 on hesse-web d=8, one-reflection control, A+ half-target, deflated-LM runs); 4 CPUs
cd /home/kejace/devel/sextics-work && timeout 170 python3 -m hodge.rjob pull --lane one-reflection-skew >/dev/null 2>&1
cd /tmp/hodge/results/one-reflection-skew && python3 w2_report.py --merge >/dev/null 2>&1
cd /home/kejace/devel/sextics-work && python3 -m hodge.rjob submit --lane one-reflection-skew --name w2-chain5 --cpus 4 "$@" -- "bash chains/chain5.sh"
