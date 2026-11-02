#!/bin/bash
# Sequential local fallback queue v3 (remote saturated): one nice'd process at a time, priorities re-ordered.
cd /tmp/hodge/results/polarized-k3-remainder/wave2-excl
export PYTHONPATH=/home/kejace/devel/sextics-work
while pgrep -fx "python3 relax.py 7 65521 --target=k3 --exact-ranks --squares --out=relax-k3-d7-p65521-local" >/dev/null; do sleep 30; done
echo "=== [chain3] E1 tables second prime d<=9 $(date)"
nice -n 19 python3 e1table.py 2147483629 9 e1table-p2147483629-local.json 2>&1 | grep -v "engine v2"
echo "=== [chain3] full-target d=7 replication p=65521 (modular; L-ranks over Q come from the K3 run, same L) $(date)"
nice -n 19 python3 relax.py 7 65521 --squares --out=relax-full-d7-p65521-local 2>&1 | grep -v "engine v2" | grep -v "^\[result\]"
echo "=== [chain3] K3-target d=8 relaxation p=65521 (modular) $(date)"
nice -n 19 python3 relax.py 8 65521 --target=k3 --squares --out=relax-k3-d8-p65521-local 2>&1 | grep -v "engine v2" | grep -v "^\[result\]"
echo "=== [chain3] full-target d=8 relaxation regeneration p=65521 (modular) $(date)"
nice -n 19 python3 relax.py 8 65521 --squares --out=relax-full-d8-p65521-local 2>&1 | grep -v "engine v2" | grep -v "^\[result\]"
echo "=== [chain3] fU*h family exact tests d=8 $(date)"
nice -n 19 bash job_fUh_d8.sh 2>&1 | grep -v "engine v2"
echo "=== [chain3] DONE $(date)"
