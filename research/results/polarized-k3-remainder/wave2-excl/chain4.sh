#!/bin/bash
# After chain3 finishes: full-target degree-8 membership tests of the thread's spurious family (were killed by a 600 s wrapper).
cd /tmp/hodge/results/polarized-k3-remainder/wave2-excl
export PYTHONPATH=/home/kejace/devel/sextics-work
until grep -q "\[chain3\] DONE" chain3.log 2>/dev/null; do sleep 60; done
for spec in "U0sq:mul(U[0],U[0])" "U0V0:mul(U[0],V[0])" "U0sq+U0V0:add(mul(U[0],U[0]),mul(U[0],V[0]))"; do
  name=${spec%%:*}; h=${spec#*:}
  echo "=== [chain4] fU*$name full $(date)"
  timeout 5400 nice -n 15 python3 relax_member.py "fU*$name" 8 "[scale(mul(U[k],$h),Fr(1,42)) for k in range(3)] + [{},{},{}]" --target=full --modular-only --out=relax-member-d8.jsonl 2>&1 | grep "^\[member\]" | cut -c1-300
done
echo "=== [chain4] DONE $(date)"
