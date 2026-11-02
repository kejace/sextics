#!/bin/bash
# Job excl-sub7: degree-7 single-involution-pure subfamilies (60 unknowns each) by msolve -g 1 (16 threads, 1.5 h cap each),
# then block-Macaulay k=4 for any not closed; then the full 120-unknown analysis (I_3) and a capped msolve run.
cd /tmp/hodge/results/quadric-contact-residual/wave2-excl
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-excl/excl-sub7.log.jsonl
J=d7mixed.quads.local.json
for pr in y0 y1 b0=0 b0=1 b2=0 b2=1; do
  tag=d7_$(echo "$pr" | tr -d '=&')
  MSOLVE_TIMEOUT=5400 python3 sub_closure.py $J $tag "$pr" 16 2>&1 | grep -v "engine v2" | tail -3
done
for pr in y0 y1 b0=0 b0=1 b2=0 b2=1; do
  tag=d7_$(echo "$pr" | tr -d '=&')
  if ! grep -q "\"tag\": \"$tag\".*EXCLUDED" $HODGE_W2LOG; then
    python3 sub_closure.py $J ${tag}_mac "$pr" 0 4 12 2>&1 | grep -v "engine v2" | tail -3
  fi
done
MSOLVE_THREADS=16 MSOLVE_TIMEOUT=10800 python3 d7mixed.py 2>&1 | grep -v "engine v2" | tail -5
echo SUB7 DONE
