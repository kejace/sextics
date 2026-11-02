#!/bin/bash
# Chain C: degree-8 Y-even family (177-dim) multistart search in batches.
cd /tmp/hodge/results/rank-one-biquartic
export PYTHONPATH=/home/kejace/devel/sextics-work
O=wave2/out
echo "chainC start $(date)"
for b in 0 1 2 3 4 5 6 7 8 9; do
  python3 wave2/e2search.py --d 8 --pD 0 --restarts 80 --seed $((b*80)) --workers 8 --tag d8e_b$b > $O/d8e_b$b.out 2> $O/d8e_b$b.err
  echo "d8e batch $b done $(date)"
done
echo "chainC end $(date)"
