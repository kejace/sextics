#!/bin/bash
# Chain B: degree-7 Y-even family (66-dim) multistart search in batches, then Y-odd degree 8.
cd /tmp/hodge/results/rank-one-biquartic
export PYTHONPATH=/home/kejace/devel/sextics-work
O=wave2/out
echo "chainB start $(date)"
for b in 0 1 2 3 4 5 6 7; do
  python3 wave2/e2search.py --d 7 --pD 0 --restarts 200 --seed $((b*200)) --workers 8 --tag d7e_b$b > $O/d7e_b$b.out 2> $O/d7e_b$b.err
  echo "d7e batch $b done $(date)"
done
for b in 0 1 2 3; do
  python3 wave2/e2search.py --d 8 --pD 1 --restarts 100 --seed $((b*100)) --workers 8 --tag d8o_b$b > $O/d8o_b$b.out 2> $O/d8o_b$b.err
  echo "d8o batch $b done $(date)"
done
echo "chainB end $(date)"
