#!/bin/bash
# Chain D: degree-9 Y-even family (405-dim) then degree-10 (>=1044-dim) multistart searches.
cd /tmp/hodge/results/rank-one-biquartic
export PYTHONPATH=/home/kejace/devel/sextics-work
O=wave2/out
echo "chainD start $(date)"
for b in 0 1 2 3 4 5; do
  python3 wave2/e2search.py --d 9 --pD 0 --restarts 32 --seed $((b*32)) --workers 8 --maxiter 300 --tag d9e_b$b > $O/d9e_b$b.out 2> $O/d9e_b$b.err
  echo "d9e batch $b done $(date)"
done
python3 wave2/e2search.py --d 10 --pD 0 --restarts 16 --seed 0 --workers 8 --maxiter 200 --tag d10e_b0 > $O/d10e_b0.out 2> $O/d10e_b0.err
echo "chainD end $(date)"
