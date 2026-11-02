#!/bin/bash
# Chain D2: pure E1-kernel probes (complement of J_d in the admissible space) at d = 7, 8, 9, then small d=9 / d=10 full-family batches.
cd /tmp/hodge/results/rank-one-biquartic
export PYTHONPATH=/home/kejace/devel/sextics-work
O=wave2/out
echo "chainD2 start $(date)"
python3 wave2/e2search.py --d 7 --pD 0 --family kernel --restarts 40 --seed 5000 --workers 8 --tag d7k_a > $O/d7k_a.out 2> $O/d7k_a.err
echo "d7 kernel probe done $(date)"
python3 wave2/e2search.py --d 8 --pD 0 --family kernel --restarts 40 --seed 5000 --workers 8 --tag d8k_a > $O/d8k_a.out 2> $O/d8k_a.err
echo "d8 kernel probe done $(date)"
python3 wave2/e2search.py --d 8 --pD 1 --family kernel --restarts 24 --seed 5000 --workers 8 --tag d8ok_a > $O/d8ok_a.out 2> $O/d8ok_a.err
python3 wave2/e2search.py --d 9 --pD 0 --family kernel --restarts 32 --seed 5000 --workers 8 --maxiter 300 --tag d9k_a > $O/d9k_a.out 2> $O/d9k_a.err
echo "d9 kernel probe done $(date)"
for b in 0 1 2 3 4 5 6 7; do
  python3 wave2/e2search.py --d 9 --pD 0 --restarts 8 --seed $((6000+b*8)) --workers 8 --maxiter 300 --tag d9e_c$b > $O/d9e_c$b.out 2> $O/d9e_c$b.err
  echo "d9e batch c$b done $(date)"
done
python3 wave2/e2search.py --d 10 --pD 0 --restarts 8 --seed 7000 --workers 8 --maxiter 200 --tag d10e_c0 > $O/d10e_c0.out 2> $O/d10e_c0.err
echo "chainD2 end $(date)"
