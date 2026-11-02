#!/bin/bash
# Chain E: relaxation-first searches (cheap) with full-E2 refinement of relaxation zeros; degrees 7, 8, 9 (Y-even),
# 8 (Y-odd), pure-kernel subfamilies; then d=8 pure-kernel full-E2 probe redo.
cd /tmp/hodge/results/rank-one-biquartic
export PYTHONPATH=/home/kejace/devel/sextics-work
O=wave2/out
echo "chainE start $(date)"
# 0. decisive exact step first: msolve on the SLICED degree-7 Y-even relaxation (432 quadrics, 63 vars), 8 threads
if [ -s $O/relax_d7_p0_q_slice.ms ] && [ ! -s $O/relax_d7_p0_q_slice.msout ]; then
  timeout 14400 msolve -t 8 -v 2 -g 2 -f $O/relax_d7_p0_q_slice.ms -o $O/relax_d7_p0_q_slice.msout > $O/msolve_d7_slice_E.out 2>&1
  echo "msolve d7 slice (chain E) exit $? $(date)"
fi
python3 wave2/relaxsearch.py --d 7 --pD 0 --restarts 400 --seed 0 --workers 8 --tag r7e_a > $O/r7e_a.out 2> $O/r7e_a.err
echo "r7e_a done $(date)"
python3 wave2/relaxsearch.py --d 8 --pD 0 --family kernel --restarts 200 --seed 0 --workers 8 --tag r8k_a > $O/r8k_a.out 2> $O/r8k_a.err
python3 wave2/relaxsearch.py --d 8 --pD 0 --restarts 1500 --seed 0 --workers 8 --tag r8e_a > $O/r8e_a.out 2> $O/r8e_a.err
echo "r8e_a done $(date)"
python3 wave2/relaxsearch.py --d 8 --pD 1 --restarts 600 --seed 0 --workers 8 --tag r8o_a > $O/r8o_a.out 2> $O/r8o_a.err
echo "r8o_a done $(date)"
python3 wave2/relaxsearch.py --d 9 --pD 0 --family kernel --restarts 100 --seed 0 --workers 8 --refine-maxiter 200 --tag r9k_a > $O/r9k_a.out 2> $O/r9k_a.err
python3 wave2/relaxsearch.py --d 9 --pD 0 --restarts 300 --seed 0 --workers 8 --refine-maxiter 200 --tag r9e_a > $O/r9e_a.out 2> $O/r9e_a.err
echo "r9e_a done $(date)"
python3 wave2/e2search.py --d 8 --pD 0 --family kernel --restarts 48 --seed 5000 --workers 8 --tag d8k_b > $O/d8k_b.out 2> $O/d8k_b.err
python3 wave2/relaxsearch.py --d 8 --pD 0 --restarts 1500 --seed 1500 --workers 8 --tag r8e_b > $O/r8e_b.out 2> $O/r8e_b.err
python3 wave2/relaxsearch.py --d 10 --pD 0 --restarts 100 --seed 0 --workers 8 --refine-maxiter 150 --tag r10e_a > $O/r10e_a.out 2> $O/r10e_a.err
# exact: Y-odd degree-7 relaxation at the msolve prime (sliced export) + msolve
python3 wave2/relax.py --d 7 --pD 1 --prime 1073741789 --tag relax_q --msolve-out $O/relax_d7_p1_q.ms > $O/relax_d7_p1_q.out 2> $O/relax_d7_p1_q.err
if [ -s $O/relax_d7_p1_q.ms ]; then timeout 14400 msolve -t 8 -v 1 -g 2 -f $O/relax_d7_p1_q.ms -o $O/relax_d7_p1_q.msout > $O/msolve_d7_p1.out 2>&1; echo "msolve d7 odd exit $?"; fi
echo "chainE end $(date)"
