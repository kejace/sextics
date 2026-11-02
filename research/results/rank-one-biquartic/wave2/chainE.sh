#!/bin/bash
# Chain E: focused searches on the first admissible denominators outside the packet's obstruction class
# (degree-9 "mixed" kernel directions), then pure-kernel probes, relaxation-first searches, msolve (degree 7),
# exact tests.  Each step logs incrementally; STOP_<tag> files end a step cooperatively.
cd /tmp/hodge/results/rank-one-biquartic
export PYTHONPATH=/home/kejace/devel/sextics-work
O=wave2/out
echo "chainE start $(date)"
python3 wave2/e2search.py --d 9 --pD 0 --family kernel-mixed --restarts 48 --seed 8000 --workers 8 --lm-iters 60 --tag d9m_a > $O/d9m_a.out 2> $O/d9m_a.err
echo "d9 mixed-kernel full-E2 (LM) done $(date)"
python3 wave2/relaxsearch.py --d 9 --pD 0 --family kernel-mixed --restarts 300 --seed 0 --workers 8 --refine-maxiter 200 --tag r9m_a > $O/r9m_a.out 2> $O/r9m_a.err
echo "d9 mixed-kernel relaxation-first done $(date)"
python3 wave2/e2search.py --d 9 --pD 0 --family kernel --restarts 32 --seed 8100 --workers 8 --lm-iters 60 --tag d9k_b > $O/d9k_b.out 2> $O/d9k_b.err
python3 wave2/e2search.py --d 8 --pD 0 --family kernel --restarts 40 --seed 8200 --workers 8 --lm-iters 60 --tag d8k_b > $O/d8k_b.out 2> $O/d8k_b.err
python3 wave2/e2search.py --d 8 --pD 1 --family kernel --restarts 24 --seed 8300 --workers 8 --lm-iters 60 --tag d8ok_b > $O/d8ok_b.out 2> $O/d8ok_b.err
echo "pure-kernel probes done $(date)"
python3 wave2/relaxsearch.py --d 8 --pD 0 --restarts 1200 --seed 0 --workers 8 --tag r8e_a > $O/r8e_a.out 2> $O/r8e_a.err
python3 wave2/relaxsearch.py --d 7 --pD 0 --restarts 400 --seed 0 --workers 8 --tag r7e_a > $O/r7e_a.out 2> $O/r7e_a.err
python3 wave2/relaxsearch.py --d 9 --pD 0 --restarts 300 --seed 0 --workers 8 --refine-maxiter 200 --tag r9e_a > $O/r9e_a.out 2> $O/r9e_a.err
python3 wave2/relaxsearch.py --d 8 --pD 1 --restarts 600 --seed 0 --workers 8 --tag r8o_a > $O/r8o_a.out 2> $O/r8o_a.err
echo "relaxation-first searches done $(date)"
for i in 6 7 8; do python3 wave2/e2exact.py --kernel-json $O/e1kernel_d8_p1.json --index $i --tag e2exact >> $O/e2exact_d8o.out 2>> $O/e2exact_d8o.err; done
for s in 0 1 2; do python3 wave2/e2exact.py --mixed-member --d 9 --parity 0 --seed $s --tag e2exact >> $O/e2exact_mixed.out 2>> $O/e2exact_mixed.err; done
for s in 3 4; do python3 wave2/e2exact.py --mixed-member --with-J --d 9 --parity 0 --seed $s --tag e2exact >> $O/e2exact_mixed.out 2>> $O/e2exact_mixed.err; done
echo "exact mixed-member tests done $(date)"
# exact E2 of structured members outside J found by structured.py (degree 10): bidegree (5,5) Y-odd (15-dim, all new),
# 5-variable cones (3 new each), Y-base-type (15 new)
python3 wave2/e2exact.py --members-from $O/structured.jsonl --family "bideg(5,5)" --d 10 --parity 1 --max-members 2 --tag e2exact >> $O/e2exact_struct10.out 2>> $O/e2exact_struct10.err
for v in "vars[X0X1Y0Y1Y2]" "vars[X0X2Y0Y1Y2]" "vars[X1X2Y0Y1Y2]"; do
  python3 wave2/e2exact.py --members-from $O/structured.jsonl --family "$v" --d 10 --parity 0 --max-members 1 --tag e2exact >> $O/e2exact_struct10.out 2>> $O/e2exact_struct10.err
  python3 wave2/e2exact.py --members-from $O/structured.jsonl --family "$v" --d 10 --parity 1 --max-members 1 --tag e2exact >> $O/e2exact_struct10.out 2>> $O/e2exact_struct10.err
done
python3 wave2/e2exact.py --members-from $O/structured.jsonl --family "Ybase P(X)M^a g^b" --d 10 --parity 0 --max-members 1 --tag e2exact >> $O/e2exact_struct10.out 2>> $O/e2exact_struct10.err
echo "exact structured-member tests done $(date)"
if [ -s $O/relax_d7_p0_q_slice.ms ] && [ ! -s $O/relax_d7_p0_q_slice.msout ] && [ ! -s $O/relax_d7_p0_q.msout ]; then
  timeout 14400 msolve -t 8 -v 2 -g 2 -f $O/relax_d7_p0_q_slice.ms -o $O/relax_d7_p0_q_slice.msout > $O/msolve_d7_slice_E.out 2>&1
  echo "msolve d7 slice (chain E) exit $? $(date)"
fi
python3 wave2/relax.py --d 7 --pD 1 --prime 1073741789 --tag relax_q --msolve-out $O/relax_d7_p1_q.ms > $O/relax_d7_p1_q.out 2> $O/relax_d7_p1_q.err
if [ -s $O/relax_d7_p1_q.ms ]; then timeout 14400 msolve -t 8 -v 1 -g 2 -f $O/relax_d7_p1_q.ms -o $O/relax_d7_p1_q.msout > $O/msolve_d7_p1.out 2>&1; echo "msolve d7 odd exit $?"; fi
python3 wave2/relaxsearch.py --d 10 --pD 0 --restarts 100 --seed 0 --workers 8 --refine-maxiter 150 --tag r10e_a > $O/r10e_a.out 2> $O/r10e_a.err
echo "chainE end $(date)"
