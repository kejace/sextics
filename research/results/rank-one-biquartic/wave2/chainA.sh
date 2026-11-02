#!/bin/bash
# Chain A: exact structure (E1 spaces, relaxation ranks), light searches, structured families, msolve, high-degree E1 spaces.
cd /tmp/hodge/results/rank-one-biquartic
export PYTHONPATH=/home/kejace/devel/sextics-work
O=wave2/out
echo "chainA start $(date)"
for pD in 0 1; do for d in 5 6; do
  python3 wave2/relax.py --d $d --pD $pD --tag relax --msolve-out $O/relax_d${d}_p${pD}.ms > $O/relax_d${d}_p${pD}.out 2> $O/relax_d${d}_p${pD}.err
done; done
echo "relax 5,6 done $(date)"
python3 wave2/e1spaces.py --dmin 4 --dmax 10 --tag e1spaces_d4_10 > $O/e1spaces_d4_10.out 2> $O/e1spaces_d4_10.err
echo "e1spaces 4-10 done $(date)"
for pD in 0 1; do for d in 7; do
  python3 wave2/relax.py --d $d --pD $pD --tag relax --msolve-out $O/relax_d${d}_p${pD}.ms > $O/relax_d${d}_p${pD}.out 2> $O/relax_d${d}_p${pD}.err
done; done
echo "relax 5-7 done $(date)"
python3 wave2/e2search.py --d 5 --pD 0 --restarts 200 --seed 0 --workers 4 --tag d5e_a > $O/d5e_a.out 2> $O/d5e_a.err
python3 wave2/e2search.py --d 5 --pD 1 --restarts 200 --seed 0 --workers 4 --tag d5o_a > $O/d5o_a.out 2> $O/d5o_a.err
python3 wave2/e2search.py --d 6 --pD 0 --restarts 2000 --seed 0 --workers 8 --tag d6e_a > $O/d6e_a.out 2> $O/d6e_a.err
python3 wave2/e2search.py --d 6 --pD 1 --restarts 1000 --seed 0 --workers 8 --tag d6o_a > $O/d6o_a.out 2> $O/d6o_a.err
echo "d5,d6 searches done $(date)"
python3 wave2/relax.py --d 8 --pD 0 --tag relax > $O/relax_d8_p0.out 2> $O/relax_d8_p0.err
python3 wave2/structured.py --dmin 6 --dmax 10 --tag structured > $O/structured_6_10.out 2> $O/structured_6_10.err
echo "structured 6-10 done $(date)"
python3 wave2/e2search.py --d 7 --pD 1 --restarts 400 --seed 0 --workers 8 --tag d7o_a > $O/d7o_a.out 2> $O/d7o_a.err
python3 wave2/e1spaces.py --dmin 7 --dmax 7 --parities 0 --prime 1073741789 --tag e1spaces_d7_q > $O/e1spaces_d7_q.out 2> $O/e1spaces_d7_q.err
python3 wave2/relax.py --d 7 --pD 0 --prime 1073741789 --tag relax_q --msolve-out $O/relax_d7_p0_q.ms > $O/relax_d7_p0_q.out 2> $O/relax_d7_p0_q.err
python3 wave2/relax.py --d 6 --pD 0 --prime 1073741789 --tag relax_q --msolve-out $O/relax_d6_p0_q.ms > $O/relax_d6_p0_q.out 2> $O/relax_d6_p0_q.err
if [ -s $O/relax_d6_p0_q.ms ]; then timeout 1800 msolve -t 4 -v 1 -g 2 -f $O/relax_d6_p0_q.ms -o $O/relax_d6_p0_q.msout > $O/msolve_d6_p0.out 2>&1; echo "msolve d6 exit $?"; fi
if [ -s $O/relax_d7_p0_q.ms ]; then timeout 14400 msolve -t 8 -v 1 -g 2 -f $O/relax_d7_p0_q.ms -o $O/relax_d7_p0_q.msout > $O/msolve_d7_p0.out 2>&1; echo "msolve d7 exit $?"; fi
echo "msolve done $(date)"
python3 wave2/e1spaces.py --dmin 11 --dmax 11 --tag e1spaces_d11 > $O/e1spaces_d11.out 2> $O/e1spaces_d11.err
python3 wave2/e1spaces.py --dmin 12 --dmax 12 --tag e1spaces_d12 > $O/e1spaces_d12.out 2> $O/e1spaces_d12.err
python3 wave2/structured.py --dmin 11 --dmax 12 --tag structured_11_12 > $O/structured_11_12.out 2> $O/structured_11_12.err
for s in 0 1 2; do python3 wave2/e2exact.py --random-member --d 7 --parity 0 --seed $s --tag e2exact >> $O/e2exact_random.out 2>> $O/e2exact_random.err; done
for s in 0 1 2; do python3 wave2/e2exact.py --random-member --d 8 --parity 0 --seed $s --tag e2exact >> $O/e2exact_random.out 2>> $O/e2exact_random.err; done
for i in $(seq 0 26); do python3 wave2/e2exact.py --kernel-json $O/e1kernel_d9_p0.json --index $i --tag e2exact >> $O/e2exact_d9.out 2>> $O/e2exact_d9.err; done
for s in 0 1 2; do python3 wave2/e2exact.py --random-member --d 9 --parity 0 --seed $s --tag e2exact >> $O/e2exact_random.out 2>> $O/e2exact_random.err; done
echo "e2exact done $(date)"
echo "chainA end $(date)"
