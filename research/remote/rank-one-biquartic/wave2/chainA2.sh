#!/bin/bash
# Chain A2 (replaces the tail of chain A): msolve on the sliced degree-7 relaxation (background, 6 threads) while
# single-core exact steps run: structured families, relax d=8 rank, E1 spaces d=11,12, structured 11-12, exact E2 of
# random members / degree-9 kernel reps.  (d=5/6 numerical searches dropped: those degrees are excluded exactly.)
cd /tmp/hodge/results/rank-one-biquartic
export PYTHONPATH=/home/kejace/devel/sextics-work
O=wave2/out
echo "chainA2 start $(date)"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 FLINT_NUM_THREADS=1
if [ ! -f $O/e1kernel_d10_p0.json ]; then python3 wave2/e1spaces.py --dmin 10 --dmax 10 --tag e1spaces_d10 > $O/e1spaces_d10.out 2> $O/e1spaces_d10.err; fi
echo "e1spaces d10 present $(date)"
( timeout 21600 msolve -t 6 -v 2 -g 2 -f $O/relax_d7_p0_q_slice.ms -o $O/relax_d7_p0_q_slice.msout > $O/msolve_d7_slice.out 2>&1; echo "msolve d7 slice exit $? $(date)" >> $O/msolve_d7_slice.out ) &
MSPID=$!
python3 wave2/relax.py --d 7 --pD 0 --tag relax --msolve-out $O/relax_d7_p0.ms > $O/relax_d7_p0.out 2> $O/relax_d7_p0.err
python3 wave2/relax.py --d 7 --pD 1 --tag relax --msolve-out $O/relax_d7_p1.ms > $O/relax_d7_p1.out 2> $O/relax_d7_p1.err
echo "relax d7 (main prime) done $(date)"
python3 wave2/structured.py --dmin 6 --dmax 10 --tag structured > $O/structured_6_10.out 2> $O/structured_6_10.err
echo "structured 6-10 done $(date)"
python3 wave2/relax.py --d 8 --pD 0 --tag relax > $O/relax_d8_p0.out 2> $O/relax_d8_p0.err
echo "relax d8 done $(date)"
for s in 0 1 2; do python3 wave2/e2exact.py --random-member --d 7 --parity 0 --seed $s --tag e2exact >> $O/e2exact_random.out 2>> $O/e2exact_random.err; done
for s in 0 1 2; do python3 wave2/e2exact.py --random-member --d 8 --parity 0 --seed $s --tag e2exact >> $O/e2exact_random.out 2>> $O/e2exact_random.err; done
echo "e2exact random d7,d8 done $(date)"
python3 wave2/e1spaces.py --dmin 11 --dmax 11 --tag e1spaces_d11 > $O/e1spaces_d11.out 2> $O/e1spaces_d11.err
echo "e1spaces d11 done $(date)"
for i in $(seq 0 26); do python3 wave2/e2exact.py --kernel-json $O/e1kernel_d9_p0.json --index $i --tag e2exact >> $O/e2exact_d9.out 2>> $O/e2exact_d9.err; done
for s in 0 1 2; do python3 wave2/e2exact.py --random-member --d 9 --parity 0 --seed $s --tag e2exact >> $O/e2exact_random.out 2>> $O/e2exact_random.err; done
echo "e2exact d9 done $(date)"
python3 wave2/e1spaces.py --dmin 12 --dmax 12 --tag e1spaces_d12 > $O/e1spaces_d12.out 2> $O/e1spaces_d12.err
echo "e1spaces d12 done $(date)"
python3 wave2/structured.py --dmin 11 --dmax 12 --tag structured_11_12 > $O/structured_11_12.out 2> $O/structured_11_12.err
echo "structured 11-12 done $(date)"
wait $MSPID
echo "chainA2 end $(date)"
