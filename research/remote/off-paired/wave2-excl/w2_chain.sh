#!/bin/bash
# Combined wave-2 remote chain for lane off-paired (one job; several processes).  Logs: per-process W2LOG files.
cd /tmp/hodge/results/off-paired
PK=wave2-excl
export PYTHONPATH=/home/kejace/devel/sextics-work
# (1) degree-18 discriminant full solve (4 threads)
( export W2LOG=/tmp/hodge/results/off-paired/wave2-excl.chain-deltat.log.jsonl; python3 $PK/w2_deltat.py excl-off-paired-chain-deltat > $PK/stdout_chain_deltat.txt 2>&1 ) &
# (2) exact rigor pipeline, degrees 7 and 8, all 36 sectors, six parallel groups, Singular cap 2 h per call
for g in "00,01,02,03,04,05" "10,11,12,13,14,15" "20,21,22,23,24,25" "30,31,32,33,34,35" "40,41,42,43,44,45" "50,51,52,53,54,55"; do
  ( export W2LOG=/tmp/hodge/results/off-paired/wave2-excl.chain-rigor-${g%%,*}.log.jsonl
    python3 $PK/w2_elim2.py --d 8 --chars $g --job excl-off-paired-chain-rigor-d8 --timeout 7200 > $PK/stdout_chain_rigor_d8_${g%%,*}.txt 2>&1
    python3 $PK/w2_elim2.py --d 7 --chars $g --job excl-off-paired-chain-rigor-d7 --timeout 3600 > $PK/stdout_chain_rigor_d7_${g%%,*}.txt 2>&1 ) &
done
# (3) degree-9 chart-complete modular sweep, six parallel groups x 6 threads
for g in "00,01,02,03,04,05" "10,11,12,13,14,15" "20,21,22,23,24,25" "30,31,32,33,34,35" "40,41,42,43,44,45" "50,51,52,53,54,55"; do
  ( export W2LOG=/tmp/hodge/results/off-paired/wave2-excl.chain-d9-${g%%,*}.log.jsonl
    python3 $PK/w2_bilinear.py --d 9 --chars $g --mode sector --charts rev --threads 6 --timeout 2400 --job excl-off-paired-chain-bilin-d9 --addJ --fast > $PK/stdout_chain_d9_${g%%,*}.txt 2>&1 ) &
done
# (4) multi-prime re-certification of the three kernel-lift-incomplete receipts; kernel dims d=12
( export W2LOG=/tmp/hodge/results/off-paired/wave2-excl.chain-misc.log.jsonl
  python3 $PK/w2_recert.py --labels 'X1*M*H,X1*S*FY0,X1*M*FY0' --job excl-off-paired-chain-recert > $PK/stdout_chain_recert.txt 2>&1
  python3 $PK/w2_kernels.py 12 12 > $PK/stdout_chain_kernels12.txt 2>&1 ) &
wait
echo CHAIN DONE
