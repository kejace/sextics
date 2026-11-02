#!/bin/bash
cd /tmp/hodge/results/cycle-q2
for d in 7 8 9 10; do
  python3 w2_search.py --d $d --family multi --subfamilies cone4,mult4,cone3 --restarts 15 --maxiter 60 --nproc 32 --tag multi$d > search-cycle-q2-multi$d.out 2>&1
done
bash w2_high.sh
echo CHAINB-DONE >> search-cycle-q2-k12.out
