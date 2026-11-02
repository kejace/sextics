#!/bin/bash
cd /tmp/hodge/results/cycle-q2
python3 w2_search.py --d 7 --family kernel --restarts 10 --batch 10 --maxiter 60 --nproc 48 --tag k7 --exclude-keys  > search-cycle-q2-k7b.out 2>&1
python3 w2_search.py --d 8 --family kernel --restarts 5 --batch 5 --maxiter 60 --nproc 48 --tag k8 > search-cycle-q2-k8.out 2>&1
python3 w2_search.py --d 9 --family kernel --restarts 3 --batch 3 --maxiter 60 --nproc 48 --tag k9 > search-cycle-q2-k9.out 2>&1
echo CHAINA-DONE >> search-cycle-q2-k9.out
