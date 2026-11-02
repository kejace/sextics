#!/bin/bash
cd /tmp/hodge/results/cycle-q2
python3 w2_search.py --d 6 --family kernel --restarts 100 --batch 25 --maxiter 200 --nproc 16 --tag k6 > search-cycle-q2-k6.out 2>&1
python3 w2_search.py --d 6 --family kernel --keys 0 --restarts 1000 --split 16 --batch 25 --maxiter 200 --nproc 16 --tag k6triv > search-cycle-q2-k6triv.out 2>&1
python3 w2_search.py --d 7 --family kernel --restarts 60 --batch 20 --maxiter 200 --nproc 16 --tag k7 > search-cycle-q2-k7.out 2>&1
python3 w2_search.py --d 8 --family kernel --restarts 30 --batch 10 --maxiter 200 --nproc 16 --tag k8 > search-cycle-q2-k8.out 2>&1
python3 w2_search.py --d 9 --family kernel --restarts 12 --batch 6 --maxiter 150 --nproc 16 --tag k9 > search-cycle-q2-k9.out 2>&1
echo CHAINA-DONE >> search-cycle-q2-k9.out
