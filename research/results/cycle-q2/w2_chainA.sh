#!/bin/bash
cd /tmp/hodge/results/cycle-q2
python3 w2_search.py --d 6 --family kernel --restarts 40 --batch 10 --maxiter 60 --nproc 32 --tag k6 > search-cycle-q2-k6.out 2>&1
python3 w2_search.py --d 6 --family kernel --keys 0 --restarts 400 --split 32 --batch 13 --maxiter 60 --nproc 32 --tag k6triv > search-cycle-q2-k6triv.out 2>&1
python3 w2_search.py --d 7 --family kernel --restarts 30 --batch 10 --maxiter 60 --nproc 32 --tag k7 > search-cycle-q2-k7.out 2>&1
python3 w2_search.py --d 8 --family kernel --restarts 15 --batch 5 --maxiter 60 --nproc 32 --tag k8 > search-cycle-q2-k8.out 2>&1
python3 w2_search.py --d 9 --family kernel --restarts 8 --batch 4 --maxiter 60 --nproc 32 --tag k9 > search-cycle-q2-k9.out 2>&1
echo CHAINA-DONE >> search-cycle-q2-k9.out
