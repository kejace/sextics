#!/bin/bash
# high-degree q=2 searches: kernel d=10 (all sectors), product/square relaxation at d=12, kernel d=11/12 (Y-even sectors only, few restarts)
cd /tmp/hodge/results/cycle-q2
python3 w2_search.py --d 10 --family kernel --restarts 12 --batch 6 --maxiter 60 --nproc 32 --tag k10 > search-cycle-q2-k10.out 2>&1
python3 w2_search.py --d 12 --family prod6 --restarts 12 --batch 6 --maxiter 60 --nproc 32 --tag prod6_12 > search-cycle-q2-prod6.out 2>&1
python3 w2_search.py --d 11 --family kernel --restarts 8 --batch 8 --maxiter 60 --nproc 32 --tag k11 > search-cycle-q2-k11.out 2>&1
python3 w2_search.py --d 12 --family kernel --restarts 6 --batch 6 --maxiter 60 --nproc 32 --tag k12 > search-cycle-q2-k12.out 2>&1
echo HIGH-DONE >> search-cycle-q2-k12.out
