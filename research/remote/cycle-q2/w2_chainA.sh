#!/bin/bash
cd /tmp/hodge/results/cycle-q2
python3 w2_search.py --d 7 --family kernel --restarts 10 --batch 10 --maxiter 60 --nproc 48 --tag k7 --exclude-keys 6,27,36,57,87,102,117,132,147,162,192,213,216,237,267,282,297,312,327,342,372,393,402,423,447,462,477,492,507,522,552,573 > search-cycle-q2-k7b.out 2>&1
python3 w2_search.py --d 8 --family kernel --restarts 5 --batch 5 --maxiter 60 --nproc 48 --tag k8 > search-cycle-q2-k8.out 2>&1
python3 w2_search.py --d 9 --family kernel --restarts 3 --batch 3 --maxiter 60 --nproc 48 --tag k9 > search-cycle-q2-k9.out 2>&1
echo CHAINA-DONE >> search-cycle-q2-k9.out
