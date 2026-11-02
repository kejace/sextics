#!/bin/bash
# q=1 sibling (F, B = X0*X1**2*Y0*Y1**2): a rational certificate family in lambda for F_lambda = F + (lambda-1/2) X0*X2**3*Y0*Y2
# differentiates to a q=2 certificate for A = M1*B with denominator D^2 (dossier P2 / A.7). Search the sibling at d = 8..12.
cd /tmp/hodge/results/cycle-q2
T='X0*X1**2*Y0*Y1**2'
python3 w2_search.py --d 8  --family kernel --target "$T" --q 1 --restarts 150 --batch 30 --maxiter 200 --nproc 8 --tag sibB8  > search-cycle-q2-sibB8.out 2>&1
python3 w2_search.py --d 9  --family kernel --target "$T" --q 1 --restarts 100 --batch 25 --maxiter 200 --nproc 8 --tag sibB9  > search-cycle-q2-sibB9.out 2>&1
python3 w2_search.py --d 10 --family kernel --target "$T" --q 1 --restarts 60  --batch 20 --maxiter 200 --nproc 8 --tag sibB10 > search-cycle-q2-sibB10.out 2>&1
python3 w2_search.py --d 11 --family kernel --target "$T" --q 1 --restarts 40  --batch 20 --maxiter 200 --nproc 8 --tag sibB11 > search-cycle-q2-sibB11.out 2>&1
python3 w2_search.py --d 12 --family kernel --target "$T" --q 1 --restarts 30  --batch 15 --maxiter 200 --nproc 8 --tag sibB12 > search-cycle-q2-sibB12.out 2>&1
echo SIBLING-DONE >> search-cycle-q2-sibB12.out
