#!/bin/bash
# chain 2nd attempt: lean reflection search (k=1 conditions, 2 charts) then Gauss-Manin rank, then Fermat continuation base Y
cd /tmp/hodge/results/reflection-coupled-residual/construct
export MSOLVE_THREADS=6
echo "=== refl_all (lean) start $(date)"
REFL_CHARTS=2 REFL_K1_ONLY=1 bash refl_all.sh > refl_all_lean.out 2>&1
echo "=== gm_rank start $(date)"
python3 gm_rank.py gmrank_K70 67108859 70 2 4 > gmrank_K70.out 2>&1
echo "=== fermatcont Y start $(date)"
python3 fermat_cont.py fermatcont_Y 2147483647 Y 5 > fermatcont_Y.out 2>&1
echo "=== chain_all2 done $(date)"
