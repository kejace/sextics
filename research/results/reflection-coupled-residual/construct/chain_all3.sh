#!/bin/bash
# chain 3: v-only (lambda eliminated) lean reflection search, then Fermat continuation base Y
cd /tmp/hodge/results/reflection-coupled-residual/construct
export MSOLVE_THREADS=6
echo "=== refl_all (elim) start $(date)"
REFL_CHARTS=2 REFL_K1_ONLY=1 REFL_ELIM=1 bash refl_all.sh > refl_all_elim.out 2>&1
echo "=== fermatcont Y start $(date)"
python3 fermat_cont.py fermatcont_Y 2147483647 Y 5 > fermatcont_Y.out 2>&1
echo "=== chain_all3 done $(date)"
