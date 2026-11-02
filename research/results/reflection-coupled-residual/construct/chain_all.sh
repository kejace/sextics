#!/bin/bash
# merged chain (one job, 6 cpus): cheap experiments first, then the reflection search, then the Gauss-Manin rank
cd /tmp/hodge/results/reflection-coupled-residual/construct
export MSOLVE_THREADS=6
echo "=== fermatcont X start $(date)"
python3 fermat_cont.py fermatcont_X 2147483647 X 5 > fermatcont_X.out 2>&1
echo "=== critpts start $(date)"
python3 critpts.py critpts 2147483647 > critpts.out 2>&1
echo "=== refl_all start $(date)"
sed -i 's/^export MSOLVE_THREADS=16/export MSOLVE_THREADS=6/; s/-t 16/-t 6/' refl_all.sh critpts.py
bash refl_all.sh > refl_all.out 2>&1
echo "=== gm_rank start $(date)"
python3 gm_rank.py gmrank_K70 67108859 70 2 4 > gmrank_K70.out 2>&1
echo "=== fermatcont Y start $(date)"
python3 fermat_cont.py fermatcont_Y 2147483647 Y 5 > fermatcont_Y.out 2>&1
echo "=== chain_all done $(date)"
