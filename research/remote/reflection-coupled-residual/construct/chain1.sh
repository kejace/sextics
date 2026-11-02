#!/bin/bash
# chain 1: reflection searches (hidden exchange on X / pencil fibres), then Gauss-Manin rank (consistency check: predicted 186)
cd /tmp/hodge/results/reflection-coupled-residual/construct
echo "=== refl_all start $(date)"
bash refl_all.sh > refl_all.out 2>&1
echo "=== gm_rank start $(date)"
python3 gm_rank.py gmrank_K70 67108859 70 2 4 > gmrank_K70.out 2>&1
echo "=== chain1 done $(date)"
