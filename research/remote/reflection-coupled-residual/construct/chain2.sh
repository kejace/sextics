#!/bin/bash
# chain 2: Fermat-pencil continuation (bases X and Y), then critical-point classification of the (Y0:Y1) pencil
cd /tmp/hodge/results/reflection-coupled-residual/construct
echo "=== fermatcont X start $(date)"
python3 fermat_cont.py fermatcont_X 2147483647 X 5 > fermatcont_X.out 2>&1
echo "=== fermatcont Y start $(date)"
python3 fermat_cont.py fermatcont_Y 2147483647 Y 5 > fermatcont_Y.out 2>&1
echo "=== critpts start $(date)"
python3 critpts.py critpts 2147483647 > critpts.out 2>&1
echo "=== chain2 done $(date)"
