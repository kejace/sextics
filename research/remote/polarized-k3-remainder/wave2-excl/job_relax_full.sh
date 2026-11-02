#!/bin/bash
# Re-verification of the thread's degree 6/7 exclusions with an independent implementation at a second prime.
cd /tmp/hodge/results/polarized-k3-remainder/wave2-excl
export PYTHONPATH=/home/kejace/devel/sextics-work
P=65521
echo "=== d=6 p=$P exact constraints + squares (full target)"; date
python3 relax.py 6 $P --exact-constraints --squares --out=relax-full-d6-p$P 2>&1 | grep -v "engine v2"
echo "=== d=6 p=101 modular cross-check of block ranks against the thread"; date
python3 relax.py 6 101 --exact-ranks --squares --out=relax-full-d6-p101 2>&1 | grep -v "engine v2"
echo "=== d=7 p=$P with exact L-ranks (flatness) + modular squares"; date
python3 relax.py 7 $P --exact-ranks --squares --out=relax-full-d7-p$P 2>&1 | grep -v "engine v2"
echo "=== DONE"; date
