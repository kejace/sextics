#!/bin/bash
# Degree-8 L-relaxation: (i) K3-sector target (new), (ii) full target regenerated at a second prime (adversarial check of
# the thread's rank 1050 / nonempty claim).  Each run ~10-30 min, ~5-10 GB.
cd /tmp/hodge/results/polarized-k3-remainder/wave2-excl
export PYTHONPATH=/home/kejace/devel/sextics-work
echo "=== K3 d=8 p=65521"; date
python3 relax.py 8 65521 --target=k3 --exact-ranks --squares --out=relax-k3-d8-p65521 2>&1 | grep -v "engine v2"
echo "=== full d=8 p=65521"; date
python3 relax.py 8 65521 --exact-ranks --squares --out=relax-full-d8-p65521 2>&1 | grep -v "engine v2"
echo "=== DONE"; date
