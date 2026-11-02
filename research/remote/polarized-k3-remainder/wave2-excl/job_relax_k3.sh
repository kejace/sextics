#!/bin/bash
# NEW: the same relaxation for the K3-sector target T_K3 = U0^2U1^2U2^2 - V0^2V1^2V2^2 at d = 5, 6, 7.
cd /tmp/hodge/results/polarized-k3-remainder/wave2-excl
export PYTHONPATH=/home/kejace/devel/sextics-work
P=65521
for d in 5 6; do
  echo "=== K3 target d=$d p=$P exact constraints + squares"; date
  python3 relax.py $d $P --target=k3 --exact-constraints --squares --out=relax-k3-d$d-p$P 2>&1 | grep -v "engine v2"
done
echo "=== K3 target d=7 p=$P exact L-ranks + modular squares"; date
python3 relax.py 7 $P --target=k3 --exact-ranks --squares --out=relax-k3-d7-p$P 2>&1 | grep -v "engine v2"
echo "=== DONE"; date
