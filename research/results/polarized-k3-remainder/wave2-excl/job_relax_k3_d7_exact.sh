#!/bin/bash
# Exact rational constraints for the degree-7 K3-target relaxation (closes the Q-vs-F_p gap of the modular run):
# augmented fmpq rref per character block (1926 x (3870 + 7260)), exact square elimination, residual ideal.
cd /tmp/hodge/results/polarized-k3-remainder/wave2-excl
export PYTHONPATH=/home/kejace/devel/sextics-work
echo "=== K3 d=7 exact constraints $(date)"
python3 relax.py 7 65521 --target=k3 --exact-constraints --squares --out=relax-k3-d7-p65521-exact 2>&1 | grep -v "engine v2" | grep -v "^\[result\]"
echo "=== survivors (modular residual ideal, Singular) $(date)"
python3 survivors.py relax-k3-d7-p65521-exact --primes 2>&1 | tail -40
echo "=== DONE $(date)"
