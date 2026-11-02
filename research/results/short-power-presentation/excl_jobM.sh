#!/bin/bash
# wave-2 excl lens (instance B): mixed-parity exclusion via annihilator functionals (d=6 exact, d=6 second prime, d=7)
cd /tmp/hodge/results/short-power-presentation
export W2JOB=exclB-spp-relax-mixed
echo "=== exact D7 (excl_kernel.py 7) $(date)"
python3 excl_kernel.py 7 2>&1 | grep -v "hodge engine" | tail -5
echo "=== mixed d=6 exact $(date)"
python3 excl_mixed.py 6 --threads 8 2>&1 | grep -v "hodge engine" | grep -v "^LOG"
echo "=== mixed d=6 pidx 1 $(date)"
python3 excl_mixed.py 6 --threads 8 --pidx 1 --no-exact 2>&1 | grep -v "hodge engine" | grep -v "^LOG"
echo "=== mixed d=7 $(date)"
python3 excl_mixed.py 7 --threads 8 2>&1 | grep -v "hodge engine" | grep -v "^LOG"
echo "=== done $(date)"
