#!/bin/bash
# Job excl-parity8: the 8 parity classes at degree 8 (all six h1-characters mixed inside a class), msolve -g 1 with 16 threads.
cd /tmp/hodge/results/quadric-contact-residual/wave2-excl
export HODGE_W2LOG=/tmp/hodge/results/quadric-contact-residual/wave2-excl/excl-parity8.log.jsonl
MSOLVE_TIMEOUT=7200 python3 parity_closure.py 8 all 16 2>&1 | grep -v "engine v2" | tail -12
echo PARITY8 DONE
