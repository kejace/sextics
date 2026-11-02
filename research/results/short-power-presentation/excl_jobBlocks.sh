#!/bin/bash
# wave-2 excl lens (instance B): block elimination (R5) at d=7 (e,o,m; with exact certificates) and d=8 (e,o)
cd /tmp/hodge/results/short-power-presentation
export W2JOB=exclB-spp-blocks
for spec in "7 e --exact" "7 o --exact" "7 m --exact" "8 e" "8 o"; do
  echo "=== blocks $spec $(date)"
  python3 excl_blocks.py $spec --scheme comp --threads 8 2>&1 | grep -v "hodge engine" | grep -v "^LOG"
done
echo "=== done $(date)"
