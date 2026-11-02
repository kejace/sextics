#!/bin/bash
cd /tmp/hodge/results/cyclic-triple-residual/wave2-excl
for c in "0,0" "0,2" "0,4" "2,0" "2,2" "2,4" "4,0" "4,2" "4,4"; do
  tag="d8_sector_${c/,/}"
  python3 b2_hessian.py 8 --cD sector:$c --kernel kernel_d8.json --msolve --threads 2 --nocheck3 --tag $tag > b2_$tag.out 2>&1
done
echo ALLDONE
