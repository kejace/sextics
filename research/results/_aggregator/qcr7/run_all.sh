#!/bin/bash
cd /tmp/hodge/results/_aggregator/qcr7
cls="000 001 010 011 100 101 110 111"
for a in $cls; do for b in $cls; do
  if [[ "$a" < "$b" ]]; then python3 qcr7.py $a $b 998244353 300; fi
done; done
echo ALLDONE
