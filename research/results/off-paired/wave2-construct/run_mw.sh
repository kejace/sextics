#!/bin/bash
cd /tmp/hodge/results/off-paired/wave2-construct
for e in -1 +1; do
  timeout 1500 msolve -t 2 -v 1 -f mw_sections_eps${e}.ms -o mw_sections_eps${e}.out > mw_sections_eps${e}.log 2>&1
  echo "eps=$e exit=$? out=$(head -c 300 mw_sections_eps${e}.out | tr '\n' ' ')"
done
