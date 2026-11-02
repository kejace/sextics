#!/bin/bash
while true; do
  for pid in $(pgrep -f "Singular -q /tmp/hodge/results/symmetric-sextic-residual/wave2-excl/relax_(sector|full)_.*_d(8|9|1[0-9])_" 2>/dev/null); do kill $pid 2>/dev/null; done
  sleep 10
done
