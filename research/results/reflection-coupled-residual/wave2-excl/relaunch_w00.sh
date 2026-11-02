#!/bin/bash
cd /tmp/hodge/results/reflection-coupled-residual/wave2-excl
while ! grep -q 'hess_d7_W00_p2147483587.ms' gbcheck.jsonl 2>/dev/null; do sleep 15; done
if grep 'hess_d7_W00_p2147483587.ms' gbcheck.jsonl | grep -q TIMEOUT; then
  nice -n 10 timeout 4000 python3 gbcheck.py hess_d7_W00_p2147483587.ms 4 3900 12 >> gbcheck_d7_W00_long.out 2>&1
fi
echo RELAUNCH DONE >> relaunch_w00.log
