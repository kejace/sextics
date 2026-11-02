#!/bin/bash
cd /tmp/hodge/results/reflection-coupled-residual/wave2-excl
while pgrep -f "msolve -f hess_d7_W00_p214748358[7]" > /dev/null; do sleep 20; done
if ! grep 'hess_d7_W00_p2147483587.ms' gbcheck.jsonl | grep -q zero_dim; then
  nice -n 10 timeout 28800 python3 gbcheck.py hess_d7_W00_p2147483587.ms 6 28000 12 >> gbcheck_d7_W00_verylong.out 2>&1
fi
echo RELAUNCH2 DONE >> relaunch_w00.log
