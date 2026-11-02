#!/bin/bash
cd /tmp/hodge/results/symmetric-sextic-residual
while pgrep -f "local_chain1.sh" > /dev/null; do sleep 30; done
bash wave2-excl/local_chain2.sh
bash wave2-excl/local_chain3.sh
echo ALL LOCAL CHAINS DONE
