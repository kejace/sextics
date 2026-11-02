#!/bin/bash
# pull + merge + report every 20 minutes for up to 14 hours (keeps the lane log and report current after the session)
end=$(( $(date +%s) + 14*3600 ))
while [ $(date +%s) -lt $end ]; do
  bash /tmp/hodge/results/symmetric-sextic-residual/w2b/sync.sh >> /tmp/hodge/results/symmetric-sextic-residual/w2b/autosync.log 2>&1
  python3 /tmp/hodge/results/symmetric-sextic-residual/w2b/final_summary.py > /tmp/hodge/results/symmetric-sextic-residual/w2b/summary-latest.txt 2>/dev/null
  sleep 1200
done
