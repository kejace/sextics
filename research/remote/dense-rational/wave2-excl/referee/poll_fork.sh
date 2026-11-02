#!/bin/bash
# poll the fork's research head until README-referee.md is present there (max ~25 min)
for i in $(seq 1 50); do
  head=$(curl -s -m 30 https://api.github.com/repos/kejace/sextics/branches/research | python3 -I -c 'import json,sys; print(json.load(sys.stdin).get("commit",{}).get("sha",""))' 2>/dev/null)
  code=$(curl -s -m 30 -o /dev/null -w "%{http_code}" "https://api.github.com/repos/kejace/sextics/contents/research/results/dense-rational/wave2-excl/referee/referee-report.md?ref=$head")
  echo "$(date -u +%H:%M:%SZ) head=$head referee-report: HTTP $code"
  [ "$code" = "200" ] && { echo "FOUND $head"; exit 0; }
  sleep 30
done
echo "TIMEOUT"; exit 1
