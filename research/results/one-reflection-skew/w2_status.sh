#!/bin/bash
# one-shot remote status for the wave-2 jobs of lane one-reflection-skew (light ssh)
cd /home/kejace/devel/sextics-work && timeout 120 python3 -m hodge.rjob sh -- 'cd /tmp/hodge/results/one-reflection-skew; echo "load $(uptime | sed "s/.*average: //" | cut -d, -f1) $(date +%H:%M) UTC"; for j in w2-chain1b w2-chain2 w2-kernels-exact w2-chain3 w2-chain4 w2-chain5; do pid=$(cat /tmp/hodge/jobs/$j.pid 2>/dev/null); if [ -n "$pid" ] && kill -0 $pid 2>/dev/null; then echo -n "$j:RUN "; fi; done; echo; for f in *.status.json; do python3 -c "
import json,sys
s=json.load(open(\"$f\"))
if not s[\"final\"] or s[\"job\"].startswith((\"calib\",\"ker-\",\"aplus\",\"ctrl\")):
    print(s[\"job\"], \"K\",s[\"K\"], s[\"done\"],\"/\",s[\"planned\"], \"final\" if s[\"final\"] else \"RUN\", \"best %.2e\" % s[\"best_rel\"] if s[\"best_rel\"] is not None else None, \"hits\", s[\"hist\"].get(\"<1e-08\"), \"conf\", len(s[\"confirmed\"]))
" 2>/dev/null; done; echo "kernels $(grep -c exact-kernel w2-kernels-exact.out) | cones $(ls cone-*.status.json 2>/dev/null | wc -l) | struct lines $(wc -l < w2-struct.struct.jsonl 2>/dev/null)"' 2>&1 | grep -v "^$" | tail -25
