#!/bin/bash
# wait until chain A finished e1spaces d=10 (both parities), then stop chain A and submit chain A2
cd /home/kejace/devel/sextics-work
until timeout 90 python3 -m hodge.rjob sh -- 'grep -c "\"d\": 10, \"parity\": 1" /tmp/hodge/results/rank-one-biquartic/wave2/out/e1spaces_d4_10.out' 2>/dev/null | grep -q "^1"; do sleep 120; done
echo "e1spaces d10 finished at $(date); stopping chain A"
timeout 90 python3 -m hodge.rjob sh -- 'pg=$(ps -eo pgid,args | grep "bash wave2/chainA.sh" | grep -v grep | awk "{print \$1}" | head -1); echo "chainA pgid $pg"; [ -n "$pg" ] && kill -TERM -- -$pg; sleep 5; [ -n "$pg" ] && kill -KILL -- -$pg 2>/dev/null; ps -eo pgid,args | grep -c "[c]hainA.sh"'
sleep 5
python3 -m hodge.rjob submit --lane rank-one-biquartic --name search-rank-one-biquartic-chainA2 --cpus 8 --force -- 'bash wave2/chainA2.sh > wave2/out/chainA2.log 2>&1'
echo "chainA2 submitted at $(date)"
