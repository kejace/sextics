#!/usr/bin/env bash
# Usage: queue.sh <cpuset> <label1> [label2 ...]   -- runs run_module2.sh Cap.S_<label> sequentially on the cpuset
CPUSET="$1"; shift
for L in "$@"; do
  echo "$(date +%T) start $L" >> /tmp/hodge/lean/capacity/logs/queue-$CPUSET.log
  /tmp/hodge/lean/capacity/run_module2.sh "Cap.S_$L" "$L" "$CPUSET" >> /tmp/hodge/lean/capacity/logs/queue-$CPUSET.log 2>&1
  python3 /tmp/hodge/lean/capacity/gen/caplog.py job_done "label=$L cpuset=$CPUSET $(tail -1 /tmp/hodge/lean/capacity/results/timings.jsonl | python3 -c "import sys,json; r=json.loads(sys.stdin.read()); print('rc=',r.get('compile_rc'),'wall=',r.get('compile',{}).get('wall'),'user=',r.get('compile',{}).get('user'),'rss=',r.get('compile',{}).get('maxrss_gb'),'krc=',r.get('kernel_rc'),'kwall=',r.get('kernel',{}).get('wall'),'kuser=',r.get('kernel',{}).get('user'),'load=',r.get('load0'))" 2>/dev/null)"
done
