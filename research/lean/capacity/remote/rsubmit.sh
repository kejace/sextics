#!/usr/bin/env bash
# rsubmit.sh <jobname> <command...>  — run a detached job on the remote box (log: /tmp/hodge/jobs/<jobname>.log)
NAME="$1"; shift
CMD="$*"
ssh -n -i "$HOME/.ssh/google_compute_engine" -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20 kejace@35.229.16.59 \
  "mkdir -p /tmp/hodge/jobs && export PYTHONPATH=/home/kejace/devel/sextics-work && nohup setsid bash -c $(printf '%q' "$CMD") > /tmp/hodge/jobs/$NAME.log 2>&1 < /dev/null & echo \$! > /tmp/hodge/jobs/$NAME.pid; echo submitted $NAME pid \$!"
echo "{\"ts\":\"$(date +%FT%T)\",\"step\":\"rsubmit\",\"job\":\"$NAME\",\"cmd\":$(python3 -c 'import json,sys;print(json.dumps(sys.argv[1]))' "$CMD")}" >> /tmp/hodge/lean/capacity/log.jsonl
