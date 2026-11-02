#!/bin/bash
# Aggregator-only job submission that does NOT rsync /tmp/hodge/notes to the remote (avoids the rjob.up() clobber).
# usage: rsubmit.sh <name> <cpus> <command...>     (command runs in /tmp/hodge/results/_aggregator on the remote)
set -e
NAME="$1"; CPUS="$2"; shift 2
CMD="$*"
HOST=kejace@35.229.16.59
KEY=/home/kejace/.ssh/google_compute_engine
SSH="ssh -i $KEY -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20"
rsync -az --exclude '__pycache__' --exclude 'venv' -e "$SSH" /home/kejace/devel/sextics-work/hodge/ $HOST:/home/kejace/devel/sextics-work/hodge/
rsync -az --exclude '*.npz' --max-size=50m -e "$SSH" /tmp/hodge/results/_aggregator/ $HOST:/tmp/hodge/results/_aggregator/
$SSH $HOST "mkdir -p /tmp/hodge/jobs && cd /tmp/hodge/results/_aggregator && export PYTHONPATH=/home/kejace/devel/sextics-work OMP_NUM_THREADS=$CPUS OPENBLAS_NUM_THREADS=$CPUS FLINT_NUM_THREADS=$CPUS && nohup setsid bash -c '$CMD' > /tmp/hodge/jobs/$NAME.log 2>&1 < /dev/null & echo \$! > /tmp/hodge/jobs/$NAME.pid; echo submitted $NAME pid \$(cat /tmp/hodge/jobs/$NAME.pid)"
echo "{\"name\": \"$NAME\", \"lane\": \"_aggregator\", \"cmd\": \"$CMD\", \"ts\": \"$(date +%FT%T)\", \"via\": \"rsubmit.sh\"}" >> /tmp/hodge/jobs-registry.jsonl
