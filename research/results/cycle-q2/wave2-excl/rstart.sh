#!/bin/bash
# start a detached remote job WITHOUT rsyncing the lane dir (avoids clobbering live remote files).
# usage: rstart.sh <jobname> <cpus> '<command run inside /tmp/hodge/results/cycle-q2>'
name="$1"; cpus="$2"; cmd="$3"
cd /home/kejace/devel/sextics-work && python3 - "$name" "$cpus" "$cmd" <<'PY'
import sys, shlex, json, time
sys.path.insert(0, '/home/kejace/devel/sextics-work')
from hodge import rjob
name, cpus, cmd = sys.argv[1], sys.argv[2], sys.argv[3]
remote_cmd = (f'mkdir -p /tmp/hodge/jobs && cd /tmp/hodge/results/cycle-q2 && '
              f'export PYTHONPATH=/home/kejace/devel/sextics-work OMP_NUM_THREADS={cpus} OPENBLAS_NUM_THREADS={cpus} FLINT_NUM_THREADS={cpus} && '
              f'nohup setsid bash -c {shlex.quote(cmd)} > /tmp/hodge/jobs/{name}.log 2>&1 < /dev/null & '
              f'echo $! > /tmp/hodge/jobs/{name}.pid; echo $!')
r = rjob.ssh(remote_cmd)
print('started', name, 'pid', r.stdout.strip().splitlines()[-1] if r.stdout.strip() else '?', r.stderr[-200:])
with open('/tmp/hodge/jobs-registry.jsonl', 'a') as f:
    f.write(json.dumps({'name': name, 'lane': 'cycle-q2', 'cmd': cmd, 'ts': time.strftime('%Y-%m-%dT%H:%M:%S'), 'via': 'rstart'}) + '\n')
PY
