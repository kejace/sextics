#!/usr/bin/env python3
"""Remote job runner: offload heavy computations to the 128-core gcloud box `hodge-sextics`.

The remote mirrors this machine's layout: engine at ~/devel/sextics-work, lane dirs at /tmp/hodge/results/<lane>,
repo at ~/devel/sextics (Lean built). Jobs run detached with nohup; outputs stay on the remote until pulled.

CLI (run from anywhere):
  python3 -m hodge.rjob up                       # push engine + all lane dirs to the remote (rsync, fast)
  python3 -m hodge.rjob up --lane <lane>          # push engine + one lane dir
  python3 -m hodge.rjob submit --lane <lane> --name <job> [--cpus N] [--wait SECS] [--force] -- <command ...>
        # refuses (exit 3) when the remote is near capacity (allocated cores > 220, < 90 GB free, or load > 700; jobs are core-pinned so load overstates contention); --wait polls for capacity
        # runs `cd /tmp/hodge/results/<lane> && <command>` on the remote, log in /tmp/hodge/jobs/<job>.log
        # the command may use python3 with the engine on sys.path (/home/kejace/devel/sextics-work), msolve, Singular
  python3 -m hodge.rjob status [--name <job>]     # running jobs, load average, free memory
  python3 -m hodge.rjob log --name <job> [--tail 50]
  python3 -m hodge.rjob wait --name <job> [--timeout 3600]
  python3 -m hodge.rjob pull --lane <lane>        # rsync the remote lane dir into /tmp/hodge/remote/<lane>/ and
                                                  # append new JSONL lines of its log.jsonl to the local lane log
  python3 -m hodge.rjob sh -- <command>           # run a short command on the remote and print its output

Conventions for agents: every job is pinned to its --cpus cores (taskset) and capped at max(64 GB, 8 GB per cpu) address space; request --cpus honestly (multiprocessing pools must not exceed it). Keep local processes light (<1 CPU-minute); anything heavier goes through `submit`.
Write per-job outputs to files named after the job inside the lane dir so pulls never collide.
"""
from __future__ import annotations
import argparse, json, os, shlex, subprocess, sys, time
from pathlib import Path

HOST = os.environ.get('HODGE_REMOTE_HOST', '35.229.16.59')
USER = os.environ.get('HODGE_REMOTE_USER', 'kejace')
KEY = os.path.expanduser(os.environ.get('HODGE_REMOTE_KEY', '~/.ssh/google_compute_engine'))
SSH_OPTS = ['-i', KEY, '-o', 'StrictHostKeyChecking=accept-new', '-o', 'ConnectTimeout=20', '-o', 'ServerAliveInterval=30']
ENGINE = Path('/home/kejace/devel/sextics-work')
RESULTS = Path('/tmp/hodge/results')
MIRROR = Path('/tmp/hodge/remote')
REGISTRY = Path('/tmp/hodge/jobs-registry.jsonl')


def ssh(cmd: str, check=False, capture=True, timeout=None):
    full = ['ssh', *SSH_OPTS, f'{USER}@{HOST}', cmd]
    return subprocess.run(full, check=check, capture_output=capture, text=True, timeout=timeout)


def rsync(src, dst, extra=()):
    cmd = ['rsync', '-az', '--exclude', '__pycache__', '--exclude', 'venv', '--exclude', '*.olean', '-e',
           'ssh ' + ' '.join(shlex.quote(o) for o in SSH_OPTS), *extra, src, dst]
    return subprocess.run(cmd, capture_output=True, text=True)


def up(lane=None):
    ssh('mkdir -p /tmp/hodge/results /tmp/hodge/jobs /tmp/hodge/notes ~/devel/sextics-work')
    r = rsync(str(ENGINE) + '/', f'{USER}@{HOST}:{ENGINE}/')
    if r.returncode: print(r.stderr[-500:], file=sys.stderr)
    lanes = [lane] if lane else [p.name for p in RESULTS.iterdir() if p.is_dir()]
    for ln in lanes:
        src = RESULTS / ln
        if src.is_dir():
            r = rsync(str(src) + '/', f'{USER}@{HOST}:{RESULTS}/{ln}/', extra=('--exclude', '*.npz', '--max-size=50m'))
            if r.returncode: print(ln, r.stderr[-300:], file=sys.stderr)
    # notes too (dossiers etc.) so remote scripts can read them
    rsync('/tmp/hodge/notes/', f'{USER}@{HOST}:/tmp/hodge/notes/')
    print('up: done', lanes if lane else f'{len(lanes)} lanes')



def remote_capacity():
    """(load1, avail_gb, running_jobs, allocated_cpus) on the remote."""
    r = ssh('cut -d" " -f1 /proc/loadavg; free -g | awk "/Mem/{print \\$7}"; n=0; c=0; for p in /tmp/hodge/jobs/*.pid; do if kill -0 $(cat $p) 2>/dev/null; then n=$((n+1)); f=${p%.pid}.cpus; c=$((c + $(cat $f 2>/dev/null || echo 1))); fi; done; echo $n $c')
    try:
        load, avail, n, c = r.stdout.split()[:4]
        return float(load), float(avail), int(n), int(c)
    except Exception:
        return 0.0, 1e9, 0, 0


def submit(lane, name, command, cpus=1, force=False, max_load=700.0, min_avail_gb=90.0, wait_s=0, max_alloc_cpus=220):
    """Submit, but refuse (exit code 3) when the remote is near capacity unless --force.
    With wait_s > 0, poll every 60 s up to wait_s seconds for capacity instead of refusing."""
    t0 = time.time()
    while not force:
        load, avail, n, alloc = remote_capacity()
        if load <= max_load and avail >= min_avail_gb and alloc + cpus <= max_alloc_cpus:
            break
        msg = f'REMOTE BUSY: load {load:.0f} (max {max_load:.0f}), available memory {avail:.0f} GB (min {min_avail_gb:.0f}), {n} jobs running with {alloc} cpus allocated (max {max_alloc_cpus})'
        if wait_s and time.time() - t0 < wait_s:
            print(msg + '; waiting 60 s', file=sys.stderr); time.sleep(60); continue
        print(msg + '; NOT submitted. Retry later (or pass --wait 1800 to wait for capacity).', file=sys.stderr)
        sys.exit(3)
    up(lane)
    # Hard caps enforced on the remote: CPU affinity to `cpus` cores (children inherit), address-space limit
    # 64 GB (or 8 GB per requested cpu if larger), so a runaway multiprocessing job cannot take the whole box.
    import random as _rnd
    start = _rnd.randint(0, max(0, 128 - cpus))
    cpuset = f'{start}-{start + cpus - 1}' if cpus > 1 else f'{start}'
    mem_kb = max(64, 8 * cpus) * 1024 * 1024
    remote_cmd = (f'mkdir -p /tmp/hodge/jobs {RESULTS}/{lane} && cd {RESULTS}/{lane} && '
                  f'export PYTHONPATH={ENGINE} OMP_NUM_THREADS={cpus} OPENBLAS_NUM_THREADS={cpus} FLINT_NUM_THREADS={cpus} HODGE_JOB_CPUS={cpus} && '
                  f'echo {cpus} > /tmp/hodge/jobs/{name}.cpus && '
                  f'nohup setsid taskset -c {cpuset} bash -c {shlex.quote("ulimit -v " + str(mem_kb) + "; " + command)} > /tmp/hodge/jobs/{name}.log 2>&1 < /dev/null & '
                  f'echo $! > /tmp/hodge/jobs/{name}.pid; echo $!')
    r = ssh(remote_cmd)
    pid = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else '?'
    with REGISTRY.open('a') as f:
        f.write(json.dumps({'name': name, 'lane': lane, 'pid': pid, 'cmd': command, 'ts': time.strftime('%Y-%m-%dT%H:%M:%S')}) + '\n')
    print(f'submitted {name} (lane {lane}) pid {pid}; log: /tmp/hodge/jobs/{name}.log on remote')
    if r.stderr.strip(): print(r.stderr[-300:], file=sys.stderr)
    return pid


def status(name=None):
    if name:
        r = ssh(f'pid=$(cat /tmp/hodge/jobs/{name}.pid 2>/dev/null); if [ -n "$pid" ] && kill -0 $pid 2>/dev/null; then echo RUNNING; else echo DONE; fi; tail -3 /tmp/hodge/jobs/{name}.log 2>/dev/null')
    else:
        r = ssh('uptime; free -g | head -2; echo "--- jobs (name cpus)"; for p in /tmp/hodge/jobs/*.pid; do n=$(basename $p .pid); pid=$(cat $p); if kill -0 $pid 2>/dev/null; then echo "RUNNING $n $(cat ${p%.pid}.cpus 2>/dev/null)"; fi; done; echo "--- top"; ps aux --sort=-%cpu | awk "NR>1 && NR<12 {print \\$3, \\$4, \\$11, \\$12, \\$13}"')
    print(r.stdout, r.stderr[-300:] if r.returncode else '')


def log(name, tail=50):
    r = ssh(f'tail -n {tail} /tmp/hodge/jobs/{name}.log 2>/dev/null')
    print(r.stdout)


def wait(name, timeout=3600, poll=30):
    t0 = time.time()
    while time.time() - t0 < timeout:
        r = ssh(f'pid=$(cat /tmp/hodge/jobs/{name}.pid 2>/dev/null); if [ -n "$pid" ] && kill -0 $pid 2>/dev/null; then echo RUNNING; else echo DONE; fi')
        if 'DONE' in r.stdout:
            print(f'{name}: done after {time.time()-t0:.0f}s'); return True
        time.sleep(poll)
    print(f'{name}: still running after {timeout}s'); return False


def pull(lane):
    dst = MIRROR / lane
    dst.mkdir(parents=True, exist_ok=True)
    r = rsync(f'{USER}@{HOST}:{RESULTS}/{lane}/', str(dst) + '/', extra=('--exclude', '*.npz', '--max-size=50m'))
    if r.returncode: print(r.stderr[-300:], file=sys.stderr)
    rsync(f'{USER}@{HOST}:/tmp/hodge/jobs/', str(MIRROR / '_jobs') + '/')
    # merge JSONL logs: append lines of the remote log not present locally
    rl, ll = dst / 'log.jsonl', RESULTS / lane / 'log.jsonl'
    added = 0
    if rl.exists():
        have = set(ll.read_text().splitlines()) if ll.exists() else set()
        with ll.open('a') as f:
            for line in rl.read_text().splitlines():
                if line and line not in have:
                    f.write(line + '\n'); have.add(line); added += 1
    print(f'pull {lane}: mirror at {dst}; {added} new log lines merged into {ll}; job logs in {MIRROR}/_jobs/')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('up'); p.add_argument('--lane')
    p = sub.add_parser('submit'); p.add_argument('--lane', required=True); p.add_argument('--name', required=True); p.add_argument('--cpus', type=int, default=1); p.add_argument('--force', action='store_true'); p.add_argument('--wait', type=int, default=0, help='seconds to wait for remote capacity'); p.add_argument('command', nargs=argparse.REMAINDER)
    p = sub.add_parser('status'); p.add_argument('--name')
    p = sub.add_parser('log'); p.add_argument('--name', required=True); p.add_argument('--tail', type=int, default=50)
    p = sub.add_parser('wait'); p.add_argument('--name', required=True); p.add_argument('--timeout', type=int, default=3600)
    p = sub.add_parser('pull'); p.add_argument('--lane', required=True)
    p = sub.add_parser('sh'); p.add_argument('command', nargs=argparse.REMAINDER)
    a = ap.parse_args()
    if a.cmd == 'up': up(a.lane)
    elif a.cmd == 'submit':
        cmd = a.command[1:] if a.command and a.command[0] == '--' else a.command
        submit(a.lane, a.name, ' '.join(cmd), a.cpus, force=a.force, wait_s=a.wait)
    elif a.cmd == 'status': status(a.name)
    elif a.cmd == 'log': log(a.name, a.tail)
    elif a.cmd == 'wait': wait(a.name, a.timeout)
    elif a.cmd == 'pull': pull(a.lane)
    elif a.cmd == 'sh':
        cmd = a.command[1:] if a.command and a.command[0] == '--' else a.command
        r = ssh(' '.join(cmd)); print(r.stdout, r.stderr[-500:] if r.returncode else '')


if __name__ == '__main__':
    main()
