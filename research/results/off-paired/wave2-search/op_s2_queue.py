#!/usr/bin/env python3
"""Local queue runner for lane off-paired / lens search (unique name op_s2_queue.py; never kill by generic pattern).
queue.txt lines: <name>|<cpus>|<command>. Keeps <= MAXRUN jobs with prefix search-off-paired-s2- running.
A job is considered submitted when its remote pid file exists (robust to hanging ssh)."""
import os, sys, time, subprocess
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
PREFIX = 'search-off-paired-s2-'
MAXRUN = int(os.environ.get('W2_MAXRUN', '4'))
ENGINE = '/home/kejace/devel/sextics-work'

def rjob(*args, timeout=600):
    r = subprocess.run([sys.executable, '-m', 'hodge.rjob', *args], cwd=ENGINE, capture_output=True, text=True, timeout=timeout)
    return r.stdout + r.stderr

def my_running():
    out = rjob('status', timeout=180)
    return [l.split()[1] for l in out.splitlines() if l.startswith('RUNNING ') and PREFIX in l]

def pid_exists(name):
    out = rjob('sh', '--', f'test -f /tmp/hodge/jobs/{PREFIX}{name}.pid && echo PIDEXISTS', timeout=120)
    return 'PIDEXISTS' in out

def pop_line(line):
    rest = [l.rstrip('\n') for l in open('queue.txt')]
    rest = [l for l in rest if l.strip() != line.strip()]
    open('queue.txt', 'w').write('\n'.join(rest) + ('\n' if rest else ''))
    open('queue.done', 'a').write(time.strftime('%Y-%m-%dT%H:%M:%S') + ' ' + line.strip() + '\n')

while True:
    q = [l.strip() for l in open('queue.txt')] if os.path.exists('queue.txt') else []
    q = [l for l in q if l and not l.startswith('#')]
    if not q:
        print(time.strftime('%H:%M:%S'), 'queue empty; sleeping', flush=True); time.sleep(180); continue
    try:
        running = my_running()
    except Exception as ex:
        print(time.strftime('%H:%M:%S'), 'status failed', ex, flush=True); time.sleep(120); continue
    if len(running) >= MAXRUN:
        print(time.strftime('%H:%M:%S'), f'{len(running)} running: {running}; waiting', flush=True); time.sleep(150); continue
    line = q[0]
    name, cpus, cmd = line.split('|', 2)
    try:
        if pid_exists(name):
            print(time.strftime('%H:%M:%S'), f'{name}: pid file exists on remote -> already submitted; popping', flush=True)
            pop_line(line); continue
    except Exception as ex:
        print('pid check failed', ex, flush=True); time.sleep(60); continue
    print(time.strftime('%H:%M:%S'), f'submitting {name} (cpus {cpus}); running now: {running}', flush=True)
    try:
        out = rjob('submit', '--lane', 'off-paired', '--name', PREFIX + name, '--cpus', cpus, '--force', '--', cmd, timeout=600)
        print(out.strip()[-300:], flush=True)
    except subprocess.TimeoutExpired:
        print(time.strftime('%H:%M:%S'), 'submit call timed out (ssh hang); verifying via pid file', flush=True)
    time.sleep(20)
    try:
        ok = pid_exists(name)
    except Exception as ex:
        print('verification failed', ex, flush=True); ok = False
    if ok:
        pop_line(line); print(time.strftime('%H:%M:%S'), f'{name} confirmed submitted', flush=True); time.sleep(30)
    else:
        print(time.strftime('%H:%M:%S'), f'{name} NOT submitted; retry after 120 s', flush=True); time.sleep(120)
